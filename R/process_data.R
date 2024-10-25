# Functions for pre-processing data before conduct MAIC

# Functions to be exported ---------------------------------------

#' Pre-process aggregate data
#'
#' This function checks the format of the aggregate data.
#' Data is required to have three columns: STUDY, ARM, and N.
#' Column names that do not have legal suffixes (MEAN, MEDIAN, SD, COUNT, or PROP) are dropped.
#' If a variable is a count variable, it is converted to proportions by dividing the sample size (N).
#' Note, when the count is specified, proportion is always calculated based on the count, that is,
#' specified proportion will be ignored if applicable.
#' If the aggregated data comes from multiple sources (i.e. different analysis population) and
#' sample size differs for each variable, one option is to specify proportion directly instead of count by using suffix
#' _PROP.
#'
#' @param raw_agd raw aggregate data should contain STUDY, ARM, and N. Variable names should be followed
#' by legal suffixes (i.e. MEAN, MEDIAN, SD, COUNT, or PROP).
#'
#' @examples
#' data(agd)
#' agd <- process_agd(agd)
#'
#' @return pre-processed aggregate level data
#' @export

process_agd <- function(raw_agd) {
  raw_agd <- as.data.frame(raw_agd)
  # make all column names to be capital letters to avoid different style
  names(raw_agd) <- toupper(names(raw_agd))

  # define column name patterns[-]
  must_exist <- c("STUDY", "ARM", "N")
  legal_suffix <- c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")

  # swap "TREATMENT" column to "ARM", if applicable
  if ("TREATMENT" %in% names(raw_agd) && (!"ARM" %in% names(raw_agd))) {
    raw_agd$ARM <- raw_agd$TREATMENT
    raw_agd <- raw_agd[, names(raw_agd) != "TREATMENT"]
    warning("'TREATMENT' column is renamed as 'ARM'")
  }

  # check: must exist
  if (!all(must_exist %in% names(raw_agd))) {
    stop("At least 1 of the must-exists columns (STUDY, ARM, N) cannot be found in raw_agd!")
  }

  # check: legal suffix
  other_colnames <- setdiff(names(raw_agd), must_exist)
  ind1 <- grepl("_", other_colnames, fixed = TRUE)
  ind2 <- sapply(other_colnames, function(xx) {
    tmp <- unlist(strsplit(xx, split = "_"))
    tmp[length(tmp)] # this deployment is robust to the cases that there are multiple _ in the column name
  })
  ind2 <- (ind2 %in% legal_suffix)

  use_cols <- other_colnames[ind1 & ind2]
  use_agd <- raw_agd[, c(must_exist, use_cols), drop = FALSE]
  if (!all(other_colnames %in% use_cols)) {
    warning(paste0(
      "following columns are ignored since it does not follow the naming conventions:",
      paste(setdiff(other_colnames, use_cols), collapse = ",")
    ))
  }

  # If the aggregate data is divided by different arms, calculate pooled arm statistics using
  # complete_agd function; complete statistics is specified by ARM=="Total"
  if (!"total" %in% tolower(use_agd$ARM)) {
    use_agd <- complete_agd(use_agd)
  }

  # calculate percentage columns
  ind <- grepl("_COUNT$", names(use_agd))
  if (any(ind)) {
    for (i in which(ind)) {
      tmp_prop <- use_agd[[i]] / use_agd$N
      # in case some count are not specified, but proportion are specified, copy over those proportions
      # this also means, in case count is specified, proportion is ignored even it is specified
      prop_name_i <- gsub("_COUNT$", "_PROP", names(use_agd)[i])
      if (prop_name_i %in% names(use_agd)) {
        tmp_prop[is.na(tmp_prop)] <- use_agd[is.na(tmp_prop), prop_name_i]
        names(use_agd)[names(use_agd) == prop_name_i] <- paste0(prop_name_i, "_redundant")
      }
      use_agd[[i]] <- tmp_prop
    }
    names(use_agd) <- gsub("_COUNT$", "_PROP", names(use_agd))
  }
  use_agd <- use_agd[, !grepl("_redundant$", names(use_agd))]

  # output
  with(use_agd, use_agd[tolower(ARM) == "total", , drop = FALSE])
}


#' Create dummy variables from categorical variables in an individual patient data (ipd)
#'
#' This is a convenient function to convert categorical variables into dummy binary variables.
#' This would be especially useful if the variable has more than two factors.
#' Note that the original variable is kept after a variable is dummized.
#'
#' @param raw_ipd ipd data that contains variable to dummize
#' @param dummize_cols vector of column names to binarize
#' @param dummize_ref_level vector of reference level of the variables to binarize
#'
#' @examples
#' data(adsl_twt)
#' dummize_ipd(adsl_twt, dummize_cols = c("SEX"), dummize_ref_level = c("Male"))
#'
#' @return ipd with dummized columns
#' @export

dummize_ipd <- function(raw_ipd, dummize_cols, dummize_ref_level) {
  for (i in seq_along(dummize_cols)) {
    yy <- raw_ipd[[dummize_cols[i]]]
    yy_levels <- na.omit(unique(yy))
    yy <- factor(as.character(yy), levels = c(dummize_ref_level[i], setdiff(yy_levels, dummize_ref_level[i])))
    new_yy <- sapply(levels(yy)[-1], function(j) {
      as.numeric(yy == j)
    })
    new_yy <- as.data.frame(new_yy)
    names(new_yy) <- toupper(paste(dummize_cols[i], levels(yy)[-1], sep = "_"))
    raw_ipd <- cbind(raw_ipd, new_yy)
  }
  raw_ipd
}


#' Center individual patient data (IPD) variables using aggregate data averages
#'
#' This function subtracts IPD variables (prognostic variables and/or effect modifiers)
#' by the aggregate data averages. This centering is needed in order to calculate weights.
#' IPD and aggregate data variable names should match.
#'
#' @param ipd IPD variable names should match the aggregate data names without the suffix.
#' This would involve either changing the aggregate data name or the ipd name.
#' For instance, if we binarize SEX variable with MALE as a reference using [dummize_ipd],
#' function names the new variable as SEX_MALE.
#' In this case, SEX_MALE should also be available in the aggregate data.
#' @param agd pre-processed aggregate data which contain STUDY, ARM, and N. Variable names
#' should be followed by legal suffixes (i.e. MEAN, MEDIAN, SD, or PROP). Note that COUNT
#' suffix is no longer accepted.
#' @examples
#' data(adsl_sat)
#' data(agd)
#' agd <- process_agd(agd)
#' ipd_centered <- center_ipd(ipd = adsl_sat, agd = agd)
#' @return centered ipd using aggregate level data averages
#' @export

center_ipd <- function(ipd, agd) {
  # regularized column name patterns
  must_exist <- c("STUDY", "ARM", "N")
  legal_suffix <- c("MEAN", "MEDIAN", "SD", "PROP")
  suffix_pat <- paste(paste0("_", legal_suffix, "$"), collapse = "|")

  for (i in seq_len(nrow(agd))) { # study i
    study_id <- agd$STUDY[i]
    use_agd <- agd[i, !names(agd) %in% must_exist, drop = FALSE]
    param_id <- gsub(suffix_pat, "", names(use_agd))

    for (j in seq_len(ncol(use_agd))) { # effect modifier j
      if (is.na(use_agd[[j]])) next

      ipd_param <- param_id[j]

      if (grepl("_MEAN$|_PROP$", names(use_agd)[j])) {
        ipd[[paste0(ipd_param, "_", "CENTERED")]] <- ipd[[ipd_param]] - use_agd[[j]]
      } else if (grepl("_MEDIAN$", names(use_agd)[j])) {
        ipd[[paste0(ipd_param, "_MEDIAN_", "CENTERED")]] <- ipd[[ipd_param]] > use_agd[[j]]
        ipd[[paste0(ipd_param, "_MEDIAN_", "CENTERED")]] <- ipd[[paste0(ipd_param, "_MEDIAN_", "CENTERED")]] - 0.5
      } else if (grepl("_SD$", names(use_agd)[j])) {
        ipd[[paste0(ipd_param, "_SQUARED_", "CENTERED")]] <- ipd[[ipd_param]]^2
        tmp_aim <- use_agd[[j]]^2 + (use_agd[[paste0(ipd_param, "_MEAN")]]^2)
        ipd[[paste0(ipd_param, "_SQUARED_", "CENTERED")]] <- ipd[[paste0(ipd_param, "_SQUARED_", "CENTERED")]] - tmp_aim
      }
    } # end of j
  } # end of i

  # output
  ipd
}


#' Calculate pooled arm statistics in Aggregated Data (AgD) based on arm-specific statistics
#'
#' This is a convenient function to pool arm statistics. This function is called
#' within process_agd and when the ARM is not equal to "Total". Note pooled
#' median can't be calculated and it is only an approximation.
#'
#' @param use_agd aggregated data that is processed within process_agd
#' @noRd
#' @return Complete N, count, mean, sd, and median for the pooled arm

complete_agd <- function(use_agd) {
  use_agd <- as.data.frame(use_agd)
  use_agd <- with(use_agd, {
    use_agd[tolower(ARM) != "total", , drop = FALSE]
  })

  if (nrow(use_agd) < 2) stop("error in call complete_agd: need to have at least 2 rows that ARM!='total' ")

  rowId <- nrow(use_agd) + 1
  use_agd[rowId, ] <- NA
  use_agd$STUDY[rowId] <- use_agd$STUDY[1]
  use_agd$ARM[rowId] <- "total"

  # complete N and count
  NN <- use_agd[["N"]][rowId] <- sum(use_agd[["N"]], na.rm = TRUE)
  nn <- use_agd[["N"]][-rowId]
  for (i in grep("_COUNT$", names(use_agd), value = TRUE)) {
    use_agd[[i]][rowId] <- sum(use_agd[[i]][-rowId], na.rm = TRUE)
  }

  # complete MEAN
  for (i in grep("_MEAN$", names(use_agd), value = TRUE)) {
    use_agd[[i]][rowId] <- sum(use_agd[[i]][-rowId] * nn) / NN
  }

  # complete SD
  for (i in grep("_SD$", names(use_agd), value = TRUE)) {
    use_agd[[i]][rowId] <- sqrt(sum(use_agd[[i]][-rowId]^2 * (nn - 1)) / (NN - 1))
  }

  # complete MEDIAN, approximately!!
  for (i in grep("_MEDIAN$", names(use_agd), value = TRUE)) {
    use_agd[[i]][rowId] <- mean(use_agd[[i]][-rowId])
  }

  # output
  rownames(use_agd) <- NULL
  use_agd
}


#' helper function: transform TTE ADaM data to suitable input for survival R package
#'
#' @param dd data frame, ADTTE read via `haven::read_sas`
#' @param time_scale a character string, 'years', 'months', 'weeks' or 'days', time unit of median survival time
#' @param trt values to include in treatment column
#'
#' @return a data frame that can be used as input to `survival::Surv`
#' @keywords internal

ext_tte_transfer <- function(dd, time_scale = "months", trt = NULL) {
  time_scale <- match.arg(time_scale, choices = c("years", "months", "weeks", "days"))
  time_units <- get_time_conversion(time_scale)

  if ("CENSOR" %in% names(dd)) {
    dd <- dd[!is.na(dd$CENSOR), ]
    dd$status <- as.numeric(dd$CENSOR == 0)
  }
  if ("EVENT" %in% names(dd)) {
    dd$status <- as.numeric(as.character(dd$EVENT))
  }
  if ("TIME" %in% names(dd)) {
    dd$AVAL <- as.numeric(as.character(dd$TIME))
  }

  dd$time <- dd$AVAL * time_units
  if (!is.null(trt)) dd$treatment <- trt
  as.data.frame(dd)
}
