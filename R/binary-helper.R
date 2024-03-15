#' Helper funciton: create pseudo IPD given AgD binary
#'
#' @param binary_agd
#' @param format
#'
#' @details
#' Additional details...
#'
#' @return a data.frame of pseudo binary IPD, with columns USUBJID, ARM, RESPONSE
#' @export
#'
#' @examples
get_pseudo_ipd_binary <- function(binary_agd, format = c("stacked", "unstacked", "brief")) {
  # pre check
  if (format == "stacked") {
    if (!is.data.frame(binary_agd)) stop("stacked binary_agd should be data.frame with columns 'ARM', 'RESPONSE', 'COUNT'")
    names(binary_agd) <- toupper(names(binary_agd))
    if (!all(c("ARM", "RESPONSE", "COUNT") %in% names(binary_agd))) stop("stacked binary_agd should be data.frame with columns 'ARM', 'Response', 'Count'")
    if (!is.logical(binary_agd$RESPONSE) & all(toupper(binary_agd$RESPONSE) %in% c("YES", "NO"))) stop("'RESPONSE' column in stacked binary_agd should be either logical vector or character vector of 'Yes' and 'No'")
    if (nrow(binary_agd) %% 2 != 0) stop("nrow(binary_agd) is not even number, you may miss to provide 1 level of binary response to certain arm")
  } else if (format == "unstacked") {
    if (!(is.data.frame(binary_agd) | is.matrix(binary_agd))) stop("unstacked binary_agd should be either a 2x2 data frame or matrix")
    if (!identical(dim(binary_agd), c(2, 2))) stop("unstacked binary_agd should be either a 2x2 data frame or matrix")
    bin_res <- toupper(colnames(binary_agd))
    bin_res <- sort(bin_res)
    if (!(identical(bin_res, c("FALSE", "TRUE")) | identical(bin_res, c("NO", "YES")))) stop("column names of unstacked binary_agd should be either TRUE/FALSE or Yes/No")
  } else if (format == "brief") {
    if (!is.list(binary_agd)) stop("brief binary_agd should be a R list including 'total' sample size")
    if (!"total" %in% names(binary_agd)) stop("brief binary_agd should be a R list including 'total' sample size")
  }

  # pre process binary_agd, depending on format
  use_binary_agd <- switch(format,
    "stacked" = {
      names(binary_agd) <- toupper(names(binary_agd))
      if (!is.logical(binary_agd$RESPONSE)) {
        binary_agd$RESPONSE <- toupper(binary_agd$RESPONSE)
        binary_agd$RESPONSE <- binary_agd$RESPONSE == "YES"
      }
      binary_agd
    },
    "unstacked" = {
      trt_names <- rownames(binary_agd)
      bin_res <- toupper(colnames(binary_agd))
      if ("YES" %in% bin_res) {
        bin_res <- ifelse(bin_res == "YES", "TRUE", "FALSE")
        colnames(binary_agd) <- bin_res
      }
      tmpout <- utils::stack(binary_agd)
      tmpout <- cbind(ARM = rep(trt_names, each = 2), tmpout)
      names(tmpout) <- c("ARM", "COUNT", "RESPONSE")
      rownames(tmpout) <- NULL
      tmpout$RESPONSE <- as.logical(tmpout$RESPONSE)
      tmpout
    },
    "brief" = {
      nr_arms <- length(brief) - 1
      tmpout <- data.frame(
        ARM = rep(setdiff(names(brief), "total"), each = 2),
        RESPONSE = rep(c(TRUE, FALSE), nr_arms),
        COUNT = NA
      )
      tmpout$COUNT <- unlist(brief[tmpout$ARM])
      tmpout
    }
  )

  # create IPD
  use_binary_agd$ARM <- factor(use_binary_agd$ARM, levels = unique(use_binary_agd$ARM))
  n_per_arm <- tapply(use_binary_agd$COUNT, use_binary_agd$ARM, sum)
  n_yes_per_arm <- use_binary_agd$COUNT[use_binary_agd$RESPONSE] # use_binary_agd is already ordered as per factor ARM

  tmpipd <- data.frame(
    USUBJID = NA,
    ARM = mapply(rep, x = levels(use_binary_agd$ARM), each = n_per_arm),
    RESPONSE = unlist(
      lapply(seq(n_per_arm), function(ii) {
        c(rep(TRUE, n_yes_per_arm[ii]), rep(FALSE, n_per_arm[ii] - n_yes_per_arm[ii]))
      })
    )
  )
  tmpipd$USUBJID <- paste0("pseudo_binary_subj_", 1:nrow(tmpipd))

  # output
  tmpipd
}
