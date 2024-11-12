#' Kaplan-Meier (KM) plot function for anchored and unanchored cases using ggplot
#'
#' This is wrapper function of \code{basic_kmplot2}.
#' The argument setting is similar to \code{maic_anchored} and \code{maic_unanchored},
#' and it is used in those two functions.
#'
#' @param weights_object an object returned by \code{estimate_weight}
#' @param tte_ipd a data frame of individual patient data (IPD) of internal trial, contain at least `"USUBJID"`,
#'   `"EVENT"`, `"TIME"` columns and a column indicating treatment assignment
#' @param tte_pseudo_ipd a data frame of pseudo IPD by digitized KM curves of external trial (for time-to-event
#'   endpoint), contain at least `"EVENT"`, `"TIME"`
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial, by default is NULL,
#'   indicating unanchored case
#' @param trt_var_ipd a string, column name in \code{tte_ipd} that contains the treatment assignment
#' @param trt_var_agd a string, column name in \code{tte_pseudo_ipd} that contains the treatment assignment
#' @param normalize_weights logical, default is \code{FALSE}. If \code{TRUE},
#'   \code{scaled_weights} (normalized weights) in \code{weights_object$data} will be used.
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param km_layout a string, only applicable for unanchored case (\code{trt_common = NULL}), indicated the
#'   desired layout of output KM curve.
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months',
#'   weeks' or 'days'
#' @param ... other arguments in \code{basic_kmplot2}
#'
#' @return In unanchored case, a KM plot with risk set table. In anchored case, depending on \code{km_layout},
#' \itemize{
#'   \item if "by_trial", 2 by 1 plot, first all KM curves (incl. weighted) in IPD trial, and then KM curves in AgD
#'   trial, with risk set table.
#'   \item if "by_arm", 2 by 1 plot, first KM curves of \code{trt_agd} and  \code{trt_ipd} (with and without weights),
#'    and then KM curves of \code{trt_common} in AgD trial and IPD trial (with and without weights). Risk set table is
#'     appended.
#'   \item if "all", 2 by 2 plot, all plots in "by_trial" and "by_arm" without risk set table appended.
#' }
#' @example inst/examples/kmplot2_unanchored_ex.R
#' @example inst/examples/kmplot2_anchored_ex.R
#' @export

kmplot2 <- function(weights_object,
                    tte_ipd,
                    tte_pseudo_ipd,
                    trt_ipd,
                    trt_agd,
                    trt_common = NULL,
                    normalize_weights = FALSE,
                    trt_var_ipd = "ARM",
                    trt_var_agd = "ARM",
                    km_conf_type = "log-log",
                    km_layout = c("all", "by_trial", "by_arm"),
                    time_scale,
                    ...) {
  if (!requireNamespace("survminer", quietly = TRUE)) stop("survminer package is required for this function")

  names(tte_ipd) <- toupper(names(tte_ipd))
  names(tte_pseudo_ipd) <- toupper(names(tte_pseudo_ipd))
  trt_var_ipd <- toupper(trt_var_ipd)
  trt_var_agd <- toupper(trt_var_agd)

  # pre check
  if (!"maicplus_estimate_weights" %in% class(weights_object)) {
    stop("weights_object should be an object returned by estimate_weights")
  }
  if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(tte_ipd))) {
    stop(paste("tte_ipd needs to include at least USUBJID, TIME, EVENT,", trt_var_ipd))
  }
  if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(tte_pseudo_ipd))) {
    stop(paste("tte_pseudo_ipd needs to include at least TIME, EVENT,", trt_var_agd))
  }
  km_layout <- match.arg(km_layout, choices = c("all", "by_trial", "by_arm"), several.ok = FALSE)

  # preparing data
  is_anchored <- !is.null(trt_common)
  tte_ipd <- tte_ipd[tte_ipd[[trt_var_ipd]] %in% c(trt_ipd, trt_common), , drop = FALSE]
  tte_pseudo_ipd <- tte_pseudo_ipd[tte_pseudo_ipd[[trt_var_agd]] %in% c(trt_agd, trt_common), , drop = FALSE]
  if (normalize_weights) {
    tte_ipd$weights <- weights_object$data$scaled_weights[match(weights_object$data$USUBJID, tte_ipd$USUBJID)]
  } else {
    tte_ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, tte_ipd$USUBJID)]
  }
  tte_pseudo_ipd$weights <- 1

  tte_ipd$TIME <- get_time_as(tte_ipd$TIME, as = time_scale)
  tte_pseudo_ipd$TIME <- get_time_as(tte_pseudo_ipd$TIME, as = time_scale)
  my_survfit <- function(data, weighted = FALSE) {
    if (weighted) {
      survfit(Surv(TIME, EVENT) ~ 1, data = data, conf.type = km_conf_type, weights = data$weights)
    } else {
      survfit(Surv(TIME, EVENT) ~ 1, data = data, conf.type = km_conf_type)
    }
  }

  if (!is_anchored) {
    kmlist <- list(
      kmobj_B = my_survfit(data = tte_pseudo_ipd),
      kmobj_A = my_survfit(data = tte_ipd),
      kmobj_A_adj = my_survfit(data = tte_ipd, weighted = TRUE)
    )
    kmlist_name <- c(trt_agd, trt_ipd, paste0(trt_ipd, " (weighted)"))
    basic_kmplot2(kmlist, kmlist_name, ...)
  } else if (is_anchored) {
    all_km <- list(
      kmobj_A = my_survfit(data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_ipd, ]),
      kmobj_B = my_survfit(data = tte_pseudo_ipd[tte_pseudo_ipd[, trt_var_agd] == trt_agd, ]),
      kmobj_A_adj = my_survfit(data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_ipd, ], weighted = TRUE),
      kmobj_C = my_survfit(data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_common, ]),
      kmobj_C_adj = my_survfit(data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_common, ], weighted = TRUE),
      kmobj_C_agd = my_survfit(data = tte_pseudo_ipd[tte_pseudo_ipd[, trt_var_agd] == trt_common, ])
    )

    kmlist_combined <- list()
    if (km_layout %in% c("by_trial", "all")) {
      kmlist_1_2 <- list(
        setNames(
          all_km[c(4, 1, 3, 5)],
          c(trt_common, trt_ipd, paste0(trt_ipd, " (weighted)"), paste0(trt_common, " (weighted)"))
        ),
        setNames(all_km[c(6, 2)], c(trt_common, trt_agd))
      )
      names(kmlist_1_2) <- c(
        paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_common, ") in the IPD trial"),
        paste0("Kaplan-Meier Curves \n(", trt_agd, " vs ", trt_common, ") in the AgD trial")
      )
      kmlist_combined <- c(kmlist_combined, kmlist_1_2)
    }
    if (km_layout %in% c("by_arm", "all")) {
      kmlist_3_4 <- list(
        setNames(all_km[c(2, 1, 3)], c(trt_agd, trt_ipd, paste0(trt_ipd, " (weighted)"))),
        setNames(all_km[c(6, 4, 5)], paste(trt_common, c("(AgD)", "(IPD)", "(IPD,weighted)")))
      )
      names(kmlist_3_4) <- c(
        paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_agd, ")"),
        paste0("Kaplan-Meier Curves of Common Comparator \n", trt_common, "(IPD vs AgD Trial)")
      )
      kmlist_combined <- c(kmlist_combined, kmlist_3_4)
    }
    if (km_layout == "all") {
      kmlist_combined <- kmlist_combined[c(1, 3, 2, 4)]
    }

    splots <- mapply(
      FUN = basic_kmplot2,
      kmlist = kmlist_combined,
      kmlist_name = lapply(kmlist_combined, names),
      main_title = names(kmlist_combined),
      MoreArgs = list(...),
      SIMPLIFY = FALSE
    )
    survminer::arrange_ggsurvplots(splots, nrow = 1 + (km_layout == "all"))
  }
}

#' Basic Kaplan Meier (KM) plot function using ggplot
#'
#' This function generates a basic KM plot using ggplot.
#'
#' @param kmlist a list of \code{survfit} object
#' @param kmlist_name a vector indicating the treatment names of each \code{survfit} object
#' @param endpoint_name a string, name of time to event endpoint, to be show in the
#'   last line of title
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param break_x_by bin parameter for \code{survminer}
#' @param censor indicator to include censor information
#' @param xlab label name for x-axis of the plot
#' @param xlim x limit for the x-axis of the plot
#' @param use_colors a character vector of length up to 4, colors to the KM curves,
#'   it will be passed to 'col' of \code{lines()}
#' @param use_line_types a numeric vector of length up to 4, line type to the KM curves,
#'   it will be passed to \code{lty} of \code{lines()}
#' @example inst/examples/basic_kmplot2_ex.R
#' @returns A Kaplan-Meier plot object created with `survminer::ggsurvplot()`.
#' @export

basic_kmplot2 <- function(kmlist,
                          kmlist_name,
                          endpoint_name = "Time to Event Endpoint",
                          show_risk_set = TRUE,
                          main_title = "Kaplan-Meier Curves",
                          break_x_by = NULL,
                          censor = TRUE,
                          xlab = "Time",
                          xlim = NULL,
                          use_colors = NULL,
                          use_line_types = NULL) {
  if (!requireNamespace("survminer", quietly = TRUE)) stop("survminer package is required for this function")

  if (is.null(use_line_types)) {
    use_line_types <- c(1, 1, 2, 2)
  }

  if (is.null(use_colors)) {
    use_colors <- c("#5450E4", "#00857C", "#6ECEB2", "#7B68EE")
  }

  # Produce the Kaplan-Meier plot
  survminer_plot <- survminer::ggsurvplot(kmlist,
    linetype = use_line_types,
    palette = use_colors,
    size = 0.2,
    combine = TRUE,
    risk.table = show_risk_set,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    break.x.by = break_x_by,
    censor = censor,
    censor.size = 2,
    xlab = xlab,
    ylab = endpoint_name,
    legend.title = "Treatment",
    legend = c(0.85, 0.82),
    title = paste0(main_title, "\nEndpoint: ", endpoint_name),
    legend.labs = kmlist_name,
    tables.theme = survminer::theme_cleantable(),
    ggtheme = ggplot2::theme_classic(base_size = 10),
    fontsize = 3,
    conf.int = FALSE,
    xlim = xlim
  )
  survminer_plot
}
