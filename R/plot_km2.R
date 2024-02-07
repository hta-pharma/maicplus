#' Kaplan-Meier (KM) plot function for anchored and unanchored cases using ggplot
#'
#' This is wrapper function of \code{basic_kmplot2}.
#' The argument setting is similar to \code{maic_anchored} and \code{maic_unanchored},
#' and it is used in those two functions.
#'
#' @param weights_object an object returned by \code{estimate_weight}
#' @param tte_ipd a data frame of individual patient data (IPD) of internal trial, contain at least `"USUBJID"`, `"EVENT"`, `"TIME"` columns and a column indicating treatment assignment
#' @param tte_pseudo_ipd a data frame of pseudo IPD by digitized KM curves of external trial (for time-to-event endpoint), contain at least `"EVENT"`, `"TIME"`
#' @param trt_ipd a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial, by default is NULL, indicating unanchored case
#' @param trt_var_ipd a string, column name in \code{tte_ipd} that contains the treatment assignment
#' @param trt_var_agd a string, column name in \code{tte_pseudo_ipd} that contains the treatment assignment
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param km_layout a string, only applicable for unanchored case (\code{trt_common = NULL}), indicated the desired layout of output KM curve.
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or 'days'
#' @param ... other arguments in \code{basic_kmplot2}
#'
#' @return
#' In unanchored case, a KM plot with risk set table. In anchored case, depending on \code{km_layout},
#' \itemize{
#'   \item if "by_trial", 2 by 1 plot, first all KM curves (incl. weighted) in IPD trial, and then KM curves in AgD trial, with risk set table.
#'   \item if "by_arm", 2 by 1 plot, first KM curves of \code{trt_agd} and  \code{trt_ipd} (with and without weights), and then KM cuvers of \code{trt_common} in AgD trial and IPD trial (with and without weights). Risk set table is appended.
#'   \item if "all", 2 by 2 plot, all plots in "by_trial" and "by_arm" without risk set table appended.
#' }
#' @example inst/examples/kmplot2_anchored_ex.R
#' @example inst/examples/kmplot2_unanchored_ex.R
#' @export

kmplot2 <- function(weights_object,
                    tte_ipd,
                    tte_pseudo_ipd,
                    trt_ipd,
                    trt_agd,
                    trt_common = NULL,
                    trt_var_ipd = "ARM",
                    trt_var_agd = "ARM",
                    km_conf_type = "log-log",
                    km_layout = c("all", "by_trial", "by_arm"),
                    time_scale,
                    ...) {
  names(tte_ipd) <- toupper(names(tte_ipd))
  names(tte_pseudo_ipd) <- toupper(names(tte_pseudo_ipd))
  trt_var_ipd <- toupper(trt_var_ipd)
  trt_var_agd <- toupper(trt_var_agd)

  # pre check
  if (!"maicplus_estimate_weights" %in% class(weights_object)) stop("weights_object should be an object returned by estimate_weights")
  if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(tte_ipd))) stop(paste("tte_ipd needs to include at least USUBJID, TIME, EVENT,", trt_var_ipd))
  if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(tte_pseudo_ipd))) stop(paste("tte_pseudo_ipd needs to include at least TIME, EVENT,", trt_var_agd))
  km_layout <- match.arg(km_layout, choices = c("all", "by_trial", "by_arm"), several.ok = FALSE)

  # preparing data
  is_anchored <- !is.null(trt_common)
  tte_ipd <- tte_ipd[tte_ipd[[trt_var_ipd]] %in% c(trt_ipd, trt_common), , drop = TRUE]
  tte_pseudo_ipd <- tte_pseudo_ipd[tte_pseudo_ipd[[trt_var_agd]] %in% c(trt_agd, trt_common), , drop = TRUE]
  tte_ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, tte_ipd$USUBJID)]
  tte_pseudo_ipd$weights <- 1

  tte_ipd$TIME <- get_time_as(tte_ipd$TIME, as = time_scale)
  tte_pseudo_ipd$TIME <- get_time_as(tte_pseudo_ipd$TIME, as = time_scale)

  if (!is_anchored) {
    kmobj_B <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_pseudo_ipd,
      conf.type = km_conf_type
    )

    kmobj_A <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd,
      conf.type = km_conf_type
    )

    kmobj_A_adj <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd,
      conf.type = km_conf_type,
      weights = weights
    )

    kmlist <- list(
      kmobj_B = kmobj_B,
      kmobj_A = kmobj_A,
      kmobj_A_adj = kmobj_A_adj
    )
    kmlist_name <- c(trt_agd, trt_ipd, paste0(trt_ipd, " (weighted)"))
    basic_kmplot2(
      kmlist,
      kmlist_name,
      ...
    )
  } else if (is_anchored) {
    kmobj_A <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_ipd, ],
      conf.type = km_conf_type
    )

    kmobj_B <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_pseudo_ipd[tte_pseudo_ipd[, trt_var_agd] == trt_agd, ],
      conf.type = km_conf_type
    )

    kmobj_A_adj <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_ipd, ],
      conf.type = km_conf_type,
      weights = weights
    )

    kmobj_C <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_common, ],
      conf.type = km_conf_type
    )

    kmobj_C_adj <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_ipd[tte_ipd[, trt_var_ipd] == trt_common, ],
      conf.type = km_conf_type,
      weights = weights
    )

    kmobj_C_agd <- survfit(Surv(TIME, EVENT) ~ 1,
      data = tte_pseudo_ipd[tte_pseudo_ipd[, trt_var_agd] == trt_common, ],
      conf.type = km_conf_type
    )

    kmlist <- list(
      kmobj_C = kmobj_C,
      kmobj_A = kmobj_A,
      kmobj_A_adj = kmobj_A_adj,
      kmobj_C_adj = kmobj_C_adj
    )
    kmlist_name <- c(trt_ipd, trt_common, paste0(trt_ipd, " (weighted)"), paste0(trt_common, " (weighted)"))

    kmlist2 <- list(
      kmobj_C_agd = kmobj_C_agd,
      kmobj_B = kmobj_B
    )
    kmlist2_name <- c(trt_common, trt_agd)

    kmlist3 <- list(
      kmobj_B = kmobj_B,
      kmobj_A = kmobj_A,
      kmobj_A_adj = kmobj_A_adj
    )
    kmlist3_name <- c(trt_agd, trt_ipd, paste0(trt_ipd, " (weighted)"))

    kmlist4 <- list(
      kmobj_C_agd = kmobj_C_agd,
      kmobj_C = kmobj_C,
      kmobj_C_adj = kmobj_C_adj
    )
    kmlist4_name <- c(paste0(trt_common, " (AgD)"), paste0(trt_common, " (IPD)"), paste0(trt_common, " (IPD,weighted)"))

    main_title <- paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_common, ") in the IPD trial")
    main_title2 <- paste0("Kaplan-Meier Curves \n(", trt_agd, " vs ", trt_common, ") in the AgD trial")
    main_title3 <- paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_agd, ")")
    main_title4 <- paste0("Kaplan-Meier Curves of Common Comparator \n", trt_common, "(IPD vs AgD Trial)")

    kmlist_combined <- list(kmlist, kmlist2, kmlist3, kmlist4)
    kmlist_name_combined <- list(kmlist_name, kmlist2_name, kmlist3_name, kmlist4_name)
    main_title_combined <- list(main_title, main_title2, main_title3, main_title4)

    plot_basic_kmplot2 <- function(ii) {
      basic_kmplot2(kmlist_combined[[ii]],
        kmlist_name_combined[[ii]],
        main_title = main_title_combined[[ii]],
        ...
      )
    }

    splots <- list()

    if (km_layout == "by_trial") {
      splots <- lapply(c(1, 2), plot_basic_kmplot2)
      arrange_ggsurvplots(splots, ncol = 2, nrow = 1)
    } else if (km_layout == "by_arm") {
      splots <- lapply(c(3, 4), plot_basic_kmplot2)
      arrange_ggsurvplots(splots, ncol = 2, nrow = 1)
    } else if (km_layout == "all") {
      splots <- lapply(c(1, 3, 2, 4), plot_basic_kmplot2)
      arrange_ggsurvplots(splots, ncol = 2, nrow = 2)
    }
  }
}

#' Basic Kaplan Meier (KM) plot function using ggplot
#'
#' This function generates a basic KM plot using ggplot.
#'
#' @param kmlist a list of survival::survfit object
#' @param kmlist_name a vector indicating the treatment names of each survfit object
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param break_x_by bin parameter for survminer
#' @param censor indicator to include censor information
#' @param xlim x limit for the x-axis of the plot
#' @param use_colors a character vector of length up to 4, colors to the KM curves, it will be passed to 'col' of \code{lines()}
#' @param use_line_types a numeric vector of length up to 4, line type to the KM curves, it will be passed to 'lty' of \code{lines()}

basic_kmplot2 <- function(kmlist,
                          kmlist_name,
                          endpoint_name = "Time to Event Endpoint",
                          show_risk_set = TRUE,
                          main_title = "Kaplan-Meier Curves",
                          break_x_by = NULL,
                          censor = T,
                          xlim = NULL,
                          use_colors = NULL,
                          use_line_types = NULL) {
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
    xlab = "Time",
    ylab = endpoint_name,
    legend.title = "Treatment",
    legend = c(0.85, 0.82),
    title = paste0(main_title, "\nEndpoint:", endpoint_name),
    legend.labs = kmlist_name,
    tables.theme = theme_cleantable(),
    ggtheme = ggplot2::theme_classic(base_size = 10),
    fontsize = 3,
    conf.int = FALSE,
    xlim = xlim
  )
  survminer_plot
}
