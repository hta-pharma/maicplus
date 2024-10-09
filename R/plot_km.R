#' Kaplan Meier (KM) plot function for anchored and unanchored cases
#'
#' It is wrapper function of \code{basic_kmplot}. The argument setting is similar to \code{maic_anchored} and
#' \code{maic_unanchored}, and it is used in those two functions.
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
#' @param km_layout a string, only applicable for unanchored case (\code{trt_common = NULL}), indicated the desired
#'   layout of output KM curve.
#' @param ... other arguments in \code{basic_kmplot}
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
#' @example inst/examples/kmplot_unanchored_ex.R
#' @example inst/examples/kmplot_anchored_ex.R
#' @export

kmplot <- function(weights_object,
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
                   ...) {
  names(tte_ipd) <- toupper(names(tte_ipd))
  names(tte_pseudo_ipd) <- toupper(names(tte_pseudo_ipd))
  trt_var_ipd <- toupper(trt_var_ipd)
  trt_var_agd <- toupper(trt_var_agd)

  # pre check
  if (!"maicplus_estimate_weights" %in% class(weights_object)) {
    stop("weights_object should be an object returned by estimate_weights")
  }
  if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(tte_ipd))) {
    stop("tte_ipd needs to include at least USUBJID, TIME, EVENT, ", trt_var_ipd)
  }
  if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(tte_pseudo_ipd))) {
    stop("tte_pseudo_ipd needs to include at least TIME, EVENT, ", trt_var_agd)
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

  # generate plot
  if (!is_anchored) {
    ## unanchored case
    kmobj_B <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_agd)),
      data = tte_pseudo_ipd,
      conf.type = km_conf_type
    )
    kmobj_A <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type
    )
    kmobj_A_adj <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type,
      weights = weights
    )

    kmdat <- do.call(
      rbind,
      c(
        survfit_makeup(kmobj_B, trt_agd),
        survfit_makeup(kmobj_A, trt_ipd),
        survfit_makeup(kmobj_A_adj, paste(trt_ipd, "(weighted)"))
      )
    )
    kmdat$treatment <- factor(kmdat$treatment, levels = unique(kmdat$treatment))

    basic_kmplot(kmdat,
      show_risk_set = TRUE,
      main_title = "Kaplan-Meier Curves",
      subplot_heights = NULL,
      suppress_plot_layout = FALSE,
      ...
    )
  } else {
    # anchored case
    # - agd trial km data
    kmobj_C_S1 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_agd)),
      data = tte_pseudo_ipd,
      conf.type = km_conf_type,
      subset = eval(parse(text = paste0("(tte_pseudo_ipd$", trt_var_agd, " == '", trt_common, "')")))
    )
    kmobj_B_S1 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_agd)),
      data = tte_pseudo_ipd,
      conf.type = km_conf_type,
      subset = eval(parse(text = paste0("(tte_pseudo_ipd$", trt_var_agd, " == '", trt_agd, "')")))
    )
    # - ipd trial km data
    kmobj_C_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type,
      subset = eval(parse(text = paste0("(tte_ipd$", trt_var_ipd, " == '", trt_common, "')")))
    )
    kmobj_A_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type,
      subset = eval(parse(text = paste0("(tte_ipd$", trt_var_ipd, " == '", trt_ipd, "')")))
    )
    # - ipd trial km data with weights
    kmobj_Cadj_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type,
      weights = weights,
      subset = eval(parse(text = paste0("(tte_ipd$", trt_var_ipd, " == '", trt_common, "')")))
    )
    kmobj_Aadj_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", trt_var_ipd)),
      data = tte_ipd,
      conf.type = km_conf_type,
      weights = weights,
      subset = eval(parse(text = paste0("(tte_ipd$", trt_var_ipd, " == '", trt_ipd, "')")))
    )
    # - plotdat for layout by trial
    kmdat_s2 <- do.call(
      rbind,
      c(
        survfit_makeup(kmobj_C_S2, trt_common),
        survfit_makeup(kmobj_A_S2, trt_ipd),
        survfit_makeup(kmobj_Aadj_S2, paste(trt_ipd, "(weighted)")),
        survfit_makeup(kmobj_Cadj_S2, paste(trt_common, "(weighted)"))
      )
    )
    kmdat_s2$treatment <- factor(kmdat_s2$treatment, levels = unique(kmdat_s2$treatment))
    kmdat_s1 <- do.call(
      rbind,
      c(
        survfit_makeup(kmobj_C_S1, trt_common),
        survfit_makeup(kmobj_B_S1, trt_agd)
      )
    )
    kmdat_s1$treatment <- factor(kmdat_s1$treatment, levels = unique(kmdat_s1$treatment))
    # - plotdat for layout by arm
    kmdat_a2 <- do.call(
      rbind,
      c(
        survfit_makeup(kmobj_B_S1, trt_agd),
        survfit_makeup(kmobj_A_S2, trt_ipd),
        survfit_makeup(kmobj_Aadj_S2, paste(trt_ipd, "(weighted)"))
      )
    )
    kmdat_a2$treatment <- factor(kmdat_a2$treatment, levels = unique(kmdat_a2$treatment))
    kmdat_a1 <- do.call(
      rbind,
      c(
        survfit_makeup(kmobj_C_S1, paste(trt_common, "(AgD)")),
        survfit_makeup(kmobj_C_S2, paste(trt_common, "(IPD)")),
        survfit_makeup(kmobj_Cadj_S2, paste(trt_common, "(IPD,weighted)"))
      )
    )
    kmdat_a1$treatment <- factor(kmdat_a1$treatment, levels = unique(kmdat_a1$treatment))


    # make plot depending on the layout
    if (km_layout == "by_trial") {
      # 1 by 2 plot, each plot is per trial
      subplot_heights <- c(7, 0.7 + 2 * 0.7, 0.8)
      layout_mat <- matrix(1:4, ncol = 2)
      layout(layout_mat, heights = subplot_heights)

      basic_kmplot(kmdat_s2,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_common, ") in the IPD trial"),
        suppress_plot_layout = TRUE, ...
      )

      basic_kmplot(kmdat_s1,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_agd, " vs ", trt_common, ") in the AgD trial"),
        suppress_plot_layout = TRUE, ...
      )
    } else if (km_layout == "by_arm") {
      # 1 by 2 plot, by 1 is for investigational arm, the other is for common comparator
      subplot_heights <- c(7, 0.7 + 2 * 0.7, 0.8)
      layout_mat <- matrix(1:4, ncol = 2)
      layout(layout_mat, heights = subplot_heights)

      basic_kmplot(kmdat_a2,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_agd, ")"),
        suppress_plot_layout = TRUE, ...
      )

      basic_kmplot(kmdat_a1,
        main_title = paste0("Kaplan-Meier Curves of Common Comparator \n", trt_common, "(IPD vs AgD Trial)"),
        suppress_plot_layout = TRUE, ...
      )
    } else {
      # 2 by 2 plot, combine by trial and by arm
      layout_mat <- matrix(1:4, ncol = 2, byrow = TRUE)
      layout(layout_mat)

      basic_kmplot(kmdat_s2,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_common, ") in the IPD trial"),
        show_risk_set = FALSE,
        suppress_plot_layout = TRUE, ...
      )

      basic_kmplot(kmdat_s1,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_agd, " vs ", trt_common, ") in the AgD trial"),
        show_risk_set = FALSE,
        suppress_plot_layout = TRUE, ...
      )

      basic_kmplot(kmdat_a2,
        main_title = paste0("Kaplan-Meier Curves \n(", trt_ipd, " vs ", trt_agd, ")"),
        show_risk_set = FALSE,
        suppress_plot_layout = TRUE, ...
      )

      basic_kmplot(kmdat_a1,
        main_title = paste0("Kaplan-Meier Curves of Common Comparator \n", trt_common, "(IPD vs AgD Trial)"),
        show_risk_set = FALSE,
        suppress_plot_layout = TRUE, ...
      )
    }
  }
  invisible(NULL)
}


#' Basic Kaplan Meier (KM) plot function
#'
#' This function can generate a basic KM plot with or without risk set table appended at the bottom. In a single plot,
#' it can include up to 4 KM curves. This depends on number of levels in 'treatment' column in the input data.frame
#' \code{kmdat}
#'
#' @param kmdat a `data.frame`, must consist `treatment`, `time` (unit in days), `n.risk`, `censor`, `surv`, similar to
#'   an output from \code{maicplus:::survfit_makeup}
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or 'days'
#' @param time_grid a numeric vector in the unit of \code{time_scale}, risk set table and x axis of the km plot will be
#'   defined based on this time grid
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param subplot_heights a numeric vector, heights argument to \code{graphic::layout()},NULL by default which means
#'   user will use the default setting
#' @param suppress_plot_layout logical, suppress the layout setting in this function so that user can specify layout
#'   outside of the function, FALSE by default
#' @param use_colors a character vector of length up to 4, colors to the KM curves, it will be passed to `col` of
#'   \code{lines()}
#' @param use_line_types a numeric vector of length up to 4, line type to the KM curves, it will be passed to `lty` of
#'   \code{lines()}
#' @param use_pch_cex a scalar between 0 and 1, point size to indicate censored individuals on the KM curves, it will be
#'   passed to `cex` of \code{points()}
#' @param use_pch_alpha a scalar between 0 and 255, degree of color transparency of points to indicate censored
#'   individuals on the KM curves, it will be passed to `cex` of \code{points()}
#'
#' @example inst/examples/basic_kmplot_ex.R
#'
#' @return a KM plot with or without risk set table appended at the bottom, with up to 4 KM curves
#' @export

basic_kmplot <- function(kmdat,
                         endpoint_name = "Time to Event Endpoint",
                         time_scale = NULL,
                         time_grid = NULL,
                         show_risk_set = TRUE,
                         main_title = "Kaplan-Meier Curves",
                         subplot_heights = NULL,
                         suppress_plot_layout = FALSE,
                         use_colors = NULL,
                         use_line_types = NULL,
                         use_pch_cex = 0.65,
                         use_pch_alpha = 100) {
  original_par <- par("bty", "tcl", "mgp", "cex.lab", "cex.axis", "cex.main", "mar")
  on.exit(par(original_par))
  # precheck
  if (!length(subplot_heights) %in% c(0, (1 + show_risk_set))) {
    stop("length of subplot_heights should be ", (1 + show_risk_set))
  }
  if (!is.factor(kmdat$treatment)) {
    stop("kmdat$treatment needs to be a factor, its levels will be used in legend and title, first level is comparator")
  }
  if (nlevels(kmdat$treatment) > 4) stop("kmdat$treatment cannot have more than 4 levels")

  # set up x axis (time)
  if (is.null(time_grid)) {
    max_t <- max(kmdat$time)
    t_range <- c(0, get_time_as(max_t, time_scale) * 1.07)
    time_grid <- pretty(t_range)
  } else {
    t_range <- c(0, max(time_grid))
  }



  # plat layout in par
  if (!suppress_plot_layout) {
    nr_subplot <- (1 + show_risk_set)
    if (is.null(subplot_heights)) subplot_heights <- c(7, 0.7 + nlevels(kmdat$treatment) * 0.7, 0.8)
    layout_mat <- matrix(1:nr_subplot, ncol = 1)
    layout(layout_mat, heights = subplot_heights)
  }

  # plot cosmetic setup
  if (is.null(use_line_types)) {
    use_lty <- c(1, 1, 2, 2)
    use_lwd <- c(1.5, 1.5, 1.2, 1.2)
  } else {
    use_lty <- use_line_types
    use_lwd <- c(1.5, 1.5, 1.2, 1.2)
  }

  if (is.null(use_colors)) {
    use_col <- c("#5450E4", "#00857C", "#6ECEB2", "#7B68EE")
  } else {
    use_col <- use_colors
  }
  use_col2 <- col2rgb(use_col) # preparing semi-transparent colors
  use_col2 <- rgb(use_col2[1, ], use_col2[2, ], use_col2[3, ], alpha = use_pch_alpha, maxColorValue = 255)

  ## : first subplot: KM curve -------------------
  # base plot
  par(bty = "n", tcl = -0.15, mgp = c(1.8, 0.4, 0), cex.lab = 0.85, cex.axis = 0.8, cex.main = 0.9, mar = c(3, 4, 5, 1))
  plot(0, 0,
    type = "n", xlab = paste0("Time in ", time_scale), ylab = "Survival Probability",
    ylim = c(0, 1), xlim = t_range, yaxt = "n", xaxt = "n",
    main = paste0(main_title, "\nEndpoint:", endpoint_name)
  )
  axis(2, las = 1)
  if (!is.null(time_grid)) {
    axis(1, at = time_grid)
  } else {
    axis(1)
  }

  # add km line
  for (ii in 1:nlevels(kmdat$treatment)) {
    tmpkmdat <- kmdat[as.numeric(kmdat$treatment) == ii, , drop = FALSE]
    lines(
      y = tmpkmdat$surv,
      x = get_time_as(tmpkmdat$time, time_scale),
      col = use_col[ii],
      lty = use_lty[ii],
      lwd = use_lwd[ii],
      type = "s"
    )
    tmpid <- (tmpkmdat$censor != 0) # cannot just ==1, anticipating weighted case
    points(
      y = tmpkmdat$surv[tmpid],
      x = get_time_as(tmpkmdat$time[tmpid], time_scale),
      col = use_col2[ii],
      pch = 3,
      cex = use_pch_cex
    )
  }

  ## : second subplot: risk set table (if applicable) -------------------
  if (show_risk_set) {
    # add legend, with treatment index
    legend("topright",
      bty = "n",
      cex = 0.8,
      lty = use_lty[1:nlevels(kmdat$treatment)],
      lwd = use_lwd[1:nlevels(kmdat$treatment)],
      col = use_col[1:nlevels(kmdat$treatment)],
      legend = paste0("(T", 1:nlevels(kmdat$treatment), ") ", levels(kmdat$treatment))
    )

    # add risk set table
    par(bty = "n", tcl = -0.15, mgp = c(1.8, 0.4, 0), mar = c(1, 4, 0, 1))
    plot(0, 0,
      type = "n", xlab = "", ylab = "", main = NULL,
      ylim = c(nlevels(kmdat$treatment) + 1.2, -0.5),
      xlim = t_range,
      yaxt = "n", xaxt = "n"
    )
    axis(2,
      at = 1:nlevels(kmdat$treatment), labels = paste0("T", 1:nlevels(kmdat$treatment)),
      line = NA, lty = "blank", las = 1
    )

    for (ii in 1:nlevels(kmdat$treatment)) {
      tmpkmdat <- kmdat[as.numeric(kmdat$treatment) == ii, , drop = FALSE]
      tmptime <- get_time_as(tmpkmdat$time, time_scale)
      tmpnr <- sapply(time_grid, function(kk) {
        tmpid <- which(tmptime > kk)
        if (length(tmpid) == 0) {
          if (min(tmpkmdat$n.risk) == 0) {
            tmpid <- which.min(tmpkmdat$n.risk)[1]
          } else {
            tmpid <- NULL
          }
        }
        tout <- ifelse(is.null(tmpid), "n/a", round(tmpkmdat$n.risk[tmpid], 1))
      })
      text(0, 0, labels = "Number at risk", pos = 4, cex = 0.8, offset = -0.8)
      text(
        y = rep(ii, length(time_grid)),
        x = time_grid,
        labels = tmpnr,
        col = use_col[ii],
        cex = 0.75
      )
      text(0, nlevels(kmdat$treatment) + 1,
        pos = 4, cex = 0.7, offset = -0.8, col = "gray30",
        labels = "Note: Number at risk for adjusted/weighted treament arm is the sum of individual weight at risk."
      )
    }
  } else {
    # add simple
    legend("topright",
      bty = "n",
      cex = 0.8,
      lty = use_lty[1:nlevels(kmdat$treatment)],
      lwd = use_lwd[1:nlevels(kmdat$treatment)],
      col = use_col[1:nlevels(kmdat$treatment)],
      legend = levels(kmdat$treatment)
    )
  }
  invisible(NULL)
}


#' Diagnosis plot of proportional hazard assumption for anchored and unanchored
#'
#' @param weights_object an object returned by \code{estimate_weight}
#' @param tte_ipd a data frame of individual patient data (IPD) of internal trial, contain at least "USUBJID", "EVENT",
#'   "TIME" columns and a column indicating treatment assignment
#' @param tte_pseudo_ipd a data frame of pseudo IPD by digitized KM curves of external trial (for time-to-event
#'   endpoint), contain at least "EVENT", "TIME"
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{tte_ipd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial
#'   \code{tte_pseudo_ipd} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial, by default is NULL,
#'   indicating unanchored case
#' @param trt_var_ipd a string, column name in \code{tte_ipd} that contains the treatment assignment
#' @param trt_var_agd a string, column name in \code{tte_pseudo_ipd} that contains the treatment assignment
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or 'days'
#' @param zph_transform a string, pass to \code{survival::cox.zph}, default is "log"
#' @param zph_log_hazard a logical, if TRUE (default), y axis of the time dependent hazard function is log-hazard,
#'   otherwise, hazard.
#'
#' @return a 3 by 2 plot, include log-cumulative hazard plot, time dependent hazard function and unscaled Schoenfeld
#'   residual plot, before and after matching
#'
#' @example inst/examples/ph_diagplot_unanchored_ex.R
#' @example inst/examples/ph_diagplot_anchored_ex.R
#' @export
ph_diagplot <- function(weights_object,
                        tte_ipd,
                        tte_pseudo_ipd,
                        trt_ipd,
                        trt_agd,
                        trt_common = NULL,
                        trt_var_ipd = "ARM",
                        trt_var_agd = "ARM",
                        endpoint_name = "Time to Event Endpoint",
                        time_scale,
                        zph_transform = "log",
                        zph_log_hazard = TRUE) {
  names(tte_ipd) <- toupper(names(tte_ipd))
  names(tte_pseudo_ipd) <- toupper(names(tte_pseudo_ipd))
  trt_var_ipd <- toupper(trt_var_ipd)
  trt_var_agd <- toupper(trt_var_agd)

  # pre check
  if (!"maicplus_estimate_weights" %in% class(weights_object)) {
    stop("weights_object should be an object returned by estimate_weights")
  }
  if (!all(c("USUBJID", "TIME", "EVENT", trt_var_ipd) %in% names(tte_ipd))) {
    stop("tte_ipd needs to include at least USUBJID, TIME, EVENT, ", trt_var_ipd)
  }
  if (!all(c("TIME", "EVENT", trt_var_agd) %in% names(tte_pseudo_ipd))) {
    stop("tte_ipd needs to include at least TIME, EVENT, ", trt_var_agd)
  }

  # preparing analysis data
  is_anchored <- ifelse(is.null(trt_common), FALSE, TRUE)
  tte_ipd <- tte_ipd[tte_ipd[[trt_var_ipd]] %in% c(trt_ipd, trt_common), , drop = TRUE]
  tte_pseudo_ipd <- tte_pseudo_ipd[tte_pseudo_ipd[[trt_var_agd]] %in% c(trt_agd, trt_common), , drop = TRUE]
  tte_ipd$weights <- weights_object$data$weights[match(weights_object$data$USUBJID, tte_ipd$USUBJID)]
  tte_pseudo_ipd$weights <- 1
  tte_ipd$TIME2 <- get_time_as(tte_ipd$TIME, as = time_scale) # for cox.zph
  tte_pseudo_ipd$TIME2 <- get_time_as(tte_pseudo_ipd$TIME, as = time_scale) # for cox.zph
  if (!"USUBJID" %in% names(tte_pseudo_ipd)) tte_pseudo_ipd$USUBJID <- paste0("ID", seq_len(nrow(tte_pseudo_ipd)))
  if (trt_var_ipd != "ARM") tte_ipd$ARM <- tte_ipd[[trt_var_ipd]]
  if (trt_var_agd != "ARM") tte_pseudo_ipd$ARM <- tte_pseudo_ipd[[trt_var_agd]]

  # prepare plot data
  retain_cols <- c("USUBJID", "TIME", "TIME2", "EVENT", "ARM", "weights")
  if (!is_anchored) {
    # unanchored case
    tte_dat <- rbind(
      tte_ipd[, retain_cols, drop = FALSE],
      tte_pseudo_ipd[, retain_cols, drop = FALSE]
    )
  } else {
    tte_dat <- tte_ipd[, retain_cols, drop = FALSE]
  }
  kmobj <- survival::survfit(Surv(TIME, EVENT) ~ ARM, tte_dat, conf.type = "log-log")
  kmobj_adj <- survival::survfit(Surv(TIME, EVENT) ~ ARM, tte_dat, conf.type = "log-log", weights = weights)
  coxobj <- survival::coxph(Surv(TIME, EVENT) ~ ARM, data = tte_dat)
  coxobj2 <- survival::coxph(Surv(TIME2, EVENT) ~ ARM, data = tte_dat)
  zphobj <- survival::cox.zph(coxobj2, transform = zph_transform, global = FALSE)
  coxobj_adj <- survival::coxph(Surv(TIME, EVENT) ~ ARM, data = tte_dat, weights = weights)
  coxobj_adj2 <- survival::coxph(Surv(TIME2, EVENT) ~ ARM, data = tte_dat, weights = weights)
  zphobj_adj <- survival::cox.zph(coxobj_adj2, transform = zph_transform, global = FALSE)

  # making the plot
  original_par <- par(mfrow = c(3, 2), cex.lab = 0.85, cex.axis = 0.8, cex.main = 0.9)
  on.exit(par(original_par))
  # log-cum-hazard plot
  ph_diagplot_lch(kmobj,
    time_scale = time_scale,
    log_time = TRUE,
    endpoint_name = endpoint_name,
    subtitle = "(Before Matching)"
  )

  ph_diagplot_lch(kmobj_adj,
    time_scale = time_scale,
    log_time = TRUE,
    endpoint_name = endpoint_name,
    subtitle = "(After Matching)"
  )
  # time dependent hazard plot
  plot(zphobj,
    main = paste0(
      "Time-dependent Hazard function (scaled Schoenfeld residual)\n",
      "Endpoint:", endpoint_name, "\n(Before Matching)"
    ),
    resid = FALSE, se = TRUE, df = 4, nsmo = 40,
    # xlim = range(0,zphobj$time),
    ylab = ifelse(zph_log_hazard, "Log Hazard", "Hazard"),
    xlab = paste("Time in", time_scale),
    lty = 1:2, lwd = 2, pch = 16, cex = 0.8,
    col = rgb(0, 0, 128, alpha = 120, maxColorValue = 255),
    hr = (!zph_log_hazard), yaxt = "n"
  )
  axis(2, las = 1)
  pv <- as.data.frame(zphobj$table)$p
  pv <- ifelse(round(pv, 4) < 0.0001, "<0.0001", format(round(pv, 4), nsmall = 4))
  legend("bottomright",
    cex = 0.75, bty = "n", text.col = "dodgerblue3",
    legend = c(paste0("p-value: ", pv), paste0("time-transform: ", zph_transform)),
    title = "PH test (survival::cox.zph)"
  )

  plot(zphobj_adj,
    main = paste0(
      "Time-dependent Hazard function (scaled Schoenfeld residual)\n",
      "Endpoint:", endpoint_name, "\n(After Matching)"
    ),
    resid = FALSE, se = TRUE, df = 4, nsmo = 40,
    # xlim = range(0,zphobj$time),
    ylab = ifelse(zph_log_hazard, "Log Hazard", "Hazard"),
    xlab = paste("Time in", time_scale),
    lty = 1:2, lwd = 2, pch = 16, cex = 0.8,
    col = rgb(0, 0, 128, alpha = 120, maxColorValue = 255),
    hr = (!zph_log_hazard), yaxt = "n"
  )
  axis(2, las = 1)
  pv <- as.data.frame(zphobj_adj$table)$p
  pv <- ifelse(round(pv, 4) < 0.0001, "<0.0001", format(round(pv, 4), nsmall = 4))
  legend("bottomright",
    cex = 0.75, bty = "n", text.col = "dodgerblue3",
    legend = c(paste0("p-value: ", pv), paste0("time-transform: ", zph_transform)),
    title = "PH test (survival::cox.zph)"
  )

  # unscaled schoenfeld residual
  ph_diagplot_schoenfeld(coxobj,
    time_scale = time_scale,
    log_time = FALSE,
    endpoint_name = endpoint_name,
    subtitle = "(Before Matching)"
  )

  ph_diagplot_schoenfeld(coxobj_adj,
    time_scale = time_scale,
    log_time = FALSE,
    endpoint_name = endpoint_name,
    subtitle = "(After Matching)"
  )
}

#' PH Diagnosis Plot of Log Cumulative Hazard Rate versus time or log-time
#'
#' This plot is also known as log negative log survival rate.
#'
#' a diagnosis plot for proportional hazard assumption, versus log-time (default) or time
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param time_scale a character string, 'years', 'months', 'weeks' or 'days', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @param exclude_censor logical, should censored data point be plotted
#' @examples
#' library(survival)
#' data(adtte_sat)
#' data(pseudo_ipd_sat)
#' combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data, conf.type = "log-log")
#' ph_diagplot_lch(kmobj,
#'   time_scale = "month", log_time = TRUE,
#'   endpoint_name = "OS", subtitle = "(Before Matching)"
#' )
#' @return a plot of log cumulative hazard rate
#' @export

ph_diagplot_lch <- function(km_fit,
                            time_scale,
                            log_time = TRUE,
                            endpoint_name = "",
                            subtitle = "",
                            exclude_censor = TRUE) {
  time_scale <- match.arg(arg = time_scale, choices = c("days", "weeks", "months", "years"))

  clldat <- survfit_makeup(km_fit)

  if (exclude_censor) {
    clldat <- lapply(clldat, function(xxt) xxt[xxt$censor == 0, , drop = FALSE])
  }

  all.times <- get_time_as(do.call(rbind, clldat)$time, time_scale)
  if (log_time) all.times <- log(all.times)
  t_range <- range(all.times)
  y_range <- range(log(do.call(rbind, clldat)$cumhaz))

  original_par <- par("mar", "bty", "tcl", "mgp")
  par(mar = c(4, 4, 4, 1), bty = "n", tcl = -0.15, mgp = c(1.5, 0.3, 0))
  on.exit(par(original_par))
  plot(0, 0,
    type = "n", xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
    ylab = "Log-Cumulative Hazard Rate",
    ylim = y_range, xlim = t_range, yaxt = "n",
    main = paste0(
      "Log Cumulative Hazard versus Log Time\nEndpoint: ", endpoint_name,
      ifelse(subtitle == "", "", "\n"), subtitle
    )
  )
  axis(2, las = 1)

  trts <- names(clldat)
  cols <- c("dodgerblue3", "firebrick3")
  pchs <- c(1, 4)
  for (i in seq_along(clldat)) {
    use_x <- get_time_as(clldat[[i]]$time, time_scale)
    if (log_time) use_x <- log(use_x)

    lines(
      y = log(clldat[[i]]$cumhaz),
      x = use_x, col = cols[i],
      type = "s",
    )
    points(
      y = log(clldat[[i]]$cumhaz),
      x = use_x,
      col = cols[i], pch = pchs[i], cex = 0.7
    )
  }
  legend("bottomright",
    bty = "n", lty = c(1, 1, 2), cex = 0.8,
    col = cols, pch = pchs, legend = paste0("Treatment: ", trts)
  )
}


#' PH Diagnosis Plot of Schoenfeld residuals for a Cox model fit
#'
#' @param coxobj object returned from \code{\link[survival]{coxph}}
#' @param time_scale a character string, 'years', 'months', 'weeks' or 'days', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @examples
#' library(survival)
#' data(adtte_sat)
#' data(pseudo_ipd_sat)
#' combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
#' unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = combined_data)
#' ph_diagplot_schoenfeld(unweighted_cox,
#'   time_scale = "month", log_time = TRUE,
#'   endpoint_name = "OS", subtitle = "(Before Matching)"
#' )
#' @return a plot of Schoenfeld residuals
#' @export

ph_diagplot_schoenfeld <- function(coxobj,
                                   time_scale = "months",
                                   log_time = TRUE,
                                   endpoint_name = "",
                                   subtitle = "") {
  # pre-check
  time_scale <- match.arg(arg = time_scale, choices = c("days", "weeks", "months", "years"))

  # prepare data
  schresid <- residuals(coxobj, type = "schoenfeld")
  plot_x <- get_time_as(as.numeric(names(schresid)), time_scale)
  if (log_time) plot_x <- log(plot_x)

  # loewss fit
  fit0 <- predict(loess(schresid ~ plot_x), se = TRUE)
  uppband <- fit0$fit + qt(0.975, fit0$df) * fit0$se
  lowband <- fit0$fit - qt(0.975, fit0$df) * fit0$se
  use_yrange <- range(schresid, uppband, lowband)

  # making the plot
  original_par <- par(bty = "n", mar = c(4, 4, 4, 1), tcl = -0.15, mgp = c(1.5, 0.3, 0))
  on.exit(par(original_par))
  plot(schresid ~ plot_x,
    type = "n",
    yaxt = "n", ylim = use_yrange,
    ylab = "Unscaled Schoenfeld Residual",
    xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
    main = paste0(
      "Unscaled Schoenfeld Residual\nEndpoint: ", endpoint_name,
      ifelse(subtitle == "", "", "\n"), subtitle
    )
  )
  axis(2, las = 1)
  lines(fit0$fit ~ plot_x, lty = 2, lwd = 1, col = rgb(0, 0, 128, 150, maxColorValue = 255))
  # lines(uppband ~ plot_x, lty =2, lwd=1, col = rgb(0,0,128,120,maxColorValue = 255))
  # lines(lowband ~ plot_x, lty =2, lwd=1, col = rgb(0,0,128,120,maxColorValue = 255))
  polygon(
    x = c(plot_x, rev(plot_x)),
    y = c(uppband, rev(lowband)),
    col = rgb(0, 0, 128, 60, maxColorValue = 255),
    border = NA
  )
  abline(h = 0, lty = 1, lwd = 1, col = "deeppink")
  points(schresid ~ plot_x,
    pch = 16,
    cex = 0.85,
    col = rgb(169, 169, 169, 120, maxColorValue = 255)
  )
}
