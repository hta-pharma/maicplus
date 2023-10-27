#' Kaplan Meier (KM) plot function for anchored and unanchored cases
#'
#' It is wrapper function of \code{basic_kmplot}.
#' The argument setting is similar to \code{maic_anchored} and \code{maic_unanchored},
#' and it is used in those two functions.
#'
#' @param ipd_weights an object returned by \code{estimate_weight}
#' @param tte_dat_ipd a data frame of individual patient data (IPD) of internal trial, contain at least "USUBJID", "EVENT", "TIME" columns and a column indicating treatment assignment
#' @param ipd_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param tte_dat_pseudo a data frame of pseudo IPD by digitized KM curves of external trial (for time-to-event endpoint), contain at least "EVENT", "TIME"
#' @param pseudo_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_agd a string, name of the interested investigation arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial, by default is NULL, indicating unanchored case
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param km_layout a string, only applicable for unachored case (\code{trt_common = NULL}), indicated the desired layout of output KM curve.
#' @param ... other arguments in \code{basic_kmplot}
#'
#' @return
#' In unanchored case, a KM plot with risk set table. In anchored case, depending on \code{km_layout},
#' \itemize{
#'   \item if "by_trial", 2 by 1 plot, first all KM curves (incl. weighted) in IPD trial, and then KM curves in AgD trial, with risk set table.
#'   \item if "by_arm", 2 by 1 plot, first KM curves of \code{trt_agd} and  \code{trt_ipd} (with and without weights), and then KM cuvers of \code{trt_common} in AgD trial and IPD trial (with and without weights). Risk set table is appended.
#'   \item if "all", 2 by 2 plot, all plots in "by_trial" and "by_arm" without risk set table appended.
#' }
#' @example kmplot_anchored_ex.R. kmplot_unanchored_ex.R
#' @export

kmplot <- function(ipd_weights,
                   tte_dat_ipd,
                   ipd_trt_var = "ARM",
                   tte_dat_pseudo,
                   pseudo_trt_var = "ARM",
                   trt_ipd,
                   trt_agd,
                   trt_common = NULL,
                   km_conf_type = "log-log",
                   km_layout = c("all","by_trial","by_arm"),
                   ...){

  names(tte_dat_ipd) <- toupper(names(tte_dat_ipd))
  names(tte_dat_pseudo) <- toupper(names(tte_dat_pseudo))
  ipd_trt_var <- toupper(ipd_trt_var)
  pseudo_trt_var <- toupper(pseudo_trt_var)

  # pre check
  if(!all(c("USUBJID","TIME", "EVENT", ipd_trt_var) %in% names(tte_dat_ipd))) stop(paste("tte_dat_ipd needs to include at least USUBJID, TIME, EVENT,", ipd_trt_var))
  if(!all(c("TIME", "EVENT", pseudo_trt_var) %in% names(tte_dat_pseudo))) stop(paste("tte_dat_ipd needs to include at least TIME, EVENT,", pseudo_trt_var))
  km_layout <- match.arg(km_layout, choices = c("all","by_trial","by_arm"), several.ok = FALSE)

  # preparing data
  is_anchored <- ifelse(is.null(trt_common), FALSE, TRUE)
  tte_dat_ipd <- tte_dat_ipd[tte_dat_ipd[[ipd_trt_var]]%in%c(trt_ipd,trt_common),,drop=TRUE]
  tte_dat_pseudo <- tte_dat_pseudo[tte_dat_pseudo[[pseudo_trt_var]]%in%c(trt_agd,trt_common),,drop=TRUE]
  tte_dat_ipd$weights <- ipd_weights$data$weights[match(ipd_weights$data$USUBJID,tte_dat_ipd$USUBJID)]
  tte_dat_pseudo$weights <- 1

  # generate plot
  if(!is_anchored){

    ## unanchored case
    kmobj_B <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", pseudo_trt_var)),
                       data = tte_dat_pseudo,
                       conf.type = km_conf_type)
    kmobj_A <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                       data = tte_dat_ipd,
                       conf.type = km_conf_type)
    kmobj_A_adj <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                           data = tte_dat_ipd,
                           conf.type = km_conf_type,
                           weights = weights)

    kmdat <- do.call(rbind,
                     c(survfit_makeup(kmobj_B, trt_agd),
                       survfit_makeup(kmobj_A, trt_ipd),
                       survfit_makeup(kmobj_A_adj, paste(trt_ipd,"(weighted)")))
    )
    kmdat$treatment <- factor(kmdat$treatment, levels = unique(kmdat$treatment))

    basic_kmplot(kmdat,
                 show_risk_set = TRUE,
                 main_title = "Kaplan-Meier Curves",
                 subplot_heights = NULL,
                 suppress_plot_layout = FALSE,
                 ...)

  }else{

    # anchored case
    # - agd trial km data
    kmobj_C_S1 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", pseudo_trt_var)),
                          data = tte_dat_pseudo,
                          conf.type = km_conf_type,
                          subset = eval(parse( text = paste0("(tte_dat_pseudo$",pseudo_trt_var," == '",trt_common,"')") ))
                  )
    kmobj_B_S1 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", pseudo_trt_var)),
                       data = tte_dat_pseudo,
                       conf.type = km_conf_type,
                       subset = eval(parse( text = paste0("(tte_dat_pseudo$",pseudo_trt_var," == '",trt_agd,"')") ))
                  )
    # - ipd trial km data
    kmobj_C_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                          data = tte_dat_ipd,
                          conf.type = km_conf_type,
                          subset = eval(parse( text = paste0("(tte_dat_ipd$",ipd_trt_var," == '",trt_common,"')") ))
                  )
    kmobj_A_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                          data = tte_dat_ipd,
                          conf.type = km_conf_type,
                          subset = eval(parse( text = paste0("(tte_dat_ipd$",ipd_trt_var," == '",trt_ipd,"')") ))
                  )
    # - ipd trial km data with weights
    kmobj_Cadj_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                             data = tte_dat_ipd,
                             conf.type = km_conf_type,
                             weights = weights,
                             subset = eval(parse( text = paste0("(tte_dat_ipd$",ipd_trt_var," == '",trt_common,"')") ))
                     )
    kmobj_Aadj_S2 <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
                             data = tte_dat_ipd,
                             conf.type = km_conf_type,
                             weights = weights,
                             subset = eval(parse( text = paste0("(tte_dat_ipd$",ipd_trt_var," == '",trt_ipd,"')") ))
                     )
    # - plotdat for layout by trial
    kmdat_s2 <- do.call(rbind,
                        c(survfit_makeup(kmobj_C_S2, trt_common),
                          survfit_makeup(kmobj_A_S2, trt_ipd),
                          survfit_makeup(kmobj_Aadj_S2, paste(trt_ipd,"(weighted)")),
                          survfit_makeup(kmobj_Cadj_S2, paste(trt_common,"(weighted)")))
    )
    kmdat_s2$treatment <- factor(kmdat_s2$treatment, levels = unique(kmdat_s2$treatment))
    kmdat_s1 <- do.call(rbind,
                        c(survfit_makeup(kmobj_C_S1, trt_common),
                          survfit_makeup(kmobj_B_S1, trt_agd))
    )
    kmdat_s1$treatment <- factor(kmdat_s1$treatment, levels = unique(kmdat_s1$treatment))
    # - plotdat for layout by arm
    kmdat_a2 <- do.call(rbind,
                        c(survfit_makeup(kmobj_B_S1, trt_agd),
                          survfit_makeup(kmobj_A_S2, trt_ipd),
                          survfit_makeup(kmobj_Aadj_S2, paste(trt_ipd,"(weighted)")))
    )
    kmdat_a2$treatment <- factor(kmdat_a2$treatment, levels = unique(kmdat_a2$treatment))
    kmdat_a1 <- do.call(rbind,
                        c(survfit_makeup(kmobj_C_S1, paste(trt_common, "(AgD)")),
                          survfit_makeup(kmobj_C_S2, paste(trt_common, "(IPD)")),
                          survfit_makeup(kmobj_Cadj_S2, paste(trt_common,"(IPD,weighted)")))
    )
    kmdat_a1$treatment <- factor(kmdat_a1$treatment, levels = unique(kmdat_a1$treatment))


    # make plot depending on the layout
    if(km_layout == "by_trial"){

      # 1 by 2 plot, each plot is per trial
      subplot_heights <- c(7, 0.7 + 2 * 0.7, 0.8)
      layout_mat <- matrix(1:4, ncol = 2)
      layout(layout_mat, heights = subplot_heights)

      basic_kmplot(kmdat_s2,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_ipd," vs ", trt_common, ") in the IPD trial"),
                   suppress_plot_layout = TRUE, ...)

      basic_kmplot(kmdat_s1,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_agd," vs ", trt_common, ") in the AgD trial"),
                   suppress_plot_layout = TRUE, ...)

    }else if(km_layout == "by_arm"){

      # 1 by 2 plot, by 1 is for investigational arm, the other is for common comparator
      subplot_heights <- c(7, 0.7 + 2 * 0.7, 0.8)
      layout_mat <- matrix(1:4, ncol = 2)
      layout(layout_mat, heights = subplot_heights)

      basic_kmplot(kmdat_a2,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_ipd," vs ", trt_agd, ")"),
                   suppress_plot_layout = TRUE, ...)

      basic_kmplot(kmdat_a1,
                   main_title = paste0("Kaplan-Meier Curves of Common Comparator \n",trt_common,"(IPD vs AgD Trial)"),
                   suppress_plot_layout = TRUE, ...)

    }else{

      # 2 by 2 plot, combine by trial and by arm
      layout_mat <- matrix(1:4, ncol = 2, byrow = TRUE)
      layout(layout_mat)

      basic_kmplot(kmdat_s2,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_ipd," vs ", trt_common, ") in the IPD trial"),
                   show_risk_set = FALSE,
                   suppress_plot_layout = TRUE, ...)

      basic_kmplot(kmdat_s1,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_agd," vs ", trt_common, ") in the AgD trial"),
                   show_risk_set = FALSE,
                   suppress_plot_layout = TRUE, ...)

      basic_kmplot(kmdat_a2,
                   main_title = paste0("Kaplan-Meier Curves \n(",trt_ipd," vs ", trt_agd, ")"),
                   show_risk_set = FALSE,
                   suppress_plot_layout = TRUE, ...)

      basic_kmplot(kmdat_a1,
                   main_title = paste0("Kaplan-Meier Curves of Common Comparator \n",trt_common,"(IPD vs AgD Trial)"),
                   show_risk_set = FALSE,
                   suppress_plot_layout = TRUE, ...)

    }
  }
}


#' Basic Kaplan Meier (KM) plot function
#'
#' This function can generate a basic KM plot with or without risk set table appended at the bottom.
#' In a single plot, it can include up to 4 KM curves. This depends on number of levels in 'treatment' column in the input data.frame \code{kmdat}
#'
#' @param kmdat a data.frame, must consist 'treatment', 'time' (unit in days), 'n.risk', 'censor', 'surv', similar to an output from \code{maicplus:::survfit_makeup}
#' @param time_scale a string, time unit of median survival time, taking a value of 'year', 'month', 'week' or 'day'
#' @param time_grid a numeric vector in the unit of \code{time_scale}, risk set table and x axis of the km plot will be defined based on this time grid
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param subplot_heights a numeric vector, heights argument to \code{graphic::layout()},NULL by default which means user will use the default setting
#' @param suppress_plot_layout logical, suppress the layout setting in this function so that user can specify layout outside of the function, FALSE by default
#' @param use_colors a character vector of length up to 4, colors to the KM curves, it will be passed to 'col' of \code{lines()}
#' @param use_line_types a numeric vector of length up to 4, line type to the KM curves, it will be passed to 'lty' of \code{lines()}
#' @param use_pch_cex a scalar between 0 and 1, point size to indicate censored individuals on the KM curves, it will be passed to 'cex' of \code{points()}
#' @param use_pch_alpha a scalar between 0 and 255, degree of color transparency of points to indicate censored individuals on the KM curves, it will be passed to 'cex' of \code{points()}
#'
#' @example basic_kmplot_ex.R
#'
#' @return a KM plot with or without risk set table appended at the bottom, with up to 4 KM curves
#' @export

basic_kmplot <- function(kmdat, time_scale,
                         time_grid,
                         show_risk_set = TRUE,
                         main_title = "Kaplan-Meier Curves",
                         subplot_heights = NULL,
                         suppress_plot_layout = FALSE,
                         use_colors = NULL,
                         use_line_types = NULL,
                         use_pch_cex = 0.65,
                         use_pch_alpha = 100) {

  time_unit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  # precheck
  if (!length(subplot_heights) %in% c(0, (1 + show_risk_set))) stop(paste("length of subplot_heights should be", (1 + show_risk_set)))
  if (!is.factor(kmdat$treatment)) stop("kmdat$treatment needs to be a factor, its levels will be used in legend and title, first level is comparator")
  if (nlevels(kmdat$treatment) > 4) stop("kmdat$treatment cannot have more than 4 levels")
  if (is.null(time_grid) & show_risk_set) stop("please provide a numeric vector as time_grid to show risk set table")

  # set up x axis (time)
  if (is.null(time_grid)) {
    max_t <- max(kmdat$time)
    t_range <- c(0, (max_t / time_unit[[time_scale]]) * 1.07)
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
  use_col2 <- rgb(use_col2[1,], use_col2[2,], use_col2[3,], alpha = use_pch_alpha, maxColorValue = 255)

  ## : first subplot: KM curve -------------------
  # base plot
  par(bty = "n", tcl = -0.15, mgp = c(1.8, 0.4, 0), cex.lab = 0.85, cex.axis = 0.8, cex.main = 0.9, mar = c(3, 4, 5, 1))
  plot(0, 0,
       type = "n", xlab = paste0("Time in ", time_scale), ylab = "Survival Probability",
       ylim = c(0, 1), xlim = t_range, yaxt = "n", xaxt = "n",
       main = main_title
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
      x = (tmpkmdat$time / time_unit[[time_scale]]),
      col = use_col[ii],
      lty = use_lty[ii],
      lwd = use_lwd[ii],
      type = "s"
    )
    tmpid <- (tmpkmdat$censor != 0) # cannot just ==1, anticipating weighted case
    points(
      y = tmpkmdat$surv[tmpid],
      x = (tmpkmdat$time[tmpid] / time_unit[[time_scale]]),
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
      tmptime <- (tmpkmdat$time / time_unit[[time_scale]])
      tmpnr <- sapply(time_grid, function(kk) {
        tmpid <- which(tmptime > kk)
        if (length(tmpid) == 0) {
          if (min(tmpkmdat$n.risk) == 0) {
            tmpid <- which.min(tmpkmdat$n.risk)[1]
          } else {
            tmpid <- NULL
          }
        }
        tout <- ifelse(is.null(tmpid), "n/a", round(tmpkmdat$n.risk[tmpid],1))
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
}

