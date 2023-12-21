#' Kaplan-Meier (KM) plot function for anchored and unanchored cases
#'
#' This is wrapper function of \code{basic_kmplot}.
#' The argument setting is similar to \code{maic_anchored} and \code{maic_unanchored},
#' and it is used in those two functions.
#'
#' @param ipd_weights an object returned by \code{estimate_weight}
#' @param tte_dat_ipd a data frame of individual patient data (IPD) of internal trial, contain at least "USUBJID", "EVENT", "TIME" columns and a column indicating treatment assignment
#' @param ipd_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param tte_dat_pseudo a data frame of pseudo IPD by digitized KM curves of external trial (for time-to-event endpoint), contain at least "EVENT", "TIME"
#' @param pseudo_trt_var a string, column name in \code{dat_ipd} that contains the treatment assignment
#' @param trt_ipd  a string, name of the interested investigation arm in internal trial \code{dat_igd} (real IPD)
#' @param trt_pseudo a string, name of the interested investigation arm in external trial \code{dat_pseudo} (pseudo IPD)
#' @param trt_common a string, name of the common comparator in internal and external trial, by default is NULL, indicating unanchored case
#' @param km_conf_type a string, pass to \code{conf.type} of \code{survfit}
#' @param km_layout a string, only applicable for unachored case (\code{trt_common = NULL}), indicated the desired layout of output KM curve.
#' @param ... other arguments in \code{basic_kmplot}
#'
#' @return
#' In unanchored case, a KM plot with risk set table. In anchored case, depending on \code{km_layout},
#' \itemize{
#'   \item if "by_trial", 2 by 1 plot, first all KM curves (incl. weighted) in IPD trial, and then KM curves in AgD trial, with risk set table.
#'   \item if "by_arm", 2 by 1 plot, first KM curves of \code{trt_pseudo} and  \code{trt_ipd} (with and without weights), and then KM cuvers of \code{trt_common} in AgD trial and IPD trial (with and without weights). Risk set table is appended.
#'   \item if "all", 2 by 2 plot, all plots in "by_trial" and "by_arm" without risk set table appended.
#' }
#' @export

kmplot2 <- function(ipd_weights,
                   tte_dat_ipd,
                   tte_dat_pseudo,
                   trt_ipd,
                   trt_pseudo,
                   trt_common = NULL,
                   ipd_trt_var = "ARM",
                   pseudo_trt_var = "ARM",
                   km_conf_type = "log-log",
                   km_layout = c("all", "by_trial", "by_arm"),
                   ...) {
  
  names(tte_dat_ipd) <- toupper(names(tte_dat_ipd))
  names(tte_dat_pseudo) <- toupper(names(tte_dat_pseudo))
  ipd_trt_var <- toupper(ipd_trt_var)
  pseudo_trt_var <- toupper(pseudo_trt_var)
  
  # pre check
  if(!"maicplus_estimate_weights" %in% class(ipd_weights)) stop("ipd_weights should be an object returned by estimate_weights")
  if (!all(c("USUBJID", "TIME", "EVENT", ipd_trt_var) %in% names(tte_dat_ipd))) stop(paste("tte_dat_ipd needs to include at least USUBJID, TIME, EVENT,", ipd_trt_var))
  if (!all(c("TIME", "EVENT", pseudo_trt_var) %in% names(tte_dat_pseudo))) stop(paste("tte_dat_pseudo needs to include at least TIME, EVENT,", pseudo_trt_var))
  km_layout <- match.arg(km_layout, choices = c("all", "by_trial", "by_arm"), several.ok = FALSE)
  
  # preparing data
  is_anchored <- ifelse(is.null(trt_common), FALSE, TRUE)
  tte_dat_ipd <- tte_dat_ipd[tte_dat_ipd[[ipd_trt_var]] %in% c(trt_ipd, trt_common), , drop = TRUE]
  tte_dat_pseudo <- tte_dat_pseudo[tte_dat_pseudo[[pseudo_trt_var]] %in% c(trt_pseudo, trt_common), , drop = TRUE]
  tte_dat_ipd$weights <- ipd_weights$data$weights[match(ipd_weights$data$USUBJID, tte_dat_ipd$USUBJID)]
  tte_dat_pseudo$weights <- 1
  
  is_anchored <- ifelse(is.null(trt_common), FALSE, TRUE)
  
  if (!is_anchored) {
    ## unanchored case
    kmobj_B <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", pseudo_trt_var)),
      data = tte_dat_pseudo,
      conf.type = km_conf_type
    )
    kmobj_A <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
      data = tte_dat_ipd,
      conf.type = km_conf_type
    )
    kmobj_A_adj <- survfit(as.formula(paste("Surv(TIME, EVENT) ~", ipd_trt_var)),
      data = tte_dat_ipd,
      conf.type = km_conf_type,
      weights = weights
    )

    kmlist <- list(
      kmobj_B = kmobj_B,
      kmobj_A = kmobj_A,
      kmobj_A_adj = kmobj_A_adj
    )
#    kmdat$treatment <- factor(kmdat$treatment, levels = unique(kmdat$treatment))

    basic_kmplot2(kmlist)
  } 
}

#' Basic Kaplan Meier (KM) plot function using ggplot
#'
#' This function generates a basic KM plot using ggplot.
#' 
#' @param kmdat a data.frame, must consist 'treatment', 'time' (unit in days), 'n.risk', 'censor', 'surv', similar to an output from \code{maicplus:::survfit_makeup}
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param show_risk_set logical, show risk set table or not, TRUE by default
#' @param main_title a string, main title of the KM plot
#' @param break.x.by bin parameter for survminer
#' @param censor indicator to include censor information

basic_kmplot2 <- function(kmdat,
                          endpoint_name = "Time to Event Endpoint",
                          show_risk_set = TRUE,
                          main_title = "Kaplan-Meier Curves",
                          break.x.by = 60,
                          censor = T,
                          xlim = NULL
                          ){

  # Produce the Kaplan-Meier plot
  survminer_plot <- survminer::ggsurvplot(kmlist,
    linetype = c(1, 1, 2),
    size = 0.2,
    combine = TRUE,
    risk.table = T,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    break.x.by = break.x.by,
    censor = censor,
    censor.size = 2,
    xlab = "Time",
    ylab = endpoint_name,
    legend.title = "Treatment",
    legend = c(0.85, 0.82),
    title = paste0(main_title, "\nEndpoint:", endpoint_name),
    legend.labs = c("Internal IPD", "Internal IPD weighted", "External comparator"),
    tables.theme = theme_cleantable(),
    ggtheme = theme_classic(base_size = 10),
    fontsize = 3,
    conf.int = FALSE,
    xlim = xlim,
  )
  survminer_plot

}

#' Basic Kaplan Meier (KM) plot function
#'
#' This function can generate a basic KM plot with or without risk set table appended at the bottom.
#' In a single plot, it can include up to 4 KM curves. This depends on number of levels in 'treatment' column in the input data.frame \code{kmdat}
#'
#' @param kmlist a list of survival::survfit object
#' @param endpoint_name a string, name of time to event endpoint, to be show in the last line of title
#' @param time_scale a string, time unit of median survival time, taking a value of 'years', 'months', 'weeks' or 'days'
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
#' @return a KM plot with or without risk set table appended at the bottom, with up to 4 KM curves
#' @export

# basic_kmplot2 <- function(kmlist,
#                          endpoint_name = "Time to Event Endpoint",
#                          time_scale = NULL,
#                          time_grid = NULL,
#                          show_risk_set = TRUE,
#                          main_title = "Kaplan-Meier Curves",
#                          subplot_heights = NULL,
#                          suppress_plot_layout = FALSE,
#                          use_colors = NULL,
#                          use_line_types = NULL,
#                          use_pch_cex = 0.65,
#                          use_pch_alpha = 100) {
#   
#   
#  
# }





# gg_default <- survfit2(Surv(TIME, EVENT) ~ treatment, weights = weights, data = kmdat) %>%
#   ggsurvfit() +
#   add_censor_mark()
# 
# add_risktable(
#   risktable_stats = c("{format(round(n.risk), nsmall = 0)}",
#                       "{format(round(n.event), nsmall = 0)}"),
#   stats_label = c("patients at risk",
#                   "events")) +
#   gg_default

# tte_dat_ipd_weighted_selected <- tte_dat_ipd[,c("TIME", "EVENT", ipd_trt_var, "weights")]
# colnames(tte_dat_ipd_weighted_selected) <- c("TIME", "EVENT", "treatment", "weights")
# 
# tte_dat_ipd_selected <- tte_dat_ipd_weighted_selected
# tte_dat_ipd_selected$weights <- 1
# 
# tte_dat_ipd_weighted_selected$treatment <- paste(tte_dat_ipd_weighted_selected$treatment,
#                                                  "(weighted)")
# 
# tte_dat_pseudo_selected <- tte_dat_pseudo[,c("TIME", "EVENT", pseudo_trt_var, "weights")]
# colnames(tte_dat_pseudo_selected) <- c("TIME", "EVENT", "treatment", "weights")
# 
# kmdat <- rbind(tte_dat_ipd_selected,
#                tte_dat_ipd_weighted_selected,
#                tte_dat_pseudo_selected)