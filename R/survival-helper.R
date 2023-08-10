#' helper function: makeup to get median survival time from a `survival::survfit` object
#'
#' extract and display median survival time with confidence interval
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param legend a character string, name used in 'type' column in returned data frame
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#'
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' medSurv <- medSurv_makeup(kmobj, legend = "before matching", time_scale = "day")
#'
#' @return a data frame with a index column 'type', median survival time and confidence interval
#' @export

medSurv_makeup <- function(km_fit, legend = "before matching", time_scale) {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  # km_fit is the returned object from survival::survfit
  km_fit <- summary(km_fit)$table
  km_fit[, 5:ncol(km_fit)] <- km_fit[, 5:ncol(km_fit)] / timeUnit[[time_scale]]

  toyadd <- data.frame(
    treatment = gsub("treatment=", "", rownames(km_fit)),
    type = rep(legend, 2)
  )

  km_fit <- cbind(toyadd, km_fit)
  rownames(km_fit) <- NULL

  km_fit
}


#' Helper function to select set of variables used for Kaplan-Meier plot
#'
#' @param km_fit returned object from \code{survival::survfit}
#'
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' survobj <- survfit_makeup(kmobj)
#'
#' @return a list of data frames of variables from survfit. Data frame is divided by treatment.
#' @export

survfit_makeup <- function(km_fit) {
  kmdat <- data.frame(
    time = km_fit$time,
    treatment = unlist(mapply(rep, 1:2, each = km_fit$strata)),
    n.risk = km_fit$n.risk,
    n.event = km_fit$n.event,
    censor = km_fit$n.censor,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper,
    cumhaz = km_fit$cumhaz
  )
  kmdat$treatment <- sapply(strsplit(names(km_fit$strata), "="), "[[", 2)[kmdat$treatment]
  split(kmdat, f = kmdat$treatment)
}


#' helper function: KM plot with unadjusted and adjusted KM
#'
#' @param km_fit_before returned object from \code{survival::survfit} before adjustment
#' @param km_fit_after returned object from \code{survival::survfit} after adjustment
#' @param time_scale time unit of median survival time, taking a value of 'year', 'month', 'week' or 'day'
#' @param trt internal trial treatment
#' @param trt_ext external trial treatment
#' @param endpoint_name name of the endpoint
#' @param line_col color of the line curves with the order of external, internal unadjusted, and internal adjusted
#'
#' @return a Kaplan-Meier plot
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' kmobj_adj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, weights = combined_data_tte$weight, conf.type = "log-log")
#' par(cex.main=0.85)
#' km_plot(kmobj, kmobj_adj, time_scale = "month", trt = "A", trt_ext = "B", endpoint_name = "OS")
#' @export

km_plot <- function(km_fit_before, km_fit_after = NULL, time_scale, trt, 
                    trt_ext, endpoint_name = "", line_col = c("#5450E4","#00857C","#6ECEB2")) {
  
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  # prepare data for plot
  pd_be <- survfit_makeup(km_fit_before)
  if (!is.null(km_fit_after)) pd_af <- survfit_makeup(km_fit_after)

  # set up x axis (time)
  if (!is.null(km_fit_after)) {
    max_t <- max(km_fit_before$time, km_fit_after$time)
  } else {
    max_t <- max(km_fit_before$time)
  }
  t_range <- c(0, (max_t / timeUnit[[time_scale]]) * 1.07)

  # base plot
  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(0, 0,
    type = "n", xlab = paste0("Time in ", time_scale), ylab = "Survival Probability",
    ylim = c(0, 1), xlim = t_range, yaxt = "n",
    main = paste0(
      "Kaplan-Meier Curves of Comparator ", ifelse(!is.null(km_fit_after), "(AgD) ", ""),
      "and Treatment ", ifelse(!is.null(km_fit_after), "(IPD)", ""),
      "\nEndpoint: ", endpoint_name
    )
  )
  axis(2, las = 1)

  # add km lines from external trial
  lines(
    y = pd_be[[trt_ext]]$surv,
    x = (pd_be[[trt_ext]]$time / timeUnit[[time_scale]]), col = line_col[1],
    type = "s"
  )
  tmpid <- pd_be[[trt_ext]]$censor == 1
  points(
    y = pd_be[[trt_ext]]$surv[tmpid],
    x = (pd_be[[trt_ext]]$time[tmpid] / timeUnit[[time_scale]]),
    col = line_col[1], pch = 3, cex = 0.7
  )

  # add km lines from internal trial before adjustment
  lines(
    y = pd_be[[trt]]$surv,
    x = (pd_be[[trt]]$time / timeUnit[[time_scale]]), col = line_col[2],
    type = "s"
  )
  tmpid <- pd_be[[trt]]$censor == 1
  points(
    y = pd_be[[trt]]$surv[tmpid],
    x = (pd_be[[trt]]$time[tmpid] / timeUnit[[time_scale]]),
    col = line_col[2], pch = 3, cex = 0.7
  )

  # add km lines from internal trial after adjustment
  if (!is.null(km_fit_after)) {
    lines(
      y = pd_af[[trt]]$surv,
      x = (pd_af[[trt]]$time / timeUnit[[time_scale]]), col = line_col[3], lty = 2,
      type = "s"
    )
    tmpid <- pd_af[[trt]]$censor == 1
    points(
      y = pd_af[[trt]]$surv[tmpid],
      x = (pd_af[[trt]]$time[tmpid] / timeUnit[[time_scale]]),
      col = line_col[3], pch = 3, cex = 0.7
    )
  }

  use_leg <- 1:ifelse(is.null(km_fit_after), 2, 3)
  # add legend
  legend("topright",
    bty = "n", lty = c(1, 1, 2)[use_leg], cex = 0.8, col = c(line_col[1], line_col[2], line_col[3])[use_leg],
    legend = c(
      paste0("Comparator: ", trt_ext),
      paste0("Treatment: ", trt),
      paste0("Treatment: ", trt, " (with weights)")
    )[use_leg]
  )
}

#' Function to plot Kaplan-Meier curves using survminer package
#'
#' @param combined_data combined data of internal IPD and psuedo IPD from aggregate data. 
#' @param trt internal trial treatment
#' @param trt_ext external trial treatment
#' @param trt_common treatment that is shared between internal and external trial. Only applies in anchored comparison.
#' @param break.x.by bin parameter for survminer
#' @param endpoint_name a character name of the endpoint
#' @param censor indicator to include censor information 
#' @param risk.table indicator to include risk table
#'
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' km_plot2(combined_data_tte, trt = "A", trt_ext = "B", censor = TRUE, risk.table = TRUE)
#' @return a Kaplan-Meier plot

km_plot2 <- function(combined_data, trt, trt_ext, trt_common = NULL,
                    break.x.by = 60, endpoint_name = "Overall survival", 
                    censor = TRUE, risk.table = TRUE) {
  
  # check if survminer package is available
  if(!requireNamespace("survminer", quietly = TRUE)){
    stop("survminer package is needed to run this function")
  }
  
  colnames(combined_data) <- c("TIME", "EVENT", "ARM", "weights")
  
  internal <- combined_data[combined_data$ARM == trt,]
  external <- combined_data[combined_data$ARM == trt_ext,]
  
  # Unweighted internal data
  km_internal <- survfit(Surv(TIME, EVENT==1) ~ 1,
                         data = internal,
                         conf.type = "log-log")
  # Weighted internal data
  km_internal_weighted <- survfit(Surv(TIME, EVENT==1) ~ 1,
                                  data = internal,
                                  weights = internal$weights,
                                  conf.type = "log-log")
  # Comparator data
  km_external <- survfit(Surv(TIME, EVENT==1) ~ 1,
                         data = external,
                         conf.type = "log-log")

  # Combine the survfit objects ready for ggsurvplot
  km_list <- list(internal = km_internal,
                  internal_weighted = km_internal_weighted,
                  external = km_external)
  
  #Produce the Kaplan-Meier plot
  survminer_plot <- survminer::ggsurvplot(km_list,
                        linetype = c(1,1,2),
                        size = 0.2,
                        combine = TRUE,
                        risk.table = risk.table,
                        risk.table.y.text.col = T,
                        risk.table.y.text = FALSE,
                        break.x.by = break.x.by, 
                        censor = censor,
                        censor.size = 2,
                        xlab = "Time",
                        ylab = endpoint_name,
                        legend.title = "Treatment",
                        legend = c(0.85,0.82),
                        title = paste0("Kaplan-Meier plot of ",  tolower(endpoint_name)),
                        legend.labs = c("Internal IPD", "Internal IPD weighted", "External comparator"),
                        tables.theme = theme_cleantable(),
                        ggtheme = theme_classic(base_size = 10),
                        fontsize = 3,
                        conf.int = FALSE)
  survminer_plot
}


#' Plot Log Cumulative Hazard Rate
#'
#' a diagnosis plot for proportional hazard assumption, versus log-time (default) or time
#'
#' @param km_fit returned object from \code{survival::survfit}
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @param exclude_censor logical, should censored data point be plotted
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
#' log_cum_haz_plot(kmobj, time_scale = "month", log_time = TRUE, endpoint_name = "OS", subtitle = "(Before Matching)")
#'
#' @return a plot
#' @export
log_cum_haz_plot <- function(km_fit, time_scale, log_time = TRUE, endpoint_name = "", subtitle = "", exclude_censor = TRUE) {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)
  
  clldat <- survfit_makeup(km_fit)
  
  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")
  
  if (exclude_censor) {
    clldat <- lapply(clldat, function(xxt) xxt[xxt$censor == 0, , drop = FALSE])
  }

  all.times <- do.call(rbind, clldat)$time / timeUnit[[time_scale]]
  if (log_time) all.times <- log(all.times)
  t_range <- range(all.times)
  y_range <- range(log(do.call(rbind, clldat)$cumhaz))

  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(0, 0,
    type = "n", xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
    ylab = "Log-Cumulative Hazard Rate",
    ylim = y_range, xlim = t_range, yaxt = "n",
    main = paste0(
      "Diagnostic plot for Proportional Hazard assumption\nEndpoint: ", endpoint_name,
      ifelse(subtitle == "", "", "\n"), subtitle
    )
  )
  axis(2, las = 1)

  trts <- names(clldat)
  cols <- c("dodgerblue3", "firebrick3")
  pchs <- c(1, 4)
  for (i in seq_along(clldat)) {
    use_x <- (clldat[[i]]$time / timeUnit[[time_scale]])
    if (log_time) use_x <- log(use_x)

    lines(
      y = log(clldat[[i]]$cumhaz),
      x = use_x, col = cols[i]
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


#' Plot Schoenfeld residuals for a Cox model fit
#'
#' @param coxobj object returned from \code{\link[survival]{coxph}}
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#' @param log_time logical, TRUE (default) or FALSE
#' @param endpoint_name a character string, name of the endpoint
#' @param subtitle a character string, subtitle of the plot
#' @examples
#' load(system.file("extdata","combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
#' unweighted_cox <- coxph(Surv(TIME, EVENT==1) ~ ARM, data = combined_data_tte)
#' resid_plot(unweighted_cox, time_scale = "month", log_time = TRUE, endpoint_name = "OS", subtitle = "(Before Matching)")
#'
#' @return a plot
#' @export

resid_plot <- function(coxobj, time_scale = "month", log_time = TRUE, endpoint_name = "", subtitle = "") {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)
  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  schresid <- residuals(coxobj, type = "schoenfeld")
  plotX <- as.numeric(names(schresid)) / timeUnit[[time_scale]]
  if (log_time) plotX <- log(plotX)
  par(mfrow = c(1, 1), bty = "n", tcl = -0.15, mgp = c(2.3, 0.5, 0))
  plot(schresid ~ plotX,
    cex = 0.9, col = "navyblue", yaxt = "n",
    ylab = "Unscaled Schoenfeld Residual", xlab = paste0(ifelse(log_time, "Log-", ""), "Time in ", time_scale),
    main = paste0(
      "Diagnosis Plot: Unscaled Schoenfeld Residual\nEndpoint: ", endpoint_name,
      ifelse(subtitle == "", "", "\n"), subtitle
    )
  )
  axis(2, las = 1)
}
