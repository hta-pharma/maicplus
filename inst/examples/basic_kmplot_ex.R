library(survival)
load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
kmdat <- do.call(rbind, survfit_makeup(kmobj))
kmdat$treatment <- factor(kmdat$treatment)

# without risk set table
basic_kmplot(kmdat, time_scale="month",
             time_grid = seq(0, 20, by =2),
             show_risk_set = FALSE,
             main_title = "Kaplan-Meier Curves",
             subplot_heights = NULL,
             suppress_plot_layout = FALSE,
             use_colors = NULL,
             use_line_types = NULL)

# with risk set table
basic_kmplot(kmdat, time_scale="month",
             time_grid = seq(0, 20, by =2),
             show_risk_set = TRUE,
             main_title = "Kaplan-Meier Curves",
             subplot_heights = NULL,
             suppress_plot_layout = FALSE,
             use_colors = NULL,
             use_line_types = NULL)
