## Example: weighted data

# library(maicplus)
devtools::load_all()
load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))

kmobj_B <- survfit(Surv(TIME, EVENT) ~ ARM,
                   data = combined_data_tte,
                   conf.type = "log-log",
                   subset = (ARM == "B"))
kmobj_A <- survfit(Surv(TIME, EVENT) ~ ARM,
                   data = combined_data_tte,
                   conf.type = "log-log",
                   subset = (ARM == "A"))
kmobj_A_adj <- survfit(Surv(TIME, EVENT) ~ ARM,
                   data = combined_data_tte,
                   conf.type = "log-log",
                   weights = weights,
                   subset = (ARM == "A" & weights !=  1))

kmdat <- do.call(rbind,
                 c(survfit_makeup(kmobj_B, "B"),
                   survfit_makeup(kmobj_A, "A"),
                   survfit_makeup(kmobj_A_adj, "A (weighted)"))
         )
kmdat$treatment <- factor(kmdat$treatment, levels = unique(kmdat$treatment))


basic_kmplot(kmdat, time_scale="month",
             time_grid = seq(0, 20, by =2),
             show_risk_set = TRUE,
             main_title = "Kaplan-Meier Curves",
             subplot_heights = NULL,
             suppress_plot_layout = FALSE,
             use_colors = NULL,
             use_line_types = NULL,
             use_pch_cex = 0.65,
             use_pch_alpha = 100)




