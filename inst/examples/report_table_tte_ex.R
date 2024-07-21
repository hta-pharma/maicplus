data(adtte_sat)
data(pseudo_ipd_sat)
combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = combined_data)

# Derive median survival time
kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data, conf.type = "log-log")
medSurv <- medSurv_makeup(kmobj, legend = "before matching", time_scale = "day")

report_table_tte(unweighted_cox, medSurv)