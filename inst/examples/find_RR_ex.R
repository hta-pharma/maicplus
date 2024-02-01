
set.seed(1)
ipd <- generate_survival_data() # fix column names and time
colnames(ipd) <- c("USUBJID", "TIME", "EVENT", "ARM", "X1", "X2", "X3", "X4")
ipd$TIME <- ipd$TIME * 365
ipd$ARM <- ifelse(ipd$ARM == 1, "A", "C")
ipd$ARM <- stats::relevel(as.factor(ipd$ARM), ref = "C")

adtte <- ipd[,c("USUBJID", "TIME", "EVENT", "ARM")]

unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = ipd)
RR <- find_RR(unweighted_cox)