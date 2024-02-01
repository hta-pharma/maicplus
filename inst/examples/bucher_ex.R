set.seed(1)
ipd <- generate_survival_data() # fix column names and time
colnames(ipd) <- c("USUBJID", "TIME", "EVENT", "ARM", "X1", "X2", "X3", "X4")
ipd$TIME <- ipd$TIME * 365
ipd$ARM <- ifelse(ipd$ARM == 1, "A", "C")
ipd$ARM <- stats::relevel(as.factor(ipd$ARM), ref = "C")

adtte <- ipd[,c("USUBJID", "TIME", "EVENT", "ARM")]

agd <- data.frame(
  STUDY = "Simulation study",
  ARM = "Total",
  N = 150,
  X1_PROP = 0.5,
  X2_PROP = 0.5,
  X3_PROP = 0.5,
  X4_PROP = 0.5)

BC_reported <- list(RR = 0.4354202, RR_SE = 0.0939458)

ipd_centered <- center_ipd(ipd = ipd, agd = agd)
head(ipd_centered)

centered_colnames <- paste0("X", 1:4, "_CENTERED")

weights_object <- estimate_weights(
  data = ipd_centered,
  centered_colnames = centered_colnames
)

ipd_matched <- weights_object$data

# find weighted relative effect of AC
weighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = ipd_matched, weights = ipd_matched$weights)
AC_weighted <- find_RR(weighted_cox)

weighted_bucher <- bucher(trt = AC_weighted, com = BC_reported)
print(weighted_bucher)