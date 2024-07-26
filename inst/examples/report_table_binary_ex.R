data(adrs_sat)
testdat <- data.frame(Yes = 280, No = 120)
rownames(testdat) <- "B"
pseudo_ipd_binary_sat <- get_pseudo_ipd_binary(
  binary_agd = testdat,
  format = "unstacked"
)
combined_data <- rbind(adrs_sat[, c("USUBJID", "RESPONSE", "ARM")], pseudo_ipd_binary_sat)
combined_data$ARM <- as.factor(combined_data$ARM)

binobj_dat <- glm(RESPONSE ~ ARM, combined_data, family = binomial(link = "logit"))
report_table_binary(binobj_dat, eff_measure = "OR")
