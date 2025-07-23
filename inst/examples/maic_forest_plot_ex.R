if (requireNamespace("ggplot2") && requireNamespace("patchwork")) {
  data(centered_ipd_sat)
  data(adtte_sat)
  data(pseudo_ipd_sat)

  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(
    data = centered_ipd_sat,
    centered_colnames = centered_colnames
  )

  result <- maic_unanchored(
    weights_object = weighted_data,
    ipd = adtte_sat,
    pseudo_ipd = pseudo_ipd_sat,
    trt_ipd = "A",
    trt_agd = "B",
    normalize_weight = FALSE,
    endpoint_name = "Overall Survival",
    endpoint_type = "tte",
    eff_measure = "HR",
    time_scale = "month",
    km_conf_type = "log-log"
  )

  maic_forest_plot(result)
}
