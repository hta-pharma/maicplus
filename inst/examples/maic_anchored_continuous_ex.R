library(maicplus)

# Load data
data(centered_ipd_twt)
data(adlb_twt)

# Estimate weights
centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")

weighted_data <- estimate_weights(
  data = centered_ipd_twt,
  centered_colnames = centered_colnames
)

# Create pseudo IPD for external comparator (Arm B and C)
pseudo_adlb <- get_pseudo_ipd_continuous(
  continuous_agd = data.frame(
    ARM = c("B", "C"),
    MEAN = c(-12, -9),
    SD = c(9, 9.5),
    N = c(120, 150)
  )
)

#### Anchored MAIC for Mean Difference (MD) ####
maic_anchored_results_MD <- maic_anchored(
  weights_object = weighted_data,
  ipd = adlb_twt,
  pseudo_ipd = pseudo_adlb,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Change from Baseline",
  endpoint_type = "continuous",
  eff_measure = "MD"
)

# View results
maic_anchored_results_MD$inferential$summary

#### Anchored MAIC for Standardized Mean Difference (SMD) ####
maic_anchored_results_SMD <- maic_anchored(
  weights_object = weighted_data,
  ipd = adlb_twt,
  pseudo_ipd = pseudo_adlb,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Change from Baseline",
  endpoint_type = "continuous",
  eff_measure = "SMD",
  smd_type = "cohen"
)

maic_anchored_results_SMD$inferential$summary

