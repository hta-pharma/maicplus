library(maicplus)

# Load data
data(centered_ipd_sat)
data(adlb_sat)

# Estimate weights
centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")

weighted_data <- estimate_weights(
  data = centered_ipd_sat,
  centered_colnames = centered_colnames
)

# Create pseudo IPD for external comparator (Arm B)
pseudo_adlb <- get_pseudo_ipd_continuous(
  continuous_agd = data.frame(
    ARM = "B",
    MEAN = -12,
    SD = 9,
    N = 120
  )
)

#### Unanchored MAIC for Mean Difference (MD) ####
maic_unanchored_results_MD <- maic_unanchored(
  weights_object = weighted_data,
  ipd = adlb_sat,
  pseudo_ipd = pseudo_adlb,
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Change from Baseline",
  endpoint_type = "continuous",
  eff_measure = "MD"
)

# View results
maic_unanchored_results_MD$inferential$summary

#### Unanchored MAIC for Standardized Mean Difference (SMD) ####
maic_unanchored_results_SMD <- maic_unanchored(
  weights_object = weighted_data,
  ipd = adlb_sat,
  pseudo_ipd = pseudo_adlb,
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Change from Baseline",
  endpoint_type = "continuous",
  eff_measure = "SMD",
  smd_type = "hedges"
)

maic_unanchored_results_SMD$inferential$summary

