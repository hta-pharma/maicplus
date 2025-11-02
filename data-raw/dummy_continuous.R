#### create continuous example datasets ####

# devtools::load_all()
library(maicplus)
set.seed(1234)

# Load existing adsl datasets to get matching USUBJIDs
data("adsl_sat", package = "maicplus")
data("adsl_twt", package = "maicplus")

# create adlb_sat (single-arm continuous outcome data)
# Use SAME USUBJIDs as adsl_sat for consistency
adlb_sat <- adsl_sat[, c("USUBJID", "ARM")]
adlb_sat$AVAL <- rnorm(nrow(adlb_sat), mean = -15, sd = 8) # negative mean indicates improvement
adlb_sat$PARAM <- "Change from Baseline in Lab Value"
adlb_sat$PARAMCD <- "CHG"

# create adlb_twt (two-arm continuous outcome data)
# Use SAME USUBJIDs as adsl_twt for consistency
adlb_twt <- adsl_twt[, c("USUBJID", "ARM")]

# Generate different means per arm
adlb_twt$AVAL <- NA
# Note: adsl_twt only has arms A and C
for (arm_name in unique(adlb_twt$ARM)) {
  arm_idx <- adlb_twt$ARM == arm_name
  if (arm_name == "A") {
    adlb_twt$AVAL[arm_idx] <- rnorm(sum(arm_idx), mean = -15, sd = 8)  # Arm A: better outcome
  } else if (arm_name == "C") {
    adlb_twt$AVAL[arm_idx] <- rnorm(sum(arm_idx), mean = -9, sd = 9.5)  # Arm C: less improvement
  }
}

adlb_twt$PARAM <- "Change from Baseline in Lab Value"
adlb_twt$PARAMCD <- "CHG"

### Output
usethis::use_data(adlb_sat, adlb_twt,
  internal = FALSE, overwrite = TRUE
)

