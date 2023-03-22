
# This example code estimates weights for individual patient data from a single
# arm study of 'intervention' based on aggregate baseline characteristics from
# the comparator trial

# Note that matching step is the same regardless the endpoint type (binary, TTE, cont., count)
# Matching by arms are not supported, due to lack of merit as from current literature (or lack of research)

devtools::load_all()
library(magrittr)

#### Prepare the data ----------------------------------------------------------

### Intervention data

# Read in relevant ADaM data and rename variables of interest
adsl <- read.csv(system.file("extdata", "adsl.csv", package = "maicplus",
                             mustWork = TRUE))
adrs <- read.csv(system.file("extdata", "adrs.csv", package = "maicplus",
                             mustWork = TRUE))
adtte <- read.csv(system.file("extdata", "adtte.csv", package = "maicplus",
                              mustWork = TRUE))

# Data containing the matching variables
adsl <- within(adsl, SEX <- ifelse(SEX=="Male", 1, 0)) # Coded 1 for males and 0 for females

# Response data
adrs <- with(adrs,{
  tmp <- adrs[PARAM=="Response", ]
  tmp$response <- tmp$AVAL
  tmp[,c("USUBJID","ARM","response")]
})

# Time to event data (overall survival)
adtte <- with(adtte,{
  tmp <- adtte[PARAMCD=="OS", ]
  tmp$Event <- 1-tmp$CNSR
  tmp$Time <- tmp$AVAL
  tmp[,c("USUBJID","ARM","Time","Event")]
})

# Combine all intervention data
intervention_input <- merge(adsl,adrs, by=c("USUBJID", "ARM"), all=TRUE)
intervention_input <- merge(intervention_input,adtte, by=c("USUBJID", "ARM"), all=TRUE)

# List out the variables in the intervention data that have been identified as
# prognostic factors or treatment effect modifiers and will be used in the
# matching
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

## Baseline data from the comparator trial
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "aggregate_data.csv",
                                   package = "maicplus", mustWork = TRUE))
#**!! change of the csv file, to follow our standard naming convention
#**!!


# Rename target population cols to be consistent with match_cov
target_pop <- within(target_pop,{
  N = N
  Treatment=ARM
  AGE=age.mean
  SEX=prop.male
  SMOKE=prop.smoke
  ECOG0=prop.ecog0
})
target_pop


#### Estimate weights ----------------------------------------------------------

### Center baseline characteristics
# (subtract the aggregate comparator data from the corresponding column of
# intervention PLD)
intervention_data <- within(intervention_input,{
  Age_centered = AGE - target_pop$age.mean
  # matching on both mean and standard deviation for continuous variables (optional)
  Age_squared_centered = (AGE^2) - (target_pop$age.mean^2 + target_pop$age.sd^2)
  Sex_centered = SEX - target_pop$prop.male
  Smoke_centered = SMOKE - target_pop$prop.smoke
  ECOG0_centered = ECOG0 - target_pop$prop.ecog0
})


########################################### ROCHE
## Define the matching covariates
cent_match_cov <- c("Age_centered",
                    "Age_squared_centered",
                    "Sex_centered",
                    "Smoke_centered",
                    "ECOG0_centered")

## Optimization procedure
# Following the centering of the baseline characteristics of the intervention
# study, patient weights can be estimated using estimate_weights
# The function output is a list containing (1) a data set of the individual
# patient data with the assigned weights "analysis_data" and (2) a vector
# containing the matching variables "matching_vars"
est_weights <- estimate_weights(intervention_data = intervention_data,
                                matching_vars = cent_match_cov)

# Are the weights sensible? ----------------------------------------------------

# The wt_diagnostics function requires the output from the estimate_weights
# function and will output:
# - the effective sample size (ESS)
# - a summary of the weights and rescaled weights (mean, standard deviation,
#   median, minimum and maximum)
# - a unique set of weights with the corresponding patient profile based on the
#   matching variables

diagnostics <- wt_diagnostics(est_weights$analysis_data,
                              vars = match_cov)
diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles

# Each of the wt_diagnostics outputs can also be estimated individually
ESS <- estimate_ess(est_weights$analysis_data)
weight_summ <- summarize_wts(est_weights$analysis_data)
wts_profile <- profile_wts(est_weights$analysis_data, vars = match_cov)

# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
library(ggplot2)
histogram <- hist_wts(est_weights$analysis_data, bin = 50)
histogram


# Has the optimization worked? -------------------------------------------------

# The following code produces a summary table of the intervention baseline
# characteristics before and after matching compared with the comparator
# baseline characteristics:

check_weights(analysis_data = est_weights$analysis_data, matching_vars = match_cov,
              target_pop_standard = target_pop_standard)


########################################### MSD



















