# This example code estimates weights for individual patient data from a single
# arm study of 'intervention' based on aggregate baseline characteristics from
# the comparator trial

# Note that matching step is the same regardless the endpoint type (binary, TTE, cont., count)
# Matching by arms are not supported, due to lack of merit as from current literature (or lack of research)

devtools::load_all()
#### load data ----------------------------------------------------------

### IPD
# Read in relevant ADaM data and rename variables of interest
adsl <- read.csv(system.file("extdata", "adsl.csv", package = "maicplus",
                             mustWork = TRUE))
adrs <- read.csv(system.file("extdata", "adrs.csv", package = "maicplus",
                             mustWork = TRUE))
adtte <- read.csv(system.file("extdata", "adtte.csv", package = "maicplus",
                              mustWork = TRUE))

### AgD
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "aggregate_data_example_1.csv",
                                   package = "maicplus", mustWork = TRUE))
# target_pop2 <- read.csv(system.file("extdata", "aggregate_data_example_2.csv",
#                                     package = "maicplus", mustWork = TRUE))
# target_pop3 <- read.csv(system.file("extdata", "aggregate_data_example_3.csv",
#                                     package = "maicplus", mustWork = TRUE))

# for time-to-event endpoints, pseudo IPD from digitalized KM
pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv", package = "maicplus",
                                   mustWork = TRUE))

#### prepare data ----------------------------------------------------------
target_pop <- process_agd(target_pop)
# target_pop2 <- process_agd(target_pop2) # demo of process_agd in different scenarios
# target_pop3 <- process_agd(target_pop3) # demo of process_agd in different scenarios
adsl <- dummize_ipd(adsl,dummize_cols=c("SEX"),dummize_ref_level=c("Female"))
use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

match_res <-  estimate_weights(data=use_adsl,
                               centered_colnames = grep("_CENTERED$",names(use_adsl)),
                               startVal = 0,
                               method = "BFGS")

plot_weights(wt = match_res$data$weights, main_title = "Unscaled Individual Weigths")




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
use_data <- intervention_data[,grepl("_centered$",names(intervention_data))]
use_weigths <- cal_weights(EM = as.matrix(use_data))

plot_weights(use_weigths$wt)
plot_weights(use_weigths$wt.rs, main.title = "Scaled Individual Weights")

### unanchored tte

# emulate IPD
ipd_dat <- adtte
ipd_dat_ext <- pseudo_ipd
ipd_dat_ext$treatment <- "B"
ipd_dat_ext$time <- ipd_dat_ext$Time
ipd_dat_ext$status <- ipd_dat_ext$Event
ipd_dat$treatment <- ipd_dat$ARM
ipd_dat$time <- ipd_dat$Time
ipd_dat$status <- ipd_dat$Event

library(survival)
fit_unanchored <- maic_tte_unanchor(useWt=use_weigths$wt,
                                    dat=ipd_dat,
                                    dat_ext=ipd_dat_ext,
                                    trt="A",
                                    trt_ext="B",
                                    time_scale = "month",
                                    endpoint_name = "OS",
                                    transform = "log")
fit_unanchored$report_median_surv


#**!! write a print method for output object from est_weights














