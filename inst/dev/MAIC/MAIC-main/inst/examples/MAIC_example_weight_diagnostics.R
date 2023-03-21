
# This example code uses the weighted individual patient data, outputted from
# the estimate_weights function to perform weight diagnostics. The weighted data
# is saved within est_weights. To check the weighted aggregate baseline
# characteristics for 'intervention' match those in the comparator data,
# standardized data "target_pop_standard" is used. Please see the package
# vignette for more information on how to use the estimate_weights function and
# derive the "target_pop_standard" data.

library(dplyr)
library(MAIC)

# load est_weights
data(est_weights, package = "MAIC")

# load target_pop_standard
data(target_pop_standard, package = "MAIC")

# List out the uncentered variables used in the matching
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

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
histogram <- hist_wts(est_weights$analysis_data, bin = 50)
histogram


# Has the optimization worked? -------------------------------------------------

# The following code produces a summary table of the intervention baseline
# characteristics before and after matching compared with the comparator
# baseline characteristics:

check_weights(analysis_data = est_weights$analysis_data, matching_vars = match_cov, 
                          target_pop_standard = target_pop_standard)
