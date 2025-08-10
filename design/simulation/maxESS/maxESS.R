
# R packages
library(maicChecks)
library(maicplus)
library(dplyr)
library(limSolve) # downloaded tar.gz because it was archived https://cran.r-project.org/web/packages/limSolve/index.html
library(stringr)

# load in function file maxESS_functions from Jackson et al paper
setwd("C:/Users/swj88/Documents/Github/maicplus/design/simulation/maxESS")
source("maxESS_functions.R") 
source("mew_functions.R")

# function to use in pipeline
maxESS <- function(ipd, adtte, agd, match_mean, id = "USUBJID"){
  
  # calculate jackson's maxESS weights
  agg_means <- agd
  agg_means$N <- NULL
  agg_colnames <- stringr::str_split(colnames(agg_means), "_")
  agg_means <- unlist(agg_means)
  names(agg_means) <- sapply(agg_colnames, function(x) { paste(x[-length(x)], collapse = "_")})
  
  jackson_result <- jackson_MAIC(the_data = ipd, id = "USUBJID", match_mean = match_mean, agg_means = agg_means, agg_sds = NA, 0.000000000001)
  
  weights_jackson <- jackson_result$the_data$weight.opt
  ipd$weights_jackson <- weights_jackson
  
  ipd_trimmed <- ipd %>% 
    filter(weights_jackson != 0)
  
  adtte_trimmed <- adtte %>%
    filter(weights_jackson != 0)
  
  ipd_centered <- center_ipd(ipd = ipd_trimmed, agd = agd)
  centered_colnames <- paste0(match_mean, "_CENTERED")
  
  weighted_sat <- estimate_weights(
    data = ipd_centered,
    centered_colnames = centered_colnames
  )
  weighted_sat$data$weights <- ipd_trimmed$weights_jackson # replace with jackson's weights
  
  #calculate ess
  #sum(jackson$the_data$weight.opt)^2 / sum(jackson$the_data$weight.opt^2)
  
  result <- maic_unanchored(
    weights_object = weighted_sat,
    ipd = adtte_trimmed,
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
  
  return(result)
}


compare_jackson_result <- function(ipd, agd, match_mean, id = "USUBJID"){
  
  # calculate jackson's maxESS weights
  agg_means <- agd
  agg_means$N <- NULL
  agg_colnames <- stringr::str_split(colnames(agg_means), "_")
  agg_means <- unlist(agg_means)
  names(agg_means) <- sapply(agg_colnames, function(x) { paste(x[-length(x)], collapse = "_")})

  jackson_result <- jackson_MAIC(the_data = ipd, id = "USUBJID", match_mean = match_mean, agg_means = agg_means, agg_sds = NA, 0.000000000001)
  
  weight.maic <- jackson_result$the_data$weight.maic
  weight.opt <- jackson_result$the_data$weight.opt
  mew_result <- calculate_mew_weights(X_internal = as.matrix(ipd[,match_mean]), X_external_means = agg_means)
    
  return(cbind(weight.maic = weight.maic, weight.opt = weight.opt, mew_result = as.numeric(mew_result)))
}
  

########## input

data("adsl_sat")
data("adtte_sat")

adsl_sat <- adsl_sat %>%
  mutate(SEX_MALE = ifelse(SEX == "Male", 1, 0)) %>%
  mutate(AGE_SQUARED = AGE^2)

ipd <- adsl_sat
adtte <- adtte_sat

agd <- data.frame(
  N = 300,
  AGE_MEAN = 51,
  SEX_MALE_PROP = 147 / 300,
  ECOG0_PROP = 0.40,
  SMOKE_PROP = 58 / (300 - 5)
)

match_mean = c("AGE","SEX_MALE", "ECOG0", "SMOKE")


## call the function
result <- maxESS(ipd, adtte, agd, match_mean, id = "USUBJID")
result$inferential$summary # just look at adjusted value. Original value is wrong because data is trimmed

weight_comparison <- data.frame(compare_jackson_result(ipd, agd, match_mean))

weight_comparison <- weight_comparison %>%
  mutate(mew_sw1 = normalize_weights(mew_result, method = "SW1"))

# jackson_result (normalizes to add up to 1; SW1) constrains weights to be 
# positive in the opitmization
# mew estimates ess weights and then forces weights to be positive
# if there are no negative weights, two methods are equal