###############################################################################
# example code for anchored MAIC using polyMAIC for weight generation and ESS calculation
# for time-to-event endpoint
#
# Data came from maicplus package: https://github.com/hta-pharma/maicplus 
# Weights and ESS came from polyMAIC package: https://github.com/Numerus-Ltd/polyMAIC/blob/main/R/polymaic.R
# polyMAIC method: https://pubmed.ncbi.nlm.nih.gov/35506464/
#
# Author: Jinjie Liu
# 
# Date created: 1 July 2025
###############################################################################
###### install required packages
devtools::install_github("Numerus-Ltd/polyMAIC")
library(polyMAIC)
if(!require(maicplus)) {install.packages("maicplus"); library(maicplus)}
if(!require(survminer)) {install.packages("survminer"); library(survminer)}
if(!require(survival)) {install.packages("survival"); library(survival)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if(!require(tibble)) {install.packages("tibble"); library(tibble)}
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
#data(centered_ipd_sat)
data(adtte_sat)
data(pseudo_ipd_sat)
data(adtte_twt)
data(pseudo_ipd_sat)
data(adrs_sat)
data(adsl_sat)
data(adrs_twt)
data(adsl_twt)
data(agd)

###### polymaic example from github: https://github.com/Numerus-Ltd/polyMAIC/blob/main/R/polymaic.R
#' ## One continuous and one binary variable
# TARGS_A    <- c(AGEMEAN = 55, AGESD = 6, SEXFPROP = 0.65)          ## targets
# TOLS_AN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005) ## tolerances (narrow)
# TOLS_AW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01)         ## tolerances (wide)
# VARSTAT_A  <- c(AGE = "MEAN AND SD", SEXF = "MEAN")                ## Variables and measures to match on
# TYPE_A     <- c(AGE = "CON", SEXF = "BIN")                         ## Type of variables to match on
# 
# polymaic_an <- polymaic(IPD=ScenarioA, VARSTAT=VARSTAT_A, TARGS=TARGS_A, TOLS=TOLS_AN, TYPE=TYPE_A, ID="Rnum")
# summary(polymaic_an)
# 
# ## Three continuous and one binary variable
# TARGS_C    <- c(AGEMEAN = 49.5, AGESD = 7.2, SEXFPROP = 0.683, PAINMEAN = 7.5, PAINSD = 1.5, TIMEDMEAN = 6.1, TIMEDSD=3.8)
# TOLS_CN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005, PAINMEAN = 0.005, PAINSD = 0.005, TIMEDMEAN = 0.01, TIMEDSD = 0.01)
# TOLS_CW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01, PAINMEAN = 0.2, PAINSD = 0.1, TIMEDMEAN = 0.2, TIMEDSD = 0.2)
# VARSTAT_C  <- c(AGE = "MEAN AND SD", SEXF = "MEAN", PAIN = "MEAN AND SD", TIMED = "MEAN AND SD")
# TYPE_C     <- c(AGE = "CON", SEXF = "BIN", PAIN = "CON", TIMED="CON")
# data(ScenarioC)
# polymaic_cn <- polymaic(IPD=ScenarioC, VARSTAT=VARSTAT_C, TARGS=TARGS_C, TOLS=TOLS_CN, TYPE=TYPE_C, ID="Rnum")
# summary(polymaic_an)

######## Use maicplus data #########
##                                ##
##                                ##
########                   #########
## ipd data
ipd <-  adsl_twt

## data processing

run_polymaic <- function(ipd, agd, match_config, id_col = "USUBJID", default_tolerance = 0.005) {
  agd <- process_agd(agd)
  
  # Handle any median-based variables by creating binary columns
  median_vars <- match_config %>% filter(match_method == "MEDIAN")
  
  if (nrow(median_vars) > 0) {
    for (i in 1:nrow(median_vars)) {
      var_info <- median_vars[i, ]
      ipd_var_name <- var_info$ipd_var
      agd_median_col <- var_info$agd_col
      new_col_name <- paste0(ipd_var_name, "_ABOVE_MEDIAN")
      
      ipd[[new_col_name]] <- ifelse(ipd[[ipd_var_name]] > agd[[agd_median_col]], 1, 0)
    }
  }
  
  # data processing for polyMAIC input
  VARSTAT <- c()
  TYPE <- c()
  TARGS <- c()
  TOLS <- c() 
  
  for (i in 1:nrow(match_config)) {
    var_info <- match_config[i, ]
    ipd_var <- var_info$ipd_var
    
    # Use the specific tolerance if provided, otherwise use the default
    current_tol <- if ("tolerance" %in% names(var_info) && !is.na(var_info$tolerance)) {
      var_info$tolerance
    } else {
      default_tolerance
    }
    
    if (var_info$match_method == "MEAN_SD") {
      VARSTAT[ipd_var] <- "MEAN AND SD"
      TYPE[ipd_var] <- "CON"
      
      mean_name <- paste0(ipd_var, "MEAN")
      sd_name <- paste0(ipd_var, "SD")
      mean_name1 <- paste0(ipd_var, "_MEAN")
      sd_name1 <- paste0(ipd_var, "_SD")
      
      TARGS[mean_name] <- agd[[mean_name1]]
      TARGS[sd_name] <- agd[[sd_name1]]
      
      # Assign tolerance for both mean and sd
      TOLS[mean_name] <- current_tol
      TOLS[sd_name] <- current_tol
      
    } else if (var_info$match_method == "MEAN") {
      # For binary proportions
      VARSTAT[ipd_var] <- "MEAN"
      TYPE[ipd_var] <- "BIN"
      
      mean_name <- paste0(ipd_var, "MEAN")
      TARGS[mean_name] <- agd[[var_info$agd_col]]
      TOLS[mean_name] <- current_tol
      
    } else if (var_info$match_method == "MEDIAN") {
      new_var_name <- paste0(ipd_var, "_ABOVE_MEDIAN")
      VARSTAT[new_var_name] <- "MEAN"
      TYPE[new_var_name] <- "BIN"
      
      mean_name <- paste0(new_var_name, "MEAN")
      TARGS[mean_name] <- 0.5 
      TOLS[mean_name] <- current_tol
    }
  }
  
  # run polyMAIC
  polymaic_results <- polymaic(
    IPD = ipd,
    VARSTAT = VARSTAT,
    TARGS = TARGS,
    TOLS = TOLS,
    TYPE = TYPE,
    ID = "USUBJID"
  )
  
  return(polymaic_results)
}

# Define the matching variables and other attributes. for mean+sd, leave agd_col to be empty
match_config_table <- tribble(
  ~ipd_var,    ~match_method, ~agd_col,             ~tolerance,
  "AGE",       "MEAN_SD",     NA,                   0.005,       
  "SEX_MALE",  "MEAN",        "SEX_MALE_PROP",      0.005,     
  "SMOKE",     "MEAN",        "SMOKE_PROP",         0.005,         
  "ECOG0",     "MEAN",        "ECOG0_PROP",         0.005,        
  "N_PR_THER", "MEDIAN",      "N_PR_THER_MEDIAN",   0.005      
)


maic_results <- run_polymaic(
  ipd = ipd,
  agd = agd,
  match_config = match_config_table
)

summary(maic_results)

# run maicplus and replace the weight with polyMAIC weights
ipd_centered <- center_ipd(ipd = ipd, agd = agd)
centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")

weighted_data <- estimate_weights(
  data = centered_ipd_twt,
  centered_colnames = centered_colnames
)

#update ESS and weights
weighted_data$data$weights <- maic_results$weights$WEIGHTS
weighted_data$ess <- maic_results$ESS

# inferential result
result <- maic_anchored(
  weights_object = weighted_data,
  ipd = adtte_twt,
  pseudo_ipd = pseudo_ipd_twt,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  normalize_weight = FALSE,
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)

result$inferential$summary

