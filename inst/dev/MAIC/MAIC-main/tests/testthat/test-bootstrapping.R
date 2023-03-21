
test_that("Test bootstrapping functions", {
  
  require(dplyr)
  require(survival)
  
  # load intervention data with weights saved in est_weights
  data(est_weights, package = "MAIC")
  
  # Combine data -----------------------------------------------------------------
  
  # Combine the the comparator pseudo data with the analysis data, outputted from
  # the estimate_weights function
  
  # Read in digitised pseudo survival data for the comparator
  comparator_surv <- read.csv(system.file("extdata", "psuedo_IPD.csv",
                                          package = "MAIC", mustWork = TRUE)) %>%
    rename(Time=Time, Event=Event)
  
  # Simulate response data based on the known proportion of responders
  comparator_n <- nrow(comparator_surv) # total number of patients in the comparator data
  comparator_prop_events <- 0.4 # proportion of responders
  # Calculate number with event
  # Use round() to ensure we end up with a whole number of people
  # number without an event = Total N - number with event to ensure we keep the
  # same number of patients
  n_with_event <- round(comparator_n*comparator_prop_events, digits = 0)
  comparator_binary <- data.frame("response"= c(rep(1, n_with_event), rep(0, comparator_n - n_with_event)))
  
  
  # Join survival and response comparator data
  # Note the rows do not represent observations from a particular patient
  # i.e. there is no relationship between the survival time and response status
  # for a given row since this is simulated data
  # In a real data set this relationship would be present
  comparator_input <- cbind(comparator_surv, comparator_binary) %>%
    mutate(wt=1, wt_rs=1, ARM="Comparator") # All patients have weight = 1
  
  # Join comparator data with the intervention data (naming of variables should be
  # consistent between both datasets)
  # Set factor levels to ensure "Comparator" is the reference treatment
  combined_data <-  bind_rows(est_weights$analysis_data, comparator_input)
  combined_data$ARM <- relevel(as.factor(combined_data$ARM), ref="Comparator")
  
  
  unweighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = combined_data)
  
  # Fit a Cox model with weights to estimate the weighted HR
  weighted_cox <- coxph(Surv(Time, Event==1) ~ ARM,
                        data = combined_data,
                        weights = wt)
  
  # check the bootstrap function returns same value on full data
  boot_test_HR <- bootstrap_HR(intervention_data =  est_weights$analysis_data, 
                               i =c(1:500),
                               comparator_data = comparator_input, # comparator pseudo data
                               matching = est_weights$matching_vars, # matching variables
                               model = Surv(Time, Event==1) ~ ARM # model to fit
  )
  
  # check matches to normal weighted cox on full data
  expect_equal(exp(weighted_cox$coefficients), boot_test_HR, tolerance = 0.01)
  
  
  # Fit a logistic regression model with weights to estimate the weighted OR
  # this is a nonsense model on censoring to test function only
  weighted_OR <- suppressWarnings(glm(formula = Event~ARM,
                                      family = binomial(link="logit"),
                                      data = combined_data,
                                      weight = wt))
  
  boot_test_OR <- bootstrap_OR(intervention_data =  est_weights$analysis_data, 
                               i =c(1:500),
                               comparator_data = comparator_input, # comparator pseudo data
                               matching = est_weights$matching_vars, # matching variables
                               model = 'Event ~ ARM' )
  
  # check matches to normal weighted logistic regression on full data
  expect_equal(exp(as.numeric(coef(weighted_OR)[2])), boot_test_OR, tolerance = 0.01)
  
})

