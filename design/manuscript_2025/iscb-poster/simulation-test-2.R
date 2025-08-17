# Basic execution
source("maic_simulation.R")
results <- main()

# Custom configuration
sim_params$n_cores <- 8
sim_params$n_sim <- 1000
results <- main()

# Testing subset
test_params <- sim_params
test_params$n_sizes <- c(50, 200)
test_params$n_sim <- 100
test_results <- run_simulation_study(test_params)
