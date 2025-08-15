# Test with reduced parameters
test_params <- sim_params
test_params$n_sizes <- c(50, 200)
test_params$true_hr <- 0.7
test_params$n_sim <- 100
test_params$n_boot <- 100
test_params$n_cores <- 2
test_params$distributions <- "weibull"
test_params$ess_ratios <- c(0.25, 0.50)

# Run test
test_results <- main()
