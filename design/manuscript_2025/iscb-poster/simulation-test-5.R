# Updated simulation parameters
sim_params <- list(
  n_sizes = c(50, 200),
  target_marginal_hr = c(0.75, 0.85),  # Specify desired marginal HRs
  conditional_hr = c(0.70, 0.85),      # For reference
  distributions = c("weibull"),
  median_surv = c(24),
  surv_48m = c(0.10),
  dropout_rates = c(0.3),
  ess_ratios = c(0.25, 0.40),
  n_sim = 1000,
  max_followup = 48,
  n_boot = 200,
  max_weight = 100,
  n_mc_truth = 5000,
  n_cores = 4,
  results_prefix = "maic_fixed",
  timestamp_files = TRUE
)

# Run simulation
results <- main(sim_params)