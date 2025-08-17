# Source the modified script
setwd("/cloud/project/design/manuscript_2025/iscb-poster")
source("simulation-4.R")

# =============================================================================
# EXAMPLE CODE TO RUN FULL SIMULATION
# =============================================================================

# 1. QUICK TEST RUN (small scale)
# --------------------------------
test_params <- list(
  n_sizes = c(50, 200),
  conditional_hr = c(0.70),
  distributions = c("weibull"),
  median_surv = 24,
  surv_48m = 0.15,
  dropout_rates = 0.2,
  event_rates = c(0.10, 0.25),
  ess_ratios = c(0.25, 0.50),
  n_sim = 100,  # Small number for testing
  max_followup = 48,
  n_boot = 100,  # Fewer bootstrap iterations
  n_mc_truth = 5000,  # Smaller MC sample for truth
  max_weight = 100,
  n_cores = 2,
  results_prefix = "test_run",
  timestamp_files = TRUE
)

# Run test
test_results <- main(test_params, run_analysis = TRUE)

# View results
print(test_results$analysis$method_comparison)
head(test_results$results)

# 2. FULL SIMULATION RUN
# ----------------------
full_params <- list(
  # Sample sizes per arm
  n_sizes = c(20, 50, 200, 500),

  # Conditional HR values (exp(gamma))
  conditional_hr = c(0.70, 0.85),

  # Distributions
  distributions = c("weibull", "lognormal"),

  # Survival parameters
  median_surv = c(24, 36),
  surv_48m = c(0.05, 0.15),

  # Censoring
  dropout_rates = c(0.2, 0.4),
  event_rates = c(0.05, 0.10, 0.25, 0.40),

  # Population overlap
  ess_ratios = c(0.15, 0.25, 0.33, 0.50),

  # Simulation settings
  n_sim = 10000,
  max_followup = 48,
  n_boot = 500,
  n_mc_truth = 10000,
  max_weight = 100,

  # Computing
  n_cores = parallel::detectCores() - 1,

  # Output
  results_prefix = "maic_iptw_atut_full",
  timestamp_files = TRUE
)

# Run full simulation
full_results <- main(full_params, run_analysis = TRUE)

# 3. CUSTOM ANALYSIS OF EXISTING RESULTS
# ---------------------------------------
# Load saved results
saved_results <- readRDS("simulation_results_20241201_143022/
                         maic_iptw_atut_full_n20-500_chr0.70-0.85_sim10000.rds")

# Run specific analyses
weight_analysis <- check_weight_distributions(saved_results)
print(weight_analysis)

# Check specific scenario
scenario_50_correct <- saved_results[
  saved_results$n == 50 &
  saved_results$weighting_scenario == "correct", ]

# Compare methods for this scenario
scenario_comparison <- compare_methods(scenario_50_correct)

# Check convergence
convergence_check <- summarize_convergence(saved_results, n_sim = 10000)
problematic <- convergence_check[convergence_check$flag_high_mc_error, ]

# 4. GENERATE PUBLICATION TABLES
# ------------------------------
# Create LaTeX tables
latex_tables <- create_summary_tables(saved_results, output_format = "latex")

# Create comparison plots (requires ggplot2)
library(ggplot2)

# Bias plot
bias_summary <- aggregate(
  cbind(bias = hr_est - true_marginal_hr) ~ method + n,
  data = saved_results,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                     se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
)

ggplot(bias_summary, aes(x = n, y = bias[,"mean"],
                         color = method, group = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = bias[,"mean"] - 1.96*bias[,"se"],
                    ymax = bias[,"mean"] + 1.96*bias[,"se"]),
                width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Bias by Sample Size and Method",
       x = "Sample Size per Arm",
       y = "Bias (HR scale)",
       color = "Method") +
  theme_minimal()

# 5. EXTRACT KEY RESULTS FOR MANUSCRIPT
# --------------------------------------
# Non-collapsibility assessment
nc_summary <- unique(saved_results[, c("conditional_hr", "true_marginal_hr")])
nc_summary$ratio <- nc_summary$true_marginal_hr / nc_summary$conditional_hr
write.csv(nc_summary, "non_collapsibility_summary.csv")

# Best performing methods
performance <- calculate_performance_metrics(saved_results)
best_methods <- performance[order(performance$abs_bias[,"mean"]), ]
head(best_methods, 10)

# Coverage by CI type
coverage_summary <- aggregate(
  cbind(robust = covered_robust,
        boot_perc = covered_boot_perc,
        boot_norm = covered_boot_norm,
        boot_bca = covered_boot_bca) ~ method,
  data = saved_results,
  FUN = function(x) mean(x, na.rm = TRUE)
)
print(coverage_summary)

# =============================================================================

# # Load and analyze existing results
# saved_results <- readRDS("path/to/results.rds")
# weight_analysis <- check_weight_distributions(saved_results)
# scenario_comparison <- compare_methods(scenario_50_correct)
# convergence_check <- summarize_convergence(saved_results, n_sim = 10000)
#
# # Create formatted tables for manuscript
# latex_tables <- create_summary_tables(saved_results, output_format = "latex")
# # Generate bias plots with ggplot2
# # Coverage summaries by CI type
#
# # Non-collapsibility assessment
# nc_summary <- unique(saved_results[, c("conditional_hr", "true_marginal_hr")])
# nc_summary$ratio <- nc_summary$true_marginal_hr / nc_summary$conditional_hr
#
# # Best performing methods
# performance <- calculate_performance_metrics(saved_results)
# best_methods <- performance[order(performance$abs_bias[,"mean"]), ]
