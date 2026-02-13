# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate performance metrics
calculate_performance_metrics <- function(results) {
  cat("Calculating performance metrics using true marginal ATUT HR...\n")
  
  # Calculate bias using true marginal HR
  results$bias <- results$hr_est - results$true_marginal_hr
  results$rel_bias <- (results$hr_est - results$true_marginal_hr) / results$true_marginal_hr
  results$abs_bias <- abs(results$bias)
  results$abs_rel_bias <- abs(results$rel_bias)
  
  # Coverage indicators using true marginal HR
  results$covered_robust <- (results$true_marginal_hr >= results$ci_robust_lower) &
    (results$true_marginal_hr <= results$ci_robust_upper)
  results$covered_boot_perc <- (results$true_marginal_hr >= results$ci_boot_perc_lower) &
    (results$true_marginal_hr <= results$ci_boot_perc_upper)
  results$covered_boot_norm <- (results$true_marginal_hr >= results$ci_boot_norm_lower) &
    (results$true_marginal_hr <= results$ci_boot_norm_upper)
  results$covered_boot_bca <- (results$true_marginal_hr >= results$ci_boot_bca_lower) &
    (results$true_marginal_hr <= results$ci_boot_bca_upper)
  
  # CI widths
  results$ci_width_robust <- results$ci_robust_upper - results$ci_robust_lower
  results$ci_width_boot_perc <- results$ci_boot_perc_upper - results$ci_boot_perc_lower
  results$ci_width_boot_norm <- results$ci_boot_norm_upper - results$ci_boot_norm_lower
  results$ci_width_boot_bca <- results$ci_boot_bca_upper - results$ci_boot_bca_lower
  
  return(results)
}

#' Create comprehensive summary tables
#'
#' @description Creates detailed summary tables with all required metrics
#' @param results Simulation results dataframe
#' @param output_format Either "data.frame" or "latex"
#' @return Summary table with performance metrics
create_summary_tables <- function(results, output_format = "data.frame") {
  cat("\nCreating comprehensive summary tables...\n")
  
  # Debug: Check input
  if (is.null(results) || nrow(results) == 0) {
    stop("Input results is NULL or has no rows")
  }
  
  cat(sprintf("  Input data: %d rows, %d columns\n", nrow(results), ncol(results)))
  cat(sprintf("  Columns: %s\n", paste(head(names(results), 20), collapse = ", ")))
  
  # First ensure performance metrics are calculated
  if (!"bias" %in% names(results)) {
    results <- calculate_performance_metrics(results)
  }
  
  # Add main_method column
  results$main_method <- NA
  results$main_method[grepl("^IPTW", results$method)] <- "IPTW"
  results$main_method[grepl("^MEW", results$method)] <- "MEW" 
  results$main_method[results$method %in% c("OW", "SW1", "SWN", "SWESS")] <- "MAIC"
  
  # Extract normalization method more carefully
  results$normalize_method <- NA
  
  # MAIC methods
  maic_idx <- results$method %in% c("OW", "SW1", "SWN", "SWESS")
  results$normalize_method[maic_idx] <- results$method[maic_idx]
  
  # MEW methods
  mew_idx <- results$method %in% c("MEW", "MEW1", "MEWN", "MEWESS")
  results$normalize_method[results$method == "MEW"] <- "OW"
  results$normalize_method[results$method == "MEW1"] <- "SW1"
  results$normalize_method[results$method == "MEWN"] <- "SWN"
  results$normalize_method[results$method == "MEWESS"] <- "SWESS"
  
  # IPTW methods
  iptw_idx <- grepl("^IPTW_", results$method)
  results$normalize_method[iptw_idx] <- substring(results$method[iptw_idx], 6)
  
  # Use d instead of ess_ratio since it's available
  # Round d to 2 decimal places for grouping
  if ("d" %in% names(results)) {
    results$d_group <- round(results$d, 2)
  }
  
  # Check which columns are available for grouping
  available_grouping <- c("method", "main_method", "normalize_method", "n", "weighting_scenario")
  if ("conditional_hr" %in% names(results)) available_grouping <- c(available_grouping, "conditional_hr")
  if ("true_marginal_hr" %in% names(results)) available_grouping <- c(available_grouping, "true_marginal_hr")
  if ("distribution" %in% names(results)) available_grouping <- c(available_grouping, "distribution")
  if ("d_group" %in% names(results)) {
    available_grouping <- c(available_grouping, "d_group")
  } else if ("d" %in% names(results)) {
    available_grouping <- c(available_grouping, "d")
  }
  
  cat(sprintf("  Grouping by: %s\n", paste(available_grouping, collapse = ", ")))
  
  # Check which metrics are available
  available_metrics <- c()
  if ("bias" %in% names(results)) available_metrics <- c(available_metrics, "bias")
  if ("rel_bias" %in% names(results)) available_metrics <- c(available_metrics, "rel_bias") 
  if ("abs_bias" %in% names(results)) available_metrics <- c(available_metrics, "abs_bias")
  if ("covered_robust" %in% names(results)) available_metrics <- c(available_metrics, "covered_robust")
  if ("covered_boot_perc" %in% names(results)) available_metrics <- c(available_metrics, "covered_boot_perc")
  if ("covered_boot_norm" %in% names(results)) available_metrics <- c(available_metrics, "covered_boot_norm")
  if ("ess" %in% names(results)) available_metrics <- c(available_metrics, "ess")
  if ("max_weight" %in% names(results)) available_metrics <- c(available_metrics, "max_weight")
  
  cat(sprintf("  Available metrics: %s\n", paste(available_metrics, collapse = ", ")))
  
  if (length(available_metrics) == 0) {
    stop("No metrics available for aggregation")
  }
  
  # Build formula dynamically
  formula_str <- sprintf("cbind(%s) ~ %s", 
                         paste(available_metrics, collapse = ", "),
                         paste(available_grouping, collapse = " + "))
  
  cat(sprintf("  Formula: %s\n", formula_str))
  
  # Remove rows with NA in grouping variables
  n_original <- nrow(results)
  for (var in available_grouping) {
    if (var %in% names(results)) {
      n_before <- nrow(results)
      results <- results[!is.na(results[[var]]), ]
      n_after <- nrow(results)
      if (n_before > n_after) {
        cat(sprintf("  Removed %d rows with NA in %s\n", n_before - n_after, var))
      }
    }
  }
  
  # Check if we still have data
  if (nrow(results) == 0) {
    stop(sprintf("No data left after removing NA values in grouping variables (started with %d rows)", n_original))
  }
  
  cat(sprintf("  Data for aggregation: %d rows\n", nrow(results)))
  
  # Aggregate metrics
  summary_table <- tryCatch({
    aggregate(
      as.formula(formula_str),
      data = results,
      FUN = function(x) {
        if (length(x) == 0 || all(is.na(x))) {
          c(mean = NA, sd = NA, median = NA, q25 = NA, q75 = NA, n_valid = 0)
        } else {
          c(mean = mean(x, na.rm = TRUE),
            sd = sd(x, na.rm = TRUE),
            median = median(x, na.rm = TRUE),
            q25 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
            q75 = as.numeric(quantile(x, 0.75, na.rm = TRUE)),
            n_valid = sum(!is.na(x)))
        }
      }
    )
  }, error = function(e) {
    cat(sprintf("  Error in aggregation: %s\n", e$message))
    cat("  Trying simpler aggregation...\n")
    
    # Try a simpler aggregation with fewer grouping variables
    simple_formula <- sprintf("cbind(%s) ~ method + n + weighting_scenario",
                              paste(available_metrics[1:min(3, length(available_metrics))], collapse = ", "))
    
    aggregate(
      as.formula(simple_formula),
      data = results,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
  })
  
  if (is.null(summary_table) || nrow(summary_table) == 0) {
    stop("Aggregation failed - no summary table created")
  }
  
  cat(sprintf("  Summary table created: %d rows\n", nrow(summary_table)))
  
  # Create a cleaner data frame
  n_grouping <- length(available_grouping)
  clean_summary <- summary_table[, 1:n_grouping, drop = FALSE]
  
  # Add statistics for each metric
  for (metric in available_metrics) {
    if (metric %in% names(summary_table)) {
      metric_data <- summary_table[[metric]]
      
      # Handle both matrix and vector outputs from aggregate
      if (is.matrix(metric_data) || is.data.frame(metric_data)) {
        if (metric %in% c("bias", "rel_bias", "abs_bias")) {
          clean_summary[[paste0(metric, "_mean")]] <- metric_data[, "mean"]
          clean_summary[[paste0(metric, "_sd")]] <- metric_data[, "sd"]
        } else if (metric %in% c("covered_robust", "covered_boot_perc", "covered_boot_norm")) {
          clean_summary[[paste0(metric, "_rate")]] <- metric_data[, "mean"]
        } else if (metric == "ess") {
          clean_summary[["ess_median"]] <- metric_data[, "median"]
          clean_summary[["ess_q25"]] <- metric_data[, "q25"]
          clean_summary[["ess_q75"]] <- metric_data[, "q75"]
        } else if (metric == "max_weight") {
          clean_summary[["max_weight_median"]] <- metric_data[, "median"]
        }
      } else {
        # Handle simple vector output
        if (metric %in% c("bias", "rel_bias", "abs_bias")) {
          clean_summary[[paste0(metric, "_mean")]] <- metric_data
        } else if (metric %in% c("covered_robust", "covered_boot_perc", "covered_boot_norm")) {
          clean_summary[[paste0(metric, "_rate")]] <- metric_data
        } else if (metric == "ess") {
          clean_summary[["ess_median"]] <- metric_data
        } else if (metric == "max_weight") {
          clean_summary[["max_weight_median"]] <- metric_data
        }
      }
    }
  }
  
  # Sort for better readability
  sort_cols <- intersect(c("n", "conditional_hr", "distribution", "d", "d_group", 
                           "weighting_scenario", "main_method", "normalize_method"),
                         names(clean_summary))
  
  if (length(sort_cols) > 0) {
    clean_summary <- clean_summary[do.call(order, clean_summary[sort_cols]), ]
  }
  
  # Print summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("Total unique scenarios: %d\n", nrow(clean_summary)))
  
  if ("main_method" %in% names(clean_summary)) {
    methods <- unique(clean_summary$main_method)
    cat(sprintf("Methods compared: %s\n", paste(methods[!is.na(methods)], collapse = ", ")))
  }
  
  if ("n" %in% names(clean_summary)) {
    cat(sprintf("Sample sizes: %s\n", paste(sort(unique(clean_summary$n)), collapse = ", ")))
  }
  
  if ("d_group" %in% names(clean_summary)) {
    cat(sprintf("d values: %s\n", paste(sort(unique(clean_summary$d_group)), collapse = ", ")))
  } else if ("d" %in% names(clean_summary)) {
    cat(sprintf("d values: %s\n", paste(sort(unique(round(clean_summary$d, 2))), collapse = ", ")))
  }
  
  # Map d values to approximate ESS ratios for display
  if ("d" %in% names(clean_summary) || "d_group" %in% names(clean_summary)) {
    d_col <- if ("d_group" %in% names(clean_summary)) "d_group" else "d"
    unique_d <- sort(unique(clean_summary[[d_col]]))
    cat("\nApproximate ESS ratios by d:\n")
    for (d_val in unique_d) {
      if (abs(d_val - 1.00) < 0.05) ess_approx <- "~0.50"
      else if (abs(d_val - 1.15) < 0.05) ess_approx <- "~0.33"
      else if (abs(d_val - 1.25) < 0.05) ess_approx <- "~0.25"
      else if (abs(d_val - 1.40) < 0.05) ess_approx <- "~0.15"
      else ess_approx <- "unknown"
      cat(sprintf("  d = %.2f: ESS ratio %s\n", d_val, ess_approx))
    }
  }
  
  # Output format
  if (output_format == "latex") {
    cat("\n\\begin{table}[h]\n")
    cat("\\caption{Performance Metrics by Method and Scenario}\n")
    cat("\\begin{tabular}{llrrrrrrr}\n")
    cat("\\toprule\n")
    cat("Method & n & d & Bias (SD) & Coverage & ESS (Med) & Max Weight \\\\\n")
    cat("\\midrule\n")
    
    # Select subset for LaTeX (too many rows otherwise)
    subset_idx <- rep(TRUE, nrow(clean_summary))
    if ("weighting_scenario" %in% names(clean_summary)) {
      subset_idx <- subset_idx & clean_summary$weighting_scenario == "correct"
    }
    if ("conditional_hr" %in% names(clean_summary)) {
      subset_idx <- subset_idx & clean_summary$conditional_hr == clean_summary$conditional_hr[1]
    }
    
    latex_subset <- clean_summary[subset_idx, ]
    
    for (i in 1:min(20, nrow(latex_subset))) {
      row <- latex_subset[i, ]
      d_val <- if ("d_group" %in% names(row)) row$d_group else if ("d" %in% names(row)) round(row$d, 2) else NA
      
      cat(sprintf("%s & %d & %.2f & %.3f (%.3f) & %.3f & %.1f & %.1f \\\\\n",
                  row$method, row$n, d_val,
                  if ("bias_mean" %in% names(row)) row$bias_mean else NA,
                  if ("bias_sd" %in% names(row)) row$bias_sd else NA,
                  if ("covered_robust_rate" %in% names(row)) row$covered_robust_rate else NA,
                  if ("ess_median" %in% names(row)) row$ess_median else NA,
                  if ("max_weight_median" %in% names(row)) row$max_weight_median else NA))
    }
    
    cat("\\bottomrule\n")
    cat("\\end{tabular}\n")
    cat("\\end{table}\n")
  }
  
  return(clean_summary)
}


#' Create comparison tables for MAIC vs IPTW
create_comparison_tables <- function(results) {
  cat("\nCreating MAIC/MEW vs IPTW ATUT comparison tables...\n")
  
  # Ensure metrics are calculated
  if (!"bias" %in% names(results)) {
    results <- calculate_performance_metrics(results)
  }
  
  # Add method type
  results$method_type <- ifelse(grepl("^IPTW", results$method), "IPTW_ATUT",
                                ifelse(grepl("^MEW", results$method), "MEW", "MAIC"))
  
  results$mse <- (results$hr_est - results$true_marginal_hr)^2
  
  # Overall comparison
  overall_comparison <- aggregate(
    cbind(bias, rel_bias, mse, ess) ~ method_type + n + weighting_scenario,
    data = results,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE))
    }
  )
  
  # Best method by scenario
  best_by_scenario <- aggregate(
    mse ~ method + n + conditional_hr + weighting_scenario,
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  best_method <- aggregate(
    mse ~ n + conditional_hr + weighting_scenario,
    data = best_by_scenario,
    FUN = function(x) {
      idx <- which.min(x)
      return(idx)
    }
  )
  
  # Relative efficiency (IPTW as reference)
  iptw_reference <- results[results$method == "IPTW_SW1", ]
  
  efficiency_comparison <- list()
  for (m in unique(results$method)) {
    if (m != "IPTW_SW1") {
      method_data <- results[results$method == m, ]
      merged <- merge(method_data, iptw_reference,
                      by = c("n", "conditional_hr", "distribution",
                             "weighting_scenario", "iteration"),
                      suffixes = c("", "_iptw"))
      
      merged$relative_mse <- merged$mse / merged$mse_iptw
      merged$relative_variance <- (merged$hr_est - mean(merged$hr_est, na.rm = TRUE))^2 /
        (merged$hr_est_iptw - mean(merged$hr_est_iptw, na.rm = TRUE))^2
      
      efficiency_comparison[[m]] <- data.frame(
        method = m,
        mean_relative_mse = mean(merged$relative_mse, na.rm = TRUE),
        median_relative_mse = median(merged$relative_mse, na.rm = TRUE),
        mean_relative_var = mean(merged$relative_variance, na.rm = TRUE)
      )
    }
  }
  
  efficiency_df <- do.call(rbind, efficiency_comparison)
  
  return(list(
    overall = overall_comparison,
    best_method = best_method,
    efficiency = efficiency_df
  ))
}



# After running the simulation
results <- saved_results

# Calculate performance metrics
results_with_metrics <- calculate_performance_metrics(results)

# Create comprehensive summary
summary_table <- create_summary_tables(results_with_metrics, output_format = "data.frame")

# View a subset of the summary
head(summary_table[, c("method", "n", "ess_ratio", "bias_mean", "bias_sd", 
                       "covered_robust_rate", "ess_median", "max_weight_median")])

# Create comparison tables
comparison <- create_comparison_tables(results_with_metrics)

# View overall comparison
print(comparison$overall)

# View efficiency comparison vs IPTW_SW1
print(comparison$efficiency)



# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

#' Create bias bar plot with facets
#'
#' @description Creates a bar plot showing mean bias with +/- 1 SD error bars
#' @param summary_data Output from create_summary_tables()
#' @param n_filter Sample size to filter for
#' @param d_filter d_group value to filter for
#' @param title Optional plot title
#' @return ggplot object
plot_bias_bars <- function(summary_data, n_filter, d_filter, title = NULL) {
  
  # Load required library
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  library(ggplot2)
  
  # Filter data
  plot_data <- summary_data[
    summary_data$n == n_filter & 
      abs(summary_data$d_group - d_filter) < 0.01,  # Use tolerance for numeric comparison
  ]
  
  if (nrow(plot_data) == 0) {
    stop(sprintf("No data found for n=%d and d_group=%.2f", n_filter, d_filter))
  }
  
  # Ensure required columns exist
  required_cols <- c("main_method", "normalize_method", "bias_mean", "bias_sd", 
                     "true_marginal_hr", "weighting_scenario")
  missing_cols <- setdiff(required_cols, names(plot_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Create factors for proper ordering
  plot_data$main_method <- factor(plot_data$main_method, 
                                  levels = c("MAIC", "MEW", "IPTW"))
  plot_data$normalize_method <- factor(plot_data$normalize_method,
                                       levels = c("OW", "SW1", "SWN", "SWESS"))
  
  # Round true_marginal_hr for cleaner facet labels
  plot_data$true_marginal_hr_label <- sprintf("True HR = %.2f", plot_data$true_marginal_hr)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = main_method, y = bias_mean, fill = normalize_method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = bias_mean - bias_sd, ymax = bias_mean + bias_sd),
                  position = position_dodge(width = 0.8), width = 0.25, 
                  color = "black", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_grid(true_marginal_hr_label ~ weighting_scenario, 
               scales = "free_y", 
               labeller = labeller(weighting_scenario = c(
                 "correct" = "Correct Specification",
                 "under" = "Under Specification", 
                 "over" = "Over Specification"
               ))) +
    scale_fill_manual(values = c("OW" = "#E41A1C", "SW1" = "#377EB8", 
                                 "SWN" = "#4DAF4A", "SWESS" = "#984EA3"),
                      name = "Normalization\nMethod") +
    labs(x = "Main Method",
         y = "Mean Bias ± 1 SD",
         title = if (!is.null(title)) title else 
           sprintf("Bias Comparison (n=%d, d=%.2f)", n_filter, d_filter)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

#' Create coverage line plot with facets
#'
#' @description Creates a line plot showing all three CI coverage rates
#' @param summary_data Output from create_summary_tables()
#' @param n_filter Sample size to filter for
#' @param d_filter d_group value to filter for
#' @param title Optional plot title
#' @return ggplot object
plot_coverage_lines <- function(summary_data, n_filter, d_filter, title = NULL) {
  
  # Load required library
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  library(ggplot2)
  
  # Filter data
  plot_data <- summary_data[
    summary_data$n == n_filter & 
      abs(summary_data$d_group - d_filter) < 0.01,  # Use tolerance for numeric comparison
  ]
  
  if (nrow(plot_data) == 0) {
    stop(sprintf("No data found for n=%d and d_group=%.2f", n_filter, d_filter))
  }
  
  # Check which coverage columns exist
  coverage_types <- c()
  coverage_labels <- c()
  
  if ("covered_robust_rate" %in% names(plot_data)) {
    coverage_types <- c(coverage_types, "covered_robust_rate")
    coverage_labels <- c(coverage_labels, "Robust")
  }
  if ("covered_boot_perc_rate" %in% names(plot_data)) {
    coverage_types <- c(coverage_types, "covered_boot_perc_rate")
    coverage_labels <- c(coverage_labels, "Bootstrap Percentile")
  }
  if ("covered_boot_norm_rate" %in% names(plot_data)) {
    coverage_types <- c(coverage_types, "covered_boot_norm_rate")
    coverage_labels <- c(coverage_labels, "Bootstrap Normal")
  }
  
  if (length(coverage_types) == 0) {
    stop("No coverage columns found in data")
  }
  
  # Reshape data to long format for multiple coverage types
  plot_data_long <- NULL
  for (i in seq_along(coverage_types)) {
    temp_data <- plot_data
    temp_data$coverage <- temp_data[[coverage_types[i]]]
    temp_data$coverage_type <- coverage_labels[i]
    
    # Select only needed columns
    keep_cols <- c("main_method", "normalize_method", "coverage", "coverage_type",
                   "true_marginal_hr", "weighting_scenario")
    temp_data <- temp_data[, keep_cols]
    
    if (is.null(plot_data_long)) {
      plot_data_long <- temp_data
    } else {
      plot_data_long <- rbind(plot_data_long, temp_data)
    }
  }
  
  # Create factors for proper ordering
  plot_data_long$main_method <- factor(plot_data_long$main_method, 
                                       levels = c("MAIC", "MEW", "IPTW"))
  plot_data_long$normalize_method <- factor(plot_data_long$normalize_method,
                                            levels = c("OW", "SW1", "SWN", "SWESS"))
  plot_data_long$coverage_type <- factor(plot_data_long$coverage_type,
                                         levels = coverage_labels)
  
  # Create numeric x-axis positions for each normalize method within main method
  # Main methods at positions 1, 2, 3
  # Normalize methods offset within each main method
  plot_data_long$x_position <- as.numeric(plot_data_long$main_method) + 
    (as.numeric(plot_data_long$normalize_method) - 2.5) * 0.15
  
  # Round true_marginal_hr for cleaner facet labels
  plot_data_long$true_marginal_hr_label <- sprintf("True HR = %.2f", plot_data_long$true_marginal_hr)
  
  # Create the plot
  p <- ggplot(plot_data_long, aes(x = x_position, y = coverage)) +
    # Plot lines for each coverage type
    geom_line(aes(color = coverage_type, group = interaction(coverage_type, main_method)),
              size = 1.2, alpha = 0.9) +
    # Add smaller points
    geom_point(aes(shape = normalize_method), size = 1.5, alpha = 0.7) +
    # Reference line at nominal coverage
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", alpha = 0.5) +
    # Facets
    facet_grid(true_marginal_hr_label ~ weighting_scenario,
               labeller = labeller(weighting_scenario = c(
                 "correct" = "Correct Specification",
                 "under" = "Under Specification", 
                 "over" = "Over Specification"
               ))) +
    # X-axis settings
    scale_x_continuous(breaks = 1:3, 
                       labels = levels(plot_data_long$main_method),
                       limits = c(0.5, 3.5)) +
    # Color scale for coverage types
    scale_color_manual(values = c("Robust" = "#E41A1C", 
                                  "Bootstrap Percentile" = "#377EB8", 
                                  "Bootstrap Normal" = "#4DAF4A"),
                       name = "CI Type") +
    # Shape scale for normalization methods
    scale_shape_manual(values = c("OW" = 16, "SW1" = 17, "SWN" = 15, "SWESS" = 18),
                       name = "Normalization\nMethod") +
    # Y-axis settings
    scale_y_continuous(limits = c(min(0.8, min(plot_data_long$coverage, na.rm = TRUE) - 0.02), 
                                  max(1.0, max(plot_data_long$coverage, na.rm = TRUE) + 0.02))) +
    # Labels
    labs(x = "Main Method",
         y = "Coverage Rate",
         title = if (!is.null(title)) title else 
           sprintf("Coverage Comparison (n=%d, d=%.2f)", n_filter, d_filter)) +
    # Theme
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.box = "vertical"
    ) +
    # Add annotation for nominal coverage
    annotate("text", x = 0.55, y = 0.95, label = "0.95", 
             hjust = 0, vjust = -0.3, color = "black", size = 3, alpha = 0.5)
  
  return(p)
}
#' Example usage function
plot_simulation_results <- function(summary_data, n_filter, d_filter, save_plots = FALSE) {
  # Create both plots
  p_bias <- plot_bias_bars(summary_data, n_filter, d_filter)
  p_coverage <- plot_coverage_lines(summary_data, n_filter, d_filter, coverage_type = "robust")
  
  # Display plots
  print(p_bias)
  print(p_coverage)
  
  # Optionally save plots
  if (save_plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for saving plots")
    }
    
    filename_base <- sprintf("simulation_plots_n%d_d%.2f", n_filter, d_filter)
    
    ggplot2::ggsave(paste0(filename_base, "_bias.png"), p_bias, 
                    width = 12, height = 8, dpi = 300)
    ggplot2::ggsave(paste0(filename_base, "_coverage.png"), p_coverage, 
                    width = 12, height = 8, dpi = 300)
    
    cat(sprintf("Plots saved as %s_bias.png and %s_coverage.png\n", 
                filename_base, filename_base))
  }
  
  return(list(bias_plot = p_bias, coverage_plot = p_coverage))
}




# # After creating summary table
# summary_table <- create_summary_tables(results_with_metrics)

# Create bias plot for n=200, d=1.25
library(ggplot2)
p_bias20 <- plot_bias_bars(summary_table, n_filter = 80, d_filter = 1.58, title = "Bias Comparison (n=80, True ESS Ratio=20%)")
print(p_bias20)

p_bias50 <- plot_bias_bars(summary_table, n_filter = 80, d_filter = 1.76, title = "Bias Comparison (n=80, True ESS Ratio=50%)")
print(p_bias50)

p_bias20_350 <- plot_bias_bars(summary_table, n_filter = 350, d_filter = 1.58, title = "Bias Comparison (n=350, True ESS Ratio=20%)")
print(p_bias20_350)

p_bias50_350 <- plot_bias_bars(summary_table, n_filter = 350, d_filter = 1.76, title = "Bias Comparison (n=350, True ESS Ratio=50%)")
print(p_bias50_350)

# Create coverage plot for the same parameters
p_coverage20 <- plot_coverage_lines(summary_table, n_filter = 80, d_filter = 1.58, title = "Coverage Comparison (n=80, True ESS Ratio=50%)")
print(p_coverage20)


p_coverage50 <- plot_coverage_lines(summary_table, n_filter = 80, d_filter = 1.76, title = "Coverage Comparison (n=80, True ESS Ratio=20%)")
print(p_coverage50)

p_coverage20_350 <- plot_coverage_lines(summary_table, n_filter = 350, d_filter = 1.58, title = "Coverage Comparison (n=350, True ESS Ratio=50%)")
print(p_coverage20_350)

p_coverage50_350 <- plot_coverage_lines(summary_table, n_filter = 350, d_filter = 1.76, title = "Coverage Comparison (n=350, True ESS Ratio=20%)")
print(p_coverage50_350)


# # Or use the convenience function to create both plots
# plots <- plot_simulation_results(summary_table, n_filter = 200, d_filter = 1.25, 
#                                  save_plots = TRUE)
# 
# # You can also plot different coverage types
# p_coverage_boot <- plot_coverage_lines(summary_table, n_filter = 200, d_filter = 1.25,
#                                        coverage_type = "boot_perc")
# 
# # Custom titles
# p_bias_custom <- plot_bias_bars(summary_table, n_filter = 50, d_filter = 1.40,
#                                 title = "Small Sample Performance (n=50, High Covariate Shift)")
# 



#' Create bias line plot with facets
#'
#' @description Creates a line plot showing mean bias across methods
#' with different lines for combinations of sample size (n) and d_group
#' @param summary_data Output from create_summary_tables()
#' @param n_values Vector of sample sizes to include (default: all available)
#' @param d_values Vector of d_group values to include (default: all available)
#' @param d_labels Named vector for d_group legend labels (e.g., c("1.00" = "ESS~0.50", "1.15" = "ESS~0.33"))
#' @param title Optional plot title
#' @return ggplot object
plot_bias_lines <- function(summary_data, n_values = NULL, d_values = NULL, 
                            d_labels = NULL, title = NULL) {
  
  # Load required library
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  library(ggplot2)
  
  # Get available values if not specified
  if (is.null(n_values)) {
    n_values <- sort(unique(summary_data$n))
  }
  if (is.null(d_values)) {
    d_values <- sort(unique(summary_data$d_group))
  }
  
  # Filter data for specified n and d_group values
  plot_data <- summary_data[
    summary_data$n %in% n_values & 
      summary_data$d_group %in% d_values,
  ]
  
  if (nrow(plot_data) == 0) {
    stop("No data found for specified n_values and d_values")
  }
  
  # Ensure required columns exist
  required_cols <- c("main_method", "normalize_method", "bias_mean", 
                     "true_marginal_hr", "weighting_scenario", "n", "d_group")
  missing_cols <- setdiff(required_cols, names(plot_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Create factors for proper ordering
  plot_data$main_method <- factor(plot_data$main_method, 
                                  levels = c("MAIC", "MEW", "IPTW"))
  plot_data$normalize_method <- factor(plot_data$normalize_method,
                                       levels = c("OW", "SW1", "SWN", "SWESS"))
  
  # Create combination variable for n and d_group
  plot_data$n_d_combo <- paste0("n=", plot_data$n, ", d=", sprintf("%.2f", plot_data$d_group))
  
  # Apply custom d_group labels if provided
  if (!is.null(d_labels)) {
    # Create a mapping for the combination labels
    d_mapping <- d_labels[as.character(sprintf("%.2f", plot_data$d_group))]
    # Replace d values in the combo labels
    for (i in 1:nrow(plot_data)) {
      d_key <- sprintf("%.2f", plot_data$d_group[i])
      if (!is.na(d_mapping[d_key])) {
        plot_data$n_d_combo[i] <- paste0("n=", plot_data$n[i], ", ", d_mapping[d_key])
      }
    }
  }
  
  # Create numeric x-axis positions for each normalize method within main method
  # Main methods at positions 1, 2, 3
  # Normalize methods offset within each main method
  plot_data$x_position <- as.numeric(plot_data$main_method) + 
    (as.numeric(plot_data$normalize_method) - 2.5) * 0.15
  
  # Round true_marginal_hr for cleaner facet labels
  plot_data$true_marginal_hr_label <- sprintf("True HR = %.2f", plot_data$true_marginal_hr)
  
  # Create color palette for n_d_combo
  n_combos <- length(unique(plot_data$n_d_combo))
  
  # Use a combination of colors and line types for better distinction
  if (n_combos <= 8) {
    color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[1:n_combos]
  } else {
    # Generate more colors if needed
    color_palette <- rainbow(n_combos)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x_position, y = bias_mean)) +
    # Plot lines for each n_d_combo grouped by main_method
    geom_line(aes(color = n_d_combo, 
                  group = interaction(n_d_combo, main_method)),
              size = 1.2, alpha = 0.9) +
    # Add points for better visibility
    geom_point(aes(color = n_d_combo, shape = normalize_method), 
               size = 2, alpha = 0.8) +
    # Reference line at zero bias
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    # Facets
    facet_grid(true_marginal_hr_label ~ weighting_scenario,
               scales = "free_y",
               labeller = labeller(weighting_scenario = c(
                 "correct" = "Correct Specification",
                 "under" = "Under Specification", 
                 "over" = "Over Specification"
               ))) +
    # X-axis settings
    scale_x_continuous(breaks = 1:3, 
                       labels = levels(plot_data$main_method),
                       limits = c(0.5, 3.5)) +
    # Color scale for n_d_combo
    scale_color_manual(values = color_palette,
                       name = "Sample Size &\nPopulation Difference") +
    # Shape scale for normalization methods
    scale_shape_manual(values = c("OW" = 16, "SW1" = 17, "SWN" = 15, "SWESS" = 18),
                       name = "Normalization\nMethod") +
    # Labels
    labs(x = "Main Method",
         y = "Mean Bias",
         title = if (!is.null(title)) title else 
           "Bias Comparison Across Methods and Scenarios") +
    # Theme
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.box = "vertical"
    ) +
    # Add annotation for zero bias reference
    annotate("text", x = 0.55, y = 0, label = "0", 
             hjust = 0, vjust = -0.3, color = "black", size = 3, alpha = 0.7)
  
  return(p)
}

#' Example usage function with custom d_group labels
plot_simulation_bias_results <- function(summary_data, n_values = NULL, d_values = NULL, 
                                         save_plots = FALSE) {
  
  # Example custom labels for d_group values (mapping d to approximate ESS ratios)
  d_labels <- c(
    "1.00" = "ESS~0.50",
    "1.15" = "ESS~0.33", 
    "1.25" = "ESS~0.25",
    "1.40" = "ESS~0.15"
  )
  
  # Create the bias line plot
  p_bias <- plot_bias_lines(summary_data, 
                            n_values = n_values,
                            d_values = d_values,
                            d_labels = d_labels)
  
  # Display plot
  print(p_bias)
  
  # Optionally save plot
  if (save_plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for saving plots")
    }
    
    # Create filename based on included values
    n_str <- if (is.null(n_values)) "all" else paste(n_values, collapse = "-")
    d_str <- if (is.null(d_values)) "all" else paste(sprintf("%.2f", d_values), collapse = "-")
    filename <- sprintf("bias_lines_n%s_d%s.png", n_str, d_str)
    
    ggplot2::ggsave(filename, p_bias, 
                    width = 14, height = 10, dpi = 300)
    
    cat(sprintf("Plot saved as %s\n", filename))
  }
  
  return(p_bias)
}

# Alternative function for specific subset plotting
plot_bias_subset <- function(summary_data, specific_n = c(50, 200), 
                             specific_d = c(1.00, 1.25), title = NULL) {
  
  # Custom labels for the specific d values
  d_labels <- c(
    "1.00" = "ESS~0.50",
    "1.25" = "ESS~0.25"
  )
  
  plot_bias_lines(summary_data, 
                  n_values = specific_n,
                  d_values = specific_d,
                  d_labels = d_labels,
                  title = title)
}



# Basic usage with all available data
plot_bias_lines(summary_table)

# With custom d_group labels
d_labels <- c("1.58" = "True-ESS-Ratio=0.5", "1.76" = "True-ESS-Ratio=0.2")
plot_bias_lines(summary_table, d_labels = d_labels)
