#' Forest Plot for One or More MAIC Objects
#'
#' This function compiles effect estimates (and their confidence intervals) from one or more
#' "MAIC objects" – typically the output from [maic_anchored()] or [maic_unanchored()] functions –
#' and creates a forest plot alongside a summary table of those effect estimates.
#'
#' @param ... One or more MAIC objects. Each object must contain an `inferential$summary` data frame
#'   with the columns `"HR"`, `"OR"`, or `"RR"` (one of these must be present), along with `"LCL"`,
#'   `"UCL"`, `"pval"`, and a `case` identifier column.
#' @param xlim A numeric vector of length two, specifying the limits of the effect-size axis
#'   in the resulting forest plot. Defaults to `c(0, 1.5)`.
#' @param reference_line A numeric value specifying where to draw the "no-effect" reference line
#'   on the forest plot. Defaults to `1`.
#'
#' @details
#' This function extracts the effect estimates (e.g., HR, OR, or RR) and their confidence intervals
#' from each provided MAIC object. It then stacks all estimates into a single data frame for plotting.
#' A forest plot is generated using **ggplot2** with vertical error bars displaying the confidence intervals.
#' The `reference_line` is drawn as a dashed line to indicate the null value (usually 1, meaning no difference).
#'
#' Below the forest plot, a table is constructed showing the point estimate and 95% confidence interval
#' for each row, along with its p-value. If the p-value is less than 0.001, it is displayed as `"< 0.001"`,
#' otherwise it is displayed to three decimal places.
#'
#' @return A [patchwork][patchwork::patchwork] object that combines:
#' \itemize{
#'   \item A forest plot of the provided effect estimates and their 95\% confidence intervals.
#'   \item A corresponding table listing each estimate in numeric form, along with the p-value.
#' }
#' Printing or plotting this returned object will display both the forest plot and the summary table.
#'
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_errorbar geom_hline coord_flip scale_y_continuous scale_x_discrete xlab ylab theme_classic theme element_text element_blank element_line element_rect
#' @importFrom dplyr mutate row_number select
#' @importFrom patchwork wrap_plots plot_layout
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose maic_obj is a MAIC object containing:
#' #  maic_obj$inferential$summary
#'
#' # Generate a forest plot with the default settings:
#' maic_forest_plot(maic_obj)
#'
#' # Specify a different x-axis limit and reference line:
#' maic_forest_plot(maic_obj, xlim = c(0, 2), reference_line = 1)
#' }
maic_forest_plot <- function(..., xlim = c(0, 1.5), reference_line = 1) {
  # 1) Gather all objects
  objs_list <- list(...)
  if (length(objs_list) == 0) {
    stop("No MAIC objects were provided. Pass at least one object with $inferential$summary.")
  }

  # 2) Extract and combine inferential summaries
  df_list <- lapply(objs_list, function(x) {
    if (!("inferential" %in% names(x)) ||
      !("summary" %in% names(x$inferential))) {
      stop("One of the objects doesn't have 'inferential$summary'. Check your inputs.")
    }
    x$inferential$summary
  })
  forest_data <- do.call(rbind, df_list)
  rownames(forest_data) <- NULL

  if ("HR" %in% names(forest_data)) {
    effect_col <- "HR"
  } else if ("OR" %in% names(forest_data)) {
    effect_col <- "OR"
  } else if ("RR" %in% names(forest_data)) {
    effect_col <- "RR"
  } else {
    stop("No recognized effect measure (HR, OR, or RR) in the summary data.")
  }


  # Convert to numeric if needed
  forest_data <- forest_data %>%
    dplyr::mutate(
      effect_est = as.numeric(.data[[effect_col]]),
      LCL = as.numeric(LCL),
      UCL = as.numeric(UCL),
      pval = as.numeric(pval),
      row_index = dplyr::row_number() # 1,2,... in the order they appear
    )

  # 2c) Make group_id a factor in reversed order so row 1 is at the TOP
  forest_data$group_id <- factor(forest_data$row_index,
    levels = forest_data$row_index
  )
  # 3) Create the forest plot
  col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

  forest <- ggplot2::ggplot(
    data = forest_data,
    ggplot2::aes(x = group_id, y = effect_est, ymin = LCL, ymax = UCL)
  ) +
    ggplot2::geom_pointrange(ggplot2::aes(color = case)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = LCL, ymax = UCL, color = case), width = 0, size = 1) +
    ggplot2::geom_hline(yintercept = reference_line, colour = "red", linetype = "dashed", alpha = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = xlim) +
    ggplot2::scale_x_discrete(
      labels = rev(forest_data$case),
      limits = rev(levels(forest_data$group_id))
    ) +
    ggplot2::xlab("Experimental vs. Comparator Treatment") +
    ggplot2::ylab(paste0(effect_col, " (95% CI)")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(colour = NA, fill = NA),
      panel.grid.major.y = ggplot2::element_line(colour = col_grid, size = 0.5),
      panel.border = ggplot2::element_rect(fill = NA, color = "black"),
      legend.position = "none",
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13)
    )

  # 4) Build a table showing [HR (LCL, UCL)] and p-value
  dat_table <- forest_data %>%
    dplyr::mutate(
      pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
      # Build a string with HR and 95% CI
      effect_est_ci_str = paste0(
        sprintf("%.2f", effect_est),
        " [", sprintf("%.2f", LCL), ", ",
        sprintf("%.2f", UCL), "]"
      )
    ) %>%
    dplyr::select(group_id, case, effect_est_ci_str, pval_str)

  df_effect <- data.frame(
    group_id = dat_table$group_id,
    case = dat_table$case,
    stat = "effect_est_ci_str",
    value = dat_table$effect_est_ci_str,
    stringsAsFactors = FALSE
  )

  df_pval <- data.frame(
    group_id = dat_table$group_id,
    case = dat_table$case,
    stat = "pval_str",
    value = dat_table$pval_str,
    stringsAsFactors = FALSE
  )


  dat_table_long <- rbind(df_effect, df_pval)


  dat_table_long$stat <- factor(dat_table_long$stat, levels = c("effect_est_ci_str", "pval_str"))


  # 5) Table plot
  table_base <- ggplot2::ggplot(dat_table_long, ggplot2::aes(x = stat, y = group_id, label = value)) +
    ggplot2::geom_text(size = 3) +
    ggplot2::scale_x_discrete(
      position = "top",
      labels = c(paste0(effect_col, " (95% CI)"), "P value")
    ) +
    ggplot2::scale_y_discrete(
      labels = forest_data$case,
      limits = rev(levels(dat_table_long$group_id))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold")
    )

  # 6) Combine forest & table
  # final_plot <- forest + table_base + patchwork::plot_layout(widths = c(10, 4))
  final_plot <- patchwork::wrap_plots(forest, table_base) + patchwork::plot_layout(widths = c(10, 4))

  return(final_plot)
}
