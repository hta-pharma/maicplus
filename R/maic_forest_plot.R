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
#' @example inst/examples/maic_forest_plot_ex.R
#' @export


maic_forest_plot <- function(...,
                             xlim = c(0, 1.5),
                             reference_line = 1) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 package is required for this function")
    }
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("patchwork package is required for this function")
    }

    # 1) Gather all objects
    objs_list <- list(...)
    if (length(objs_list) == 0) {
      stop("No MAIC objects were provided. Pass at least one object with $inferential$summary.")
    }

    change_case_name <-
      function(data0_case_col, A_name, B_name, C_name) {
        case_renamed <- data0_case_col
        for (i in 1:length(data0_case_col)) {
          if (data0_case_col[i] == "AC") {
            case_renamed[i] <- paste0(A_name, " vs. ", C_name)
          } else if (data0_case_col[i] == "adjusted_AC") {
            case_renamed[i] <- paste0("Adjusted ", A_name, " vs. ", C_name)
          } else if (data0_case_col[i] == "BC") {
            case_renamed[i] <- paste0(B_name, " vs. ", C_name)
          } else if (data0_case_col[i] == "AB") {
            case_renamed[i] <- paste0(A_name, " vs. ", B_name)
          } else if (data0_case_col[i] == "adjusted_AB") {
            case_renamed[i] <- paste0("Adjusted ", A_name, " vs. ", B_name)
          }
        }
        return(case_renamed)
      }

    # 2) Extract and combine inferential summaries and descriptive summaries
    df_list <-
      lapply(objs_list, function(x) {
        # FIX: Opening brace on same line
        if (!("inferential" %in% names(x)) ||
            !("summary" %in% names(x$inferential))) {
          stop("One of the objects doesn't have 'inferential$summary'. Check your inputs.")
        }
        inferential_df <- x$inferential$summary
        if (!("descriptive" %in% names(x)) ||
            !("summary" %in% names(x$descriptive))) {
          stop("One of the objects doesn't have 'descriptive$summary'. Check your inputs.")
        }
        descriptive_df <- x$descriptive$summary
        inferential_fit_obj <- x$inferential$fit
        safely_extract_name <- function(df, trt_char) {
          if (trt_char %in% df$trt_ind) {
            return(df[df$trt_ind == trt_char, ]$treatment[1])
          } else {
            return(NA_character_)
          }
        }
        C_name <- safely_extract_name(descriptive_df, "C")
        A_name <- safely_extract_name(descriptive_df, "A")
        B_name <- safely_extract_name(descriptive_df, "B")
        effect_measure_col_name <- NULL
        if ("HR" %in% names(inferential_df)) {
          effect_measure_col_name <- "HR"
        } else if ("OR" %in% names(inferential_df)) {
          effect_measure_col_name <- "OR"
        } else if ("RR" %in% names(inferential_df)) {
          effect_measure_col_name <- "RR"
        }
        # consider the bootstrap result if exists
        if (!is.null(effect_measure_col_name) &&
            "boot_res_AB" %in% names(inferential_fit_obj)) {
          boot_results <- inferential_fit_obj$boot_res_AB
          adjusted_AB_row_index <-
            which(inferential_df$case == "adjusted_AB")

          # If the "adjusted_AB" row exists and bootstrap results are valid
          if (length(adjusted_AB_row_index) > 0 &&
              !is.null(boot_results$est) &&
              !is.null(boot_results$ci_l) && !is.null(boot_results$ci_u)) {
            # Update the values for the 'adjusted_AB' row
            inferential_df[[effect_measure_col_name]][adjusted_AB_row_index] <-
              boot_results$est
            inferential_df$LCL[adjusted_AB_row_index] <-
              boot_results$ci_l
            inferential_df$UCL[adjusted_AB_row_index] <-
              boot_results$ci_u
          }
        }
        inferential_df$case <-
          change_case_name(inferential_df$case, A_name, B_name, C_name)
        return(inferential_df)
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
    forest_data$effect_est <- as.numeric(forest_data[[effect_col]])
    forest_data$LCL <- as.numeric(forest_data$LCL)
    forest_data$UCL <- as.numeric(forest_data$UCL)
    forest_data$pval <- as.numeric(forest_data$pval)
    forest_data$row_index <- seq_len(nrow(forest_data))

    # 2c) Make group_id a factor in reversed order so row 1 is at the TOP
    forest_data$group_id <- factor(forest_data$row_index,
                                   levels = forest_data$row_index)
    # 3) Create the forest plot
    col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

    forest <- ggplot2::ggplot(data = forest_data,
                              ggplot2::aes(
                                x = group_id,
                                y = effect_est,
                                ymin = LCL,
                                ymax = UCL
                              )) +
      ggplot2::geom_pointrange(ggplot2::aes(color = case)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = LCL,
          ymax = UCL,
          color = case
        ),
        width = 0,
        linewidth = 1
      ) +
      ggplot2::geom_hline(
        yintercept = reference_line,
        colour = "red",
        linetype = "dashed",
        alpha = 0.5
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(limits = xlim) +
      ggplot2::scale_x_discrete(labels = rev(forest_data$case),
                                limits = rev(levels(forest_data$group_id))) +
      ggplot2::xlab("Experimental vs. Comparator Treatment") +
      ggplot2::ylab(paste0(effect_col, " (95% CI)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(colour = NA, fill = NA),
        panel.grid.major.y = ggplot2::element_line(colour = col_grid, linewidth = 0.5),
        panel.border = ggplot2::element_rect(fill = NA, color = "black"),
        legend.position = "none",
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 13
        )
      )

    # 4) Build a table showing [HR (LCL, UCL)] and p-value
    dat_table <- forest_data
    dat_table$pval_str <-
      ifelse(dat_table$pval < 0.001,
             "< 0.001",
             sprintf("%.3f", dat_table$pval))
    dat_table$effect_est_ci_str <- paste0(
      sprintf("%.2f", dat_table$effect_est),
      " [",
      sprintf("%.2f", dat_table$LCL),
      ", ",
      sprintf("%.2f", dat_table$UCL),
      "]"
    )
    dat_table <-
      dat_table[, c("group_id", "case", "effect_est_ci_str", "pval_str")]

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

    dat_table_long$stat <-
      factor(dat_table_long$stat,
             levels = c("effect_est_ci_str", "pval_str"))

    # 5) Table plot
    table_base <-
      ggplot2::ggplot(dat_table_long,
                      ggplot2::aes(x = stat, y = group_id, label = value)) +
      ggplot2::geom_text(size = 3) +
      ggplot2::scale_x_discrete(position = "top",
                                labels = c(paste0(effect_col, " (95% CI)"), "P value")) +
      ggplot2::scale_y_discrete(labels = forest_data$case,
                                limits = rev(levels(dat_table_long$group_id))) +
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
    final_plot <-
      patchwork::wrap_plots(forest, table_base) + patchwork::plot_layout(widths = c(10, 4))

    return(final_plot)
  }
