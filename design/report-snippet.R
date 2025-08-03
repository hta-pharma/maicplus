# : make print-ready analysis report table
res$inferential[["report_overall_robustCI"]] <- rbind(
  report_table_tte(coxobj_ipd, medSurv_ipd, tag = paste0("IPD/", endpoint_name)),
  report_table_tte(coxobj_ipd_adj, medSurv_ipd_adj, tag = paste0("weighted IPD/", endpoint_name)),
  report_table_tte(coxobj_agd, medSurv_agd, tag = paste0("Agd/", endpoint_name)),
  c(
    paste0("Unadj.", trt_ipd, " vs ", trt_agd, "(by Bucher)"),
    rep("--", 4),
    reformat(res_AB_unadj, pval_digits = 3)
  ),
  c(
    paste0("** adj.", trt_ipd, " vs ", trt_agd),
    rep("--", 4),
    reformat(res_AB, pval_digits = 3)
  )
)

if (is.null(res$report[["boot_res"]])) {
  res$inferential[["report_overall_bootCI"]] <- NULL
} else {
  temp_boot_res <- boot_res_AB
  temp_boot_res$ci_l <- boot_res_AB$ci_l
  temp_boot_res$ci_u <- boot_res_AB$ci_u
  class(temp_boot_res) <- class(res_AB)

  res$inferential[["report_overall_bootCI"]] <- rbind(
    report_table_tte(coxobj_ipd, medSurv_ipd, tag = paste0("IPD/", endpoint_name)),
    report_table_tte(coxobj_ipd_adj, medSurv_ipd_adj, tag = paste0("weighted IPD/", endpoint_name)),
    report_table_tte(coxobj_agd, medSurv_agd, tag = paste0("AgD/", endpoint_name)),
    c(
      paste0("** adj.", trt_ipd, " vs ", trt_agd),
      rep("--", 4),
      reformat(boot_res_AB, pval_digits = 3)
    )
  )
}



# : make print-ready analysis report table
tags <- paste0(c("IPD/", "weighted IPD/", "AgD/"), endpoint_name)
res$inferential[["report_overall_robustCI"]] <- rbind(
  report_table_binary(binobj_ipd, tag = tags[1], eff_measure = eff_measure),
  report_table_binary(binobj_ipd_adj, res_AC, tag = tags[2], eff_measure = eff_measure),
  report_table_binary(binobj_agd, tag = tags[3], eff_measure = eff_measure),
  c(
    paste0("** adj.", trt_ipd, " vs ", trt_agd),
    rep("--", 3),
    reformat(res_AB, pval_digits = 3)
  )
)

if (is.null(res$inferential[["boot_est"]])) {
  res$inferential[["report_overall_bootCI"]] <- NULL
} else {
  res$inferential[["report_overall_bootCI"]] <- rbind(
    report_table_binary(binobj_ipd, tag = tags[1], eff_measure = eff_measure),
    report_table_binary(binobj_ipd_adj, boot_res_AC, tag = tags[2], eff_measure = eff_measure),
    report_table_binary(binobj_agd, tag = tags[3], eff_measure = eff_measure),
    c(
      paste0("** adj.", trt_ipd, " vs ", trt_agd),
      rep("--", 3),
      reformat(boot_res_AB, pval_digits = 3)
    )
  )
}





















