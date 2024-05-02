# simulate data to demo binary maic
devtools::load_all()
load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
ipd_centered <- center_ipd(ipd = ipd, agd = agd)

centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")
weighted_data <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames)
weighted_data2 <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames, n_boot_iteration = 100)
check <- weighted_data2$boot


adrs <- read.csv(system.file("extdata", "adrs.csv", package = "maicplus", mustWork = TRUE))
adrs$RESPONSE <- adrs$AVAL

pseudo_adrs <- get_pseudo_ipd_binary(binary_agd = data.frame(ARM = rep("B",2),
                                                             RESPONSE = c("YES","NO"),
                                                             COUNT = c(280,120)),
                                     format = "stacked")


maic_unanchored(weights_object = weighted_data,
                ipd = adrs,
                pseudo_ipd = pseudo_adrs,
                trt_ipd = "A",
                trt_agd = "B",
                trt_var_ipd = "ARM",
                trt_var_agd = "ARM",
                endpoint_type = "binary",
                endpoint_name = "Binary Endpoint",
                eff_measure = "RR",
                # binary specific args
                binary_robust_cov_type = "CR2")



















