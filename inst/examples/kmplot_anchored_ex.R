# anchored example using kmplot

### IPD
# Read in relevant ADaM data and rename variables of interest
adsl <- read.csv(system.file("extdata", "adsl.csv",
                             package = "maicplus",
                             mustWork = TRUE
))
adtte <- read.csv(system.file("extdata", "adtte.csv",
                              package = "maicplus",
                              mustWork = TRUE
))
adtte$TIME <- adtte$AVAL
adtte$EVENT <- adtte$EVNT
adtte2 <- adtte
adtte2$ARM <- "C"
adtte2$TIME <- adtte2$TIME + 7
adtte <- rbind(adtte, adtte2)

### AgD
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "aggregate_data_example_1.csv",
                                   package = "maicplus", mustWork = TRUE
))
# for time-to-event endpoints, pseudo IPD from digitalized KM
pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv",
                                   package = "maicplus",
                                   mustWork = TRUE
))
pseudo_ipd$ARM <- "B"
pseudo_ipd2 <- pseudo_ipd
pseudo_ipd2$ARM <- "C"
pseudo_ipd2$Time <- pseudo_ipd2$Time +5
pseudo_ipd <- rbind(pseudo_ipd, pseudo_ipd2)

#### prepare data
target_pop <- process_agd(target_pop)
adsl <- dummize_ipd(adsl, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

#### derive weights
match_res <- estimate_weights(
  data = use_adsl,
  centered_colnames = grep("_CENTERED$", names(use_adsl)),
  start_val = 0,
  method = "BFGS"
)

# plot by trial
kmplot( ipd_weights = match_res,
        tte_dat_ipd = adtte,
        ipd_trt_var = "ARM",
        tte_dat_pseudo = pseudo_ipd,
        pseudo_trt_var = "ARM",
        endpoint_name = "Overall Survival",
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        km_conf_type = "log-log",
        km_layout = "by_trial",
        time_scale="month",
        time_grid = seq(0, 20, by =2),
        use_colors = NULL,
        use_line_types = NULL,
        use_pch_cex = 0.65,
        use_pch_alpha = 100)


# plot by arm
kmplot( ipd_weights = match_res,
        tte_dat_ipd = adtte,
        ipd_trt_var = "ARM",
        tte_dat_pseudo = pseudo_ipd,
        pseudo_trt_var = "ARM",
        endpoint_name = "Overall Survival",
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        km_conf_type = "log-log",
        km_layout = "by_arm",
        time_scale="month",
        time_grid = seq(0, 20, by =2),
        use_colors = NULL,
        use_line_types = NULL,
        use_pch_cex = 0.65,
        use_pch_alpha = 100)

# plot all
kmplot( ipd_weights = match_res,
        tte_dat_ipd = adtte,
        ipd_trt_var = "ARM",
        tte_dat_pseudo = pseudo_ipd,
        pseudo_trt_var = "ARM",
        endpoint_name = "Overall Survival",
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        km_conf_type = "log-log",
        km_layout = "all",
        time_scale="month",
        time_grid = seq(0, 20, by =2),
        use_colors = NULL,
        use_line_types = NULL,
        use_pch_cex = 0.65,
        use_pch_alpha = 100)
