# anchored example using maic_anchored for tte
library(flexsurv)
devtools::load_all()
### IPD
set.seed(1234)
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
adtte2$TIME <- adtte2$TIME * runif(nrow(adtte2), 0.15,0.3)
fit_C <- flexsurv::flexsurvspline(formula = Surv(TIME, EVENT) ~ 1, data = adtte2, k=3)
tmp <- simulate(fit_C, nsim=1, seed=1234, newdata=adtte2, censtime = max(adtte$TIME))
adtte2$TIME <- tmp$time_1
adtte2$EVENT <- tmp$event_1
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
pseudo_ipd2 <- adtte2[,c("TIME","EVENT","ARM")]
names(pseudo_ipd2) <- c("Time", "Event", "ARM")
tmp <- simulate(fit_C, nsim=1, seed=4321, newdata=adtte2, censtime = max(pseudo_ipd$Time))
pseudo_ipd2$Time<- tmp$time_1
pseudo_ipd2$Event <- tmp$event_1
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

# inferential result
result <- maic_anchored(ipd_weights = match_res,
                        dat_ipd = adtte,
                        ipd_trt_var = "ARM",
                        dat_pseudo = pseudo_ipd,
                        pseudo_trt_var = "ARM",
                        trt_ipd = "A",
                        trt_agd = "B",
                        trt_common = "C",
                        endpoint_name = "Overall Survival",
                        endpoint_type = "tte",
                        eff_measure = "HR",
                        time_scale = "month",
                        km_conf_type = "log-log")
result$inferential$report_median_surv
result$inferential$report_overall

ph_diagplot(ipd_weights = match_res,
            tte_dat_ipd = adtte,
            ipd_trt_var = "ARM",
            tte_dat_pseudo = pseudo_ipd,
            pseudo_trt_var = "ARM",
            trt_ipd = "A",
            trt_agd = "B",
            trt_common = "C",
            endpoint_name = "Overall Survival",
            time_scale="week",
            zph_transform = "log",
            zph_log_hazard = TRUE)

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









