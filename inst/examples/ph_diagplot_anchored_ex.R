# anchored example using ph_diagplot

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


ph_diagplot(weights_object = match_res,
            tte_ipd = adtte,
            trt_var_ipd = "ARM",
            tte_pseudo_ipd = pseudo_ipd,
            trt_var_agd = "ARM",
            trt_ipd = "A",
            trt_agd = "B",
            trt_common = "C",
            endpoint_name = "Overall Survival",
            time_scale="week",
            zph_transform = "log",
            zph_log_hazard = TRUE)

