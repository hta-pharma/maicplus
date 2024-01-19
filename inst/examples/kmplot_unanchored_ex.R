# unanchored example using maicplus_kmplot

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

# plot
kmplot(
  weights_object = match_res,
  tte_ipd = adtte,
  trt_var_ipd = "ARM",
  tte_pseudo_ipd = pseudo_ipd,
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = NULL,
  km_conf_type = "log-log",
  time_scale = "month",
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)
