#### create datasets for unanchored case ####
## adsl_sat, adtte_sat, adrs_sat, agd (AgD of effect modifiers), pseudo_ipd_sat (AgD, tte data)

devtools::load_all()
# Read in relevant ADaM data and rename variables of interest
adsl_sat <- read.csv(system.file("extdata", "adsl.csv",
  package = "maicplus",
  mustWork = TRUE
))
adsl_sat$X <- NULL
adtte_sat <- read.csv(system.file("extdata", "adtte.csv",
  package = "maicplus",
  mustWork = TRUE
))
adtte_sat$TIME <- adtte_sat$AVAL
adtte_sat$EVENT <- adtte_sat$EVNT
adtte_sat <- adtte_sat[adtte_sat$ARM == "A", , drop = FALSE]
adsl_sat <- adsl_sat[adsl_sat$USUBJID %in% adtte_sat$USUBJID, , drop = FALSE]

### AgD
# Baseline aggregate data for the comparator population
agd <- read.csv(system.file("extdata", "aggregate_data_example_1.csv",
  package = "maicplus", mustWork = TRUE
))
# for time-to-event endpoints, pseudo IPD from digitalized KM
pseudo_ipd_sat <- read.csv(system.file("extdata", "psuedo_IPD.csv",
  package = "maicplus",
  mustWork = TRUE
))
pseudo_ipd_sat$ARM <- "B"

### Centered IPD
agd <- process_agd(agd)
adsl_sat <- dummize_ipd(adsl_sat, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
centered_ipd_sat <- center_ipd(ipd = adsl_sat, agd = agd_sat)

### Binary
adrs_sat <- read.csv(system.file("extdata", "adrs.csv", package = "maicplus", mustWork = TRUE))
adrs_sat$RESPONSE <- adrs_sat$AVAL

### Output
usethis::use_data(adsl_sat, adtte_sat, agd, pseudo_ipd_sat, centered_ipd_sat, adrs_sat,
  internal = FALSE, overwrite = TRUE
)
