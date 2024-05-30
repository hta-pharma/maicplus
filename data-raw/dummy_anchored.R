#### create anchored example datasets ####

devtools::load_all()
library(flexsurv)
set.seed(2024)

# create adsl_twt
adsl <- read.csv(system.file("extdata", "adsl.csv",
  package = "maicplus",
  mustWork = TRUE
))
adsl$USUBJID <- paste0("xx", adsl$USUBJID)
adsl2 <- adsl
adsl2$USUBJID <- sample(size = nrow(adsl2), paste0("yy", adsl2$USUBJID), replace = FALSE)
adsl2 <- adsl2[order(adsl2$USUBJID), ]

adsl_twt <- rbind(adsl, adsl2)

# create adtte_twt
adtte <- read.csv(system.file("extdata", "adtte.csv",
  package = "maicplus",
  mustWork = TRUE
))
adtte$TIME <- adtte$AVAL
adtte$EVENT <- adtte$EVNT
adtte$USUBJID <- paste0("xx", adtte$USUBJID)

adtte2 <- adtte
adtte2$ARM <- "C"
adtte2$TIME <- adtte2$TIME * runif(nrow(adtte2), 0.15, 0.3)
fit_C <- flexsurv::flexsurvspline(formula = Surv(TIME, EVENT) ~ 1, data = adtte2, k = 3)
tmp <- simulate(fit_C, nsim = 1, seed = 1234, newdata = adtte2, censtime = max(adtte$TIME))
adtte2$TIME <- tmp$time_1
adtte2$EVENT <- tmp$event_1
adtte2$USUBJID <- paste0("yy", adtte2$USUBJID)

adtte_twt <- rbind(adtte, adtte2)


data("agd")
# create pseudo_ipd_twt
pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv",
  package = "maicplus",
  mustWork = TRUE
))
pseudo_ipd$ARM <- "B"
pseudo_ipd2 <- adtte2[, c("TIME", "EVENT", "ARM")]
names(pseudo_ipd2) <- c("Time", "Event", "ARM")
tmp <- simulate(fit_C, nsim = 1, seed = 4321, newdata = adtte2, censtime = max(pseudo_ipd$Time))
pseudo_ipd2$Time <- tmp$time_1
pseudo_ipd2$Event <- tmp$event_1

pseudo_ipd_twt <- rbind(pseudo_ipd, pseudo_ipd2)

# create centered adsl_twt
agd <- process_agd(agd)
adsl_twt <- dummize_ipd(adsl_twt, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
centered_ipd_twt <- center_ipd(ipd = adsl_twt, agd = agd)


### Output
usethis::use_data(adsl_twt, adtte_twt, pseudo_ipd_twt, centered_ipd_twt,
  internal = FALSE, overwrite = TRUE
)
