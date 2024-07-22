library(survival)
data(adtte_sat)
data(pseudo_ipd_sat)

kmobj_A <- survfit(Surv(TIME, EVENT) ~ ARM,
  data = adtte_sat,
  conf.type = "log-log"
)

kmobj_B <- survfit(Surv(TIME, EVENT) ~ ARM,
  data = pseudo_ipd_sat,
  conf.type = "log-log"
)

kmlist <- list(kmobj_A = kmobj_A, kmobj_B = kmobj_B)
kmlist_name <- c("A", "B")

basic_kmplot2(kmlist, kmlist_name)
