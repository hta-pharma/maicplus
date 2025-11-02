# Example with single arm
continuous_agd <- data.frame(
  ARM = "B",
  MEAN = -12,
  SD = 9,
  N = 120
)

pseudo_ipd <- get_pseudo_ipd_continuous(continuous_agd)
head(pseudo_ipd)

# Example with multiple arms
continuous_agd_multi <- data.frame(
  ARM = c("B", "C"),
  MEAN = c(-12, -9),
  SD = c(9, 9.5),
  N = c(120, 150)
)

pseudo_ipd_multi <- get_pseudo_ipd_continuous(continuous_agd_multi)
head(pseudo_ipd_multi)
aggregate(AVAL ~ ARM, data = pseudo_ipd_multi, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

