# example of unstacked
testdat <- data.frame(Yes = 280, No = 120)
rownames(testdat) <- "B"
get_pseudo_ipd_binary(
  binary_agd = testdat,
  format = "unstacked"
)

# example of stacked
get_pseudo_ipd_binary(
  binary_agd = data.frame(
    ARM = rep("B", 2),
    RESPONSE = c("YES", "NO"),
    COUNT = c(280, 120)
  ),
  format = "stacked"
)
