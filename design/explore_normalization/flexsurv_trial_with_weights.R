# R Script to fit a weighted parametric survival model using the flexsurvreg function
# directly, without a separate function definition.
#
# This script performs the following steps:
# 1. Loads the required 'flexsurv' package.
# 2. Generates sample survival data.
# 3. Creates a numeric vector of case weights.
# 4. Fits a weighted Weibull model using flexsurvreg.
# 5. Prints a summary of the fitted weighted model.
# 6. For comparison, fits and summarizes an unweighted model.
#
# To run this script, you must have the 'flexsurv' package installed.
# You can install it with: install.packages("flexsurv")

# Load the flexsurv package.
library(flexsurv)

# --- Step 1: Generate some sample survival data ---
set.seed(42)
n <- 200
sample_data <- data.frame(
  time = rexp(n, rate = 0.1),
  status = rbinom(n, 1, 0.8),
  group = factor(rep(c("A", "B"), each = n / 2)),
  age = rnorm(n, 50, 10)
)

# --- Step 2: Create a weighting variable ---
# Let's say observations from group 'B' are twice as important.
weights_vector <- ifelse(sample_data$group == "A", 1, 2)
weights_vector_SW1 <- weights_vector / sum(weights_vector)

# --- Step 3: Fit the weighted Weibull model directly using flexsurvreg ---
weighted_weibull_model <- flexsurvreg(
  formula = Surv(time, status) ~ group + age,
  data = sample_data,
  weights = weights_vector,
  dist = "weibull"
)

weighted_weibull_model_SW1 <- flexsurvreg(
  formula = Surv(time, status) ~ group + age,
  data = sample_data,
  weights = weights_vector_SW1,
  dist = "weibull"
)

# --- Step 4: Print a summary of the fitted weighted model ---
print("Summary of the weighted Weibull model:")
weighted_weibull_model
weighted_weibull_model_SW1

print(summary(weighted_weibull_model)[[1]][1:10, ])
print(summary(weighted_weibull_model_SW1)[[1]][1:10, ])

## Conclusion: estimates are the same with different normalization
