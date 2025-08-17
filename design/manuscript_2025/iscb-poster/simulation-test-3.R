setwd("/cloud/project/design/manuscript_2025/iscb-poster")
source("simulation-3.R")
test_params <- sim_params  # Copy defaults
# Modify as needed...
test_results <- main(test_params)
