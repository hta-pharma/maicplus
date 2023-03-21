

test_that("Test weight functions", {
  require(dplyr)
  require(ggplot2)
  
  # test set
  # wt is uniform 1
  # wt_rs is 1 for age 50, 0 for age 60
  test_weights <- data.frame(wt = 1,
                             wt_rs = c(rep(1,5), rep(0,5)),
                             age = c(rep(50,5), rep(60,5))
  )
  
  test_diag <- wt_diagnostics(test_weights, vars = "age")
  
  #####################################################
  # check diagnostics as expected
  expect_equal(test_diag$ESS, 10)
  expect_equal(test_diag$Summary_of_weights$mean, c(1, 0.5))
  expect_equal(test_diag$Weight_profiles$age, c(50,60))
  
  
  test_profile <- profile_wts(test_weights, vars = "age")
  #####################################################
  # check diagnostics as expected
  expect_equal(test_profile$age, c(50,60))
  expect_equal(test_profile$wt, c(1,1))
  expect_equal(test_profile$wt_rs, c(1,0))
  
  test_hist <- hist_wts(test_weights)
  
  # extract plot data
  # panel 1  is the rescaled weights and should be 5 for 0 and 5 for 1
  test_hist_data <- layer_data(test_hist, 1) %>%
    filter(count > 0, PANEL == 1)
  
  #####################################################
  # check diagnostics as expected
  expect_equal(test_hist_data$x, c(0,1))
  expect_equal(test_hist_data$y, c(5,5))
  
  
})
