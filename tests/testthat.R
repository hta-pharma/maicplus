pkg_name <- "maicplus"
library(pkg_name, character.only = TRUE)
library(flexsurv)
library(clubSandwich)
testthat::test_check(pkg_name)
