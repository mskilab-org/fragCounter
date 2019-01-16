Sys.setenv("R_TESTS" = "")

library(testthat)
library(fragCounter)

test_check("fragCounter")

