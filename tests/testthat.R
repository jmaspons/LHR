library(testthat)
library(LHR)

test_check("LHR")
# test_check("LHR", filter="^(LH|Env|Sim|Model)")

# devtools::test(filter="^(LH|Env|Sim|Model)") # test main classes only
