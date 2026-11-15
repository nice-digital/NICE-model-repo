suppressPackageStartupMessages(library(data.table))
source(here::here("3_Functions", "misc", "fpnma_fns.R"))

testthat::test_that("zero FP powers use log time in either position", {
  # With log(HR) = log(t), HR must equal t at positive times.
  first <- f_HRs_by_cycle(4, 1, c(0, 1), c(0, 1, 0))[-1]
  second <- f_HRs_by_cycle(4, 1, c(1, 0), c(0, 0, 1))[-1]
  testthat::expect_equal(first, 1:4)
  testthat::expect_equal(second, 1:4)
  unchanged <- f_HRs_by_cycle(4, 1, c(1, 2), c(0, 0.1, 0.2))
  testthat::expect_equal(unchanged, exp(0.1 * (0:4) + 0.2 * (0:4)^2))
})

testthat::test_that("CODA time zero is neutral even when exp(log HR) is zero", {
  coefficients <- data.table::data.table(Population = 0, Line = 1,
    Molecule = 1, Endpoint = 0, Reference.treatment = 7, Reference.trial = 0,
    Exponent.1 = 0, Exponent.2 = 1, intd = 0, s1 = 1, s2 = 0)
  output <- f_generate_FPNMA_coda(coefficients, TH = 4, wks_per_month = 1)
  testthat::expect_equal(output$HR, c(1, 1:4))
})
