suppressPackageStartupMessages(library(data.table))
source(here::here("3_Functions", "misc", "fpnma_fns.R"))

testthat::test_that("reversing the reference inverts OS and PFS hazard ratios", {
  input <- data.table::data.table(time = rep(2:2080, 2), Molecule = 8,
    Reference.treatment = 10, Line = 2, Endpoint = rep(0:1, each = 2079),
    Population = 0, Reference.trial = 1,
    HR = c(rep(0.4, 2079), rep(0.8, 2079)))
  output <- f_rebase_for_cabo_as_ref_in_2L(input)
  reverse <- output[Molecule == 10 & Reference.treatment == 8]
  # h_A / h_B multiplied by h_B / h_A must equal one.
  testthat::expect_equal(reverse[Endpoint == 0, HR], rep(2.5, 2079))
  testthat::expect_equal(reverse[Endpoint == 1, HR], rep(1.25, 2079))
  testthat::expect_equal(output[Reference.treatment == 10, HR], input$HR)
})
