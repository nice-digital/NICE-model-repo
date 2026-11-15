suppressPackageStartupMessages(library(data.table))
source(here::here("3_Functions", "survival", "Survival_functions.R"))
source(here::here("3_Functions", "misc", "fpnma_fns.R"))

testthat::test_that("the FP network consumes the current posterior draw", {
  ratios <- data.table::data.table(Population = 0, Line = 1, Molecule = 1,
    Endpoint = 0, Reference.treatment = 7, Reference.trial = 0,
    time = 0:2, HR = 0.7)
  destinations <- unique(ratios[, c("Population", "Line", "Molecule",
    "Endpoint", "Reference.treatment", "Reference.trial"), with = FALSE])
  p <- list(basic = list(th = 2), releff = list(
    CODA = list(FP = ratios), fp_dest = destinations))
  code <- readLines(here::here("2_Scripts", "probabilistic_model.R"))
  # Evaluate the actual production call, without sourcing the full model script.
  call <- parse(text = code[654:659])
  for (sampled_hr in c(0.4, 1.6)) {
    fresh <- data.table::copy(ratios)
    fresh$HR <- sampled_hr
    p_psa <- list(releff = list(CODA = list(FP = fresh),
      fp_dest = destinations))
    network <- list()
    eval(call)
    actual <- network$pop_0$line_1$mol_1$trial_0$endpoint_0$fp$HR
    testthat::expect_equal(actual, rep(sampled_hr, 3))
  }
})

testthat::test_that("repeated reference expansion has unique lookup keys", {
  ratios <- data.table::data.table(Population = 0, Line = 2, Molecule = 8,
    Endpoint = 0, Reference.treatment = 10, Reference.trial = rep(1:2, each = 3),
    time = rep(0:2, 2), HR = 0.6)
  p_psa <- list(releff = list(CODA = list(FP = ratios)))
  i <- list(FPNMA = list(destinations = f_gen_destinations(ratios)))
  code <- readLines(here::here("2_Scripts", "probabilistic_model.R"))
  eval(parse(text = code[638:644]))
  keys <- c("Population", "Line", "Molecule", "Endpoint",
    "Reference.treatment", "Reference.trial", "time")
  testthat::expect_false(anyDuplicated(p_psa$releff$CODA$FP[, ..keys]) > 0)
  testthat::expect_false(anyDuplicated(p_psa$releff$fp_dest) > 0)
  testthat::expect_equal(unique(p_psa$releff$CODA$FP$HR), 0.6)
})
