source(here::here("3_Functions", "survival", "treatment_effect_waning.R"))

testthat::test_that("treatment effect wanes over time", {
  
  start_cycle <- 24
  finish_cycle <- 36
  
  t <- seq(0, 5, length.out = 61)
  
  surv_active <- pweibull(t, shape = 1.2, scale = 5.0, lower.tail = FALSE)
  surv_ref    <- pweibull(t, shape = 1.2, scale = 2.5, lower.tail = FALSE)
  
  res <- treatment_effect_waning(surv_active, surv_ref, start_cycle, finish_cycle)
  
  # Calculated in Excel
  target <- c(1,0.992678102,0.983258533,0.972909693,0.961955215,0.950566622,0.938853589,0.926893604,0.914744939,0.902453319,0.890055742,0.877582838,0.865060413,0.852510516,0.83995219,0.827402031,0.814874602,0.802382763,0.789937915,0.77755021,0.765228711,0.752981531,0.740815937,0.728738455,0.716754939,0.704233849,0.690574022,0.675834556,0.660081347,0.643386456,0.625827422,0.607486516,0.588449961,0.568807117,0.548649641,0.528070644,0.507163842,0.486498698,0.466570363,0.447359589,0.428847082,0.411013568,0.393839841,0.377306815,0.361395562,0.346087351,0.331363678,0.3172063,0.30359725,0.290518868,0.277953811,0.265885074,0.254295999,0.243170288,0.232492009,0.222245606,0.2124159,0.202988094,0.193947777,0.18528092,0.176973879 )
  
  testthat::expect_equal(res, target)
  
})

testthat::test_that("treatment effect bypass flag works", {
  
  start_cycle <- 24
  finish_cycle <- 36
  
  t <- seq(0, 5, length.out = 61)
  
  surv_active <- pweibull(t, shape = 1.2, scale = 5.0, lower.tail = FALSE)
  surv_ref    <- pweibull(t, shape = 1.2, scale = 2.5, lower.tail = FALSE)
  
  res <- treatment_effect_waning(surv_active, surv_ref, start_cycle, finish_cycle, apply_waning = FALSE)
  
  testthat::expect_equal(res, surv_active)
  
})

testthat::test_that("treatment effect wanes immediately", {
  
  start_cycle <- 24
  finish_cycle <- 24
  
  t <- seq(0, 5, length.out = 61)
  
  surv_active <- pweibull(t, shape = 1.2, scale = 5.0, lower.tail = FALSE)
  surv_ref    <- pweibull(t, shape = 1.2, scale = 2.5, lower.tail = FALSE)
  
  res <- treatment_effect_waning(surv_active, surv_ref, start_cycle, finish_cycle)
  
  expected <- c(
    surv_active[1:25],
    surv_active[25] * surv_ref[26:61] / surv_ref[25]
  )
  
  testthat::expect_equal(res, expected)
  
})

testthat::test_that("treatment effect wanes based on end of cycle", {
  
  start_cycle <- 24
  finish_cycle <- 36
  
  t <- seq(0, 5, length.out = 61)
  
  surv_active <- pweibull(t, shape = 1.2, scale = 5.0, lower.tail = FALSE)
  surv_ref    <- pweibull(t, shape = 1.2, scale = 2.5, lower.tail = FALSE)
  
  res <- treatment_effect_waning(surv_active, surv_ref, start_cycle, finish_cycle, wcc = 1)
  
  # Calculated in Excel
  target <- c(1,0.992678102,0.983258533,0.972909693,0.961955215,0.950566622,0.938853589,0.926893604,0.914744939,0.902453319,0.890055742,0.877582838,0.865060413,0.852510516,0.83995219,0.827402031,0.814874602,0.802382763,0.789937915,0.77755021,0.765228711,0.752981531,0.740815937,0.728738455,0.716754939,0.703597625,0.689321812,0.673989976,0.65767116,0.64044032,0.6223776,0.603567568,0.584098411,0.564061103,0.543548561,0.522654784,0.501474013,0.481040709,0.461335949,0.442340699,0.424035882,0.406402441,0.389421385,0.373073841,0.357341095,0.342204625,0.327646137,0.313647589,0.300191218,0.287259561,0.27483547,0.262902131,0.251443073,0.24044218,0.229883701,0.219752251,0.210032824,0.200710788,0.191771893,0.183202269,0.174988424 )
  
  testthat::expect_equal(res, target)
  
})

testthat::test_that("treatment_effect_waning produces a warning", {
  
  start_cycle <- 24
  finish_cycle <- 36
  
  t <- seq(0, 5, length.out = 61)
  
  surv_active <- pweibull(t, shape = 2.5, scale = 4.0, lower.tail = FALSE)
  surv_ref    <- pweibull(t, shape = 1.2, scale = 2.5, lower.tail = FALSE)
  
  testthat::expect_warning(
    res <- treatment_effect_waning(surv_active, surv_ref, start_cycle, finish_cycle, wcc = 1),
    regexp = paste(
      "Hazard rate in the active treatment is more than the hazard rate in the",
      "reference treatment in 18 cycles, the first of which is 42"
    )
  )
  
})
