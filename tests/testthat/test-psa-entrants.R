# Isolate the original entrants block; no survival extrapolation or full PSA is needed.
src <- readLines(here::here("3_Functions", "psa", "psa functions.R"), warn = FALSE)
start <- grep("^        ENTRANTS <- matrix", src)
end <- grep("^        # Discounted versions:", src)[1] - 1

testthat::test_that("entrants include flows from each occupied source state", {
  M_nam <- c("L1_on", "L1_off", "L2_on", "L2_off", "BSC", "dead")
  M <- list(rbind(c(.6,.1,.2,0,0,.1), c(0,.6,.3,0,0,.1),
    c(0,0,.5,.2,.2,.1), c(0,0,0,.4,.4,.2), c(0,0,0,0,.7,.3), c(0,0,0,0,0,1)))
  tTRACE <- matrix(c(.6,.1,.15,.05,.05,.05), ncol = 1)
  eval(parse(text = src[start:end]))
  # L2 entry=.6*.2+.1*.3=.15; death=.06+.01+.015+.01+.015=.11.
  testthat::expect_equal(as.numeric(ENTRANTS[1, ]), c(.06, .15, .03, .05, .11))
  testthat::expect_equal(as.numeric(10000 * ENTRANTS[1, "dead"]), 1100)
})
