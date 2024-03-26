source(here::here("3_Functions", "excel", "extract.R"))
source(here::here("3_Functions", "misc", "severity_modifier.R"))

testthat::test_that("get_severity_modifier works", {
  .i <- list(
    R_s_num_qalyShort_abs_LB  = 12,
    R_s_num_qalyShort_abs_UB  = 18,
    R_s_num_qalyShort_prop_LB = 0.85,
    R_s_num_qalyShort_prop_UB = 0.95,
    R_s_num_qalyShort_w_LB    = 1,
    R_s_num_qalyShort_w_bet   = 1.2,
    R_s_num_qalyShort_w_UB    = 1.7
  )
  
  sm <- severity_modifier(.i)
  
  test_cases <- data.frame(
    qalys_disease   = c(0.8, 0.4, 0.6,  0.4,  3.1,  5.4,  3.5, 16.4,  9.0,  0.5),
    qalys_nodisease = c(1.6, 4.1, 4.5, 10.3, 10.5, 11.2, 15.7, 20.5, 23.0, 28.1),
    weights         = c(1.0, 1.2, 1.2,  1.7,  1.0,  1.0,  1.2,  1.0,  1.2,  1.7)
  )
  
  testthat::expect_equal(
    mapply(
      FUN             = get_severity_modifier,
      qalys_disease   = test_cases$qalys_disease,
      qalys_nodisease = test_cases$qalys_nodisease,
      MoreArgs        = list(.severity_modifier = sm)
    ),
    test_cases$weights
  )
})

testthat::test_that("it calculates severity weights for aggregate data", {
  i <- f_excel_extract(here::here("1_Data", "PATT RCC_model inputs.xlsx"))
  
  p <- list(basic = list(cl_y = 1/52, disc_q = 0.035))
  p <- add_population_utility_params(p, .i = i)
  
  # Modifier: 1.7x (absolute shortfall > 18)
  sm1 <- calc_severity_modifier(5.0, 0.5, 2.1, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm1), 1.7)
  
  # Modifier: 1.2x (12 < abs. shortfall < 18 and 0.85 < prop. shortfall < 0.95)
  sm2 <- calc_severity_modifier(40.8, 0.5, 2.5, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm2), 1.2)
  
  # Modifier: 1.0x (abs. shortfall < 12 and prop. shortfall < 0.85)
  sm3 <- calc_severity_modifier(40.8, 0.5, 6.9, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm3), 1.0)
  
  # Modifier: 1.7x (prop. shortfall > 0.95)
  sm4 <- calc_severity_modifier(60.0, 0.5, 0.3, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm4), 1.7)
  
})

testthat::test_that("it calculates severity weights for IPD", {
  i <- f_excel_extract(here::here("1_Data", "PATT RCC_model inputs.xlsx"))
  
  p <- list(basic = list(cl_y = 1/52, disc_q = 0.035))
  p <- add_population_utility_params(p, .i = i)
  
  age <- round(rnorm(100, 8760, 1000)) / 365 # Mean age 24 years
  sex <- c("female", "male")[rbinom(100, 1, 0.5)+1]
  
  sm1 <- calc_severity_modifier(age, sex, 0.5, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm1), 1.7)
  
  sm2 <- calc_severity_modifier(age, sex, 2.5, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm2), 1.7)
  
  sm3 <- calc_severity_modifier(age, sex, 6.0, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm3), 1.2)
  
  sm4 <- calc_severity_modifier(age, sex, 11.0, .i = i, .p = p)
  testthat::expect_equal(as.vector(sm4), 1.0)
  
})
