library(data.table)
source(here::here("3_Functions", "costs_and_QALYs", "cost_processing.R"))

testthat::test_that("monthly and integer-week schedules retain their dosing periods", {
  drugs <- data.frame(Treatment.name = c("monthly", "fortnightly"), Molecule = 1:2,
    Type = "per cycle", Drug.cost.per.dose = 100,
    Admin.cost..per.administration..Mean = 0,
    Applied.every.x.cycles = c(30/7, 2), Number.of.doses.when.applied = 1,
    Time.from = 1, Time.to..cycles. = 104)
  for (line in 1:4) drugs[[paste0("RDI.Mean_Trial_", line, "L")]] <- 1
  mru <- data.frame(Type.of.cost = c(NA, "On", "Off", "On", "Off"),
    Type = c(NA, rep("per cycle", 4)), Molecule = c(NA, 1, 1, 2, 2),
    Cost = c("Mean", rep("0", 4)), X5 = c("SE", rep("0", 4)),
    Time.from..cycle. = c(NA, rep(1, 4)), Time.to..cycles. = c(NA, rep(104, 4)))
  result <- f_process_cost_data(drugs, mru, 103, "Base case", 4)
  # Weekly indices 1..104 span 721 days: starts at days 0,30,...,720 give 25 packs.
  testthat::expect_equal(sum(result$drug$line_1$mol_1), 25 * 100)
  testthat::expect_equal(sum(result$drug$line_1$mol_2), 52 * 100)
})

# Check both scheduling blocks, including PSA, without its optional progressr dependency.
src <- readLines(here::here("3_Functions", "costs_and_QALYs", "cost_processing.R"), warn = FALSE)
by_line <- grep("current_row$Applied.every.x.cycles", src, fixed = TRUE)
testthat::test_that("both costing branches assign the same calendar doses to weeks", {
  current_row <- list(Time.from = 1, Time.to..cycles. = 104, Applied.every.x.cycles = 30/7)
  calendar_weeks <- floor((0:24) * 30 / 7) + 1 # 25 starts at days 0..720.
  testthat::expect_length(by_line, 2)
  for (j in by_line) {
    eval(parse(text = src[(j - 3):(j + 1)]))
    testthat::expect_equal(dose_at_cycles, calendar_weeks)
  }
})
