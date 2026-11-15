source(here::here("3_Functions", "sequencing", "sequences.R"), local = TRUE)

testthat::test_that("a later-line treatment requires its immediate predecessor", {
  paths <- data.frame(
    L1 = c("IO", "cabo", "cabo", "pazo", "other", "IO"),
    L2 = c("pazo", "pazo", "other", "cabo", "cabo", "other"),
    L3 = c("BSC", "BSC", "pazo", "BSC", "BSC", "pazo"),
    L4 = c("", "", "BSC", "", "", "BSC")
  )
  invisible(capture.output(
    actual <- f_path_2L_only_immediate_after(paths, "pazo", "IO")
  ))
  # Row1 satisfies the rule, row4 is first-line pazo, row5 never uses pazo.
  testthat::expect_equal(actual, paths[c(1, 4, 5), , drop = FALSE])
})
