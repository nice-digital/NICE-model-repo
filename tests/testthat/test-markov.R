library(here)
source(here::here("3_Functions", "markov", "markov.R"))

testthat::test_that("f_markov_M_prep works as expected", {
  
  N_cycle <- 10
  N_state <- 2
  
  nr <- N_cycle
  nc <- 1 + 3 * (N_state) + 1
  
  tp <- matrix(rbeta(nr * nc, shape1 = 1, shape2 = 9), nrow = nr, ncol = nc)
  tp[,1] <- 1:N_cycle
  
  testthat::expect_silent( res <- f_markov_M_prep(tp, N_state) )
  
  testthat::expect_equal(
    lapply(res, names),
    list(
      L1 = c("disc", "next", "death", "stay"),
      L2 = c("disc", "next", "death", "stay"),
      NT = c("death", "stay")
    )
  )
  
  testthat::expect_equal(res$L1$disc, tp[,2])
  testthat::expect_equal(res$L2$`next`, tp[,6])
  testthat::expect_equal(res$NT$death, tp[,8])
})

testthat::test_that("f_markov_calcN works as expected", {
  testthat::expect_equal(f_markov_calcN(4, 20), 143)
})

testthat::test_that("f_markov_topleftFinder works as expected", {
  
  # `nt` is the number of non-tunnel rows, i.e., the number of rows at the top
  # of the transition matrix for health states that do not have tunnel states
  nt <- 2
  
  # `tun_n` is which of the tunnel-based health states the coordinates are being
  # requested for
  tun_n <- 1:4
  
  # `TH` is the number of cycles (time horizon)
  TH <- 5
  
  expected <- list(
    
    # 2nd line
    matrix(
      c(
        3,  4,  # On-treatment -> On-treatment (stay)
        3,  9,  # On-treatment -> Off-treatment (disc)
        3, 13,  # On-treatment -> Next line (next)
        8,  9,  # Off-treatment -> Off-treatment (stay)
        8, 13   # Off-treatment -> Next line (next)
      ),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        c("L2_stay_on", "L2_disc_on", "L2_next_on", "L2_stay_off", "L2_next_off"),
        c("row", "col")
      )
    ),
    
    # 3rd line
    matrix(
      c(
        13, 14,  # On-treatment -> On-treatment (stay)
        13, 19,  # On-treatment -> Off-treatment (disc)
        13, 23,  # On-treatment -> Next line (next)
        18, 19,  # Off-treatment -> Off-treatment (stay)
        18, 23   # Off-treatment -> Next line (next)
      ),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        c("L3_stay_on", "L3_disc_on", "L3_next_on", "L3_stay_off", "L3_next_off"),
        c("row", "col")
      )
    ),
    
    # 4th line
    matrix(
      c(
        23, 24,  # On-treatment -> On-treatment (stay)
        23, 29,  # On-treatment -> Off-treatment (disc)
        23, 33,  # On-treatment -> Next line (next)
        28, 29,  # Off-treatment -> Off-treatment (stay)
        28, 33   # Off-treatment -> Next line (next)
      ),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        c("L4_stay_on", "L4_disc_on", "L4_next_on", "L4_stay_off", "L4_next_off"),
        c("row", "col")
      )
    ),
    
    # 5th line
    matrix(
      c(
        33, 34,  # On-treatment -> On-treatment (stay)
        33, 39,  # On-treatment -> Off-treatment (disc)
        33, 43,  # On-treatment -> Next line (next)
        38, 39,  # Off-treatment -> Off-treatment (stay)
        38, 43   # Off-treatment -> Next line (next)
      ),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        c("L5_stay_on", "L5_disc_on", "L5_next_on", "L5_stay_off", "L5_next_off"),
        c("row", "col")
      )
    )
  )
  
  testthat::expect_equal(
    expected,
    lapply(tun_n, f_markov_topleftFinder, TH = TH, nt = nt)
  )
  
})

testthat::test_that("f_markov_M_compiler works as expected", {
  tp <- list(
    L1 = list(
      disc   = c(0.20, 0.22, 0.24),
      `next` = c(0.20, 0.25, 0.30),
      death  = c(0.05, 0.06, 0.07),
      stay   = c(0.55, 0.47, 0.39)
    ),
    L2 = list(
      disc   = c(0.06, 0.02, 0.12),
      `next` = c(0.05, 0.08, 0.09),
      death  = c(0.08, 0.09, 0.10),
      stay   = c(0.81, 0.81, 0.69)
    ),
    L3 = list(
      disc   = c(0.05, 0.08, 0.06),
      `next` = c(0.03, 0.04, 0.05),
      death  = c(0.12, 0.13, 0.14),
      stay   = c(0.80, 0.75, 0.75)
    ),
    NT = list(
      death  = c(0.15, 0.16, 0.17),
      stay   = c(0.85, 0.84, 0.83)
    )
  )
  
  # `n_lines` is the total number of lines of treatment, excluding BSC
  n_lines <- 3
  
  # `TH` is the time horizon (number of cycles)
  TH <- 3
  
  # `N` is the dimensionality of the transition matrix
  N <- 18
  
  # `nt` is the number of non-tunnel rows
  nt <- 2
  
  calculated <- f_markov_M_compiler(tp, n_lines, TH, N, nt)
  
  expected <- matrix(0, nrow = N, ncol = N)
  
  expected[1,1] <- tp$L1$stay[1]
  expected[1,2] <- tp$L1$disc[1]
  expected[1,3] <- tp$L1$`next`[1]
  expected[1,N] <- tp$L1$death[1]
  
  expected[2,2] <- tp$L1$disc[1] + tp$L1$stay[1]
  expected[2,3] <- tp$L1$`next`[1]
  expected[2,N] <- tp$L1$death[1]
  
  diag(expected[3:4,4:5]) <- tp$L2$stay[1:2]
  diag(expected[3:4,7:8]) <- tp$L2$disc[1:2]
  expected[3:5,9] <- tp$L2$`next`
  expected[3:5,N] <- tp$L2$death
  
  diag(expected[6:7,7:8]) <- tp$L2$disc[1:2] + tp$L2$stay[1:2]
  expected[6:8,9] <- tp$L2$`next`
  expected[6:8,N] <- tp$L2$death
  
  diag(expected[9:10,10:11]) <- tp$L3$stay[1:2]
  diag(expected[9:10,13:14]) <- tp$L3$disc[1:2]
  expected[9:11,15] <- tp$L3$`next`
  expected[9:11,N] <- tp$L3$death
  
  diag(expected[12:13,13:14]) <- tp$L3$disc[1:2] + tp$L3$stay[1:2]
  expected[12:14,15] <- tp$L3$`next`
  expected[12:14,N] <- tp$L3$death
  
  diag(expected[15:16,16:17]) <- tp$NT$stay[1:2]
  expected[15:17,N] <- tp$NT$death
  
  expected[N-1,N-1] <- 1 - expected[N-1,N]
  expected[N,N] <- 1
  
  testthat::expect_equal(as.matrix(calculated), expected)
  
})

testthat::test_that("f_markov_sequenceExtrapolator works as expected", {
  TH <- 5
  L1_tp <- list(
    seq(0.90, 0.82, length.out = 5),  # On 1st line treatment - stay
    seq(0.07, 0.11, length.out = 5),  # On 1st line treatment - discontinue
    seq(0.02, 0.06, length.out = 5),  # On 1st line treatment - jump to 2nd line treatment
    rep(0.01, 5)                      # On 1st line treatment - die
  )
  names(L1_tp) <- c("stay", "disc", "next", "death")
  N <- 18
  
  M <- Matrix::Matrix(data = 0, nrow = 18, ncol = 18, sparse = TRUE)
  M[1,1:3] <- c(0.9, 0.07, 0.02)
  M[2,2:3] <- c(0.97, 0.02)
  M[,N] <- rep(c(0.01, 0.02, 0.1, 1), c(2, 10, 5, 1))
  diag(M[3:6,4:7]) <- 0.9
  M[7,7] <- 0.9
  diag(M[3:6,9:12]) <- 0.05
  M[7,12] <- 0.05
  M[3:12,13] <- 0.03
  diag(M[8:11,9:12]) <- 0.95
  M[12,12] <- 0.95
  diag(M[13:16,14:17]) <- 0.9
  M[17,17] <- 0.9
  
  res <- f_markov_sequenceExtrapolator(L1_tp, TH, M, N)
  
  expected <- Matrix::Matrix(data = c(
    scan(
      text = "
      1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
      0.9	0.07	0.02	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.01
      0.792	0.1392	0.0291	0.018	0	0	0	0	0.001	0	0	0	0.0006	0	0	0	0	0.0201
      0.68112	0.20352	0.037248	0.02619	0.0162	0	0	0	0.001455	0.00185	0	0	0.001443	0.00054	0	0	0	0.030434
      0.5721408	0.2594208	0.044232	0.0335232	0.023571	0.01458	0	0	0.0018624	0.00269175	0.0025675	0	0.00248829	0.0012987	0.000486	0	0	0.04113756
      0.469155456	0.304196832	0.049893696	0.0398088	0.03017088	0.0212139	0.013122	0	0.0022116	0.00344544	0.003735713	0.003168125	0.003690836	0.002239461	0.00116883	0.0004374	0	0.052341032
      ",
      quiet = TRUE
    )
  ), nrow = TH + 1, ncol = N, byrow = TRUE)
  
  # At the moment we are only not including the final row of `expected` but we
  # might want to change `f_markov_sequenceExtrapolator` so that we do want to
  # include it, in which case, just use `expected` in this test instead of
  # `expected[1:TH,]`
  testthat::expect_equal(res, expected[1:TH,])
  
})

testthat::test_that("f_markov_traceConsolidator works as expected", {
  
  full_trace <- new(
    "dgCMatrix",
    i = c(0L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 2L, 3L, 4L, 3L,
      4L, 4L, 2L, 3L, 4L, 3L, 4L, 4L, 2L, 3L, 4L, 3L, 4L, 4L, 1L, 2L, 3L, 4L),
    p = c(0L, 5L, 9L, 13L, 16L, 18L, 19L, 19L, 19L, 22L, 24L, 25L, 25L, 28L,
      30L, 31L, 31L, 31L, 35L),
    Dim = c(5L, 18L),
    Dimnames = list(NULL, NULL),
    x = c(1, 0.9, 0.792, 0.68112, 0.5721408, 0.07, 0.1392, 0.20352, 0.2594208,
      0.02, 0.0291, 0.037248, 0.044232, 0.018, 0.02619, 0.0335232, 0.0162,
      0.023571, 0.01458, 0.001, 0.001455, 0.0018624, 0.00185, 0.00269175,
      0.0025675, 6e-04, 0.001443, 0.00248829, 0.00054, 0.0012987, 0.000486,
      0.01, 0.0201, 0.030434, 0.04113756),
    factors = list()
  )
  
  split_list = list(
    L1_on  = list(any_split = FALSE),
    L1_off = list(any_split = FALSE),
    L2_on  = list(any_split = TRUE, t = c(2, 4)),
    L2_off = list(any_split = FALSE),
    BSC    = list(any_split = TRUE, t = 3)
  )
  
  TH <- 5
  n_lines <- 2
  
  res <- f_markov_traceConsolidator(full_trace, split_list, TH, n_lines)
  
  expected_os <- 1 - c(0, 0.01, 0.0201, 0.030434, 0.04113756, 0.052341032)
  expected_L1 <- c(1, 0.97, 0.9312, 0.88464, 0.8315616, 0.773352288)
  expected_full_lines <- array(
    data = scan(quiet = TRUE, text = "
      1	0.9	0.792	0.68112	0.5721408	0.469155456
      0	0.07	0.1392	0.20352	0.2594208	0.304196832
      0	0.02	0.0471	0.079638	0.1159062	0.154209276
      0	0	0.001	0.003305	0.00712165	0.012560878
      0	0	0.0006	0.001983	0.00427299	0.007536527
      0	0.01	0.0201	0.030434	0.04113756	0.052341032
      "),
    dim = c(6, 6),
    dimnames = list(NULL, c(L1_on = "L1_on", L1_off = "L1_off", L2_on = "L2_on", L2_off = "L2_off", BSC = "BSC", dead = "dead"))
  )
  expected_split_pop <- list(
    L2_on = array(
      data = scan(quiet = TRUE, text = "
        0	0.02	0.0291	0.037248	0.044232 0.049893696
        0	0	0.018	0.04239	0.0716742 0.10431558
        0	0.02	0.0471	0.079638	0.1013262 0.119873376
        0	0	0	0	0.01458 0.0343359
        "),
      dim = c(6, 4),
      dimnames = list(NULL, c("L2_on_split2_before", "L2_on_split2_after", "L2_on_split4_before", "L2_on_split4_after"))
    ),
    L2_off = NULL,
    BSC = array(
      data = scan(quiet = TRUE, text = "
        0	0	0.0006	0.001983	0.00378699	0.005930297
        0	0	0	0	0.000486	0.00160623
        "),
      dim = c(6, 2),
      dimnames = list(NULL, c("BSC_split3_before", "BSC_split3_after"))
    )
  )
  
  testthat::expect_equal(res$OS, expected_os[1:TH])
  testthat::expect_equal(res$L1, expected_L1[1:TH])
  
  # When a data.frame is subset, dimnames[[1]] and dimnames[[2]] for the result
  # are no longer named vectors if they previously were
  testthat::expect_equal(res$full_lines[1:TH,], expected_full_lines[1:TH,])
  
  testthat::expect_equal(res$split_pop$L2_on,  expected_split_pop$L2_on[1:TH,])
  testthat::expect_equal(res$split_pop$L2_off, expected_split_pop$L2_off)
  testthat::expect_equal(res$split_pop$L2_BSC, expected_split_pop$L2_BSC[1:TH,])
  
})
