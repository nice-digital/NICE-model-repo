# There is no qdirichlet function in R.
# 
# Fortunately gtools has a rdirichlet function which is easily adapted.
# 
# We simply change the rgamma call to qgamma!

qdirichlet <- function (alpha, rands) {
  
  # rands must be of length l*n
  
  l  <- length(alpha)
  rl <- length(rands)
  
  if ((rl/l) %% 1 != 0) stop(paste0("The number of rands (",rl,") is not a multiple of the length of alpha (",l,")"))
  
  x <- matrix(qgamma(rands, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}

if(FALSE) {
  # Example:
  counts <- c(95,2,1)
  ndraw  <- 1000
  rands  <- runif(length(counts)*ndraw)
  
  dirich_dr <- qdirichlet(
    alpha = counts,
    rands = rands
  )
  
  all(round(rowSums(dirich_dr),12)==1)
  
  
}

#' Function to estimate a and b in a beta distribution from mean and variance
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#' Function to estimate `sd(log(x))` from mean, LB and UB
estSDlog <- function(mean, LB, UB) {
  LB_Sdlog <- (log(mean) - log(LB)) / 1.96
  UB_Sdlog <- (log(UB) - log(mean)) / 1.96
  sdlog    <- mean(LB_Sdlog, UB_Sdlog)
  return(sdlog)
}