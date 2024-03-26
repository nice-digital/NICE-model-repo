#' Simple function to compute a discount factor given time vector in years from
#' `t=0` to time horizon, the discount rate, and a flag for method.
#'
#' @param r discount rate (annual)
#' @param t_yr vector of time in years of length time horizon
#' @method flag for method
#'
#'
f_discFac <- function(r, t_yr, method = "classic cycle time") {
  stopifnot(method %in% c("classic floor","classic cycle time","continuous"))
  if(method == "classic floor") {
    # classic increment discounting using floor of year
    return(1 / ((1 + r) ^ floor(t_yr)))
  } else if (method == "classic cycle time") {
    return(1 / ((1 + r) ^ t_yr))
  } else if (method == "continuous") {
    # Continuous discounting using e^rt
    1/exp(r * t_yr)
  }
}

#' Time step converter using exponential assumption (Box 3.1, Briggs et al.)
#' Get the ratio of target time step to origin time step, Convert to r
#' in SAME time unit (i.e. / 1), THEN convert to P in target time unit
#' 
#' @param cl_orig original cycle length
#' @param cl_target target cycle length
#' @param P probability in one original cycle length
#' 
f_misc_tConvertExp <- function(cl_orig, cl_target, P) {
  t <- cl_target / cl_orig
  r <- -log(1-P)
  return(1-exp(-r*t))
}
