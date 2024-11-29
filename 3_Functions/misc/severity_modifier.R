# Ensure we have access to necessary commands from elsewhere
if (!("utility_genpop" %in% names(globalenv()))) {
  source(here::here("3_Functions", "utility", "age_related.R"))
}

if (!("get_lifetables" %in% names(globalenv()))) {
  source(here::here("3_Functions", "survival", "other_cause_mortality.R"))
}

#' Extract the severity modifier parameters from Excel and wrap them up nicely
#' in a list
#'
#' @param .i  You can inject `.i`, which is the result of extracting variables
#'            from an Excel file. If you don't inject it, the function will look
#'            for `i` in the global environment. If it can't find that it will
#'            provide suitable default values accompanied by a warning.
#'
#' @param hard_coded If TRUE, use the thresholds from NICE methods guidance circa 2022
#'
severity_modifier <- function(.i = NULL, hard_coded = FALSE) {
  # If `.i` was not provided then grab `i` from the global environment
  if (hard_coded) {
    return(list(
      breaks_absolute = c(-Inf, 12, 18, Inf),
      breaks_proportional = c(-Inf, 0.85, 0.95, Inf),
      qaly_weights = c(1, 1.2, 1.7)
    ))
  } else if (is.null(.i)) {
    if ("i" %in% names(globalenv())) {
      .i <- get("i", envir = globalenv())
    } else {
      warning(paste("`.i` was not provided and `i` is not in the global environment,", "hard-coded values are being used"))
      return(list(
        breaks_absolute = c(-Inf, 12, 18, Inf),
        breaks_proportional = c(-Inf, 0.85, 0.95, Inf),
        qaly_weights = c(1, 1.2, 1.7)
      ))
    }
  }

  # Different names may be used
  param_names <- array(
    data = c(
      "R_s_num_qalyShort_abs_LB",
      "R_s_num_qalyShort_abs_UB",
      "R_s_num_qalyShort_prop_LB",
      "R_s_num_qalyShort_prop_UB",
      "R_s_num_qalyShort_w_LB",
      "R_s_num_qalyShort_w_bet",
      "R_s_num_qalyShort_w_UB",
      "i_QALY_short_abs_LB",
      "i_QALY_short_abs_UB",
      "i_QALY_short_prop_LB",
      "i_QALY_short_prop_UB",
      "i_qaly_weight_LB",
      "i_qaly_weight_between",
      "i_qaly_weight_UB"
    ),
    dim = c(7, 2),
    dimnames = list(
      item = c(
        "Abs. shortfall moderate threshold",
        "Abs. shortfall extreme threshold",
        "Prop. shortfall extreme threshold",
        "Prop. shortfall moderate threshold",
        "Basic QALY weight",
        "Moderate severity weight",
        "Extreme severity weight"
      ),
      naming_system = list("R_s_", "i_")
    )
  )

  # Check we can at least access one per row
  variables_present <- apply(
    param_names,
    1,
    function(row) row[1] %in% names(.i) | row[2] %in% names(.i)
  )
  if (!all(variables_present)) {
    stop(
      "Could not retrieve parameters from `i` or `.i`:\n",
      paste("-", dimnames(param_names)[[1]][!variables_present], collapse = "\n")
    )
  }

  extract_parameter <- function(item) {
    possible_names <- param_names[item, ]
    values <- lapply(possible_names, function(name) .i[[name]])
    values[[which.max(!sapply(values, is.null))]]
  }

  # Grab the values from .i and return in a list
  list(
    breaks_absolute = c(
      -Inf,
      extract_parameter("Abs. shortfall moderate threshold"),
      extract_parameter("Abs. shortfall extreme threshold"),
      Inf
    ),
    breaks_proportional = c(
      -Inf,
      extract_parameter("Prop. shortfall moderate threshold"),
      extract_parameter("Prop. shortfall extreme threshold"),
      Inf
    ),
    qaly_weights = c(
      extract_parameter("Basic QALY weight"),
      extract_parameter("Moderate severity weight"),
      extract_parameter("Extreme severity weight")
    )
  )
}

#' Get the appropriate severity modifier according to the discounted QALYs for
#' those with the disease and those without the disease
#'
#' @param qalys_disease       The calculated (discounted) QALYs for people with
#'                            the disease receiving standard care.
#' @param qalys_nodisease     The calculated (discounted) QALYs for people who
#'                            do not have the disease, i.e., applying population
#'                            utility norms and general population mortality
#'                            rates.
#' @param .i                  Allows to inject `i`, otherwise it will be sourced
#'                            from the global environment.
#' @param .severity_modifier  Allows to inject the results of a call to
#'                            `severity_modifier` instead of it being called
#'                            every time.
#' @param hard_code_SM        hard-code the severity modifier using NICE methods
#'                            guidance, circa 2022
#' @param format              either `table` or `console`. Table produces a one-row table, console prints results to console
#'
get_severity_modifier <- function(qalys_disease, qalys_nodisease, .i = NULL, .severity_modifier = NULL, hard_code_SM = FALSE, format = "table") {
  stopifnot(format %in% c("table", "console"))
  if (is.null(.severity_modifier)) {
    .severity_modifier <- severity_modifier(.i = .i, hard_coded = hard_code_SM)
  }

  abs_shortfall <- qalys_nodisease - qalys_disease
  prop_shortfall <- 1 - qalys_disease / qalys_nodisease

  wt_abs <- .severity_modifier$qaly_weights[cut(abs_shortfall, .severity_modifier$breaks_absolute, labels = FALSE)]
  wt_prop <- .severity_modifier$qaly_weights[cut(prop_shortfall, .severity_modifier$breaks_proportional, labels = FALSE)]

  if (format == "table") {
    return(data.frame(
      abs_sf = abs_shortfall,
      prop_sf = prop_shortfall,
      modifer = max(wt_abs, wt_prop)
    ))
  } else {
    return(max(wt_abs, wt_prop))
  }
}

#' Calculate the Quality-Adjusted Life Expectancy (QALE) for a person, given
#' their age and sex
#'
#' @param age       Age of person
#' @param sex       Sex of person
#' @param .i        Inject `i` instead of looking in the global environment
#' @param .p        Inject `p` instead of looking in the global environment
#' @param .age_step Numerical integration step (choosing a value lower than the
#'                  default value of 0.2 does not change the results to 5
#'                  decimal places; if the function is running slowly then it
#'                  can be changed to 1 for a 5x speed-up and still be accurate
#'                  to 3 decimal places)
#' @param .use_time_horizon If this is `TRUE` then the time horizon given in
#'                  `p$basic$TH_yr` will be used (if this is less than the time
#'                  till the person reaches their 100th birthday). This is
#'                  probably not what we want to do, so the default value is
#'                  `FALSE`...
#'
#' @details
#'    In addition to being (somewhat) sensitive to the choice of which life
#'    tables are used (e.g., 2017-2019 and 2018-2020 give different results
#'    because of the Covid-19 pandemic which affected mortality in 2020) and
#'    population utility norms, there are some calculation details which
#'    distinguish this from, e.g., the shortfall calculator at
#'    https://r4scharr.shinyapps.io/shortfall/
#'
#'    First, we calculate restricted QALE up to the 100th birthday, i.e., no
#'    further QALE is accrued after reaching 100th birthday.
#'
#'    Second, we discount continuously, so that the discount factor is a smooth
#'    (exponential) function of time rather than stepping down each year.
#'
#'    Note that the default value for `.age_step` is appropriate for integer
#'    `age` but if `age` is non-integer then smaller values of `.age_step` may
#'    be necessary for numerical convergence.
baseline_qale <- function(age, sex = c("female", "male"), .i = NULL, .p = NULL,
                          .age_step = 0.2, .use_time_horizon = FALSE) {
  sex <- match.arg(sex)

  if (is.null(.i)) .i <- get("i", envir = globalenv())
  if (is.null(.p)) .p <- get("p", envir = globalenv())

  max_age <- if (.use_time_horizon) min(100, age + .p$basic$th_y) else 100

  ages <- seq(age, max_age, by = .age_step)

  if (max(ages) < max_age) ages <- c(ages, max_age)
  steps <- ages[2:length(ages)] - ages[1:(length(ages) - 1)]

  utilities <- utility_genpop(ages, sex, .p)

  lifetables <- get_lifetables(.i)

  df <- data.frame(
    age = ages[1:(length(ages) - 1)],
    t0 = ages[1:(length(ages) - 1)] - ages[1],
    t1 = ages[2:length(ages)] - ages[1]
  )
  df$q <- approx(x = lifetables$x, y = lifetables[[paste0("q_", sex)]], xout = df$age, method = "constant")$y
  df$m <- -log(1 - df$q)
  df$s <- exp(-df$m * (df$t1 - df$t0))
  df$S0 <- c(1, cumprod(df$s[1:(length(df$s) - 1)]))
  df$r <- .p$basic$discQ
  df$u0 <- utilities[1:(length(ages) - 1)]
  df$u1 <- utilities[2:length(ages)]

  df$auc <- mapply(
    auc_step,
    t0 = df$t0,
    t1 = df$t1,
    S0 = df$S0,
    u0 = df$u0,
    u1 = df$u1,
    m = df$m,
    MoreArgs = list(r = log(1 + .p$basic$discQ))
  )

  sum(df$auc)
}


#' Function to estimate area under a curve assuming exponential line between points
auc_step <- function(t0, t1, S0, u0, u1, m, r) {
  # Calculate the area under the curve
  #   S0 * exp(-m * (t - t0)) * exp(-r * t) * u(t)
  # between t = t0 and t = t1, where
  #   u(t)  = a * t + b
  #   u(t0) = u0
  #   u(t1) = u1
  a <- (u1 - u0) / (t1 - t0)
  b <- u0 - a * t0
  emrt0 <- exp(-(m + r) * t0)
  emrt1 <- exp(-(m + r) * t1)
  S0 * exp(m * t0) * ((a / ((m + r)^2)) * (emrt0 * (1 + (m + r) * t0) - emrt1 * (1 + (m + r) * t1)) + (b / (m + r)) * (emrt0 - emrt1)
  )
}


#' Calculate the severity modifier (and other outputs)
#'
#' @param age            Either a vector of ages (if supplying individual
#'                       patient data) or the mean age for a cohort.
#' @param sex            Either a vector of sexes (which can be a logical
#'                       vector with TRUE for male and FALSE for female
#'                       or a character vector with "male" and "female",
#'                       or even "m" and "f"), or the proportion of the
#'                       cohort which is male.
#' @param qalys          The discounted QALYs for somebody with the disease
#'                       receiving standard care.
#' @param .patient_level Logical. If TRUE, will treat `age` and `sex` as
#'                       vectors with individual patient data which are
#'                       aligned. If FALSE, will assume that `age` gives
#'                       the mean age, and `sex` gives the proportion of
#'                       the cohort which are male. If not provided, the
#'                       function will attempt to infer whether cohort
#'                       or individual patient data has been provided.
#' @param .i             Allows to inject `i`, otherwise it will be sourced
#'                       from the global environment.
#' @param .p             Allows for injection of `p` instead of looking
#'                       for it in the global environment.
calc_severity_modifier <- function(age, sex, qalys, .patient_level = NULL, .i = NULL, .p = NULL, format = "table") {
  # Check that sex and age are conformable
  stopifnot(format %in% c("table", "console"))
  stopifnot(length(age) == length(sex))

  # If .patient_level is not specified, infer it
  if (is.null(.patient_level)) .patient_level <- (length(sex) > 1)

  if (.patient_level) {
    # We know data.table is already a dependency, so go ahead and use it...
    dt <- data.table::data.table(age = age, sex = sex)
    dt[, n := .N, keyby = .(sex, age)]
    dt[, qale := mapply(baseline_qale, age, sex, MoreArgs = list(.i = i, .p = p))]
    qale <- dt[, sum(n * qale) / sum(n)]
  } else {
    qale_f <- baseline_qale(age = age, sex = "female", .i = .i, .p = .p)
    qale_m <- baseline_qale(age, "male", .i = .i, .p = .p)
    qale <- sex * qale_m + (1 - sex) * qale_f
  }

  val <- get_severity_modifier(qalys_disease = qalys, qalys_nodisease = qale, hard_code_SM = TRUE, format = format)
  # val <- get_severity_modifier(qalys_disease = qalys, qalys_nodisease = qale, .i = .i)

  if (format == "table") {
    data.frame(
      qaly_soc = qalys,
      qaly_gpop = qale,
      abs_sf = val$abs_sf,
      prop_sf = val$prop_sf,
      modifier = val$modifer
    )
  } else {
    attr(val, "QALE") <- qale
    attr(val, "QALY") <- qalys
    class(val) <- "severity_modifier"

    return(val)
  }
}

format.severity_modifier <- function(val) {
  stringr::str_glue(
    "{as.vector(val)}\uD7 (",
    "QALE: {format(attr(val, 'QALE'), digits = 3, nsmall = 2)}; ",
    "shortfall: {format(attr(val, 'QALE')-attr(val, 'QALY'), digits = 3, nsmall = 2)} [absolute], ",
    "{format(1 - attr(val, 'QALY') / attr(val, 'QALE'), digits = 3, nsmall = 4)} [proportional]",
    ")"
  )
}

print.severity_modifier <- function(val) {
  cat(format.severity_modifier(val), "\n")
}
