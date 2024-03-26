#' Functions to assist with creating age-related utility values (QALY weights)
#' 
#' These are based on the model fits described (in limited detail) in Hernandez
#' Alava (2022) and more detailed model description in Hernandez Alava (2012).
#' 
#' The statistical model is a mixture model but with added floor and ceiling
#' effects.
#' 
#' REFERENCES
#' 
#' Hernandez Alava M, Pudney S, Wailoo A. Estimating EQ-5D by age and sex for
#'   the UK. NICE DSU Report. 2022.
#' Hernandez Alava M, Wailoo AJ, Ara R. Tails from the Peak District: Adjusted
#'   Limited Dependent Variable Mixture Models of EQ-5D Questionnaire Health
#'   State Utility Values. Value in Health 2012; 15(3):550-561.

#' Pull the HSE 2014 models from Excel
#' 
#' @param sex Choose from "female" and "male" to load the relevant model.
#' @param psa If `psa == TRUE` then instead of returning the maximum
#'            likelihood parameter values, it will instead sample
#'            from the MLE distribution.
#' @param .i  Provides ability to inject a list `.i` which will be used
#'            to obtain the tables instead of a full extract from Excel.
#'            If this is not provided the function will take `i` from
#'            the global environment.
extract_utility_ageadjust_coefs <- function(
    sex = c("female", "male"), psa = FALSE, .i = NULL
  ) {
  
  if (is.null(.i)) .i <- get("i", envir = globalenv())
  
  sex <- match.arg(sex)
  
  tbl <- .i[[paste0("R_table_ageadjust_util_", sex)]]
  
  v_covariate_nm <- c("Age/10", "(Age/10)^2", "intercept")
  v_mix_nm <- c("Age/10", "intercept")
  
  # Extract maximum likelihood estimates ----
  
  tbl_coef <- if (sex == "female") tbl$Females.Coefficient else tbl$Males.Coefficient
  
  l_coefficients <- list(
    # Component 1
    setNames(tbl_coef[2:4], v_covariate_nm),
    
    # Component 2
    setNames(tbl_coef[5:7], v_covariate_nm),
    
    # Component 3
    setNames(tbl_coef[8:10], v_covariate_nm)
  )
  
  l_mix_coef <- list(
    # Component 1
    setNames(tbl_coef[11:12], v_mix_nm),
    
    # Component 2
    setNames(tbl_coef[13:14], v_mix_nm)
  )
  
  v_sigma <- tbl_coef[15:17]
  
  # Sample from MLE distribution if PSA ----
  
  if (psa) {
    
    stopifnot(requireNamespace("MASS"))
    
    # Note that the MLE values of sigma have been transformed into
    # the linear scale but the covariance matrix is for ln(sigma)
    #
    # As a result, we set the mean of the multivariate normal
    # distribution to zero for these elements, and then we
    # use the results as scale factors for the sigma values.
    
    mu <- c(
      unlist(l_coefficients),
      unlist(l_mix_coef),
      rep(0, length.out = length(v_sigma))
    )
    
    Sigma <- unname(
      as.matrix(
        as.data.frame(
          lapply(tbl[2:17, 4:19], as.numeric)
        )
      )
    )
    
    x <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
    
    l_coefficients[[1]][1:3] <- x[1:3]
    l_coefficients[[2]][1:3] <- x[4:6]
    l_coefficients[[3]][1:3] <- x[7:9]
    
    l_mix_coef[[1]][1:2] <- x[10:11]
    l_mix_coef[[2]][1:2] <- x[12:13]
    
    v_sigma <- v_sigma * exp(x[14:16])
    
    
  }
  
  # Return ----
  
  list(
    l_coefficients = l_coefficients,
    l_mix_coef     = l_mix_coef,
    v_sigma        = v_sigma
  )
}

#' Add parameters for calculating EQ-5D population norms to a parameter list `p`
#' 
#' @param p   An existing parameter list to which parameters will be added
#' @param psa If `FALSE` then maximum likelihood estimates for parameters will
#'            be added, while if `TRUE` values will be sampled from the MLE
#'            asymptotic distribution
#' @param .i  Provides ability to inject a list `.i` which will be used
#'            to obtain the tables instead of a full extract from Excel.
#'            If this is not provided the function will take `i` from
#'            the global environment.
#' @returns   A *new* list (because R will do copy-on-modify) with the
#'            parameters added/replaced
add_population_utility_params <- function(p, psa = FALSE, .i = NULL) {
  p$util$pop_norms <- list(
    female = extract_utility_ageadjust_coefs("female", psa, .i),
    male = extract_utility_ageadjust_coefs("male", psa, .i)
  )
  return(p)
}

#' Calculate the utility value for a single component of the ALDVMM model
#'
#' @param covariates   Either a vector of patient covariates or a matrix (or
#'                     something coercible to a matrix) where each row gives the
#'                     covariates for a single patient (e.g. each row could
#'                     correspond to different ages)
#' @param coefficients A vector of coefficients which determine the mean of the
#'                     untruncated distribution. This must have the same length
#'                     as `covariates` if `covariates` is a vector or the same
#'                     length as the number of columns of `covariates`.
#' @param sigma        The standard deviation of the untruncated distribution.
#' @param upper_limit  The upper limit at which a ceiling effect is observed.
#' @param upper_to     What happens to values greater than `upper_limit`.
#' @param lower_limit  The lower limit at which a floor effect is observed.
#' @param lower_to     What happens to values less than `lower_limit`.
#' @param type         Either "mean", in which case the mean utility is
#'                     returned, or "sampled", in which case a sample or samples
#'                     from the distribution will be returned.
#' @param .n           The number of samples (if greater than 1) when
#'                     `type == "sampled"`. This will be the number of samples
#'                     per patient if there is more than one patient given in
#'                     `covariates`, in which case different samples will be
#'                     returned in different columns.
utility_single_component <- function(
  covariates,
  coefficients,
  sigma,
  upper_limit = 0.883,
  upper_to    = 1.0,
  lower_limit = -0.594,
  lower_to    = lower_limit,
  type        = c("mean", "sampled"),
  .n          = NULL
  ) {
  
  type <- match.arg(type)
  
  ym <- if (is.null(dim(covariates))) sum(covariates * coefficients) else rowSums(as.matrix(covariates) %*% coefficients)
  
  if (type == "mean") {
    
    # Calculate mass on `upper_to`
    w_upper <- pnorm(upper_limit, ym, sigma, lower.tail = FALSE)
    
    # Calculate mass on `lower_to`
    w_lower <- pnorm(lower_limit, ym, sigma, lower.tail = TRUE)
    
    # Remaining mass
    w_rem <- 1 - w_upper - w_lower
    
    # Calculate mean of remaining mass
    # Note - If the location parameter of the truncated distribution is outside
    #        the truncation limits, this can result in significant numerical
    #        issues. We can replace with `truncnorm` package if this occurs, but
    #        for now I want to avoid requiring too many packages.
    alpha  <- (lower_limit - ym) / sigma
    beta   <- (upper_limit - ym) / sigma
    mu_rem <- ym - sigma * (dnorm(beta) - dnorm(alpha)) / (pnorm(beta) - pnorm(alpha))
    
    # Combine all components
    return(w_upper * upper_to + w_lower * lower_to + w_rem * mu_rem)
    
  } else if (type == "sampled") {
    
    n_pt <- if (is.null(dim(covariates))) 1 else nrow(covariates)
    
    if (is.null(.n)) {
      ys <- rnorm(n_pt, ym, sd = sigma)
      
      if (ys > upper_limit) ys <- upper_to
      if (ys < lower_limit) ys <- lower_to
      
      return(ys)
    } else {
      ys <- array(rnorm(.n * n_pt, ym, sd = sigma), dim = c(n_pt, .n))
      
      ys[ys > upper_limit] <- upper_to
      ys[ys < lower_limit] <- lower_to
      
      return(ys)
    }
  }
  
}

#' Calculate the utility value for a full ALDVMM model
#'
#' @param covariates     A vector of patient covariates.
#' @param l_coefficients A list of vectors of coefficients which determine the
#'                       mean of the untruncated distributions for the
#'                       components. Each vector of coefficients must have the
#'                       same length as `covariates`.
#' @param v_sigma        The standard deviations of the untruncated
#'                       distributions.
#' @param l_mix_coef     A list of vectors of coefficients which determine the
#'                       mixture weights using a multinomial logit approach.
#'                       Note that `l_mix_coef` should have one less component
#'                       than `l_coefficients` and `sigma` because the final
#'                       component is assumed to have linear predictor = 0.
#' @param upper_limit    The upper limit at which a ceiling effect is observed.
#' @param upper_to       What happens to values greater than `upper_limit`.
#' @param lower_limit    The lower limit at which a floor effect is observed.
#' @param lower_to       What happens to values less than `lower_limit`.
#' @param type           Either "mean", in which case the mean utility is
#'                       returned, or "sampled", in which case a sample or
#'                       samples from the distribution will be returned.
#' @param .n             The number of samples (if greater than 1) when
#'                       `type == "sampled"`.
utility_mixture <- function(
  covariates,
  l_coefficients,
  v_sigma,
  l_mix_coef,
  upper_limit = 0.883,
  upper_to    = 1.0,
  lower_limit = -0.594,
  lower_to    = lower_limit,
  type        = c("mean", "sampled"),
  .n          = NULL
  ) {
  
  type <- match.arg(type)
  
  stopifnot(is.list(l_coefficients))
  stopifnot(length(l_coefficients) == length(v_sigma))
  stopifnot(length(l_coefficients) == length(l_mix_coef) + 1)
  
  n_c <- length(l_coefficients)
  
  if (is.null(dim(covariates))) {
    mix_eta <- sapply(l_mix_coef, function(mix_coef) sum(covariates * mix_coef))
    mix_wt <- c(exp(mix_eta), 1)
    mix_wt <- mix_wt / sum(mix_wt)
  } else {
    mix_eta <- sapply(l_mix_coef, function(mix_coef) rowSums(as.matrix(covariates) %*% mix_coef))
    mix_wt <- cbind(exp(mix_eta), 1)
    mix_wt <- sweep(mix_wt, 1, rowSums(mix_wt), "/")
  }
  
  if (type == "mean") {
    
    component_utility <- mapply(
      utility_single_component,
      l_coefficients,
      v_sigma,
      MoreArgs = list(
        covariates  = covariates,
        upper_limit = upper_limit,
        upper_to    = upper_to,
        lower_limit = lower_limit,
        lower_to    = lower_to,
        type        = type
      )
    )
    
    return(if (is.null(dim(covariates))) sum(mix_wt * component_utility) else rowSums(mix_wt * component_utility))
    
  } else if (type == "sampled") {
    
    if (is.null(.n)) {

      cls <- which(rmultinom(1, 1, mix_wt) == 1)
      
      return(
        utility_single_component(
          covariates,
          l_coefficients[[cls]],
          v_sigma[cls],
          upper_limit,
          upper_to,
          lower_limit,
          lower_to,
          "sampled"
        )
      )
    
    } else {
      
      cls <- rmultinom(1, .n, mix_wt)
      
      u_sampled <- mapply(
        utility_single_component,
        coefficients = l_coefficients,
        sigma = v_sigma,
        .n = cls,
        MoreArgs = list(
          covariates  = covariates,
          upper_limit = upper_limit,
          upper_to    = upper_to,
          lower_limit = lower_limit,
          lower_to    = lower_to,
          type        = type
        ),
        SIMPLIFY = FALSE
      )
      
      u_sampled <- unlist(u_sampled)
      
      return(sample(u_sampled, size = length(u_sampled)))
      
    }
    
  }
  
}

#' Calculate the utility for the general population for given age(s) and sex
#' 
#' @param age The age at which to evaluate the general population utility (the
#'            function is vectorised over `age`)
#' @param sex The sex for which to evaluate general population utility (the
#'            function is *NOT* vectorised over `sex`)
#' @param .p  Allows for injection of `p` instead of looking for it in the
#'            global environment. `p` is assumed to have the parameters for the
#'            female and male ALDVMMs.
utility_genpop <- function(age, sex = c("female", "male"), .p = NULL) {
  
  if (is.null(.p)) .p <- get("p", envir = globalenv())
  
  # Fetch the appropriate model from .p
  model <- .p$util$pop_norms[[sex]]
  
  # l_mix_coef doesn't use the (age/10)^2 but we need to put it in otherwise
  # we cannot use a single covar
  
  l_mix_coef <- lapply(model$l_mix_coef, function(mc) {
    c(mc[1], `(Age/10)^2` = 0, mc[2])
  })
  
  covar <- cbind(0.1 * age, 0.01 * age ^ 2, 1.0)
  
  utility_mixture(
    covariates     = covar,
    l_coefficients = model$l_coefficients,
    v_sigma        = model$v_sigma,
    l_mix_coef     = l_mix_coef
  )
  
}

#' Adjust utilities for population norms by age and sex
#' 
#' This function operates by applying population norms according to age
#' and sex, assuming that in the first cycle no adjustment is required
#' (i.e., it is already consistent with population norms) and making
#' subsequent adjustments to continue to reflect population norms.
#' 
#' @param age            Either a vector of ages (if supplying individual
#'                       patient data) or the mean age for a cohort.
#' @param sex            Either a vector of sexes (which can be a
#'                       character vector with "male" and "female"), or
#'                       the proportion of the cohort which is male.
#' @param utilities      The utilities to be adjusted: a list whose 2nd
#'                       is the unadjusted utility value for each cycle,
#'                       or a `data.frame` whose 2nd column is the un-
#'                       adjusted utility.
#' @param .patient_level Logical. If TRUE, will treat `age` and `sex` as
#'                       vectors with individual patient data which are
#'                       aligned. If FALSE, will assume that `age` gives
#'                       the mean age, and `sex` gives the proportion of
#'                       the cohort which are male. If not provided, the
#'                       function will attempt to infer whether cohort
#'                       or individual patient data has been provided.
#' @param .p             Allows for injection of `p` instead of looking
#'                       for it in the global environment.
adjust_utility <- function(age, sex, utilities, .patient_level = NULL, .p = NULL) {
  
  # Check that sex and age are conformable
  stopifnot(length(age) == length(sex))
  
  # If .patient_level is not specified, infer it
  if (is.null(.patient_level)) .patient_level <- (length(sex) > 1)
  
  # If .p was not injected, look for p in the global environment
  if (is.null(.p)) .p <- get("p", envir = globalenv())
  
  # Extract v_utilities
  v_utilities <- if (is.list(utilities)) utilities[[2]] else utilities[,2]
  
  if (.patient_level) {
    
    ages <- lapply(age, function(age) age + .p$basic$cl_y * (seq_along(v_utilities) - 1))
    
    a_genpop <- mapply(
      FUN = utility_genpop,
      age = ages,
      sex = sex,
      MoreArgs = list(.p = .p),
      SIMPLIFY = TRUE
    )
    
    v_genpop <- apply(a_genpop, MARGIN = 1, FUN = mean)
    
    return((v_genpop / v_genpop[1]) * v_utilities)
    
  } else {
    
    ages <- age + .p$basic$cl_y * (seq_along(v_utilities) - 1)
    
    v_genpop_m <- utility_genpop(ages, "male", .p = .p)
    v_genpop_f <- utility_genpop(ages, "female", .p = .p)
    v_genpop <- sex * v_genpop_m + (1 - sex) * v_genpop_f
    
    return((v_genpop / v_genpop[1]) * v_utilities)
    
  }
}
