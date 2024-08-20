#' Functions to adjust for general population mortality
#'

#' Convenience function to pick the right table out of what was loaded from
#' Excel and rename the columns
get_lifetables <- function(i) {
  stopifnot("R_table_mort_lifeTable" %in% names(i))
  
  res <- i$R_table_mort_lifeTable
  
  names(res) <- c("x", "q_male", "q_female")
  
  res <- rbind(res, data.frame(x = 101, q_male = 1, q_female = 1))
  
  res
}

#' Fast function to get gpop OS line, optionally returns full matrix for QC
#' or just vector of age and sex adjusted OS line for gpop mortality.
f_surv_getgpopOSFast <- function(bl_mal, bl_age, t_yr, lt, full_output = FALSE) {
  
  # The cycle length is always the 2nd element of t_yr, as its always 1 time step :)
  cl_yr <- t_yr[2]
  
  # Make vector of age
  age_t <- pmin(floor(bl_age + t_yr),100)
  
  # Get vector of qx given age for m and f
  lt_aligned <- lt[match(age_t,age..x.),]
  
  # Calculate hazard rates (time units are still years)
  lt_aligned[, `:=`(h_ma = -log(1-qx..males.), h_fe = -log(1-qx..females.))]
  
  # Calculate per-cycle mortality probabilities
  lt_aligned[, `:=`(qt_ma = 1 - (1-qx..males.)^cl_yr, qt_fe = 1 - (1-qx..females.)^cl_yr)]
  
  # Calculate overall survival
  lt_aligned[, `:=`(
    os_ma = cumprod(exp(-lag(h_ma, 1, 0)*cl_yr)),
    os_fe = cumprod(exp(-lag(h_fe, 1, 0)*cl_yr))
  )]
  
  # Compute age and sex adjusted OS line
  lt_aligned[, s_t := bl_mal * os_ma + (1 - bl_mal) * os_fe]
  
  # calculate balance of male to female for weighting
  lt_aligned[, w_ma := nafill((os_ma * bl_mal) / s_t, "locf")]
  
  # Calculate weighted hazard
  lt_aligned[, h_w := (h_ma * w_ma) + (h_fe * (1-w_ma))]
  
  
  # Return matrix
  
  if (full_output) {
    return(
      as.matrix(
        lt_aligned[, .(age=age..x.,qx_m_orig=qx..males.,qx_f_orig=qx..females.,qx_m_target=qt_ma,qx_f_target=qt_fe,os_ma=os_ma,os_fe=os_fe,w_ma=w_ma,h_w=h_w,gpop_os=s_t)]
      )
    )
  } else {
    return(lt_aligned$s_t)
  }
  
}

# Testing for it:
# f_surv_getgpopOSFast(
#   bl_mal = 0.5,
#   bl_age = 50,
#   t_yr = p$basic$t_yr,
#   lt = data.table(i$R_table_mort_lifeTable),
#   full_output = FALSE
# )

#' Calculate general population survival curves for IPD data on age/sex/line
#' 
#' @param R_table_patientagesex IPD including columns Sex ("M"/"F"), Age (numeric) and Line (1-4)
#' @param R_table_mort_lifeTable the named range R_table_mort_lifeTable from the excel input sheet
#' @param t_yr vector of time in years at each cycle, stored in p$basic$t_yr
#' @param lookups list of lookups
f_surv_GenOSLines_ipd <- function(R_table_patientagesex, R_table_mort_lifeTable, t_yr, lookups) {
  
  # Calculate a matrix with survival curves across columns where every row
  # corresponds to a row in R_table_patientagesex
  all_genpop_os <- t(mapply(
    FUN      = f_surv_getgpopOSFast,
    bl_mal   = 1 * (R_table_patientagesex$Gender == "M"),
    bl_age   = R_table_patientagesex$Age,
    MoreArgs = list(t_yr = p$basic$t_yr, lt = data.table(R_table_mort_lifeTable))
  ))
  
  # Calculate the mean survival in each cycle within each line
  lines_genpop_os <- lapply(unique(R_table_patientagesex$Line), function(li)
      colMeans(all_genpop_os[R_table_patientagesex$Line == li,]))
  
  # Prepare the data structure to match what is produced by f_surv_GenOSLines_det
  pop_ref          <- paste0("pop_", lookups$ipd$pop$Number)
  names(pop_ref)   <- pop_ref
  line_nums        <- unique(R_table_patientagesex$Line)
  names(line_nums) <- paste0("line_", line_nums)

  # Create the output  
  lapply(pop_ref, function(popu) {
    lapply(line_nums, function(lin) {
      list(
        d  = list(pop = popu, line = paste0("line_", lin)),
        os = lines_genpop_os[[lin]]
      )
    })
  })
}

#' Deterministic gen pop OS line - one per row of the table in Excel (R_table_ptchar)
#' 
#' @param R_table_ptchar the named range R_table_ptchar from the excel input sheet
#' @param R_table_mort_lifeTable the named range R_table_mort_lifeTable from the excel input sheet
#' @param t_yr vector of time in years at each cycle, stored in p$basic$t_yr
#' @param lookups list of lookups
#' 
f_surv_GenOSLines_det <- function(R_table_ptchar, R_table_mort_lifeTable, t_yr, lookups) {
  
  n_pops <- nrow(R_table_ptchar)
  
  pop_ref <- paste0("pop_",lookups$ipd$pop$Number[match(unique(R_table_ptchar$Population),lookups$ipd$pop$Description)])
  names(pop_ref) <- pop_ref
  
  line_ref <- paste0("line_",unique(R_table_ptchar$Treatment.line))
  names(line_ref) <- line_ref
  
  # Make an empty structure using the labelling above (consistent with
  # labelling from excel and R model)
  empty_structure <- lapply(pop_ref, function(popu) {
    lapply(line_ref, function(lin) {
      list(
        d = list(),
        os = NULL
      )
    })
  })
  
  # Generate flat list containing our OS lines and some information on destination
  os_list <- lapply(1:n_pops, function(ptchar_row) {
    
    # Note that in 1st line it needs to divide up into 6 different populations
    # for prior IO and no prior IO
    
    dat <- as.list(R_table_ptchar[ptchar_row, ])
    
    return(list(
      d = list(
        pop  = paste0("pop_", lookups$ipd$pop[match(dat$Population, Description), ]$Number),
        line = paste0("line_", dat$Treatment.line)
      ),
      os = f_surv_getgpopOSFast(
        bl_mal      = 1 - dat$Starting...female.Mean,
        bl_age      = dat$Starting.age..years..Mean,
        t_yr        = t_yr,
        lt          = data.table(R_table_mort_lifeTable),
        full_output = FALSE
      )
    ))
    
  })
  
  # Dynamically slot the stuff from the flat list into the structured list:
  return(Reduce(
    x = 1:length(os_list),
    init = empty_structure,
    accumulate = FALSE,
    f = function(prev, excel_row) {
      
      dat                     <- os_list[[excel_row]]
      d                       <- dat$d
      pl                      <- prev[[d$pop]][[d$line]]
      pl$d                    <- d
      pl$os                   <- dat$os
      prev[[d$pop]][[d$line]] <- pl
      
      return(prev)
    }
  ))
}






#' Adjust a single survival curve for general population mortality
#' 
#' This function will ensure that the rate of mortality does not drop below the
#' rate of general population mortality.
#' 
#' @param base_age      Age when t = 0 (at model start)
#' @param cycle_length  Length of cycle (in years, e.g., 1/52)
#' @param v_t_os        Vector of times for overall survival (t = 0 at model
#'                      start)
#' @param v_p_os        Vector of overall survival probability, i.e.,
#'                      `v_p_os[i]` gives the probability of surviving to at
#'                      least `v_t_os[i]`
#' @param v_x_lifetable Vector of ages for life table
#' @param v_q_lifetable Vector of life table per-year mortality risk values,
#'                      i.e., `v_q_lifetable[i]` gives the probability that
#'                      somebody who reaches their `v_x_lifetable[i]`th birthday
#'                      will die before reaching their next birthday
#' @param .warn         If TRUE (default) will issue a warning detailing (if
#'                      applicable) when the life table gives a lower risk of
#'                      mortality than the overall survival curve
#'
#' @return A `data.frame` containing five columns: `t` (which will equal
#'         `v_t_os`), `q_genmort` (the per-cycle *probability* of general
#'         population mortality), `q_adjusted` (the per-cycle *probability* of
#'         death following adjustment to ensure the rate of mortality does not
#'         fall below the general population mortality), `s_genmort` (the
#'         *survival curve* if only general population mortality applied), and
#'         `s_adjusted` (the *survival curve* following adjustment)
adjust_single_survival <- function(
  base_age,
  cycle_length,
  v_t_os,
  v_p_os,
  v_x_lifetable,
  v_q_lifetable,
  .warn = TRUE
) {
  
  # Transform v_q_lifetable to give hazard rates
  v_r_lifetable <- -log(1 - v_q_lifetable)
  
  # Transform v_p_os to also give hazard rates (note that v_r_os will be one
  # element shorter than v_p_os)
  n <- length(v_t_os)
  v_r_os <- (log(v_p_os[1:(n-1)]) - log(v_p_os[2:n])) / cycle_length
  
  # Line v_r_lifetable up with v_r_os
  v_x_aligned <- v_t_os[1:(n-1)] + base_age
  v_r_lifetable_aligned <- approx(
    x      = v_x_lifetable,
    y      = v_r_lifetable,
    xout   = v_x_aligned,
    method = "constant"
  )$y
  
  # Check if v_r_os goes below v_r_lifetable_aligned at any point
  if (.warn & any(v_r_os < v_r_lifetable_aligned)) {
    warning_msg <- paste(
      "Mortality rate from life table exceeds extrapolated mortality at time",
      v_t_os[which.max(v_r_os < v_r_lifetable_aligned)]
    )
    warning(warning_msg)
  }
  
  # Stitch together to get combined hazard rate vector
  v_r_combined <- pmax(v_r_os, v_r_lifetable_aligned)
  
  # Convert to per-cycle probabilities of avoiding death
  v_p_genmort  <- exp(-v_r_lifetable_aligned * cycle_length)
  v_p_combined <- exp(-v_r_combined * cycle_length)
  
  # Convert to survival
  v_p_os_genmort <- Reduce(
    f          = `*`,
    x          = v_p_genmort,
    init       = 1,
    accumulate = TRUE
  )
  v_p_os_adjusted <- Reduce(
    f          = `*`,
    x          = v_p_combined,
    init       = 1,
    accumulate = TRUE
  )
  
  # Return
  data.frame(
    t          = v_t_os,
    q_genmort  = c(1 - v_p_genmort, NA_real_),
    q_adjusted = c(1 - v_p_combined, NA_real_),
    s_genmort  = v_p_os_genmort,
    s_adjusted = v_p_os_adjusted
  )

}

#' Adjust a survival curve using patient-level data on age and sex
#' 
#' @param pts  A `data.frame` with two columns, `age` (numeric) and `sex` (see
#'             Details). Has one row for each patient.
#' @param s_os The survival curve that needs adjusting (calculated for each
#'             cycle in the model)
#' @param .i   A list containing inputs which have been loaded from the Excel
#'             inputs file. If this is not provided then the function will look
#'             for `i` in the global environment and stop with an error if it is
#'             not found.
#' @param .p   A list containing "cultivated" model parameters. If this is not
#'             provided then the function will look for `p` in the global
#'             environment and stop with an error if it is not found.
#' @param .warn         If TRUE (default) will issue a warning detailing (if
#'                      applicable) when the life table gives a lower risk of
#'                      mortality than the overall survival curve
#' 
#' @details The sex of participants can be specified in a number of ways.
#'          Option 1: Character with "M" or "m" for men and "F" or "f" for women.
#'          Option 2: Factor which coerces to characters consistent with Option 1.
#'          Option 3: Integer or Boolean vector with 0/FALSE representing female
#'          and 1/TRUE representing male.
#'
#' @return A `data.frame` containing seven columns: `t`, `q_genmort` (the per-
#'         cycle *probability* of general population mortality), `q_adjusted`
#'         (the per-cycle *probability* of death following adjustment to ensure
#'         the rate of mortality does not fall below the general population
#'         mortality), `s_genmort` (the *survival curve* if only general
#'         population mortality applied), `s_adjusted` (the *survival curve*
#'         following adjustment), `prop_male.genmort` (the proportion of the
#'         remaining cohort which would be male if only general population
#'         mortality applied), and `prop_male.adjusted` (the proportion of the
#'         remaining cohort which are male after OS is adjusted for general
#'         population mortality).
adjust_survival_individuals <- function(pts, s_os, .i = NULL, .p = NULL, .warn = TRUE) {
  i <- if (is.null(.i)) get("i", envir = globalenv()) else .i
  p <- if (is.null(.p)) get("p", envir = globalenv()) else .p
  
  age <- pts$age
  sex <- pts$sex
  
  lifetables <- get_lifetables(i)
  
  if (is.factor(sex)) sex <- as.character(sex)
  
  if (is.character(sex)) {
    sex <- toupper(sex)
    stopifnot(all(sex %in% c("M", "F")))
    sex <- (sex == "M")
  }
  
  if (is.numeric(sex)) stopifnot(all(sex %in% c(0, 1)))
  
  # Generate individual survival curves
  individual_curves <- mapply(
    FUN = function(age, sex, ...) {
      adjust_single_survival(
        base_age      = age,
        v_q_lifetable = if (sex) lifetables$q_male else lifetables$q_female,
        ...
      )
    },
    age = age,
    sex = sex,
    MoreArgs = list(
      cycle_length  = p$basic$cl_y,
      v_t_os        = p$basic$cl_y * (1:length(s_os) - 1),
      v_p_os        = s_os,
      v_x_lifetable = lifetables$x,
      .warn         = .warn
    ),
    SIMPLIFY = FALSE
  )
  
  # Combine the survival curves
  #
  # Details:
  # - Combined survival is simply the mean survival across individual survival
  #   curves
  # - Combined rate and per-cycle probability of death are calculated from the
  #   individual rates and probabilities after adjusting for survival
  t_combined <- individual_curves[[1]]$t
  s_adjusted_combined <- apply(
    X      = do.call(cbind, lapply(individual_curves, `[[`, "s_adjusted")),
    MARGIN = 1,
    FUN    = mean
  )
  s_genmort_combined <- apply(
    X      = do.call(cbind, lapply(individual_curves, `[[`, "s_genmort")),
    MARGIN = 1,
    FUN    = mean
  )
  q_genmort_combined <- apply(
    X      = do.call(cbind, lapply(individual_curves, function(x) x$q_genmort * x$s_adjusted)),
    MARGIN = 1,
    FUN    = mean
  ) / s_adjusted_combined
  q_adjusted_combined <- apply(
    X      = do.call(cbind, lapply(individual_curves, function(x) x$q_adjusted * x$s_adjusted)),
    MARGIN = 1,
    FUN    = mean
  ) / s_adjusted_combined
  
  prop_male_genmort <- apply(
    X      = mapply(function(.c, .s) .s * .c$s_genmort, .c = individual_curves, .s = sex),
    MARGIN = 1,
    FUN    = mean
  ) / s_genmort_combined
  prop_male_combined <- apply(
    X = mapply(function(.c, .s) .s * .c$s_adjusted, .c = individual_curves, .s = sex),
    MARGIN = 1,
    FUN = mean
  ) / s_adjusted_combined
  
  # Return
  data.frame(
    t                  = t_combined,
    q_genmort          = q_genmort_combined,
    q_adjusted         = q_adjusted_combined,
    s_genmort          = s_genmort_combined,
    s_adjusted         = s_adjusted_combined,
    prop_male.genmort  = prop_male_genmort,
    prop_male.adjusted = prop_male_combined
  )
  
}


#' Adjust a survival curve
#' 
#' @param sex            Either a number between 0 and 1 giving the proportion
#'                       of patients who are male, or a vector of patient sex in
#'                       the form specified in `adjust_survival_individuals`
#' @param age            Either a number giving the mean age of the population,
#'                       or a vector of individual patient ages which aligns
#'                       with `sex`
#' @param survivor       The survivor function, with one observation per 1-week
#'                       model cycle. This is expected to be either a 2D array/
#'                       matrix with two columns or a `data.frame` with two
#'                       columns.
#' @param .patient_level If `TRUE` will force the use of
#'                       `adjust_survival_individuals`. If not specified (or
#'                       specified as `NULL` the function will attempt to
#'                       determine whether patient level or aggregate data has
#'                       been supplied)
#' @param .i             A list containing inputs which have been loaded from
#'                       the Excel inputs file. If this is not provided then the
#'                       function will look for `i` in the global environment
#'                       and stop with an error if it is not found.
#' @param .p             A list containing "cultivated" model parameters. If
#'                       this is not provided then the function will look for
#'                       `p` in the global environment and stop with an error if
#'                       it is not found.
#' @param .warn          If TRUE (default) will issue a warning detailing (if
#'                       applicable) when the life table gives a lower risk of
#'                       mortality than the overall survival curve
#'
#' @return A `data.frame` containing seven columns: `t`, `q_genmort` (the per-
#'         cycle *probability* of general population mortality), `q_adjusted`
#'         (the per-cycle *probability* of death following adjustment to ensure
#'         the rate of mortality does not fall below the general population
#'         mortality), `s_genmort` (the *survival curve* if only general
#'         population mortality applied), `s_adjusted` (the *survival curve*
#'         following adjustment), `prop_male.genmort` (the proportion of the
#'         remaining cohort which would be male if only general population
#'         mortality applied), and `prop_male.adjusted` (the proportion of the
#'         remaining cohort which are male after OS is adjusted for general
#'         population mortality).
adjust_survival <- function(sex, age, survivor, .patient_level = NULL,
  .i = NULL, .p = NULL, .warn = TRUE) {
  
  # Check that sex and age are conformable
  stopifnot(length(sex) == length(age))
  
  # If .patient_level is not specified, infer it
  if (is.null(.patient_level)) .patient_level <- (length(sex) > 1)
  
  # Extract v_survivor (the vector of only survival probabilities)
  v_survivor <- if (is.list(survivor)) survivor[[2]] else survivor[,2]
  
  if (.patient_level) {
    
    return(
      adjust_survival_individuals(
        pts   = data.frame(sex = sex, age = age),
        s_os  = v_survivor,
        .i    = .i,
        .p    = .p,
        .warn = .warn
      )
    )
    
  } else {
    
    m <- adjust_survival_individuals(
      pts   = data.frame(sex = "M", age = age),
      s_os  = v_survivor,
      .i    = .i,
      .p    = .p,
      .warn = .warn
    )
    f <- adjust_survival_individuals(
      pts   = data.frame(sex = "F", age = age),
      s_os  = v_survivor,
      .i    = .i,
      .p    = .p,
      .warn = .warn
    )
    
    t_combined <- m$t
    
    s_genmort_combined  <- sex * m$s_genmort + (1 - sex) * f$s_genmort
    s_adjusted_combined <- sex * m$s_adjusted + (1 - sex) * f$s_adjusted
    
    q_genmort_combined <- (
      sex * m$q_genmort * m$s_adjusted + (1 - sex) * f$q_genmort * f$s_adjusted
    ) / s_adjusted_combined
    
    q_adjusted_combined <- (
      sex * m$q_adjusted * m$s_adjusted + (1 - sex) * f$q_adjusted * f$s_adjusted
    ) / s_adjusted_combined
    
    prop_male_genmort  <- (sex * m$s_genmort) / s_genmort_combined
    prop_male_combined <- (sex * m$s_adjusted) / s_adjusted_combined
    
    return(data.frame(
      t                  = t_combined,
      q_genmort          = q_genmort_combined,
      q_adjusted         = q_adjusted_combined,
      s_genmort          = s_genmort_combined,
      s_adjusted         = s_adjusted_combined,
      prop_male.genmort  = prop_male_genmort,
      prop_male.adjusted = prop_male_combined
    ))
    
  }
}



# testing ground ----------------------------------------------------------

if (FALSE) {
  # Make plot for output:
  gp <- get_elem(p$surv$gpop,"os")
  
  gp <- rbindlist(list(
    pop_0 = data.table(
      population = p$basic$lookup$ipd$pop[match(0,p$basic$lookup$ipd$pop$Number),]$Description,
      gpop_os = gp$pop_0$line_1,
      t = p$basic$t_yr
    ),
    pop_1 = data.table(
      population = p$basic$lookup$ipd$pop[match(1,p$basic$lookup$ipd$pop$Number),]$Description,
      gpop_os = gp$pop_1,
      t = p$basic$t_yr
    ),
    pop_2 = data.table(
      population = p$basic$lookup$ipd$pop[match(2,p$basic$lookup$ipd$pop$Number),]$Description,
      gpop_os = gp$pop_2,
      t = p$basic$t_yr
    )
  ))
  gpop_pop_plot <- ggplot(gp, aes(x = t, y = gpop_os, colour = population)) + 
    geom_line() + 
    theme_classic() +
    theme(legend.position = "bottom", legend.title=element_blank()) + 
    labs(title = NULL, x = "Time (years)", y = "% Survival") + 
    scale_x_continuous(expand = expansion(mult = c(0,0.05))) + 
    scale_y_continuous(labels = scales::percent)
  
  ggsave(
    filename = file.path("./4_Output","gpop_plot_pop_level.png"),
    plot     = gpop_pop_plot,
    device = "png",
    units = "cm",
    width = 15
  )
  
  
}

