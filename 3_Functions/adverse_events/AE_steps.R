# adverse event coding

#' Calculate expected cost per patient per week cost and QALY penalty by treatment
#' 
#' @param AE_costs       Should be i$R_table_AE_costs
#' @param AE_disutil     Should be i$R_table_AE_util
#' @param AE_duration    Should be i$R_table_duration
#' @param AE_rate        Should be i$R_table_AE_rates
#' @param comparators    Should be i$lookup$trt
#' @param weeks_per_year Number of weeks per year, e.g., calculated as p$basic$cl_w / p$basic$cl_y
#' @param PSA            Sample from PSA distribution (TRUE/FALSE)
#' 
#' @details This function takes the individual inputs for adverse events and processes
#'          them into a format that fits in with the way information is organised
#'          throughout the rest of the model structure. This uses `data.table` syntax,
#'          so some practice/experience with that package is recommended before
#'          conducting a review. 
#'          
#'          To test this function: The recommended way to check functionality would
#'          be to load an input excel file within `Model_structure.R`, running that
#'          code up to the point at which `f_process_adverse_events` is called. 
#'          One can then use the values from `i` and `p` which populate the
#'          arguments to the function, assigning each one in order to then
#'          go line by line through the below code. This will enable the user
#'          to replicate the calculations from within this function side-by-side
#'          in excel to check that they are correct.
#' 
#' 
f_process_adverse_events <- function(AE_costs, AE_disutil, AE_duration, AE_rate, comparators, weeks_per_year, PSA) {
  
  AE_rate <- data.table(AE_rate)
  AE_rate[, duration_weeks := NULL] # Delete as not used in this function and risk of confusion
  
  if (PSA == FALSE) {
    AE_cost_QALYs_per_episode <- merge(
      merge(
        AE_duration[c("AE", "Duration_Mean")],
        AE_disutil[c("AE", "Disutil_mean")],
        all = TRUE,
        by  = "AE"),
      AE_costs[c("AE", "Cost.per.event_Mean")],
      all = TRUE,
      by  = "AE"
    )
    colnames(AE_cost_QALYs_per_episode) <- c("AE", "duration", "utility", "cost")
  
    # AE durations are given in weeks, but QALYs are given in years
    AE_cost_QALYs_per_episode$QALYs <- with(AE_cost_QALYs_per_episode, utility * duration / weeks_per_year)
    
    # Estimate AE rates for each line of each treatment
    all_combinations <- CJ(
      Treatment.name = names(comparators),
      Treatment.line = 1:4,
      AE             = unique(AE_costs$AE)
    )
    all_combinations[, Molecule := comparators[Treatment.name]]
    AE_rates_expanded <- AE_rate[
      all_combinations,
      on = .(AE = AE, Molecule = Molecule, Treatment.name = Treatment.name, Treatment.line = Treatment.line)
    ]
    setorder(AE_rates_expanded, Molecule, AE, Treatment.line)
    
    AE_rates_expanded[, Rate.per.patient.per.week := nafill(Rate.per.patient.per.week, "locf"), by = .(Molecule, AE)]
    AE_rates_expanded[, Rate.per.patient.per.week := nafill(Rate.per.patient.per.week, "nocb"), by = .(Molecule, AE)]
    AE_rates_expanded[, Rate.per.patient.per.week := nafill(Rate.per.patient.per.week, "const", fill = 0.0), by = .(Molecule, AE)]
    
    # Join the costs and QALYs per episode
    AE_cost_QALYs_by_treatment_and_AE <- AE_rates_expanded[
      AE_cost_QALYs_per_episode[c("AE", "cost", "QALYs")],
      on = .(AE = AE)
    ]
    
    # Weight by rate and sum for the result
    result <- AE_cost_QALYs_by_treatment_and_AE[
      ,
      .(cost = sum(Rate.per.patient.per.week * cost), QALYs = sum(Rate.per.patient.per.week * QALYs)),
      by = .(Treatment.name, Molecule, Treatment.line)
    ]
    setorder(result, Treatment.line, Molecule)
    result[, Molecule := NULL]
    setnames(result, c("trt", "line", "cost", "QALYs"))
    setDF(result)
    
    # Print a helpful message and return
    cat("Mean cost and QALY impact due to AEs per patient per week:\n")
    # print(result)
    return(result)
    
  } else {
    #code to sample iterations times from distributions of duration, cost and utility decrements
    stop("PSA sampling for AE costs and QALYs not implemented")
  }
  
}