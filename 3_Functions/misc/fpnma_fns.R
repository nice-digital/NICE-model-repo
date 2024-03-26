# code to manipulate FPNMA inputs, generate HRs by week from FPNMA coefficients etc

#' generate HR(t) from coefficients of a fractional polynomial network meta analysis (FPNMA)
#' 
#' @param TH time horizon in cycles. Note that in CE models time usually starts from 0.
#' @param wks_per_month this should be `i$i_mon_to_weeks` from the Excel model. defaults to 4.348214
#' 
#' @details `lnHR = intercept + s1*time^(power1) + s2*time^(power2), where power=0 means ln(time)` is
#' the specification of the model. rearranging that for HR gives the return in the function.
#' 
#' NOTE THAT THIS ASSUMES THE TIME UNIT OF THE ANALYSIS IS MONTHS AND THE 
#' COST EFFECTIVENESS MODEL IS IN WEEKS!
#' 
#' 
f_HRs_by_cycle <- function(TH, wks_per_month = 4.348214, exponents, coeffs) {
  exponents <- as.numeric(unlist(exponents))
  coeffs    <- as.numeric(unlist(coeffs))
  t         <- 0:TH
  t_month   <- t / wks_per_month
  return(exp(coeffs[1] + coeffs[2] * t_month ^ exponents[1] + coeffs[3] * t_month ^ exponents[2]))
}



#' Function to tidy up the input data for the FPNMA specifically for the PATT RCC model
#' 
#' @details this function should not be used anywhere else as it is highly specific
#' to this decision problem, including implicit assumptions!
#' 
f_FPNMA_tidy_and_add_exponents <- function(PSAcoefficients, exponents) {
  
  exponents$Population <- exponents$Endpoint <- NA
  exponents$Endpoint[exponents$Outcome == "OS"] <- 0
  exponents$Endpoint[exponents$Outcome == "PFS"] <- 1
  exponents$Population[exponents$Risk.group == "All"] <- 0
  exponents$Population[exponents$Risk.group == "Intermediate/poor"] <- 1
  exponents$Line[exponents$Line == "1L"] <- 1
  exponents$Line[exponents$Line == "2L+"] <- 2
  exponents$Line <- as.numeric(exponents$Line)
  
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "molecule"] <- "Molecule"
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "risk"] <- "Population"
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "endpoint"] <- "Endpoint"
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "line"] <- "Line"
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "referencetreatment"] <- "Reference.treatment"
  colnames(PSAcoefficients)[colnames(PSAcoefficients) == "referencetrial"] <- "Reference.trial"
  
  PSAcoefficients <- merge(PSAcoefficients, exponents, all.x = TRUE)
  PSAcoefficients <- PSAcoefficients[order(PSAcoefficients$run, 
                                           PSAcoefficients$Population,
                                           PSAcoefficients$Line,
                                           PSAcoefficients$Endpoint),]
  data.table(PSAcoefficients)
}

#' Function to generate FPNMA coda from coefficients. Highly specific to PATT RCC model.
f_generate_FPNMA_coda <- function(coeffs, TH, wks_per_month) {
  rbindlist(lapply(1:nrow(coeffs), function(row_in_table) {
    # Get data for this row:
    id <- as.data.table(lapply(coeffs[row_in_table,list(Population,Line,Molecule,Endpoint,Reference.treatment,Reference.trial)], function(x) rep(x,TH+1)))
    id$time <- 0:TH
    id$HR <- f_HRs_by_cycle(
      TH = TH, 
      wks_per_month = wks_per_month, 
      exponents = coeffs[row_in_table, list(Exponent.1, Exponent.2)], 
      coeffs = coeffs[row_in_table, list(intd, s1, s2)]
    )
    # correct for non-finite HRs. These are usually in time 0 where the HR is 1
    # in all contexts as no events are possible yet.
    id$HR[!is.finite(id$HR)] <- 1
    return(id)
  }))
}


#' Quick function for FPNMA: rebase the FPNMA HRs to use cabo in second line
#' 
#' @param FPNMAdata FPNMA data. Deterministic means or extrapolated "PSA" versions (one of the posterior draws / CODA sample)
#' 
#' @details rebases the hazard ratios to be in reference to molecule 8 not molecule
#' 10 in second-line. This is done for practical reasons (the reference curve available
#' and used in the base-case is molecule 8 for 2nd line therapy)
#' 
f_rebase_for_cabo_as_ref_in_2L <- function(FPNMAdata) {
  
  # Rebasing to allow use of cabo as reference treatment in 2nd line
  
  FPNMArebasecabo <- FPNMAdata[Line==2,]
  EvevscaboHROS <- FPNMArebasecabo[Molecule == 8 & Endpoint == 0]
  EvevscaboHRPFS <- FPNMArebasecabo[Molecule == 8 & Endpoint == 1]
  FPNMArebasecabo[Endpoint == 0]$HR <- FPNMArebasecabo[Endpoint == 0]$HR / EvevscaboHROS$HR
  FPNMArebasecabo[Endpoint == 1]$HR <- FPNMArebasecabo[Endpoint == 1]$HR / EvevscaboHRPFS$HR
  FPNMArebasecabo$Reference.treatment <- 8
  
  EveOSHRs <-
    cbind(
      time = seq(2, 2080),
      Molecule = rep(10, 2079),
      Reference.treatment = rep(8, 2079),
      Line = rep(2, 2079),
      Endpoint = rep(0, 2079),
      Population = rep(0, 2079),
      Reference.trial = rep(1, 2079),
      HR = EvevscaboHROS$HR
    )
  EvePFSHRs <-
    cbind(
      time = seq(2, 2080),
      Molecule = rep(10, 2079),
      Reference.treatment = rep(8, 2079),
      Line = rep(2, 2079),
      Endpoint = rep(1, 2079),
      Population = rep(0, 2079),
      Reference.trial = rep(1, 2079),
      HR = EvevscaboHRPFS$HR
    )
  
  FPNMArebasecabo <- rbind(FPNMArebasecabo, EveOSHRs, EvePFSHRs)
  
  FPNMAdata <- rbind(FPNMAdata, FPNMArebasecabo)

  FPNMAdata
}

#' Quick function for FPNMA: assume that HRs for 3L are the same as for 2L  
#' 
#' @param FPNMAdata FPNMA data. Deterministic means or extrapolated "PSA" versions (one of the posterior draws / CODA sample)
#' 
#' @details Simply applys HRs for 2L to 3L so that the relative efficacy network 
#' reflects this assumption duing compilation and propagation.
#' 
f_3L_rel_effect_same_as_2L <- function(FPNMAdata) {   

  assume3L      <- FPNMAdata[Line==2,]
  assume3L$Line <- 3
  FPNMAdata     <- rbind(FPNMAdata,assume3L)

  assumeTTD <- FPNMAdata[Endpoint==1,]
  assumeTTD$Endpoint <- 2
  FPNMAdata <- rbind(FPNMAdata,assumeTTD)

  assumeTTP <- FPNMAdata[Endpoint==1,]
  assumeTTP$Endpoint <- 3
  FPNMAdata <- rbind(FPNMAdata,assumeTTP)

  FPNMAdata
}


#' generate "destinations" object, which directs the later functions in terms of
#' where to place the FPNMA hazard ratios.
#' 
#' @details this function simply takes the first row for each grouping and then
#' drops the hazard ratio part, returning just information on "destinations". These
#' "destinations" then are fed into later functions like `f_NMA_linkFPNMA` and the
#' information therein is used to place the extrapolated hazard ratios from the
#' FPNMA to the right population, line, molecule, trial, endpoint (PLMTE) location.
#' The terms "origin" and "destination" are used instead of "intervention" and
#' "reference" because the reference curves can be at multiple degrees of separation away
#' from their ultimate destinations (e.g. curve --> FPNMA HRs --> assumed same as --> apply HR to).
#' 
#' IMPORTANT NOTE: note that a copy of the destinations for reference.trial of 2 
#' (real-world evidence) is made here. This is an assumption of the specific
#' adaptation for the RCC PATT model, which may not hold in other circumstances.
#' Consequently, this function is specific to the RCC model and not the generic
#' pathway
#' 
#' 
f_gen_destinations <- function(fp_data){
  
  destinations <- fp_data[, (head(.SD, 1)), by = list(Population,Line,Molecule,Endpoint,Reference.treatment,Reference.trial)]
  destinations <- destinations[,list(Population, Line, Molecule, Endpoint, Reference.treatment, Reference.trial)]
  
  destinations_temp <- destinations
  destinations_temp$Reference.trial <- 2
  destinations <- rbind(destinations,destinations_temp)
  
  return(destinations)
}


#' Function to add a copy of the FPNMA data which sets the reference.trial to 2,
#' which then assumes that the RWE HRs are the same as the trial-based
f_add_reference_trial_2 <- function(fp_data){
  means_temp <- fp_data
  means_temp$Reference.trial <- 2
  output <- rbind(fp_data,means_temp)
  output
}
