# utilities processing into QALYs

#' uses named range `R_table_util` from Excel to generate values applied to states
#' in the models. Populates both ST and PS models as these are simple HSUVs.
#' 
#' @param raw_utilities named range `R_table_util` from Excel front end file
#' @param PSA Flag for generating samples for probabilistic setting
#' @param samples Number of sample if `PSA==TRUE`
#' 
#' @details Simple function to tidy up the table from Excel and produce a format 
#'          which can then be further processed for the ST or PS models.
#' 
#' 
f_process_utilities <- function(raw_utilities,
                                PSA = FALSE,
                                samples = NULL) {
  
  if (PSA == TRUE & is.null(samples)) {
    stop("Warning: request to return PSA samples for utility data and number of samples not specified in call to f_process_utilities")
  }
  
  utilities <- raw_utilities[,1:9]
  colnames(utilities)[6:9] <- gsub("Mean","",colnames(utilities)[6:9])
  
  if (PSA == FALSE) {
    return(utilities)
  } else {
    sampled <- list()
    # cat("Generating",samples,"PSA utility samples")
    
    raw_utilities <- data.table(raw_utilities)
    
    out_tab <- raw_utilities[,list(Population, Population.name, Treatment.line, Molecule, Treatment.name)]
    
    # beta parameters for all rows, all with identifiers with Mean taken out:
    u_beta <- lapply(estBetaParams(raw_utilities[, 6:9], raw_utilities[, 10:13] ^ 2), function(param) {
      out <- cbind(out_tab,param)
      colnames(out) <- gsub("Mean","",colnames(out))
      out
    })
    
    # Cycle down the rows of this generating all PSA iterations all in one go and
    # producing one big table with iteration id
    s <- 1:samples
    return(rbindlist(lapply(1:nrow(u_beta$alpha), function(param_row) {
      
      id <- u_beta$alpha[param_row,list(Population, Population.name, Treatment.line, Molecule, Treatment.name)]
      
      id <- as.data.table(lapply(id, function(x) rep(x,samples)))
      
      id$iteration <- s
      
      # Get data for this HSUV
      dat <- list(
        alpha = u_beta$alpha[param_row,list(OnTxt, OffTxt, PFS, PD)],
        beta = u_beta$beta[param_row,list(OnTxt, OffTxt, PFS, PD)]
      )
      
      # Generate common random numbers
      random_gen <- runif(samples)
      
      hsuv_table <- as.data.table(lapply(structure(names(dat$alpha),.Names=names(dat$alpha)), function(state) {
        qbeta(
          random_gen,
          shape1 = .subset2(dat$alpha,state),
          shape2 = .subset2(dat$beta,state)
        )
      }))
      
      # stick the id table and the HSUV table together to produce the result - samples PSA iterations of all states for this category
      return(data.table(id,hsuv_table))
    })))
  }
}



#' NOT USED - an early attempt at applying utiltiies (including AE utilities)
#' to the trace produced by the ST model
#' 
f_apply_utilities_to_trace <- function(utilities,
                                       markov_trace,
                                       apply_AE = c("one-off", "per cycle"),
                                       AE_rate,
                                       AE_QALY_penalty,
                                       pop, 
                                       seq, 
                                       cycle_length,
                                       .p,
                                       age_adjust,
                                       starting_age,
                                       verbose) {
  
  apply_AE <- match.arg(apply_AE)
  age_adjust <- (age_adjust == "Yes")
  
  QALY_matrix <- matrix(0, nrow = nrow(markov_trace), ncol = ncol(markov_trace))
  colnames(QALY_matrix) <- colnames(markov_trace)
  for (tx_line in seq_along(seq)) {
    mol <- seq[tx_line]
    
    if(verbose) cat("Calculating QALY matrix for population", pop, "Line", tx_line, "= molecule", mol,"\n")
    
    #calculate QALYs accrued per cycle for on and off tx for this line and molecule
    QALY_on <- utilities$OnTxt[utilities$Population == pop &
                               utilities$Treatment.line == tx_line &
                               utilities$Molecule == mol] * cycle_length
    
    if(length(QALY_on) != 1) {
      stop("multiple or no entries found in R_table_util for population ", pop, 
           " Line ", tx_line, " Molecule ", mol," on treatment.\nSuggest check excel inputs file, utilities sheet.\n")
    }
    
    #add in AE QALY penalty
    if (apply_AE != "one-off") {
      if (mol != 999) {
        QALY_on <- QALY_on +
          AE_QALY_penalty$QALYs[
            AE_QALY_penalty$trt == i$List_comparators[mol+1]
            & AE_QALY_penalty$line == tx_line
          ]
      }
      #no AE decrement is applied here either for BSC or when one off-setting is chosen
    }
    
    
    # don't calculate off treatment QALYs if treatment is BSC (as there is no 'off tx' for BSC)
    if (mol != 999) {
      QALY_off <- utilities$OffTx[utilities$Population == pop &
                                  utilities$Treatment.line == tx_line &
                                  utilities$Molecule == mol] * cycle_length
      
      
      if(length(QALY_off) != 1) {
        stop("multiple or no entries found in R_table_util for population ", pop, 
             " Line ", tx_line, " Molecule ", mol," off treatment.\nSuggest check excel inputs file, utilities sheet.\n")
      }
    }    
    
    if (age_adjust) {
      #age adjust QALY_on and QALY_off here - vector of th length
      QALY_on <- rep(QALY_on, nrow(markov_trace))
      QALY_on <- adjust_utility(age = starting_age, sex = 0.5, 
                                utilities = data.frame(cycle = 1:nrow(markov_trace), QALY_on = QALY_on), 
                                .p = .p)
      
      if (mol != 999) {
        QALY_off <- rep(QALY_off, nrow(markov_trace))
        QALY_off <- adjust_utility(
          age = starting_age,
          sex = 0.5, 
          utilities = data.frame(cycle = 1:nrow(markov_trace), QALY_off = QALY_off), 
          .p = .p
        )
      }
    }
    
    #add to QALYs_matrix
    #note this works whether QALY_on is a single value or vector of length nrow(markov_trace)
    QALY_matrix[, 1 + (2 * (tx_line - 1))] <- markov_trace[, 1 + (2 * (tx_line - 1))] * QALY_on
    if (mol != 999) {
      QALY_matrix[, 2 + (2 * (tx_line - 1))] <- markov_trace[, 2 + (2 * (tx_line - 1))] * QALY_off
    }

    
    if (apply_AE == "one-off") {
      if (mol != 999) {
        AE_qaly_weight <- AE_QALY_penalty$QALYs[AE_QALY_penalty$trt == i$List_comparators[mol+1] &
                                                AE_QALY_penalty$line == tx_line]
        AE_qaly_duration <- mean(
          i$R_table_AE_rates$duration_weeks[
            i$R_table_AE_rates$Molecule == mol
            & i$R_table_AE_rates$Treatment.line == if (tx_line>2) 2 else tx_line
          ]
        )
        AE_oneoff_penalty <- AE_qaly_weight * AE_qaly_duration
        
         
        if (tx_line == 1) {
          QALY_matrix[1, 1 + (2 * (tx_line - 1))] <- QALY_matrix[1, 1 + (2 * (tx_line - 1))] + AE_oneoff_penalty
        } else {
          QALY_matrix[, 1 + (2 * (tx_line - 1))] <-  c(0,diff(QALY_matrix[, 1 + (2 * (tx_line - 1))])) * (1 + AE_oneoff_penalty)  
          #### work in progress needs to link to number entering the state
        }
      }
    }
    
  }
 
  return (QALY_matrix) 
}
