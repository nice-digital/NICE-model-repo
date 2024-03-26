#' Function to "clean" the drug cost input table from excel. This function is
#' not currently used, but may be adapted in the future.
#' 
#' @param drug_cost_table the input table from excel for drug costs.
#' 
f_clean_drug <- function(drug_cost_table) {
  
  tab <- as.data.table(drug_cost_table)
  
  # Now cycle through the individual treatments using the id of trt for each one,
  # retaining their grouping but separating different combination regimen - 
  # 
  # NOTE THAT THIS MEANS THAT THE SAME MOLECULE WILL BE IN THIS TABLE MULTIPLE
  # TIMES, ASSOCIATED WITH DIFFERENT COMBINATIONS.
  # 
  # This actually minimises repetition and allows for different formulations of 
  # the same chemical (like flat-dose vs weight-based nivo)
  
  mols <- paste0("mol_",unique(tab$Molecule))
  
  
  trt_list <- lapply(unique(tab$Molecule), function(trt) {
    
    # Filter down to this "Molecule", i.e. regimen (or reg for short)
    reg <- tab[Molecule == trt,]
    
    # Record the number of components
    n_comp         <- nrow(reg[!is.na(reg$Drug.cost.per.dose),])
    comp_id        <- reg$Treatment.name[!is.na(reg$Drug.cost.per.dose)]
    names(comp_id) <- comp_id
    
    # Each component needs its own list entry so that different regimen can
    # have different numbers of combinates and the resulting data don't need
    # any special code for individual regimen. Means that future expansions
    # require no further work
    
    component_list <- lapply(comp_id, function(combinate) {
      list(
        type      = reg[Drug.name == combinate, Per.cycle.or.one.off.cost],
        cpd       = reg[Drug.name == combinate, Drug.cost.per.dose],
        acpd      = reg[Drug.name == combinate, Admin.cost.per.admin],
        acpd_se   = reg[Drug.name == combinate, Admin.cost.SE],
        dos_freq  = reg[Drug.name == combinate, Applied.every.x.cycles],
        dos_quant = reg[Drug.name == combinate, Number.of.doses.when.applied],
        rdi       = reg[Drug.name == combinate, RDI],
        rdi_se    = reg[Drug.name == combinate, RDI.SE],
        start     = reg[Drug.name == combinate, Time.from],
        stop      = reg[Drug.name == combinate, Time.to..cycles.]
      )
    })
    return(component_list)
    
  })
  names(trt_list) <- mols
  return(trt_list)
}

#' This function calculates the components for each part of the drug and admin costs. NOT CURRENTLY USED
#' 
#' @details This function was used in a previous version of the model, and is no
#'          longer applied. 
#' 
f_drug_calcComponentCost <- function(health_state, TH, multi = FALSE, n = NULL) {
  
  # Take into account whether iterative or not:
  if (multi) {
    q <- drug_i$dos_quant[n]
    f <- drug_i$dos_freq[n]
    r <- drug_i$rdi[n]
    if (is.na(drug_i$cpd)) {
      c <- 0
    } else {
      c <- drug_i$cpd[n] * drug_i$rdi[n] # note some drug costs vary via eMIT
    }
    ac <- drug_i$acpd[n]
  } else {
    q <- drug_i$dos_quant
    f <- drug_i$dos_freq
    r <- drug_i$rdi
    if (is.na(drug_i$cpd)) {
      c <- 0
    } else {
      c <- drug_i$cpd * drug_i$rdi
    }
    ac <- drug_i$acpd
  }
  
  if (drug_i$type == "One off") {
    sched <- c(1,rep(0,TH-1))
  } else {
    # Calculate the dosing schedule over time:
    sched <- rep(c(q,rep(0,f - 1)),(TH + 1) / floor(f))[1:TH]
    
    # Add min start duration if it's active - THIS IS ABSOLUTE TIME at the moment which will need updating
    if (drug_i$start > 0) {
      sched[1:drug_i$start] <- 0
    }
    # Add max treatment duration if it's active
    if (drug_i$stop < TH) {
      sched[drug_i$stop:TH] <- 0
    }
  }
  
  # Multiply the dosing schedule by the drug cost per dose, and multiply that by RDI
  if (c == 0) {
    c_dr <- rep(0,TH)
  } else {
    c_dr <- sched * c
  }
  if (ac == 0) {
    c_ad <- rep(0,TH)
  } else {
    c_ad <- sched * ac
  }
  
  # Return the cost per cycle 
  return(list(
    schedule = sched,
    drug     = c_dr,
    admin    = c_ad
  ))
}