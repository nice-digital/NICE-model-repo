# Code to calculate costs

#' manipulate costs data pulled from inputs spreadsheet outputs vector of drug 
#' costs by week for each molecule where the column number of each matrix 
#' corresponds to the cycle px starts drug at
#' 
#' @param drug_and_admin Excel named range `R_table_drug_admin_costs` from the inputs file
#' @param per_cycle_costs Excel named range `R_table_MRU` from the inputs file
#' @param one_off_costs Excel named range `R_table_MRU_oneoff` from the inputs file
#' @param time_horizon Time horizon from `p$basic$th` from the inputs file
#' @param max_trt_lines Maximum number of lines the R model can handle from `p$basic$R_maxlines`
#' @param RDI_source Excel named range `dd_sc_RDI` from the inputs file
#' @param PSA NOT IMPLEMENTED YET
#' @param verbose Set to TRUE for more console output to track issues
#' 
#' @details Uses the inputs tables directly from the excel input file to populate
#'          drug and admin costs (note that AE costs are done separately) in a 
#'          format that is conducive to the way that relative efficacy and the
#'          overall model structure is organised. The output is organised by:
#'          
#'          - category (drug, admin, mru_on, mru_off)
#'          - molecule
#'          - treatment line
#'          
#'          Allowing for much flexibility in the application of treatment costs.
#'          Each bottom-level element is a vector of values from time zero to
#'          the time horizon (as determined by `time_horizon`). This means that
#'          in a probabilistic setting this full function would be used to generate
#'          `n_psa` "sets" of values. Similarly, in the scenario analysis a different
#'          excel file is loaded each time, which then leads to the inputs being
#'          different if affected by that scenario. 
#' 
f_process_cost_data <- function(drug_and_admin,
                                per_cycle_costs,
                                time_horizon,
                                RDI_source,
                                max_trt_lines,
                                PSA     = FALSE,
                                samples = NULL,
                                verbose = FALSE) {

  # Step one - tidying up the tables:
  
  # cleaning raw_tables
  drug_and_admin_raw <-  data.table(drug_and_admin )
  drug_and_admin_raw[is.na(drug_and_admin_raw)] <- 0
  
  # Populating RDI by line, will give an error message if the line does not exist in Excel
  if (RDI_source == "RDI per RWE") {
    rdi_table <- drug_and_admin_raw[,c("Treatment.name","Molecule","Type",paste0("RDI.Mean_RWE_",1:max_trt_lines,"L")), with=FALSE]
    setnames(rdi_table, c(paste0("RDI.Mean_RWE_",1:max_trt_lines,"L")),paste0("RDI_line_",1:max_trt_lines))
  } else if (RDI_source == "Set all RDI to 100%") {
    rdi_table    <- as.list(drug_and_admin_raw[,c("Treatment.name","Molecule","Type"), with=FALSE])
    RDI_1        <- lapply(1:max_trt_lines, function(x) rep(1,nrow(drug_and_admin)))
    names(RDI_1) <- paste0("RDI_line_",1:max_trt_lines)
    rdi_table    <- as.data.table(c(rdi_table,RDI_1))
  } else {
    rdi_table <- drug_and_admin_raw[,c("Treatment.name","Molecule","Type",paste0("RDI.Mean_Trial_",1:max_trt_lines,"L")), with=FALSE]
    setnames(rdi_table, c(paste0("RDI.Mean_Trial_",1:max_trt_lines,"L")),paste0("RDI_line_",1:max_trt_lines))
  }
  
  # RDI is used to calculate drug and admin costs and represents missed doses whilst on treatment for which no drug or admin charge applies

  if(PSA == TRUE) {
    
    # Standard errors for relative dose intensity depending on the option the user has selected
    if (RDI_source == "RDI per RWE") {
      rdi_SE_table <- drug_and_admin_raw[,c("Treatment.name","Molecule","Type",paste0("RDI.SE_RWE_",1:max_trt_lines,"L")), with=FALSE]
      setnames(rdi_SE_table, c(paste0("RDI.SE_RWE_",1:max_trt_lines,"L")),paste0("RDI_SE_line_",1:max_trt_lines))
    } else if (RDI_source == "Set all RDI to 100%") {
      rdi_SE_table <- as.list(drug_and_admin_raw[,c("Treatment.name","Molecule","Type"), with=FALSE])
      RDI_0        <- lapply(1:max_trt_lines, function(x) rep(0,nrow(drug_and_admin)))
      names(RDI_0) <- paste0("RDI_SE_line_",1:max_trt_lines)
      rdi_SE_table <- as.data.table(c(rdi_SE_table,RDI_0))
    } else {
      rdi_SE_table <- drug_and_admin_raw[,c("Treatment.name","Molecule","Type",paste0("RDI.SE_Trial_",1:max_trt_lines,"L")), with=FALSE]
      setnames(rdi_SE_table, c(paste0("RDI.SE_Trial_",1:max_trt_lines,"L")),paste0("RDI_SE_line_",1:max_trt_lines))
    }
    
    # Merge the mean and SE together into one table irrespective of the option selected
    rdi_table_se   <- merge.data.table(rdi_table,rdi_SE_table)
    
    # cycle row by row creating draws, using common RNG
    rdi_random_gen <- runif(samples)
    
    # generate samples for RDI using the set random numbers generated above
    RDI_samples <- rbindlist(lapply(1:nrow(rdi_table_se), function(row_in_table) {
      
      dat       <- rdi_table_se[row_in_table,]
      
      id        <- dat[,list(Treatment.name,Molecule,Type)]
      id        <- rbindlist(lapply(1:samples, function(x) id))
      
      mean_nams <- paste0("RDI_line_",1:max_trt_lines)
      se_nams   <- paste0("RDI_SE_line_",1:max_trt_lines)
      
      names(mean_nams) <- mean_nams
      names(se_nams)   <- se_nams
      
      # Compile a data.table which has 1 row per iteration with RDI for that iteration,
      # which also takes into account if mean is down as 1 or SE is 0
      data.table(
        id,
        as.data.table(lapply(structure(1:max_trt_lines,.Names=mean_nams), function(trt_line) {
          
          # Estimate alpha and beta parameters for a beta distribution for this treatment
          # line's RDI, then draw samples times from there using random numbs:
          mn <- dat[[mean_nams[trt_line]]]
          se <- dat[[se_nams[trt_line]]]
          
          if (mn < 1 & se > 0) {
            params <- estBetaParams(mn, se ^ 2)
            qbeta(rdi_random_gen, .subset2(params,"alpha"), .subset2(params,"beta"))
          } else if (mn == 1 & se > 0) {
            params <- estBetaParams(0.999, se ^ 2)
            qbeta(rdi_random_gen, .subset2(params,"alpha"), .subset2(params,"beta"))
          } else if (se == 0) {
            rep(mn,samples)
          } else {
            rep(1,samples)
          }
        })),
        iteration = 1:samples
      )
    }))
    
    rm(rdi_table, rdi_table_se, rdi_random_gen)
    
    # This can now be used to generate and merge with the full tables:
    
    
    # Admin costs - similarly to the rdi table, we pull the right columns and then
    # generate them so we can merge it. note that "RDI" is applied to admin costs
    admin_table <-
      drug_and_admin_raw[, c(
        "Treatment.name",
        "Molecule",
        "Type",
        "Admin.cost..per.administration..Mean",
        "Admin.cost..per.administration..SE"),
        with = FALSE]
    
    setnames(admin_table, c("Admin.cost..per.administration..Mean", "Admin.cost..per.administration..SE"), c("mean","se"))
    
    admin_samples <- rbindlist(lapply(1:nrow(admin_table), function(row_in_table) {
      
      dat <- admin_table[row_in_table,]
      
      id  <- dat[,list(Treatment.name,Molecule,Type)]
      id  <- rbindlist(lapply(1:samples, function(x) id))
      
      # Call it the same name that is used in the deterministic model:
      id$Admin.cost..per.administration..Mean <- rnorm(samples,mean = dat$mean,sd = dat$se)
      id$iteration <- 1:samples
      return(id)
    }))
    
    
    # Finally, drug costs. Some of these may have parameter uncertainty in the future 
    # because they come from eMIT or are generics, which have variance in pricing
    
    drug_table <-
      drug_and_admin_raw[, c(
        "Treatment.name",
        "Molecule",
        "Type",
        "Drug.cost.per.dose",
        "Applied.every.x.cycles",
        "Number.of.doses.when.applied",
        "Time.from",
        "Time.to..cycles."
        ),
        with = FALSE]
    
    drug_table <- rbindlist(lapply(1:samples, function(psa_iteration) {
      drug_table$iteration <- psa_iteration
      return(drug_table)
    }))
    
    
    # Now that we have all of the components, we can simply merge them all
    # since they all have id columns that are the same
    
    drug_and_admin <- merge.data.table(merge.data.table(drug_table,RDI_samples),admin_samples)
    
    # per cycle costs table is messy because the named range in excel doesn't include
    # the top header row (because there are 2. This breaks the function.)
    pcc <- as.list(as.data.table(per_cycle_costs[2:nrow(per_cycle_costs),]))
    pcc_nam <- names(pcc)
    pcc <- as.data.table(lapply(1:length(pcc), function(column_index) {
      if (column_index %in% c(1,2)) {
        return(pcc[[column_index]])
      } else {
        return(as.numeric(pcc[[column_index]]))
      }
    }))
    names(pcc) <- pcc_nam
    rm(pcc_nam)
    setnames(pcc,"Cost","mean")
    setnames(pcc,"X5","SE")
    
    # Now perform RNG on per-cycle costs. In the lines that follow using pcc
    # the table will be filtered down by type, molecule and iteration
    
    pcc <- rbindlist(lapply(1:nrow(pcc), function(row_in_table) {
      dat          <- pcc[row_in_table,]
      
      id           <- dat[,list(Type.of.cost,Type,Molecule,Time.from..cycle.,Time.to..cycles.)]
      pcc_m        <- dat$mean
      pcc_se       <- dat$SE
      id           <- rbindlist(lapply(1:samples, function(x) id))
      id$mean      <- rnorm(samples,pcc_m,pcc_se)
      id$iteration <- 1:samples
      return(id)
    }))
    
    
    
    # Produce an empty cost vector to use as a starting point for all components
    # (we can re-use it this way):
    cost <- numeric(time_horizon + 1)
    
    # create a list with $mol = molecule and $data = vector of per-cycle costs
    
    line_labs <- structure(1:max_trt_lines,.Names=paste0("line_",1:max_trt_lines))
    mol_labs  <- structure(unique(drug_and_admin$Molecule),.Names=paste0("mol_",unique(drug_and_admin$Molecule)))
    
    # For the PSA version of the function, we wrap the whole process for the deterministic
    # function in a lapply - one for each PSA iteration!
    
    # one-off costs are actually not part of this function!
    # 
    # 
    # Now, the PSA version is just like the deterministic version but we filter the
    # table by PSA iteration AND line and molecule.
    # 
    # the main model function for the probabilistic version for a specific PSA iteration 
    # then instead of pulling out p$costs$mk pulls ONE PSA ITERATION from p_psa$costs
    # thereby slotting in those inputs with the same format they have in the
    # deterministic model. This means the results are compatible with the lambda
    # approximation method AND ALSO doing the "full-blown" PSA (i.e. full expanded
    # TRACE calculations)
    # 
    # 
    # generate samples structural clones of the cost structure used for deterministic with
    # probabilistic values
    # 
    # NOTE that although this uses multicore, it DOES NOT include any of the RNG.
    # This preserves seeds and avoids issues without loss of speed.
    # 
    s <- 1:samples
    return(with_progress({
      prog <- progressr::progressor(along=s)
      # future_lapply(s, function(psa_iteration) {
      lapply(s, function(psa_iteration) {
        
        prog(paste0("PSA | cost inputs #",psa_iteration))
        
        # note that these only ever apply to on treatment states:
        da_cost_vecs  <- lapply(line_labs, function(tx_line) {
          lapply(mol_labs, function(mol) {
            
            if (mol == 999) {
              return(list(
                drug  = drug,
                admin = admin
              ))
            }
            
            # make a message optionally
            if (verbose) cat(paste0("Drug costs: line_",tx_line, "$mol_",mol,"\n"))
            
            # pull out the data and the correct RDI for this treatment line
            data          <- drug_and_admin[Molecule == mol & iteration == psa_iteration,]
            data$RDI.Mean <- data[[paste0("RDI_line_",tx_line)]]
            
            # We don't have any component identifier (i.e. a number for component)
            # We also don't have a number for type, and type is simply lumped as per cycle or one off
            # meaning that there's no way to uniquely identify the rows.
            # instead we have to cycle through the components:
            # 
            # empty vectors to populate:
            drug    <- cost
            admin   <- cost
            
            cost_vectors <- lapply(1:nrow(data), function(data_row) {
              current_row <- data[data_row,]
              stopifnot(current_row$Type %in% c("per cycle", "one-off"))
              
              # when costs are incurred:
              if (current_row$Type == "per cycle") {
                dose_at_cycles <- seq(
                  floor(current_row$Time.from),
                  ceiling(current_row$Time.to..cycles.),
                  by = floor(current_row$Applied.every.x.cycles)
                )
              } else if (current_row$Type == "one-off") {
                dose_at_cycles <- 1
              }
              
              # drug costs
              cpd   <- current_row$Drug.cost.per.dose
              ndose <- current_row$Number.of.doses.when.applied
              if (current_row$Type == "one-off") {
                rdi <- 1
              } else {
                rdi <- current_row$RDI.Mean
              }
              drug[dose_at_cycles] <- cpd*ndose*rdi
              
              
              # Admin costs
              cpa <- current_row$Admin.cost..per.administration..Mean
              naa <- current_row$Number.of.doses.when.applied
              admin[dose_at_cycles] <- cpa*naa*rdi
              
              # Return a list containing the 3 components of cost:
              return(list(
                drug = drug,
                admin = admin
              ))
            })
            
            # Cycle through our list, cumulative adding up drug costs, and cumulatively
            # adding up admin costs
            if(length(cost_vectors) == 1) {
              # If there's only one component, give me drug and admin as 2 vectors:
              return(list(
                drug  = cost_vectors[[1]]$drug,
                admin = cost_vectors[[1]]$admin
              ))
            } else {
              # If there's more than one component, add them up by drug and admin
              Reduce(
                x = 2:length(cost_vectors),
                accumulate = FALSE,
                init = cost_vectors[[1]],
                function(prev, row_n) {
                  list(
                    drug = prev$drug + cost_vectors[[row_n]]$drug,
                    admin = prev$admin + cost_vectors[[row_n]]$admin
                  )
                }
              )
            }
          })
        })
        
        # Drug costs and MRU work on different schedules, so we need a separate process
        # for MRU. it follows the same structure though:
        
        mol_labs  <- structure(unique(pcc$Molecule),.Names=paste0("mol_",unique(pcc$Molecule)))
        
        mru_vecs  <- lapply(line_labs, function(tx_line) {
          lapply(mol_labs, function(mol) {
            if(verbose) cat(paste0("Resource use: line_", tx_line, "$mol_",mol,"\n"))
            
            dat <- pcc[Molecule == mol & iteration == psa_iteration,]
            
            onoff <- structure(c("On", "Off"), .Names = c("on", "off"))
            
            if(mol != 999) {
              lapply(onoff, function(treatment_status) {
                tx <- dat[grep(treatment_status,Type.of.cost,ignore.case = FALSE),]
                
                # Work out each cost vector and then add them up (reduce with `+` will
                # add up a bunch of list elements to give you the result :-)
                return(Reduce(
                  `+`,
                  lapply(1:nrow(tx), function(tx_row) {
                    inp <- tx[tx_row,]
                    
                    dose_at_cycles <- floor(seq(inp$Time.from..cycle., inp$Time.to..cycles.))
                    out <- cost
                    out[dose_at_cycles] <-  inp$mean
                    out
                  })
                ))
              })
            } else {
              # BSC MRU, just has 1 row
              dose_at_cycles <- floor(seq(floor(dat$Time.from..cycle.), ceiling(dat$Time.to..cycles.)))
              out <- cost
              out[dose_at_cycles] <-  dat$mean
              return(list(on=out,off=cost))
            }
          })
        })
        
        # Return a reorganised list by category, line, molecule
        return(list(
          drug    = lapply(da_cost_vecs, function(li) {lapply(li, function(mol) mol[["drug"]])}),
          admin   = lapply(da_cost_vecs, function(li) {lapply(li, function(mol) mol[["admin"]])}),
          mru_on  = lapply(mru_vecs, function(li)     {lapply(li, function(mol) mol[["on"]])}),
          mru_off = lapply(mru_vecs, function(li)     {lapply(li, function(mol) mol[["off"]])})
        ))
        
      })
    }))
    
    
  } else {
    
    # DETERMINISTIC MODEL:
    
    
    # Subsetting to the required columns
    drug_and_admin <- merge.data.table(drug_and_admin_raw[, list(
        Treatment.name,
        Type,
        Molecule,
        Drug.cost.per.dose,
        Admin.cost..per.administration..Mean,
        Applied.every.x.cycles,
        Number.of.doses.when.applied,
        Time.from,
        Time.to..cycles.
      )],rdi_table)
    
    
    # per cycle costs table is messy because the named range in excel doesn't include
    # the top header row (because there are 2. This breaks the function.)
    pcc     <- as.list(as.data.table(per_cycle_costs[2:nrow(per_cycle_costs),]))
    pcc_nam <- names(pcc)
    pcc     <- as.data.table(lapply(1:length(pcc), function(column_index) {
      if (column_index %in% c(1,2)) {
        return(pcc[[column_index]])
      } else {
        return(as.numeric(pcc[[column_index]]))
      }
    }))
    names(pcc) <- pcc_nam
    rm(pcc_nam)
    setnames(pcc,"Cost","mean")
    setnames(pcc,"X5","SE")
    
    
    # Produce an empty cost vector to use as a starting point for all components
    # (we can re-use it this way):
    cost <- numeric(time_horizon + 1)
    
    # create a list with $mol = molecule and $data = vector of per-cycle costs
    
    line_labs <- structure(1:max_trt_lines,.Names=paste0("line_",1:max_trt_lines))
    mol_labs  <- structure(unique(drug_and_admin$Molecule),.Names=paste0("mol_",unique(drug_and_admin$Molecule)))
    
    # note that these only ever apply to on treatment states:
    da_cost_vecs  <- lapply(line_labs, function(tx_line) {
      lapply(mol_labs, function(mol) {
        
        # make a message optionally
        if (verbose) cat(paste0("Drug costs: line_",tx_line, "$mol_",mol,"\n"))
        
        # pull out the data and the correct RDI for this treatment line
        data <- drug_and_admin[Molecule == mol,]
        data$RDI.Mean <- data[[paste0("RDI_line_",tx_line)]]
        
        # We don't have any component identifier (i.e. a number for component)
        # We also don't have a number for type, and type is simply lumped as per cycle or one off
        # meaning that there's no way to uniquely identify the rows.
        # instead we have to cycle through the components:
        # 
        # empty vectors to populate:
        drug    <- cost
        admin   <- cost
        
        cost_vectors <- lapply(1:nrow(data), function(data_row) {
          current_row <- data[data_row,]
          stopifnot(current_row$Type %in% c("per cycle", "one-off"))
          
          # when costs are incurred:
          if (current_row$Type == "per cycle") {
            dose_at_cycles <- seq(
              floor(current_row$Time.from),
              ceiling(current_row$Time.to..cycles.),
              by = floor(current_row$Applied.every.x.cycles)
            )
          } else if (current_row$Type == "one-off") {
            dose_at_cycles <- 1
          }
          
          # drug costs
          cpd   <- current_row$Drug.cost.per.dose
          ndose <- current_row$Number.of.doses.when.applied
          if (current_row$Type == "one-off") {
            rdi <- 1
          } else {
            rdi <- current_row$RDI.Mean
          }
          drug[dose_at_cycles] <- cpd*ndose*rdi
          
          
          # Admin costs
          cpa <- current_row$Admin.cost..per.administration..Mean
          naa <- current_row$Number.of.doses.when.applied
          admin[dose_at_cycles] <- cpa*naa*rdi
          
          # Return a list containing the 3 components of cost:
          return(list(
            drug = drug,
            admin = admin
          ))
        })
        
        # Cycle through our list, cumulative adding up drug costs, and cumulatively
        # adding up admin costs
        if(length(cost_vectors) == 1) {
          # If there's only one component, give me drug and admin as 2 vectors:
          return(list(
            drug  = cost_vectors[[1]]$drug,
            admin = cost_vectors[[1]]$admin
          ))
        } else {
          # If there's more than one component, add them up by drug and admin
          Reduce(
            x = 2:length(cost_vectors),
            accumulate = FALSE,
            init = cost_vectors[[1]],
            function(prev, row_n) {
              list(
                drug = prev$drug + cost_vectors[[row_n]]$drug,
                admin = prev$admin + cost_vectors[[row_n]]$admin
              )
            }
          )
        }
      })
    })
    
    # Drug costs and MRU work on different schedules, so we need a separate process
    # for MRU. it follows the same structure though:
    
    mol_labs  <- structure(unique(pcc$Molecule),.Names=paste0("mol_",unique(pcc$Molecule)))
    
    mru_vecs  <- lapply(line_labs, function(tx_line) {
      lapply(mol_labs, function(mol) {
        if(verbose) cat(paste0("Resource use: line_", tx_line, "$mol_",mol,"\n"))
        
        dat <- pcc[Molecule == mol,]
        
        onoff <- structure(c("On", "Off"), .Names = c("on", "off"))
        
        if(mol != 999) {
          lapply(onoff, function(treatment_status) {
            tx <- dat[grep(treatment_status,Type.of.cost,ignore.case = FALSE),]
            
            # Work out each cost vector and then add them up (reduce with `+` will
            # add up a bunch of list elements to give you the result :-)
            return(Reduce(
              `+`,
              lapply(1:nrow(tx), function(tx_row) {
                inp <- tx[tx_row,]
                
                dose_at_cycles <- floor(seq(inp$Time.from..cycle., inp$Time.to..cycles.))
                out <- cost
                out[dose_at_cycles] <-  inp$mean
                out
              })
            ))
          })
        } else {
          # BSC MRU, just has 1 row
          dose_at_cycles <- floor(seq(floor(dat$Time.from..cycle.), ceiling(dat$Time.to..cycles.)))
          out <- cost
          out[dose_at_cycles] <-  dat$mean
          return(list(on=out,off=cost))
        }
      })
    })
    
    # Return a reorganised list by category, line, molecule
    return(list(
      drug    = lapply(da_cost_vecs, function(li) {lapply(li, function(mol) mol[["drug"]])}),
      admin   = lapply(da_cost_vecs, function(li) {lapply(li, function(mol) mol[["admin"]])}),
      mru_on  = lapply(mru_vecs, function(li)     {lapply(li, function(mol) mol[["on"]])}),
      mru_off = lapply(mru_vecs, function(li)     {lapply(li, function(mol) mol[["off"]])})
    ))
  }
}

#' a catch-all function for "other" costs that don't fit into the usual structure.
#' 
#' @param one_off_costs Named range `R_table_MRU_oneoff` from Excel input file
#' 
#' 
#' @details At this stage, this function only incorporates one-off costs upon
#'          entering each line, to represent radiotherapy, palliative surgery, 
#'          and end of life costs.
#'          
#' 
#' 
f_process_other_cost_data <- function(one_off_costs,
                                      PSA   = FALSE,
                                      n_psa = NULL) {
  #cleaning
  one_off_costs <- data.table(one_off_costs)
  one_off_costs <- one_off_costs[-1,]
  setnames(one_off_costs,old = c("Cost","X5"),new = c("mean", "SE"))
  one_off_costs[, (c("mean", "SE")) := lapply(.SD, as.numeric), .SDcols = c("mean", "SE")]
  
  #remove treatment initiation as this is handled in f_process_cost_data
  one_off_costs <- one_off_costs[Type.of.cost != "Treatment initiation\r\n",]

  if (PSA == FALSE) {
    #return mean
    one_off_costs[,`:=`(SE = NULL)]
    setnames(one_off_costs,old = c("mean"),new = c("cost"))
    
  } else {
    #return a set of sampled values
    stop("Code to be written")
    draws <- do.call(
      rbind,
      lapply(1:nrow(one_off_costs), function(row_n) {
        dat <- one_off_costs[row_n,]
        # Normal distribution as population level (gamma tends to normal - Prof. Stevenson)
        return(rnorm(n_psa,dat$mean,dat$SE))
      })
    )
    one_off_costs <- cbind(one_off_costs,draws)
  }
  
  return(one_off_costs)
}


#' used to clean up some subsequent treatment data for the partitioned survival model. not used.
f_process_subs_txt_data <- function(subs_txt, PSA = FALSE){
  
  #cleaning
  colnames(subs_txt)[1] <- c("Treatment")
  
    if (PSA == FALSE) {
    #return mean
      subs_txt = subs_txt[,1:6]    
     colnames(subs_txt)[3] <- "drug_cost"
     colnames(subs_txt)[4] <- "admin_cost"
     colnames(subs_txt)[5] <- "AE_cost"
     colnames(subs_txt)[6] <- "AE_QALY_impact"
  } else {
    #return a set of sampled values
    stop("Code to be written")
  }
  
  return(subs_txt)
  
}

#' Apply per cycle costs to expanded Markov trace (expanded = with tunnel states). NOT USED.
#' 
#' @details A legacy function which was used during the development of the state transition model.
#'          Not used in the final version, which incorporates these calculations
#'          within function `f_pf_computePF_mk` directly instead. The reason for this
#'          is that `f_pf_computePF_mk` is dealing with a lot of memory at once
#'          and instead of copy-pasting a potentially 500MB of RAM-using expanded trace
#'          (with up to 18,000 health states...) it uses that within the function,
#'          calculating drug costs there.
#'          
#'          The idea with this function is to find the appropriate columns
#'          corresponding to the on or off-treatment tunnel states for 2L+
#'          treatments, and then multiply those columns of the trace by
#'          the corresponding vector of costs going from time 0 to the time horizon
#'          
#' 
f_apply_costs_to_trace <- function(cost_by_cycle, 
                                   AE_costs,
                                   per_cycle_costs,
                                   one_off_cost,
                                   markov_trace, 
                                   seq) {
  
  MRU_cost_at_tx_start <- sum(one_off_cost$cost[one_off_cost$Event.to.apply == "Progression"])
  cost_death <- one_off_cost$cost[one_off_cost$Event.to.apply == "Death"]

  cost_matrix <- matrix(data = 0, 
                        nrow = nrow(markov_trace), 
                        ncol = ncol(markov_trace))
  trt_line    <- 0
  th          <- nrow(markov_trace)
  
  for (molecule in seq[1:length(seq)-1]) {
    trt_line  <- trt_line + 1
    costs     <- cost_by_cycle[[which(
                    lapply(cost_by_cycle, 
                           function(x) x$mol == molecule) == TRUE)]]

    #add in AE costs
    if (i$dd_apply_AE_options!="one-off") {costs$on_trt <- costs$on_trt + AE_costs[molecule+1,"cost"]}
    
    if (trt_line == 1) {
      cost_matrix[,1] <- markov_trace[,1] * costs$on_trt
      cost_matrix[,2] <- markov_trace[,2] * costs$off_trt
    } else {

      # On treatment
      #add in MRU costs for starting a new line of treatment
      costs$on_trt[1] <- costs$on_trt[1] + MRU_cost_at_tx_start
        
      startcol        <- (3 + (trt_line - 2) * 2 * th)
      endcol          <- startcol + th - 1
      # note the t() in the line below - apply[x,1,fun] returns a column vector 
      # even though input is a row so needs transposing
      cost_matrix[,startcol:endcol]  <- 
        t(apply(markov_trace[,startcol:endcol], 1, function(x) x * costs$on_trt))
      
      # Off treatment
      startcol <- endcol+1
      endcol   <- startcol + th - 1
      # note the t() in the line below - apply[x,1,fun] returns a column vector 
      # even though input is a row so needs transposing
      cost_matrix[,startcol:endcol]  <- 
        t(apply(markov_trace[,startcol:endcol], 1, function(x) x * costs$off_trt))
      }
  }
  
  # Add in additional cost of starting BSC
  molecule <- seq[length(seq)]
  if (molecule != 999) stop("Last molecule in sequence is not BSC! Full sequence: ", seq)
  
  startcol <- endcol + 1
  cost_matrix[,startcol] <- markov_trace[,startcol] * MRU_cost_at_tx_start
  
  #cost of death
  propn_entering_by_cycle <- markov_trace[,ncol(markov_trace)] - 
                              shift(markov_trace[,ncol(markov_trace)], n=1, fill=0)
  
  cost_matrix[,ncol(cost_matrix)] <- propn_entering_by_cycle * cost_death
  
  
  cost_matrix
}


