#' Function to compute the TPs that populate matrix M using a set of extrapolations
#' which follow the structure resulting from using function `f_seq_extrapCollector`
#' to "collect" extrapolations.
#' 
#' @param st_list the result of `f_seq_extrapCollector` (i.e. a set of extrapolated survival)
#' 
#' @details cycles through treatment lines (`f_seq_extrapCollector` reduces down to those required)
#'          pulling out the TTD, TTP (censoring for death) and PFS lines. These
#'          are then used to compute transition probabilities between health states
#'          for use within `M` later. For on-treatment states, next_on death_on
#'          disc_on and stay_on are computed, and for off-treatment states
#'          next_off death_off and stay_off are computed. This then creates
#'          a full set of time-varying transition probabilities from time 0
#'          to the time horizon. For BSC, only OS is required.
#'          
#' 
#' 
f_pf_mk_ComputeTPs <- function(st_list) {
  out <- lapply(1:length(st_list), function(trt_li) {
    
    # Figure out the time horizon:
    tl    <- st_list[[trt_li]]
    th    <- length(tl$OS$st)
    zeros <- rep(0,th)
    
    if (trt_li < length(st_list)) {
      # Pull out the endpoints for ease of reading:
      TTD <- tl$TTD$st
      TTP <- tl$TTP$st
      PFS <- tl$PFS$st
      
      # Calculate estimated transition probabilities out of those curves:
      tp_TTD <- 1-(TTD/shift(TTD,fill=1))
      tp_TTP <- 1-(TTP/shift(TTP,fill=1))
      tp_PFS <- 1-(PFS/shift(PFS,fill=1))
      
      # Apply assumptions to figure out TP from each state to each other state:
      disc <- tp_TTD
      
      # ON-treatment transition probabilities
      
      # probability of death directly from on treatment state:
      death_on <- 1 - tp_TTP - (1-tp_PFS)
      
      # Going directly onto next line from the on treatment state
      
      # The probability of discontinuation removing the probability of any other event
      disc_only <- disc - tp_TTP - death_on
      disc_only_adj <- pmax(disc_only, 0)
      
      # The probability of going to next therapy but not death directly from on-treatment
      next_on  <- tp_TTP - (disc_only_adj - disc_only)
      
      # the probability of staying on treatment is simply the inverse of the probability
      # of going off treatment for any reason
      stay_on  <- 1 - disc
      
      # The probability of discontinuation when on treatment
      disc_on  <- disc_only_adj
      
      # Just to note some equivalence here:
      # 1 - next_on - disc - stay_on == 1 - tp_TTP - (1-tp_PFS)
      
      next_off  <- next_on
      death_off <- death_on
      stay_off  <- 1 - death_on - next_on
      
    } else {
      # Pull out the endpoints for ease of reading:
      OS  <- tl$OS$st
      
      # Calculate estimated transition probabilities out of those curves:
      tp_OS  <- 1-(OS/shift(OS,fill=1))
      
      # Apply assumptions to figure out TP from each state to each other state:
      stay_on  <- 1-tp_OS
      disc_on  <- zeros
      next_on  <- zeros
      death_on <- tp_OS
      
      next_off  <- zeros
      death_off <- zeros
      stay_off  <- zeros
      
      # Testing sum to 1:
      
    }
    
    # return a list of TPs for each of the 7 required to compute the Markov model
    return(list(
      disc_on   = disc_on,
      stay_on   = stay_on,
      next_on   = next_on,
      death_on  = death_on,
      stay_off  = stay_off,
      next_off  = next_off,
      death_off = death_off
    ))
    
  })
  names(out) <- c(paste0("L",1:(length(st_list)-1)),"NT")
  return(out)
}


#' Function to compute Markov traces when using a state-transition framework. 
#' 
#' @param pops string with e.g. 'pop_0' to run risk population 0. make sure those populations exist.
#' @param basic p$basic from p
#' @param sequences p$seq from p
#' @param survival p$surv from p
#' @param costs p$costs$mk from p
#' @param util List containing hsuv (p$util$mk) and gpop (p$util$gpop)
#' @param ae p$ae$mk from p
#' @param eff_table p$releff$excel_table from p
#' @param verbose whether or not the user wants verbose output whilst computing
#' @param include_plots Whether to generate trace plots or not
#' @param just_nlines optional. lets you run just one set of nlines (e.g. just 1 active treatment line or 2)
#' 
f_pf_computePF_mk <-
  function(pops,
           basic,
           demo,
           sequences,
           survival,
           costs,
           util,
           ae,
           eff_table,
           verbose = FALSE,
           include_plots = FALSE,
           just_nlines = NULL,
           just_seq = NULL) {
  
  # Compute the maximum active treatment lines possible by going through all sequences
  # and figuring it out
  max_active <- max(sapply(sequences$qc,ncol))-1

  # pull out our lookup tables for cross referencing so it's automatic and easy
  # to read
  lookup <- basic$lookup
  pop_map <- lookup$pop_map
  
  lapply(pops, function(popu) {
    
    # Get a number for the population for use subsetting and filtering stuff.
    popu_n <- as.integer(gsub("pop_","",popu))
    
    overall_pop <- as.list(pop_map[Overall.population.number == popu_n,])
    
    f_misc_colcat(
      paste0(
        "ST model. Population ",
        popu_n,
        "\t| Sequences: '",
        overall_pop$Sequencing.population,
        "'\t | 1L Risk: '",
        overall_pop$Risk.population,
        "', further lines assumed all population"
      )
    )
    # Map the population appropriately using pop_map as derived above:
    
    # Basically, get the risk and sequence populations and use those to extract
    # the correct sequences and extrapolations, then go from there.
    rpop_n <- overall_pop$Risk.population.number
    rpop <- paste0("pop_",rpop_n)
    
    spop_n <- overall_pop$Sequencing.population.number
    spop <- paste0("pop_",spop_n)
    
    # Pull out the demographic data for our patients. We can use this later
    # (e.g. line 1 baseline age and sex for utiltiy adjustment)
    # 
    # This is based on risk population - note these are baseline population
    # and everyone is in the risk population at baseline!
    demog <- demo[[rpop]]
    
    
    # First-line included treatments for THIS risk population
    L1_inc_rpop <- sort(unique(eff_table[Treatment.line == 1 &
                            Population == rpop_n &
                            Include.in.this.analysis. == "Yes", ]$Molecule))
    
    # All further lines are risk population 0 only:
    L2p_inc_rpop <- lapply(2:max_active, function(line_n) {
      sort(unique(eff_table[Treatment.line == line_n &
                              Population == 0 & Include.in.this.analysis. == "Yes", ]$Molecule))
    })
    
    
    # Stick them together so that we can go through all the treatment sequences
    # further reducing them down to those that are allowed via efficacy, NMAs 
    # and relative efficacy assumptions / exogenous data.
    trt_filter <- c(list(L1_inc_rpop),L2p_inc_rpop)
    rm(L1_inc_rpop)
    rm(L2p_inc_rpop)
    
    # Get the treatment sequences for this population
    ateachline <- sequences$n[[spop]]
    ateachline$nlines <- ncol(ateachline)-rowSums(is.na(ateachline))
    
    # Split them up by the number of treatment lines:
    molecule_list_full <- lapply(2:(ncol(ateachline)-1), function(n_lines) {
      ateachline[nlines==n_lines,1:n_lines]
    })
    
    # Go through each of these, filtering down using the per line allowable molecules.
    # This then imposes exclusion via excel (i.e. data/assumption) as well as the rule set for providing
    # each treatment sequence (i.e. policy).
    
    sequence_list <- lapply(molecule_list_full, function(lines_amount) {
      Reduce(
        x = 1:(ncol(lines_amount)-1),
        init = lines_amount,
        accumulate = FALSE,
        f = function(prev, treatment_line) {
          
          # Find out which molecules are in that line column and filter the table
          molecules_allowed <- trt_filter[[treatment_line]]
          which_to_allow <- which(prev[[treatment_line]] %in% molecules_allowed)
          return(prev[which_to_allow,])
        }
      )
    })
    
    # User can make it just run one set of nlines. mainly for QC
    if (!is.null(just_nlines)) {
      sequence_list <- sequence_list[just_nlines]
      names(sequence_list) <- paste0("active_lines_",just_nlines)
    } else {
      names(sequence_list) <- paste0("active_lines_",1:length(sequence_list))
    }
    
    # Before we cycle down the individual treatment sequences, we need to generate
    # the appropriate utility multiplier to be applied to health state utility values
    # to account for the ageing of the population. This vector is then compiled into
    # a matrix with TH rows and the maximum possible number of columns in the 
    # consolidated trace. This matrix can then be multiplied by the dagonal 
    # matrix of HSUVs relevant to a sequence to produce the HSUV matrix accounting for
    # ageing. This matrix can be element-wise multiplied by the consolidated trace
    # to produce undiscounted QALYs. The undiscounted QALYs can be discounted 
    # using the discount factor in basic$discFacQ to produce the discounted version
    

    # Generate gpop line
    util_gpop_mat <- matrix(
      util$gpop[[rpop]], # R will recycle this automatically
      nrow = basic$th+1,
      ncol = (max_active*2) +1 + 1
    )
    
    
    # Each one of these sequences requires a Markov trace and a unique id. 
    # when qcing: amount_of_lines <- sequence_list[[1]]
    
    return(lapply(sequence_list, function(amount_of_lines) {
      
      cat(paste0("M: pathways with ", ncol(amount_of_lines)-1, " active treatment lines \n"))
      
      seq_id <- 1:nrow(amount_of_lines)
      names(seq_id) <- paste0("seq_",seq_id)
      
      if (!is.null(just_seq)) {
        seq_id <- seq_id[just_seq]
      }
      
      with_progress({
        pr <- progressr::progressor(along = seq_id, message = paste0("Pop ", popu, ", ", ncol(amount_of_lines)-1, " ATLs:"))
        future_lapply(seq_id, future.chunk.size = 1, function(seq_n) {
          
          # lapply(seq_id, function(seq_n) {
          number_of_lines <- ncol(amount_of_lines)
          
          # Pull out this specific treatment sequence:
          trt_seq <- structure(paste0("mol_",amount_of_lines[seq_n,]),.Names=names(amount_of_lines[seq_n,]))
          trts_n <- as.integer(gsub("mol_","",trt_seq))
          
          pr(paste0("sequence: ",paste(trts_n,collapse = "|")))
          
          trt_seq_plotLab <- trt_seq
          names(trt_seq_plotLab) <- gsub("line_","L",names(trt_seq))
          names(trt_seq_plotLab)[length(names(trt_seq_plotLab))] <- "NT"
          
          treatment_names <- basic$lookup$ipd$mol[match(trts_n,basic$lookup$ipd$mol$Number),Description]
          treatment_abbr  <- basic$lookup$ipd$mol[match(trts_n,basic$lookup$ipd$mol$Number),RCC_input_desc]
          
          if (verbose) f_misc_colcat(
            paste0(
              "M: ",
              ncol(amount_of_lines) - 1,
              " active line(s). seq: ",
              seq_n,
              " - ",
              paste(paste(trt_seq,treatment_names, collapse = " -> "), collapse = " ")
            ),
            31
          )
          
          # Collect the extrapolations that are available for that sequence
          # 
          # IMPORTANT NOTE: THIS FUNCTION USES POP 0 FOR ALL 2L+ LINES!!!
          # 
          # Note the pop_n argument is numeric and for RISK POPULATION
          extraps <- f_seq_extrapCollector(
            treatment_sequence = trt_seq,
            st                 = survival$st,
            lookups            = lookup,
            pop_n              = rpop_n,
            pop_0_2Lp          = TRUE
          )
          
          # Calculate what we need
          TPs <- f_pf_mk_ComputeTPs(st_list = extraps$st)
          
          # Correct for non-finite TPs:
          TPs <- lapply(TPs, function(li) {
            lapply(li, function(tp) {
              tp[!is.finite(tp)] <- 0
              if (length(tp < (basic$th+1))) {
                tp <- c(tp,rep(0,basic$th+1 - length(tp)))
              }
              tp
            })
          })
          
          # Corrected time horizon:
          th <- length(TPs$L1$disc_on)
          
          # death is certain at the end of the time horizon. impossible to get here
          # for any 2L+ state anyway since it takes 1 cycle to get there and therefore
          # the time horizon + line - 1 is the last cycle possible.
          TPs <- lapply(TPs, function(li) {
            if ("death_on" %in% names(li)) li$death_on[th] <- 1
            if ("death_off" %in% names(li)) li$death_off[th] <- 1
            return(li)
          })
          
          # For computing L1 off trt entrants:
          L1_onoff_adjustment <- TPs$L1$death_on + TPs$L1$next_on
          
          # Compile M using these TPs
          
          M <- f_markov_M_compiler_TPupdate(
            tp = TPs,
            n_lines = number_of_lines-1,
            TH = th,
            verbose = verbose
          )
          
          # Compute Markov trace
          TRACE <- f_markov_M_extrapolator(
            L1_tp   = TPs$L1,
            TH      = th,
            M       = M,
            N       = f_markov_calcN(n_lines = number_of_lines-1, TH = th),
            prog    = verbose,
            verbose = verbose
          )
          
          # discounted trace for the costs (diag.M is rowwise, M.diag is columnwise)
          DTRACE <- basic$discFacC * TRACE
          
          # column-wise recycling tested per the below:
          # DTRACE2 <- Diagonal(length(basic$discFacC),basic$discFacC) %*% TRACE
          # identical(DTRACE,DTRACE2)
          # rm(DTRACE2)
          
          # we can drop M to save memory, as we don't use it again
          rm(M)
          rm(TPs)
          
          # Generates entrant list
          
          # Return the reduced traces per the model's requirements.
          consolidated_trace <- f_markov_traceConsolidator(
            full_trace = TRACE,
            split_list = NULL,
            L1_onoff_adjustment = L1_onoff_adjustment,
            TH         = th,
            n_lines    = number_of_lines-1,
            discFacQ   = basic$discFacQ,
            discFacC   = basic$discFacC
          )
          
          # Here we will use TRACE (the expanded trace) to compute drug costs
          # and mru costs because they are too complicated to multiply by simple
          # numbers.
          # by line and treatment status (and component for HCRU + drug costs).
          # 
          # This will replace returning TRACE in the output
          
          if (verbose) cat("Computing model outcomes per cycle per cycle-in-state")
          
          # rename the colstarts to something we can refer to automatically:
          names(consolidated_trace$col_starts) <- Reduce(
            x = 1:(length(trt_seq)-1),
            init = names(consolidated_trace$col_starts),
            accumulate = FALSE,
            f = function(prev, li) {
              gsub(paste0("L",li),trt_seq[li],prev)
            }
          )
          names(consolidated_trace$col_starts)[which(names(consolidated_trace$col_starts) == "BSC")] <- "mol_999_on"
          
          # cycle through the treatment sequence, computing the different cost components:
          
          # abbreviate cost per cycle (cpc) and costs on initation (coi). we only need the
          # mols which we are using for this sequence. one off are for all.
          cpc <- costs$per_cycle
          coi  <- costs$one_off$Prog
          cod  <- costs$one_off$Death
          
          pf_costs <- lapply(1:length(trt_seq), function(line_n) {
            mol <- trt_seq[line_n]
            
            if (line_n == 1) {
              # First line, we don't need to do tunnel multiplication:
              
              tr_1l_on <- TRACE[,1]
              tr_1l_off <- TRACE[,2]
              
              drug_cost    <- cpc$drug[[names(mol)]][[mol]] * tr_1l_on
              admin_cost   <- cpc$admin[[names(mol)]][[mol]] * tr_1l_on
              ruc_on       <- cpc$mru_on[[names(mol)]][[mol]] * tr_1l_on
              ruc_off      <- cpc$mru_off[[names(mol)]][[mol]] * tr_1l_off
              
              rm(tr_1l_on)
              rm(tr_1l_off)
              
              # Discounted version:
              dtr_1l_on <- DTRACE[,1]
              dtr_1l_off <- DTRACE[,2]
              
              ddrug_cost    <- cpc$drug[[names(mol)]][[mol]] * dtr_1l_on
              dadmin_cost   <- cpc$admin[[names(mol)]][[mol]] * dtr_1l_on
              druc_on       <- cpc$mru_on[[names(mol)]][[mol]] * dtr_1l_on
              druc_off      <- cpc$mru_off[[names(mol)]][[mol]] * dtr_1l_off
              
              rm(dtr_1l_on)
              rm(dtr_1l_off)
              
            } else if (line_n < length(trt_seq)) {
              # For the central states we define diagonal matrices, multiply 
              # by the appropriate columns in TRACE, then calculate rowsums!
              
              # Add in one-off costs upon entering this treatment line:
              mru_on_tx <- cpc$mru_on[[names(mol)]][[mol]]
              mru_on_tx[1] <- mru_on_tx[1] + coi
              
              dc     <- Diagonal(n = th, x = cpc$drug[[names(mol)]][[mol]])
              ac     <- Diagonal(n = th, x = cpc$admin[[names(mol)]][[mol]])
              ru_on  <- Diagonal(n = th, x = mru_on_tx)
              ru_off <- Diagonal(n = th, x = cpc$mru_off[[names(mol)]][[mol]])
              
              start_col <- consolidated_trace$col_starts[grep(paste0(mol,"_"),names(consolidated_trace$col_starts))]
              on_cols <- start_col[grep("_on",names(start_col))]:(start_col[grep("_on",names(start_col))] + th - 1)
              off_cols <- start_col[grep("_off",names(start_col))]:(start_col[grep("_off",names(start_col))] + th - 1)
              
              tr_on <- TRACE[,on_cols]
              
              drug_cost  <- rowSums(tr_on %*% dc)
              admin_cost <- rowSums(tr_on %*% ac)
              ruc_on     <- rowSums(tr_on %*% ru_on)
              ruc_off    <- rowSums(TRACE[,off_cols] %*% ru_off)
              
              rm(tr_on)
              
              dtr_on <- DTRACE[,on_cols]
              
              ddrug_cost  <- rowSums(dtr_on %*% dc)
              dadmin_cost <- rowSums(dtr_on %*% ac)
              druc_on     <- rowSums(dtr_on %*% ru_on)
              druc_off    <- rowSums(DTRACE[,off_cols] %*% ru_off)
              
              rm(dtr_on)
              rm(dc)
              rm(ac)
              rm(ru_on)
              rm(ru_off)
              rm(start_col)
              rm(on_cols)
              rm(off_cols)
              
            } else {
              
              # Future switch? Add in one-off costs upon entering this treatment line:
              # In future, link this to the excel switch for applying this cost
              # upon starting BSC.
              if (TRUE) {
                # Sometimes BSC is coming after any active treatment line.
                # The cost should come from the last active treatment line's BSC
                # values:
                if (is.null(cpc$mru_on[[names(mol)]])) {
                  # reduce line id by 1 just for pulling out the right data
                  # from cpc:
                  namline_n <- as.numeric(gsub("line_","",names(mol)))
                  line_nam_temp <- paste0("line_",namline_n-1)
                  
                  # Now use this adjusted name to pull out the right values
                  mru_on_tx <- cpc$mru_on[[line_nam_temp]][[mol]]
                } else {
                  mru_on_tx <- cpc$mru_on[[names(mol)]][[mol]]
                }
                # Now add the cost on initiation (coi)
                mru_on_tx[1] <- mru_on_tx[1] + coi
              }
              
              # we're in BSC, we do the same as above but have to pull some of the costs
              # for BSC from previous line, and there's no drug costs
              ru_on  <- Diagonal(n = th, x = mru_on_tx)
              
              drug_cost <- rep(0,th)
              admin_cost <- rep(0,th)
              ddrug_cost <- rep(0,th)
              dadmin_cost <- rep(0,th)
              
              start_col <- consolidated_trace$col_starts[grep(mol,names(consolidated_trace$col_starts))]
              on_cols <- start_col:(start_col + th - 1)  
              
              ruc_on   <- rowSums(TRACE[,on_cols] %*% ru_on)
              druc_on  <- rowSums(DTRACE[,on_cols] %*% ru_on)
              ruc_off  <- rep(0,th)
              druc_off <- ruc_off
              
              rm(ru_on)
              rm(on_cols)
            }
            
            # Ok now no matter how many lines and their order we've done drug cost
            # with full memory for our treatment sequence -_-
            
            return(list(
              molecule = mol,
              undisc = list(
                drug = drug_cost,
                admin = admin_cost,
                mru_on = ruc_on,
                mru_off = ruc_off
              ),
              disc = list(
                drug = ddrug_cost,
                admin = dadmin_cost,
                mru_on = druc_on,
                mru_off = druc_off
              )
            ))
          })
          names(pf_costs) <- names(trt_seq)
          
          
          # We now no longer need TRACE (or discounted trace) and can drop it to free up (a lot of) memory
          rm(TRACE)
          rm(DTRACE)
          rm(cpc)
          rm(coi)
          
          # Finally, add in end of life costs by multiplying the entrants to death from consolidated
          # trace by the cost upon death
          pf_costs$eol <- list(
            undisc = consolidated_trace$entrants[,ncol(consolidated_trace$entrants)] * cod,
            disc   = consolidated_trace$disc$C$entrants[,ncol(consolidated_trace$disc$C$entrants)] * cod
          )
          rm(cod)
          
          # QALYs and AEs are calculated using the consolidated traces:
          # pf_qalys
          
          # Get the HSUVs for the active treatment line by treatment status:
          hsuv_active <- unlist(lapply(1:(number_of_lines-1), function(line_n) {
            
            # ASSUMEs POPULATION 0 IF LATER LINES!!!!
            if (line_n == 1) {
              line_hsuv <- as.list(util$hsuv[Population == rpop_n & Treatment.line == line_n & Molecule == trts_n[line_n],list(OnTxt,OffTxt)])
            } else {
              line_hsuv <- as.list(util$hsuv[Population == 0 & Treatment.line == line_n & Molecule == trts_n[line_n],list(OnTxt,OffTxt)])
            }
            
            lab_line <- paste0("L",line_n)
            names(line_hsuv) <- paste0(lab_line,c("_on","_off"))
            unlist(line_hsuv)
          }))
          # Get the HSUV for BSC:
          hsuv_bsc <- structure(
            util$hsuv[Population == 0 & Treatment.line == number_of_lines & Molecule == trts_n[number_of_lines],]$OnTxt,
            .Names = paste0("L",number_of_lines,"_on")
          )
          
          # Stick them together to form the vector of HSUVs per possible state for
          # this treatment sequence:
          hsuv_unadj <- c(hsuv_active,hsuv_bsc,"dead"=0)
          
          # Expand this to the time horizon by multiplying column-wise by the gpop line.
          # The result is U, which we can element-wise multiply by our consolidated trace
          # to get our undiscounted QALYs :)
          U <- util_gpop_mat[,1:length(hsuv_unadj)] %*% diag(hsuv_unadj)
          
          rm(hsuv_active)
          rm(hsuv_bsc)
          rm(hsuv_unadj)
          
          # Undiscounted QALYs (note AEs are separate as they're a separate thing):
          pf_qalys <- list(
            undisc = (consolidated_trace$full_lines * basic$cl_y) * U,
            disc   = (consolidated_trace$disc$Q$full_lines * basic$cl_y) * U
          )
          
          rm(U)
          
          # AEs undiscounted
          
          # pf_ae
          # Same as above, compute your AE stuff here
          
          if(ae$approach[1] == "one-off") {
            
            ae_impact_atl <- do.call(
              rbind,
              lapply(1:(number_of_lines-1), function(line_n) {
                c(cost = ae$per_cycle[ae$per_cycle$trt == treatment_abbr[line_n] & ae$per_cycle$line == line_n]$cost,
                  qaly = ae$per_cycle[ae$per_cycle$trt == treatment_abbr[line_n] & ae$per_cycle$line == line_n]$QALYs,
                  dur  = mean(
                    ae$one_off$duration_weeks[
                      ae$one_off$Molecule == trts_n[line_n]
                      & ae$one_off$Treatment.line == min(line_n,2)
                    ]))
              })
            )
            
            ae_impact_atl <- data.table(ae_impact_atl)
            ae_oo_cost <- ae_impact_atl$cost * ae_impact_atl$dur
            ae_oo_qaly <- ae_impact_atl$qaly * ae_impact_atl$dur
            
            ent <- data.table(as.matrix(consolidated_trace$entrants))
            ent <- rbindlist(list(data.frame(t(structure(rep(0,ncol(ent)),.Names=colnames(ent)))),ent))
            ent <- cbind(L1_on = c(1,rep(0,nrow(ent)-1)),ent)
            
            dent_c <- data.table(as.matrix(consolidated_trace$disc$C$entrants))
            dent_c <- rbindlist(list(data.frame(t(structure(rep(0,ncol(dent_c)),.Names=colnames(dent_c)))),dent_c))
            dent_c <- cbind(L1_on = c(1,rep(0,nrow(dent_c)-1)),dent_c)
            
            dent_q <- data.table(as.matrix(consolidated_trace$disc$Q$entrants))
            dent_q <- rbindlist(list(data.frame(t(structure(rep(0,ncol(dent_q)),.Names=colnames(dent_q)))),dent_q))
            dent_q <- cbind(L1_on = c(1,rep(0,nrow(dent_q)-1)),dent_q)
            
            nonzero_ind <- 1:number_of_lines + 0:(number_of_lines-1)
            
            nstate <- ncol(consolidated_trace$full_lines)
            n_cyc  <- nrow(consolidated_trace$full_lines)
            nam_st <- colnames(consolidated_trace$full_lines)
            
            # clever trick to expand a vector
            ae_oo_cost <- "[<-"(numeric(nstate), nonzero_ind, c(ae_oo_cost,0))
            ae_oo_qaly <- "[<-"(numeric(nstate), nonzero_ind, c(ae_oo_qaly,0))
            
            # multiply payoff by entrants
            pf_ae <- list(
              undisc = list(
                costs = matrix(
                  as.matrix(ent)[1:th,] %*% diag(ae_oo_cost),
                  ncol = nstate,
                  nrow = n_cyc,
                  dimnames = list(NULL,nam_st)
                ),
                qalys = matrix(
                  as.matrix(ent)[1:th,] %*% diag(ae_oo_qaly),
                  ncol = nstate,
                  nrow = n_cyc,
                  dimnames = list(NULL,nam_st)
                )
              ),
              disc = list(
                costs = matrix(
                  as.matrix(dent_c)[1:th,] %*% diag(ae_oo_cost),
                  ncol = nstate,
                  nrow = n_cyc,
                  dimnames = list(NULL,nam_st)
                ),
                qalys = matrix(
                  as.matrix(dent_q)[1:th,] %*% diag(ae_oo_qaly),
                  ncol = nstate,
                  nrow = n_cyc,
                  dimnames = list(NULL,nam_st)
                )
              )
            )
            
            
          } else {
            # compute cost per cycle on treatment and qalys lost per cycle on treatment:
            ae_cpc_ontrt <- unlist(lapply(1:number_of_lines, function(act_line_n) {
              as.list(ae$per_cycle[line == act_line_n & molecule == trts_n[act_line_n],])$cost
            }))
            ae_qlpc_ontrt <- unlist(lapply(1:number_of_lines, function(act_line_n) {
              as.list(ae$per_cycle[line == act_line_n & molecule == trts_n[act_line_n],])$QALYs
            }))
            # Make a correction for if we have 4 atl's:
            if (number_of_lines == max_active+1) {
              ae_cpc_ontrt[5] <- ae$per_cycle[line == max_active & molecule == trts_n[number_of_lines],]$cost
              ae_qlpc_ontrt[5] <- ae$per_cycle[line == max_active & molecule == trts_n[number_of_lines],]$QALYs
            }
            nonzero_ind <- 1:number_of_lines + 0:(number_of_lines-1)
            
            # pad the cpc and qlpc with 0s:
            ae_cpc <- "[<-"(numeric(ncol(consolidated_trace$full_lines)), nonzero_ind, ae_cpc_ontrt)
            ae_qlpc <- "[<-"(numeric(ncol(consolidated_trace$full_lines)), nonzero_ind, ae_qlpc_ontrt)
            
            names(ae_cpc) <- colnames(consolidated_trace$full_lines)
            names(ae_qlpc) <- colnames(consolidated_trace$full_lines)
            
            # Implement adverse events (yay)
            pf_ae <- list(
              undisc = list(
                costs = matrix(
                  consolidated_trace$full_lines %*% diag(ae_cpc),
                  ncol = ncol(consolidated_trace$full_lines),
                  nrow = nrow(consolidated_trace$full_lines),
                  dimnames = list(NULL,colnames(consolidated_trace$full_lines))
                ),
                qalys = matrix(
                  consolidated_trace$full_lines %*% diag(ae_qlpc),
                  ncol = ncol(consolidated_trace$full_lines),
                  nrow = nrow(consolidated_trace$full_lines),
                  dimnames = list(NULL,colnames(consolidated_trace$full_lines))
                )
              ),
              disc = list(
                costs = matrix(
                  consolidated_trace$disc$C$full_lines %*% diag(ae_cpc),
                  ncol = ncol(consolidated_trace$full_lines),
                  nrow = nrow(consolidated_trace$full_lines),
                  dimnames = list(NULL,colnames(consolidated_trace$full_lines))
                ),
                qalys = matrix(
                  consolidated_trace$disc$Q$full_lines %*% diag(ae_qlpc),
                  ncol = ncol(consolidated_trace$full_lines),
                  nrow = nrow(consolidated_trace$full_lines),
                  dimnames = list(NULL,colnames(consolidated_trace$full_lines))
                )
              )
            )
          }
          
          if(include_plots) {
            f_plot_mk_draw_consol_trace(consol_trace = consolidated_trace$full_lines,
                                        treatment_names = treatment_names,
                                        tmax = 15)
          }
          
          
          return(list(
            population   = popu,
            trt_nam      = treatment_names,
            trace_consol = consolidated_trace$full_lines,
            costs        = pf_costs,
            qalys        = pf_qalys,
            aes          = pf_ae
          ))
        })
      })
      
      
      
    }))
  })
  
}



# Summary functions -------------------------------------------------------


#' Function which takes the output of the state transition model and 
#' sums up to various different levels which are optional. For e.g. a PSA, 
#' only the top line should be returned, or if drug cost is required to be separated
#' from total costs, then breakdown should be true. In a full deterministic analysis
#' full breakdowns can inform the tables which will go into reporting.
#' 
#' @param pf_list the output from the state transition model (i.e. `pf$mk`)
#' @param disc_undisc either `disc` or `undisc` as text to tell the model to return discounted or undiscounted results
#' @param lookups the lookup tables we have used throughout the model (i.e. `p$basic$lookup`)
#' @param full_breakdown logical. if TRUE (default), return results by sequence, category and treatment line
#' @param breakdown logical. if TRUE (default), return results by sequence and category in one neat table
#' @param ypc years per cycle - used when "undisc" is entered for disc_undisc. `p$basic$cl_y`
#' 
#' 
f_pf_mk_summary <- function(pf_list, disc_undisc,lookups, full_breakdown = TRUE, breakdown = TRUE, ypc) {
  
  stopifnot(disc_undisc %in% c("disc","undisc"))
  
  if (disc_undisc == "disc") {
    items_to_get <- c("trt_nam", "disc")
  } else {
    items_to_get <- c("trt_nam", "undisc", "trace_consol")
  }
  
  # Tip: to get the first popu, do this:
  # popu <- get_elem(pf$mk,items_to_get)[[1]]
  lapply(get_elem(pf_list,items_to_get), function(popu) {
    
    index_names <- Reduce(
      x = 1:length(popu),
      init = lapply(popu,length),
      accumulate = FALSE,
      f = function(prev, this_line) {
        prev[[this_line]] <- 1:prev[[this_line]]
        names(prev[[this_line]]) <- paste0("seq_",prev[[this_line]])
        if (this_line > 1) {
          prev[[this_line]] <- max(prev[[this_line]]) + prev[[this_line - 1]]
        }
        return(prev)
      }
    )
    
    # Now collapse the list by one level
    fixed_names <- Reduce(
      x = names(popu),
      init = list(),
      accumulate = FALSE,
      f = function(prev, lines_txt) {
        out <- popu[[lines_txt]]
        names(out) <- paste0("seq_",index_names[[lines_txt]])
        prev <- c(prev,out)
        return(prev)
      }
    )
    
    n_seq <- length(fixed_names)
    
    # Costs:
    # Now we can just cycle through this flattened list summarising one at a time
    # this_sequence <- fixed_names[[1]]
    cost_full_breakdown <- lapply(fixed_names, function(this_sequence) {
      costs_less_aes <- do.call(
        rbind,lapply(this_sequence$costs, function(tx_line) {
          if(class(tx_line) == "numeric") {
            sum(tx_line)
          } else {
            unlist(lapply(tx_line,sum))
          }
        })
      )
      
      eol <- costs_less_aes["eol","mru_on"]
      costs_less_aes <- costs_less_aes[paste0("line_",1:(length(this_sequence$trt_nam))),]
      
      costs_ae <- colSums(this_sequence$aes$costs)
      costs_ae <- unlist(lapply(1:(length(this_sequence$trt_nam)), function(tx_n) {
        if (tx_n < length(this_sequence$trt_nam)) {
          sum(costs_ae[which(names(costs_ae) %in% paste0("L",tx_n,c("_on","_off")))])
        } else {
          as.numeric(costs_ae["BSC"])
        }
      }))
      
      return(cbind(
        costs_less_aes, 
        ae_cost = costs_ae, 
        eol = c(rep(0,length(this_sequence$trt_nam)-1),eol)
      ))
    })
    cost_breakdown      <- do.call(rbind,lapply(cost_full_breakdown,colSums))
    cost_total          <- rowSums(cost_breakdown)
    
    # Now do the same with QALYs
    qaly_full_breakdown <- lapply(fixed_names, function(this_sequence) {
      qaly_less_aes <- matrix(colSums(this_sequence$qalys),ncol=1,dimnames = list(colnames(this_sequence$qalys),"qaly"))
      qaly_ae <- matrix(colSums(this_sequence$aes$qalys),ncol=1,dimnames = list(colnames(this_sequence$aes$qalys),"ae_qaly"))
      return(cbind(qaly_less_aes,qaly_ae))
    })
    qaly_breakdown <- do.call(rbind,lapply(qaly_full_breakdown,colSums))
    qaly_total <- rowSums(qaly_breakdown)
    
    out <- list()
    if (full_breakdown) out$full_breakdowns <- list()
    if (breakdown) out$breakdowns <- list()
    
    # make some unique id's for each sequence:
    trt_n <- unlist(lapply(lapply(fixed_names, function(trt_seq) trt_seq$trt_nam), function(drug_names) {
      paste(lookups$ipd$mol[match(drug_names,lookups$ipd$mol$Description),]$Number,collapse = "→")
    }))
    trt_numb <- lapply(lapply(fixed_names, function(trt_seq) trt_seq$trt_nam), function(drug_names) {
      lookups$ipd$mol[match(drug_names,lookups$ipd$mol$Description),]$Number
    })
    trt_txt <- unlist(lapply(fixed_names, function(x) paste(x$trt_nam, collapse="→")),use.names = FALSE)
    
    
    # Now, if we're returning full breakdowns, we want to go by sequence, returning
    # costs and qalys as 2 objects.
    
    if (full_breakdown) out$full_breakdowns <-  lapply(structure(1:n_seq,.Names=paste0("seq_",1:n_seq)), function(this_seq) {
      list(
        n = trt_n[this_seq],
        numb = trt_numb[[this_seq]],
        txt = trt_txt[this_seq],
        cost = cost_full_breakdown[[this_seq]],
        qaly = qaly_full_breakdown[[this_seq]]
      )
    })
    if (breakdown) {
      out$breakdowns <-  data.table(
        trt_n = trt_n,
        trt = trt_txt,
        cost_breakdown,
        qaly_breakdown
      )
      
    }
    
    out$res <- data.table(
      trt_n = trt_n,
      trt   = trt_txt,
      costs = cost_total,
      qalys = qaly_total
    )
    
    if (disc_undisc == "undisc") {
      life_years <- lapply(fixed_names, function(this_sequence) {
        colSums(this_sequence$trace_consol * ypc)[1:(ncol(this_sequence$trace_consol)-1)]
      })
      if (breakdown) {
        out$ly <- list()
        out$ly$breakdown <- data.table(
          trt_n = trt_n,
          trt = trt_txt,
          rbindlist(lapply(life_years, function(x)data.table(t(data.frame(x)))),fill = TRUE)
        )
      }
      out$res$ly <- unlist(lapply(life_years,sum))
    }
    
    # Return the results table!
    return(out)
    
  })
}

