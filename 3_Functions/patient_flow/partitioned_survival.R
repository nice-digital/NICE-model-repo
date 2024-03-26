#' Function for patient flow for state residency for partitioned survival for active treatment lines
#' 
#' @param st_for_1L survival s at time t for first-line, in the list structure using names OS PFS TTD for endpoints
#' 
f_pf_sr_ps_atl <- function(st_for_1L) {
  
  s <- get_elem(st_for_1L,"st")
  
  # compute % in each possible state at each point in time for this 1L treatment,
  # creating a partitioned survival model for the 1L+ for this sequence
  
  PFS_on  <- c(1,pmin(s$PFS,  s$TTD))[1:length(s$OS)]
  PFS_off <- c(0,pmax(s$PFS - s$TTD,0))[1:length(s$OS)]
  PPS_on  <- c(0,pmax(s$TTD - s$PFS,0))[1:length(s$OS)]
  PPS_off <- c(0,s$OS - pmax(s$TTD, s$PFS))[1:length(s$OS)]
  dead    <- c(0,1-s$OS)[1:length(s$OS)]
  
  stopifnot(f_misc_all_same_length(PFS_on, PFS_off, PPS_on, PPS_off, dead))
  
  matrix(
    data = c(PFS_on, PFS_off, PPS_on, PPS_off, dead),
    ncol = 5,
    dimnames = list(NULL,c(
      "PFS_on",
      "PFS_off",
      "PPS_on",
      "PPS_off",
      "dead"
    ))
  )
}




#' Function to compute state residency for use in partitioned survival modelling
#' 
#' 
#' 
f_pf_partSA_state_res <-
  function(treatment_sequence,
           st,
           t_yr,
           lookup,
           pop_n,
           line_n,
           discFacC,
           discFacQ,
           include_plot = FALSE,
           plot_x_lim = 20) {
    
  
  # Collect the correct extrapolations
  extraps <- f_seq_extrapCollector(
    treatment_sequence = treatment_sequence,
    st                 = st,
    lookups            = lookup,
    pop_n              = pop_n
  )
  
  # Create an empty list to house the output and populate it as we go to avoid
  # repeatedly copy pasting things
  out <- list()
  
  # compute undiscounted state residency in each state at first-line ONLY
  out$sr$undisc <- as.list(as.data.table(f_pf_sr_ps_atl(st_for_1L = extraps$st[[paste0("line_", line_n)]])))
  
  # Compute discounted state residency for the two different discount rates:
  out$sr$disc <- list(
    C = lapply(out$sr$undisc, function(endp){
      endp * discFacC
    }),
    Q = lapply(out$sr$undisc, function(endp){
      endp * discFacQ
    })
  )
  
  # Calculate the proportion that will die every cycle per OS for our 1L+ model:
  out$sr$undisc$ddeath_dt <- c(0,diff(out$sr$undisc$dead))
  out$sr$disc$C$ddeath_dt <- out$sr$undisc$ddeath_dt * discFacC
  
  # Calculate the proportion that progress every cycle for our 1L+ model:
  # This assumes deaths are distributed evenly between progression free and progressed health states
  prop_deathinPFS <- p$surv$prop_deathinPFS # proportion of PFS events which are death
  progtotal <- out$sr$undisc$PFS_on + out$sr$undisc$PFS_off
  out$sr$undisc$prog_dt <- c(0, diff(-progtotal)) * (1 - prop_deathinPFS)
  out$sr$disc$C$prog_dt <- out$sr$undisc$prog_dt * discFacC
  out$sr$disc$Q$prog_dt <- out$sr$undisc$prog_dt * discFacQ
  
  # The state PPS_off is then the proportion that are in subsequent lines of treatment
  # beyond the first.
  
  if(include_plot) out$sr_plot <- f_plot_srPlot(out$sr$undisc,t_yr, plot_x_lim)
  
  return(out)
  
}


f_pf_computePF_ps <- function(pops, basic, p, cyclecosts, oneoffcosts, util, ae, substrt){
  
  # pull out our lookup tables for cross referencing so it's automatic and easy
  # to read
  lookup <- basic$lookup
  pop_map <- lookup$pop_map  
  
  
  lapply(pops, function(popu) {
    
    # Get a numeric index for population for the functions, whilst also
    # preserving the auto naming by the lapply:
    popu_n <- as.integer(gsub("pop_","",popu))
    cat(paste0("population: ", popu_n, "\n"))
    
    overall_pop <- as.list(pop_map[Overall.population.number == popu_n,])
    
    f_misc_colcat(
      paste0(
        "PartSA model. Population ",
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
    spop <- paste0("pop",spop_n)
    
    # Pull out the demographic data for our patients. We can use this later
    # (e.g. line 1 baseline age and sex for utiltiy adjustment)
    # 
    # This is based on risk population - note these are baseline population
    # and everyone is in the risk population at baseline!
    demog <- p$demo$live[[rpop]]
    
    # Generate gpop line
    util_gpop_mat <- matrix(
      util$gpop[[rpop]], # R will recycle this automatically
      nrow = basic$th+1,
      ncol = (1*2) +1 + 1 # only one line of treatment in the PartSA model
    )
    
    # make a vector to apply costs and QALYs for AEs as one off
    
    first_cycle <- rep(0,basic$th+1)
    first_cycle[1] <- 1
    
    # Make an empty list to populate with the results:
    sm <- list(
      sr = NULL,
      trt_1L = unique(p$releff$excel_table[Treatment.line == 1 & Population == rpop_n & Include.in.this.analysis. == "Yes",]$Molecule)
    )
    
    # Further filter down the list of 1st line treatments that are available using
    # the output of the sequencing 
    sm$trt_1L <- sm$trt_1L[sm$trt_1L %in% unique(p$seq$n[[paste0("pop_",spop_n)]]$line_1)]
    
    
    # generate state residency by treatment status for all of these 1L treatments:
    PartSA_pf <- lapply(1:length(sm$trt_1L), function(fl_tx) {
      
      tx_seq <- c("line_1" = paste0("mol_",sm$trt_1L[fl_tx]))
      mol_num <- sm$trt_1L[fl_tx]
      
      cat(paste0("molecule: ", tx_seq[1],"\n"))
      mol_name <- p$basic$lookup$ipd$mol[match(mol_num ,p$basic$lookup$ipd$mol$Number),RCC_input_desc]
      
      sr <- f_pf_partSA_state_res(
        treatment_sequence = tx_seq,
        st                 = p$surv$st,
        t_yr               = p$basic$t_yr,
        lookup             = p$basic$lookup,
        pop_n              = rpop_n,
        line_n             = 1,
        discFacC           = p$basic$discFacC,
        discFacQ           = p$basic$discFacQ,
        include_plot       = TRUE,
        plot_x_lim         = 20
      )
      
      # Undiscounted costs
      
      drug_cost <-  cyclecosts$drug[[1]][[tx_seq ]] * (sr$sr$undisc$PFS_on + sr$sr$undisc$PPS_on)
      admin_cost <-  cyclecosts$admin[[1]][[tx_seq ]] * (sr$sr$undisc$PFS_on + sr$sr$undisc$PPS_on) 
      
      substrt_admin_cost <- substrt$partsa[Treatment==mol_name & Population == spop, admin_cost] * sr$sr$undisc$prog_dt
      substrt_drug_cost <- substrt$partsa[Treatment==mol_name & Population == spop, drug_cost] * sr$sr$undisc$prog_dt
      substrt_AE_cost <- substrt$partsa[Treatment==mol_name & Population == spop, AE_cost] * sr$sr$undisc$prog_dt
      
      ruc_preprog       <- cyclecosts$mru_on[[1]][[tx_seq ]] * sr$sr$undisc$PFS_on + cyclecosts$mru_off[[1]][[tx_seq ]] * sr$sr$undisc$PFS_off +
        sum(oneoffcosts[Apply.to == "Init"]$cost) * first_cycle
      ruc_postprog     <- cyclecosts$mru_on[[1]][[tx_seq ]] * sr$sr$undisc$PPS_on + cyclecosts$mru_off[[1]][[tx_seq ]] * sr$sr$undisc$PPS_off
      
      EOL_cost     <- sum(oneoffcosts[Apply.to == "Death"]$cost) * sr$sr$undisc$ddeath_dt
      prog_cost    <- (sum(oneoffcosts[Apply.to == "Prog"]$cost) + sum(oneoffcosts[Apply.to == "Init"]$cost))* sr$sr$undisc$prog_dt 
      
      if(ae$approach == "one-off") {
        AE_dur <- mean(ae$one_off[Treatment.line==1 & Treatment.name==mol_name, duration_weeks])
        AE_cost <- AE_dur * ae$per_cycle[line==1 & trt==mol_name, cost] * first_cycle
      } else {
      
          AE_cost <-  ae$per_cycle[line==1 & trt==mol_name, cost]*(sr$sr$undisc$PFS_on+sr$sr$undisc$PPS_on)
      }
      
      # Discounted costs
      
      ddrug_cost <-  cyclecosts$drug[[1]][[tx_seq ]] * (sr$sr$disc$C$PFS_on + sr$sr$disc$C$PPS_on)
      dadmin_cost <-  cyclecosts$admin[[1]][[tx_seq ]] * (sr$sr$disc$C$PFS_on + sr$sr$disc$C$PPS_on)
      
      dsubstrt_admin_cost <- substrt$partsa[Treatment==mol_name & Population == spop, admin_cost] * sr$sr$disc$C$prog_dt
      dsubstrt_drug_cost <- substrt$partsa[Treatment==mol_name & Population == spop, drug_cost] * sr$sr$disc$C$prog_dt
      dsubstrt_AE_cost <- substrt$partsa[Treatment==mol_name & Population == spop, AE_cost] * sr$sr$disc$C$prog_dt
      
      druc_preprog       <- cyclecosts$mru_on[[1]][[tx_seq ]] * sr$sr$disc$C$PFS_on + cyclecosts$mru_off[[1]][[tx_seq ]] * sr$sr$disc$C$PFS_off +
        sum(oneoffcosts[Apply.to == "Init"]$cost) * first_cycle
      druc_postprog     <- cyclecosts$mru_on[[1]][[tx_seq ]] * sr$sr$disc$C$PPS_on + cyclecosts$mru_off[[1]][[tx_seq ]] * sr$sr$disc$C$PPS_off
      
      dEOL_cost     <- sum(oneoffcosts[Apply.to == "Death"]$cost) * sr$sr$disc$C$ddeath_dt
      dprog_cost    <- (sum(oneoffcosts[Apply.to == "Prog"]$cost) + sum(oneoffcosts[Apply.to == "Init"]$cost)) * sr$sr$disc$C$prog_dt
      
      if(ae$approach == "one-off") {
        AE_dur <- mean(ae$one_off[Treatment.line==1 & Treatment.name==mol_name, duration_weeks])
        dAE_cost <- AE_dur * ae$per_cycle[line==1 & trt==mol_name, cost] * first_cycle
      } else {
        
        dAE_cost <-  ae$per_cycle[line==1 & trt==mol_name, cost]*(sr$sr$disc$C$PFS_on+sr$sr$disc$C$PPS_on)
      }
      
      costs <- list(
        undisc = list(
          drug = drug_cost,
          admin = admin_cost,
          AE = AE_cost,
          substrt_drug_cost = substrt_drug_cost,
          substrt_admin_cost = substrt_admin_cost,
          substrt_AE_cost = substrt_AE_cost,
          mru_preprog = ruc_preprog,
          mru_postprog = ruc_postprog,
          EOL_cost = EOL_cost,
          prog_cost = prog_cost
        ),
        disc = list(
          drug = ddrug_cost,
          admin = dadmin_cost,
          AE = dAE_cost,
          substrt_drug_cost = dsubstrt_drug_cost,
          substrt_admin_cost = dsubstrt_admin_cost,
          substrt_AE_cost = dsubstrt_AE_cost,
          mru_preprog = druc_preprog,
          mru_postprog = druc_postprog,
          EOL_cost = dEOL_cost,
          prog_cost = dprog_cost
        )
      ) 
      
      preprog_hsuv <- as.list(util$hsuv[Population == rpop_n & Treatment.line == 1 & Molecule == mol_num,list(OnTxt,OffTxt)])
      postprog_hsuv <- as.list(util$hsuv[Population == 0 & Treatment.line == 2 & Molecule == mol_num,list(OnTxt,OffTxt)]) # always population 0 at line 2+
      
            # Expand this to the time horizon by multiplying column-wise by the gpop line.
      # The result is U, which we can element-wise multiply by our consolidated trace
      # to get our undiscounted QALYs :)
      Upreprog <- util_gpop_mat[,1:length(preprog_hsuv)] %*% diag(preprog_hsuv)
      Upostprog <- util_gpop_mat[,1:length(postprog_hsuv)] %*% diag(postprog_hsuv)
      
      qalys_PFS <- (sr$sr$undisc$PFS_on * Upreprog[,1] + sr$sr$undisc$PFS_off* Upreprog[,2])  * basic$cl_y 
      qalys_PPS <- (sr$sr$undisc$PPS_on * Upostprog[,1] + sr$sr$undisc$PPS_off* Upostprog[,2])  * basic$cl_y 
      qalys_PPS_AEs <- substrt$partsa[Treatment==mol_name & Population == spop, AE_QALY_impact] * sr$sr$undisc$prog_dt * util_gpop_mat[spop_n+1]
      
      dqalys_PFS <- (sr$sr$disc$Q$PFS_on * Upreprog[,1] + sr$sr$disc$Q$PFS_off* Upreprog[,2])  * basic$cl_y 
      dqalys_PPS <- (sr$sr$disc$Q$PPS_on * Upostprog[,1] + sr$sr$disc$Q$PPS_off* Upostprog[,2])  * basic$cl_y 
      dqalys_PPS_AEs <- substrt$partsa[Treatment==mol_name & Population == spop, AE_QALY_impact] * sr$sr$disc$Q$prog_dt * util_gpop_mat[spop_n+1]
      
      
      if(ae$approach == "one-off") {
        AE_dur <- mean(ae$one_off[Treatment.line==1 & Treatment.name==mol_name, duration_weeks])
        AE_QALYs <- AE_dur * ae$per_cycle[line==1 & trt==mol_name, QALYs] * first_cycle 
        dAE_QALYs <- AE_QALYs 
      } else {
        
        AE_QALYs <-  ae$per_cycle[line==1 & trt==mol_name, QALYs]*(sr$sr$undisc$PFS_on+sr$sr$undisc$PPS_on) * util_gpop_mat[spop_n+1]
        dAE_QALYs <- ae$per_cycle[line==1 & trt==mol_name, QALYs]*(sr$sr$disc$Q$PFS_on+sr$sr$disc$Q$PPS_on) * util_gpop_mat[spop_n+1]
      }
      
      
      
      qalys <- list(
        undisc = list(
          PFS = qalys_PFS,
          PPS = qalys_PPS,
          AE = AE_QALYs,
          AE_PPS = qalys_PPS_AEs
        ),
        disc   = list(
          PFS = dqalys_PFS,
          PPS = dqalys_PPS,
          AE = dAE_QALYs,
          AE_PPS = dqalys_PPS_AEs
        )
      )
      
      ly_PFS_on <- sr$sr$undisc$PFS_on  * basic$cl_y
      ly_PFS_off <- sr$sr$undisc$PFS_off * basic$cl_y 
      ly_PPS_on <- sr$sr$undisc$PPS_on  * basic$cl_y
      ly_PPS_off <- sr$sr$undisc$PPS_off * basic$cl_y 
 
      
      lys <- list(
        PFS_on = ly_PFS_on,
        PFS_off = ly_PFS_off,
        PPS_on = ly_PPS_on,
        PPS_off = ly_PPS_off
      )
    
      out <- list(
        sr = sr$sr,
        sr_plot = sr$sr_plot,
        costs = costs,
        qalys = qalys, 
        lys = lys
      )
        
      return(out)
      
    })
    names(PartSA_pf) <- p$basic$lookup$ipd$mol[match(sm$trt_1L,p$basic$lookup$ipd$mol$Number),RCC_input_desc]
 
    
    return(PartSA_pf)
  })
}