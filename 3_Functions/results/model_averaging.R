#' Function to generate weighted average cost and QALY breakdown per OVERALL
#' population.
#' 
#' 
#' @param res_obj Model results from the ST model, usually `res$mk`, but if loaded in via rds file, then that.
#' @param pop_oo numeric of length one. overall_pop corresponding to the lookup table `p$basic$lookup$pop_map`
#' @param subs `p$costs$settings$subsTx`, which is named range `R_table_sub_txts_prop_n_costs` from Excel (containing the weightings)
#' @param ptchar named range `R_table_ptchar` from excel (patient characteristics)
#' @param disc logical. discounted or undiscounted results
#' @param lookups usually `p$basic$lookup` containing all the lookup tables to translate numerical id's into labels
#' @param `no_active_lines` number of active lines coming from the Excel input
#' @param `max_lines` maximum number of lines including BSC, set to 5 by default
#' 
f_res_wa_model <-
  function(res_obj,
           pop_oo,
           subs,
           ptchar = NULL,
           disc,
           lookups,
           no_active_lines,
           max_lines) {
  
  stopifnot(no_active_lines %in% 2:4)
    
  # Lookup tables:
  lu_mol <- lookups$ipd$mol
  lu_pop <- lookups$pop_map
  
  # Weighting table:
  # The subsTx table only needs to be computed once, so just get it done:
  subs$L1 <- lu_mol[match(subs$Line.1,RCC_input_desc,nomatch = NA),]$Number
  subs$L2 <- lu_mol[match(subs$Line.2,RCC_input_desc,nomatch = NA),]$Number
  subs$L3 <- lu_mol[match(subs$Line.3,RCC_input_desc,nomatch = NA),]$Number
  subs$L4 <- lu_mol[match(subs$Line.4,RCC_input_desc,nomatch = NA),]$Number
  subs$L5 <- 999
  subs    <-
    subs[!is.na(Population), list(
      Population,
      L1,
      L2,
      L3,
      L4,
      L5,
      Adj.proportion.given.line.1
    )]
  
  # Checking the maximum number of active lines from the Excel input file and creating a subsetted filed which sums up the number receiving each possible pathway
  # after the maximum line and assigns it to the prior treatment
  
  if (no_active_lines ==3 ) {
    new_subs <- subs[, sum(Adj.proportion.given.line.1), by = .(Population, L1, L2, L3)]
    colnames(new_subs)[colnames(new_subs) == 'V1']  <- "Adj.proportion.given.line.1"
    new_subs <- cbind(new_subs, L4 = rep("NA", nrow(new_subs)), L5 = rep(999, nrow(new_subs)))
    
    }
  if (no_active_lines ==2 ) {
    new_subs <- subs[, sum(Adj.proportion.given.line.1), by = .(Population, L1, L2)]
    colnames(new_subs)[colnames(new_subs) == 'V1']  <- "Adj.proportion.given.line.1"
    new_subs <- cbind(new_subs, L3 = rep("NA", nrow(new_subs)),L4 = rep("NA", nrow(new_subs)), L5 = rep(999, nrow(new_subs)))
    
  }  
  
  
  # Compile trt_n to match with the breakdown table
  subs$trt_n <- do.call(paste, c(subs[,paste0("L",1:max_lines),with=FALSE], sep="→"))
  
  
  # now we have a unique identified assigned for each sequence assign the new values taking into account the maximum number of lines
  
  if (no_active_lines <4) {
    new_subs$trt_n <- do.call(paste, c(new_subs[, paste0("L", 1:max_lines), with = FALSE], sep ="→"))
    subs$Adj.proportion.given.line.1 <- new_subs$Adj.proportion.given.line.1[match(subs$trt_n, new_subs$trt_n)]
    subs$Adj.proportion.given.line.1[is.na(subs$Adj.proportion.given.line.1)] <- 0
  }
  
  subs$trt_n <- gsub("→NA","",subs$trt_n)
  subs <- subs[,list(Population,L1,trt_n,Adj.proportion.given.line.1)]
  
  # Finally, filter down to only the overall population we're doing 
  pop_s_n <- lu_pop[match(pop_oo, lu_pop$Overall.population.number)]$Sequencing.population.number
  
  subs <- subs[Population == paste0("pop",pop_s_n),]
  
  # Moving onto the results to proces for all sequences for this population:
  stopifnot(disc %in% c(TRUE,FALSE))
  
  # Determine which objects to pull out of the results object (discounted or undisc)
  if (disc) {
    d <- "disc"
  } else {
    d <- "undisc"
    res_ly  <- res_obj[[d]][[paste0("pop_",pop_oo)]]$ly$breakdown
  }
  
  # Pull out the full breakdown table for this OVERALL population, which has
  # outcomes by line and category for costs and qalys
  res_fbd <- res_obj[[d]][[paste0("pop_",pop_oo)]]$full_breakdowns
  
  
  # work out max dimensions for individual breakdown tables:
  max_rows_c <- max_lines
  max_rows_q <- (max_lines)*2
  
  
  # cycle through the full breakdown applying weightings, then pivoting the tables
  # wider. It should end up 1 row data.table with trt_n as one column.
  # (fbdt = full breakdown table)
  fbdt <- rbindlist(
    lapply(res_fbd, function(tx) {
      tx_li <- nrow(tx$cost)
      
      # expand tx$cost to all possible treatment lines if required
      if (tx_li < max_rows_c) {
        tx$cost <- rbind(tx$cost,
                         matrix(
                           rep(0, ncol(tx$cost) * (max_rows_c - tx_li)),
                           ncol = ncol(tx$cost),
                           nrow = max_rows_c - tx_li,
                           dimnames = list(paste0("line_", (nrow(
                             tx$cost
                           ) + 1):(
                             nrow(tx$cost) + (max_rows_c - tx_li)
                           )), colnames(tx$cost))
                         ))
        eol_temp <- max(tx$cost[,"eol"])
        tx$cost[,"eol"] <- c(rep(0,max_rows_c-1),eol_temp)
      }
      
      # Do the same for QALYs, which has 2 for each active treatment line
      # + BSC + dead
      tx_liq <- nrow(tx$qaly)
      if (tx_liq < max_rows_q) {
        rows_to_add <- max_rows_q - tx_liq
        n_lines <- (nrow(tx$qaly) - 2) / 2
        
        # names for extra columns
        li_txt <- unlist(lapply(paste0("L", (n_lines + 1):(max_lines-1)), function(x)
          paste0(x, c("_on", "_off"))))
        
        tx$qaly <- rbind(tx$qaly,
                         matrix(
                           data = rep(0, length(li_txt) * ncol(tx$qaly)),
                           ncol = ncol(tx$qaly),
                           nrow = length(li_txt),
                           dimnames = list(li_txt, colnames(tx$qaly))
                         ))
        
        tx$qaly <- tx$qaly[c(
          rownames(tx$qaly)[grep("L", rownames(tx$qaly))], 
          rownames(tx$qaly)[grep("L", rownames(tx$qaly), invert = TRUE)]), 
        ]
      }
      
      # Now the tables have been padded out, we can multiply by the weighting
      # and spread them wide:
      line_id <- structure(1:nrow(tx$cost),.Names=rownames(tx$cost))
      
      tx$cost <- do.call(cbind,lapply(line_id, function(tx_line) {
        dat <- t(tx$cost[paste0("line_",tx_line),])
        colnames(dat) <- paste0(colnames(dat),"_L",tx_line)
        return(dat)
      }))
      
      qaly_index <- structure(1:nrow(tx$qaly),.Names=rownames(tx$qaly))
      
      tx$qaly <- do.call(cbind,lapply(1:length(qaly_index), function(cata) {
        state_name <- names(qaly_index)[cata]
        dat <- t(tx$qaly[state_name,])
        colnames(dat) <- paste0(colnames(dat),"_",state_name)
        return(dat)
      }))
      
      return(data.table(trt_n = tx$n,trt = tx$txt,tx$cost,tx$qaly))
      
    })
  )
  
  
  # apply the weightings and sum up by first-line therapy:
  # Sum up by first line therapy:
  
  subs$Population<- NULL
  if (disc == FALSE) {
    
    # If undiscounted then add in the additional columns for state residency:
    
    res_ly <- res_ly[,-"trt"]
    ly_cols <- colnames(res_ly)
    ly_cols <- ly_cols[!ly_cols %in% c("trt_n")]
    
    res_lym <- merge.data.table(
      subs,
      res_ly,
      by="trt_n"
    )
    
    res_lym <- as.data.frame(res_lym)
    res_lym[,ly_cols] <- res_lym[,ly_cols] * res_lym$Adj.proportion.given.line.1
    res_lym <- as.data.table(res_lym)
    
    res_lym[is.na(res_lym)] <- 0
    res_lym[,(ly_cols) := lapply(.SD, sum),.SDcols = ly_cols,by="L1"]
    res_lym <- res_lym[,head(.SD,1),by="L1"]
    res_lym <- res_lym[,-c("Adj.proportion.given.line.1","trt_n")]
    
    colnames(res_lym)[2:length(colnames(res_lym))] <- paste0("ly_",colnames(res_lym)[2:length(colnames(res_lym))])
  }
  
  # subs$Adj.proportion.given.line.1<- NULL
  # Merge
  
  final_table <- merge.data.table(
    subs,
    fbdt,
    by = "trt_n"
  )
  
  cols <- colnames(final_table)
  cols <- cols[!cols %in% c("trt_n", "L1", "Adj.proportion.given.line.1", "trt")]
  
  # apply the weightings all at once:
  
  final_table <- as.data.frame(final_table)
  final_table[,cols] <- final_table[,cols] * final_table$Adj.proportion.given.line.1
  final_table <- as.data.table(final_table)
  
  # add everything up by first line treatment
  
  final_table[, (cols) := lapply(.SD, sum), .SDcols = cols, by = "L1"]
  final_table <- final_table[,head(.SD,1),by="L1"]
  
  final_table$trt_n <- NULL
  final_table$trt <- NULL
  final_table$Adj.proportion.given.line.1 <- NULL
  
  cols_to_drop <- c(paste0("eol_L",1:4),"ae_qaly_BSC","qaly_dead", "ae_qaly_dead")
  final_table[,(cols_to_drop) := NULL]
  
  if (disc) {
    return(final_table)
  } else {
    return(merge.data.table(final_table,res_lym,by = "L1"))
  }
  
}



#' WORK IN PROGRESS - NOT USED CURRENTLY. Weighting by prior adjuvant. Requires the same treatments to be avaialable 
#' in both populations as it produces a weighted average table
f_res_wam_prior_adjuvant <- function(lookups, char, wam_disc, wam_undisc) {
  
  # lookup tables:
  lu_pop <- lookups$pop_map
  lu_rpop <- lookups$ipd$pop
  
  # cycling index:
  rpop_labs <- structure(
    sort(unique(p$basic$lookup$pop_map$Risk.population.number)),
    .Names = paste0("risk_pop_",sort(unique(p$basic$lookup$pop_map$Risk.population.number)))
  )
  
  # Cycle through risk pops, calculating weighted average breakdown tables:
  lapply(rpop_labs, function(risk_pop) {
    
    char$rpop <- lu_rpop[match(char$Population,lu_rpop$Description),]$Number
    char <- char[Treatment.line == 1 & rpop==risk_pop,]
    
    # pull out our weighting for prior (N.B. this is 0 in base-case)
    w_prior <- char$Prior.IO...in.12.months.Mean
    
    # our risk population is then our lu_pop for this risk population:
    tab_rpop <- lu_pop[Risk.population.number == risk_pop,]
    
    # We then pull the OVERALL populations from the weighted table, weighting the
    # outcomes further:
    pops_oo <- tab_rpop$Overall.population.number
    
    # pull out discounted and undiscounted
    disc   <- wam_disc[pops_oo]
    undisc <- wam_undisc[pops_oo]
    
    # Determine which of the populations is the no prior adjuvant one
    # in case different ones are in a different order from each other:
    npa_which <- grep("no prior adjuvant",tab_rpop$Sequencing.population)
    
    # Make a separate number for each one so we know which is our no prior and which
    # is our prior
    pa_i  <- pops_oo[-npa_which]
    npa_i <- pops_oo[npa_which]
    
    # Discounted weighted table:
    wa_disc   <- (w_prior * disc[[paste0("pop_",pa_i)]])   + ((1-w_prior)*disc[[paste0("pop_",npa_i)]])
    wa_undisc <- (w_prior * undisc[[paste0("pop_",pa_i)]]) + ((1-w_prior)*undisc[[paste0("pop_",npa_i)]])
    
    return(list(
      disc = wa_disc,
      undisc = wa_undisc
    ))
  })
}


#' quick function to add up stuff in the weighted model tables:
f_res_sum_weighted_model <- function(rd, rud ) {
  
  # Make some flags (bits of strings inside columns)
  cflag <- c("drug", "admin", "mru", "ae_cost", "eol")
  qflag <- c("qaly_")
  lyflag <- c("ly_")
  
  rd$costs <- rowSums(dplyr::select(rd,contains(cflag)))
  rd$qalys <- rowSums(dplyr::select(rd,contains(qflag)))
  ly       <- dplyr::select(rud,contains(lyflag))
  rd$ly    <- rowSums(dplyr::select(ly,!contains(qflag,)))
  # make a results table and return that
  return(rd[,list(L1,costs,qalys,ly)])
}



#' Function to generate weighted traces by first-line therapy received according to the
#' weightings provided in the excel named range `R_table_sub_txts_prop_n_costs`.
#' 
#' @param pf_list_pop the patient flow for this markov run for this OVERALL population
#' @param oo_pop_n numeric for the overall population (1-6)
#' @param subs `p$costs$settings$subsTx`, which is named range `R_table_sub_txts_prop_n_costs` from Excel (containing the weightings)
#' @param lookups lookup list from `p$basic$lookup`
#' @param max_lines maximum number of treatment lines, used to pad out smaller traces.
#' 
f_pf_wa_sr_plots <- function(pf_list_pop, oo_pop_n, subs, lookups, max_lines, no_active_lines) {
  
  # Lookup tables:
  lu_mol <- lookups$ipd$mol
  lu_pop <- lookups$pop_map
  
  # Weighting table:
  # The subsTx table only needs to be computed once, so just get it done:
  subs$L1 <- lu_mol[match(subs$Line.1,RCC_input_desc,nomatch = NA),]$Number
  subs$L2 <- lu_mol[match(subs$Line.2,RCC_input_desc,nomatch = NA),]$Number
  subs$L3 <- lu_mol[match(subs$Line.3,RCC_input_desc,nomatch = NA),]$Number
  subs$L4 <- lu_mol[match(subs$Line.4,RCC_input_desc,nomatch = NA),]$Number
  subs$L5 <- 999
  subs    <-
    subs[!is.na(Population), list(
      Population,
      L1,
      L2,
      L3,
      L4,
      L5,
      Adj.proportion.given.line.1
    )]
  
  if (no_active_lines ==3 ) {
    new_subs <- subs[, sum(Adj.proportion.given.line.1), by = .(Population, L1, L2, L3)]
    colnames(new_subs)[colnames(new_subs) == 'V1']  <- "Adj.proportion.given.line.1"
    new_subs <- cbind(new_subs, L4 = rep("NA", nrow(new_subs)), L5 = rep(999, nrow(new_subs)))
    
  }
  if (no_active_lines ==2 ) {
    new_subs <- subs[, sum(Adj.proportion.given.line.1), by = .(Population, L1, L2)]
    colnames(new_subs)[colnames(new_subs) == 'V1']  <- "Adj.proportion.given.line.1"
    new_subs <- cbind(new_subs, L3 = rep("NA", nrow(new_subs)),L4 = rep("NA", nrow(new_subs)), L5 = rep(999, nrow(new_subs)))
    
  }  
  
  
  # Compile trt_n to match with the breakdown table
  subs$trt_n <- do.call(paste, c(subs[,paste0("L",1:max_lines),with=FALSE], sep="→"))
  
  
  # now we have a unique identified assigned for each sequence assign the new values taking into account the maximum number of lines
  
  if (no_active_lines <4 ) {
    new_subs$trt_n <- do.call(paste, c(new_subs[,paste0("L",1:max_lines),with=FALSE], sep="→"))
    subs$Adj.proportion.given.line.1 <- new_subs$Adj.proportion.given.line.1[match(subs$trt_n, new_subs$trt_n)]
    subs$Adj.proportion.given.line.1[is.na(subs$Adj.proportion.given.line.1)] <- 0
  }
  
  
  subs$trt_n <- gsub("→NA","",subs$trt_n)
  subs <- subs[,list(Population,L1,trt_n,Adj.proportion.given.line.1)]
  
  # Finally, filter down to only the overall population we're doing 
  pop_s_n <- lu_pop[match(oo_pop_n, lu_pop$Overall.population.number)]$Sequencing.population.number
  subs <- subs[Population == paste0("pop",pop_s_n),]
  
  
  # now, we need to "flatten" the list, as it's currently broken up into 
  sequences_flattened <- Reduce(
    x = 1:length(pf_list_pop),
    init = list(),
    accumulate = FALSE,
    f = function(prev, n_lines) {
      if (n_lines == 1) {
        prev <- pf_list_pop[[n_lines]]
        return(prev)
      } else {
        sequences <- pf_list_pop[[n_lines]]
        starting_n <- sum(unlist(lapply(pf_list_pop[1:(n_lines-1)],length))) + 1
        len <- length(sequences)
        names(sequences) <- paste0("seq_",starting_n:(starting_n+len-1))
        prev <- c(prev,sequences)
        return(prev)
      }
    }
  )
  
  # Now we need to translate the treatment names into numerical versions, which match
  # the subs table. we can then filter the subs table down to the row in question to
  # get the weighting and compute a weighted average consolidated trace.
  sequences_flattened <- lapply(sequences_flattened, function(this_seq) {
    this_seq$trt_n <- paste(lu_mol[match(this_seq$trt_nam,Description,nomatch = NA),]$Number,collapse = "→")
    this_seq$L1 <- lu_mol[match(this_seq$trt_nam[1],Description,nomatch = NA),]$Number
    
    tl <- length(this_seq$trt_nam)
    atl <- tl - 1
    ncol_actual <- ncol(this_seq$trace_consol)
    ncol_target <- (max_lines * 2)
    
    # if we need to pad with empty columns:
    if (ncol_actual < ncol_target) {
      lines_to_add <- (ncol_target - ncol_actual)/2
      line_txt <- paste0("L",((atl+1):(4)))
      line_txt <- unlist(lapply(line_txt, function(x) paste0(x,c("_on", "_off"))))
      empty_matrix <- matrix(
        data = rep(0,length(line_txt)*nrow(this_seq$trace_consol)),
        nrow = nrow(this_seq$trace_consol),
        ncol = length(line_txt),
        dimnames = list(NULL,line_txt)
      )
      this_seq$trace_consol <- as.data.table(cbind(this_seq$trace_consol,empty_matrix))
      setcolorder(
        this_seq$trace_consol,
        c("L1_on","L1_off","L2_on","L2_off","L3_on","L3_off","L4_on","L4_off","BSC","dead")
      )
      this_seq$trace_consol <- as.matrix(this_seq$trace_consol)
    }
    
    return(list(
      trt_nam = this_seq$trt_nam,
      trt_n = this_seq$trt_n,
      L1 = this_seq$L1,
      trace_consol = this_seq$trace_consol
    ))
  })
  
  # So now we have a "full" consolidated trace for each possible treatment pathway.
  # We can now simply weight each one according to its corresponding weighting in
  # subs, and sum the results by first line treatment!
  w_seq <- lapply(sequences_flattened, function(this_sequence) {
    
    # filter down subs:
    weighting <- subs[trt_n == this_sequence$trt_n,]
    
    stopifnot(nrow(weighting) %in% c(0,1))
    
    if (nrow(weighting) == 0) {
      out <- this_sequence$trace_consol * 0
    } else if (nrow(weighting) == 1) {
      w <- weighting$Adj.proportion.given.line.1
      out <- this_sequence$trace_consol * w
    }
    
    return(list(
      trt_n = this_sequence$trt_n,
      L1 = this_sequence$L1,
      trace_weighted = out
    ))
  })
  
  # With the weighted sequences grouping by first-line treatment (i.e. L1),
  # add up the traces.
  L1        <- unique(unlist(lapply(w_seq, function(x) x$L1),use.names = F))
  tab_L1    <- table(unlist(lapply(w_seq, function(x) x$L1),use.names = F))
  mat_dim   <- dim(w_seq[[1]]$trace_weighted)
  state_nam <- colnames(w_seq[[1]]$trace_weighted)
  
  L1_trt_pos_lab <- structure(1:length(L1), .Names = paste0("L1_",L1))
  
  # Produce weighted average consolidated trace lines per first line therapy
  weighted_traces <- lapply(L1_trt_pos_lab, function(L1_trt) {
    
    # Start from an empty matrix, and then cycle through all of the weighted
    # traces. if L1 matches L1_trt, add it to the pile, if not don't.
    empty_mat <- matrix(
      data = 0,
      nrow = mat_dim[1],
      ncol = mat_dim[2],
      dimnames = list(NULL,state_nam)
    )
    
    # Starting from an empty matrix
    Reduce(
      x = 1:length(w_seq),
      init = empty_mat,
      accumulate = FALSE,
      f = function(prev, trt_seq_n) {
        tr <- w_seq[[trt_seq_n]]
        if (tr$L1 == L1[L1_trt]) {
          prev <- prev + tr$trace_weighted
          return(prev)
        } else {
          return(prev)
        }
      }
    )
    
  })
  
  # produce plots for each of these traces
  plot_list <- lapply(weighted_traces, function(first_line_trt) {
    f_plot_mk_draw_consol_trace(consol_trace = first_line_trt,
                                treatment_names = paste0("L",1:5),
                                tmax = 15)
  })
  
  undiscounted_tis <- do.call(rbind,lapply(weighted_traces,function(x) colSums(x)/52.17857))
  
  undiscounted_tis <- data.table(rownames_to_column(as.data.frame(undiscounted_tis),"L1"))
  
  undiscounted_tis$L1 <- gsub("L1_","",undiscounted_tis$L1)
  
  # Return a list of traces, plots, and lifetime outcomes:
  return(list(
    traces = weighted_traces,
    plots  = plot_list,
    tis    = undiscounted_tis
  ))
  
}

