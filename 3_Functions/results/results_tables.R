#' function to compute pairwise ICERs against a specific treatment. expects
#' the model-averaged 1L treatment sequencing results. 
#' 
#' @param tab `res$mk$wa_summarised` for a specific population (top line table with L1 costs qalys ly columns)
#' @param int_L1 the treatment to treat as the "intervention" to calculate pairwise ICERs
#' @param lookup the lookup tables, usually `p$basic$lookup`
#' 
f_res_ICER_pairwiseVsoneTrt <- function(tab, int_L1, lookup) {
  
  int <- tab[L1 == int_L1,]
  
  comp <- tab[L1 != int_L1,]
  
  ICERs <- rbindlist(
    lapply(1:nrow(comp), function(comparator) {
      out <- list()
      
      out$ic <- int$costs - comp[comparator,]$costs
      out$iq <- int$qalys - comp[comparator,]$qalys
      out$il <- int$ly - comp[comparator,]$ly
      out$icer <- out$ic / out$iq
      
      data.table(t(unlist(out)))
      
    })
  )
  
  ICERs <- rbind(
    data.table(t(structure(rep(NA,4),.Names=colnames(ICERs)))),
    ICERs
  )
  
  results_table <- cbind(tab,ICERs)
  results_table$L1 <- lookup[match(results_table$L1, lookup$Number),]$Description
  
  return(results_table)
}


#' Function to route to lower level functions per model structure in order to
#' process patient flow (i.e., pf) and return results.
#' 
#' @param pf patient flow object resulting from `f_pf_computePF`
#' @param structure `State transition` or `Partitioned survival`
#' @param p the full p list for this model run
#' @param detail_level 1 is top line tables only, 2 is 
#' 
f_res_compute_results <- function(pf, structure, p, detail_level = 1,vs_mol=1, no_active_lines = 4) {
  stopifnot(structure %in% c("State transition","Partitioned survival"))
  if (structure == "State transition") {
    f_res_compute_results_mk(pf,p,detail_level,vs_mol, no_active_lines)
  } else {
    f_res_compute_results_ps(pf,p,detail_level)
  }
}


#' Results summarizer for the partitioned survival model
#' 
#' @param pf patient flow object for markov model
#' @param p full p list for this run
#' @param detail level 1-5
#' 
f_res_compute_results_ps <- function(pf,p,detail_level = 1) {
  
  stopifnot(detail_level %in% 1:5)
  out <- list()
  
  out$tables <- f_res_summary_ps(
    pf_ps = pf,
    lookups = p$basic$lookup,
    top_line = FALSE
  )
  
  if(detail_level == 1) return(out$tables$top_line)
  
  pop_labs <- structure(paste0("pop_",1:length(out$tables$top_line)),.Names=paste0("pop_",1:length(out$tables$top_line)))
  
  out$incremental <- lapply(pop_labs, function(popu) {
    dat <- out$tables$top_line[[popu]]
    f_res_mk_incremental(
      res_d = dat,
      res_ud = dat,
      produce_plot = TRUE,
      lu_mol = p$basic$lookup$ipd$mol,
      output_weighted = "Yes"
    )
  })
  
  
  if(detail_level == 2) return(list(
    summary = out$tables$top_line,
    incremental = lapply(out$incremental, function(x) {
      if(!"non_dominated" %in% names(x)) {
        x
      } else {
        x$non_dominated
      }
    })
  ))
  
  if(detail_level >= 3) {
    out$ly <- lapply(pf, function(popu) {
      do.call(rbind,
              lapply(popu, function(L1_treatment) {
                sapply(L1_treatment$lys, sum)
              }))
      
    })
    out$disc_qaly <- lapply(pf, function(popu) {
      do.call(rbind,
              lapply(popu, function(L1_treatment) {
                sapply(L1_treatment$qalys$disc, sum)
              }))
      
    })
    out$disc_cost <- lapply(pf, function(popu) {
      do.call(rbind,
              lapply(popu, function(L1_treatment) {
                sapply(L1_treatment$costs$disc, sum)
              }))
      
    })
  }
  
  # Sample code is provided below which allows the production of PartSA outputs to view as plots or in the console
  # Highest detail is plots:
  if(detail_level > 4){
    out$sr_plots <- lapply(pf, function(popu){
      lapply(popu, function(L1_trt) {
        L1_trt$sr_plot
      })
    })
  }
  
  return(out)
  
}
  
  
  
  



#' Results summarizer for the state transition model
#' 
#' @param pf patient flow object for markov model
#' @param p full p list for this run
#' @param detail level 1-5
#' @param vs_mol decision problem molecule for pairwise ICERs
#' 
f_res_compute_results_mk <- function(pf,p,detail_level = 1, vs_mol = 1, no_active_lines) {
  
  stopifnot(detail_level %in% 1:5)
  
  out <- list()
  
  # In any case, we need to run the main results calculator by treatment sequence
  out$undisc = f_pf_mk_summary(
    pf_list = pf,
    disc_undisc = "undisc",
    lookups = p$basic$lookup,
    full_breakdown = TRUE,
    breakdown = TRUE,
    ypc = p$basic$cl_y
    )
  out$disc = f_pf_mk_summary(
    pf_list = pf,
    disc_undisc = "disc",
    lookups = p$basic$lookup,
    full_breakdown = TRUE,
    breakdown = TRUE,
    ypc = p$basic$cl_y
   )
  
  # If detail level is 4 or above return per sequence results. per sequence
  # incrementals aren't needed for weighting analysis so optional here for higher
  # detail levels
  if (detail_level >= 4) {
    out$incremental <- lapply(structure(names(out$disc),.Names=names(out$disc)),function(popu) {
      f_res_mk_incremental(
        res_d  = out$disc[[popu]]$res,
        res_ud = out$undisc[[popu]]$res,
        produce_plot = detail_level > 4,
        no_active_lines = no_active_lines,
        output_weighted = "No")
    })
  }
  
  pop_index <- structure(1:length(out$undisc), .Names = paste0("pop_", 1:length(out$undisc)))
  
  # compute weighted models - we need to do this for ALL levels of detail
  out$weighted_model_disc <- lapply(pop_index, function(overall_pop) {
    f_res_wa_model(
      res_obj = out,
      pop_oo = overall_pop,
      subs = p$costs$settings$subsTx,
      ptchar = p$demo$table,
      lookups = p$basic$lookup,
      disc = TRUE,
      no_active_lines = no_active_lines,
      max_lines = p$basic$R_maxlines+1
    )
  })
  out$weighted_model_undisc <- lapply(pop_index, function(overall_pop) {
    f_res_wa_model(
      res_obj = out,
      pop_oo = overall_pop,
      subs = p$costs$settings$subsTx,
      ptchar = i$R_table_ptchar,
      lookups = p$basic$lookup,
      disc = FALSE,
      no_active_lines = no_active_lines,
      max_lines = p$basic$R_maxlines+1
    )
  })
  
  # only for full detail do we want the trace plots:
  if (detail_level > 4) {
    out$weighted_trace_plots <- lapply(pop_index, function(pop_n) {
      f_pf_wa_sr_plots(
        pf_list_pop = pf[[pop_n]],
        oo_pop_n = pop_n,
        subs = p$costs$settings$subsTx,
        lookups = p$basic$lookup,
        max_lines = p$basic$R_maxlines+1,
        no_active_lines = no_active_lines
      )
    })
  }
  
  # the top line result is this - the weighted average results
  out$wa_summarised <- lapply(pop_index, function(popu) {
    pop_txt <- names(pop_index)[popu]
    f_res_sum_weighted_model(
      rd  = out$weighted_model_disc[[pop_txt]],
      rud = out$weighted_model_undisc[[pop_txt]]    
    )
  })
  
  # If it's detail level 1, JUST give me these tables
  if (detail_level == 1) return(out$wa_summarised)
  
  
  out$weighted_incremental <- lapply(out$wa_summarised, function(popu) {
    popu <- popu[order(popu$costs),]
    f_res_mk_incremental(
      res_d = popu,
      res_ud = popu,
      produce_plot = detail_level > 4,
      lu_mol = p$basic$lookup$ipd$mol,
      no_active_lines = no_active_lines, 
      output_weighted = "Yes"
    )
  })
  
  
  # Detail level 2 is just top line table and non-dominated incremental results:
  if (detail_level == 2) {
    return(list(
      summary = out$wa_summarised,
      incremental = lapply(out$weighted_incremental, function(x) x$non_dominated)
    ))
  }
  
  if(detail_level >= 3) {
    out$pairwise_vs_mol <- lapply(out$wa_summarised[1:3], function(oo_pop) {
      f_res_ICER_pairwiseVsoneTrt(oo_pop,vs_mol,p$basic$lookup$ipd$mol)
    })
  }
  
  
  # For the higher detail levels we already have breakdown tables and per-sequence
  # results done, depending on detail_level
  return(out)
}


#' function to create LY and QALY output breakdowns for automated reporting
#' 
#' @param tab `cabo_nivo_outcomes` for the weighted population
#' @param tab `comparator_outcomes` for the weighted population - this should use the same comparator as defined in comparator_no
#' @param comparator_no the treatment to treat as the "comparator" to  produce the results breakdown - this should be the closest comparator on the incremental analysis
#' @param LYorQALY whether to report LYs or QALYs


ff_report_outcomes_breakdown <-
  function(cabo_nivo_outcomes,
           comparator_outcomes,
           comparator_no,
           LYorQALY) {
    
    stopifnot(LYorQALY %in% c("LY","QALY"))
    
    # Make a character vector for the names to assign depending on settings:
    if (LYorQALY == "LY") {
      nam_dat <- gsub("ly_","",names(cabo_nivo_outcomes))
    } else {
      nam_dat <- gsub("qaly_","",names(cabo_nivo_outcomes))
    }
    
    # Put suffix expansion:
    nam_dat <- gsub("_on",": on treatment",nam_dat)
    nam_dat <- gsub("_off",": off treatment",nam_dat)
    
    # Note: column names can't start with a number in R so we're doing that after
    
    names(cabo_nivo_outcomes)  <- nam_dat
    names(comparator_outcomes) <- nam_dat
    
    ft_outcomes_breakdown <-
      data.table(
        Health_state = nam_dat,
        cabo_nivo = as.numeric(cabo_nivo_outcomes),
        comp = as.numeric(comparator_outcomes)
      )
    
    
    ft_outcomes_breakdown <- ft_outcomes_breakdown[order(Health_state),]
    ft_outcomes_breakdown <- ft_outcomes_breakdown[c(2:nrow(ft_outcomes_breakdown),1),]
    ft_outcomes_breakdown$Health_state <- gsub("L","",ft_outcomes_breakdown$Health_state)
    ft_outcomes_breakdown$Health_state <- sub("(.*)(\\d)", "\\1\\2L", ft_outcomes_breakdown$Health_state)
    
    
    ft_outcomes_breakdown$inc <- ft_outcomes_breakdown$cabo_nivo - ft_outcomes_breakdown$comp
    ft_outcomes_breakdown$absinc <- abs(ft_outcomes_breakdown$inc)
    ft_outcomes_breakdown$percentabs <- ft_outcomes_breakdown$absinc / sum(ft_outcomes_breakdown$absinc) * 100
    
    # Add in death row and totals row:
    ft_outcomes_breakdown <- rbindlist(list(
      ft_outcomes_breakdown, 
      lapply(1:ncol(ft_outcomes_breakdown), function(x) {
        if (x == 1) {
          "Death"
        } else {
          0
        }
      }),
      lapply(1:ncol(ft_outcomes_breakdown), function(x) {
        if (x == 1) {
          "Total"
        } else {
          sum(ft_outcomes_breakdown[,x,with = FALSE])
        }
      })
      
    ))
    
    ft_outcomes_breakdown <- flextable(ft_outcomes_breakdown)
    
    if (LYorQALY == "LY") {
      ft_outcomes_tab <- ff_LY_table(ft_outcomes_breakdown, comparator_no)
    } else {
      ft_outcomes_tab <-
        ff_QALY_table(ft_outcomes_breakdown, comparator_no)
    }
    
    return(ft_outcomes_tab)
    
  }

#' function to create LY output breakdown for automated reporting 
#' 
#' @param tab `ft_LY_breakdown` for the weighted population based upon output of ff_report_outcomes_breakdown
#' @param comparator_no the treatment to treat as the "comparator" to  produce the results breakdown - this should be the closest comparator on the incremental analysis


ff_LY_table <- function(ft_LY_breakdown, comparator_no) {
  
  ft_LY_tab  <- ft_LY_breakdown  %>%
    theme_box() |>
    set_header_labels(
      values = list(
        Health_state = "Health state",
        cabo_nivo = "LY Cabozantinib plus nivolumab (X)",
        comp = paste0("LY ", p$basic$lookup$ipd$mol[Number == comparator_no]$Description, " (Y)"),
        inc = "Increment",
        absinc = "Absolute increment",
        percentabs = "% absolute increment"
      )
    ) %>%
    flextable::colformat_double(j = c(6), digits = 0, suffix = "%") %>%
    flextable::colformat_double(j = c(2:5), digits = 3) %>%
    add_footer_lines("Abbreviations: 1L, 1st line; 2L, 2nd line; 3L, 3rd line; 4L, 4th line; BSC, best supportive care; LY, life years; vs, versus") %>%
    add_footer_lines("Discrepancies in sums due to rounding errors: totals shown are calculated on unrounded numbers") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |>
    bold(bold = TRUE, part = "header") %>%
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9,  part = c("footer")) %>% 
    style(i=11, pr_t = fp_text_default(bold = TRUE)) %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_LY_tab)
  
}

#' function to create QALY output breakdown for automated reporting 
#' 
#' @param tab `ft_QALY_breakdown` for the weighted population based upon output of ff_report_outcomes_breakdown
#' @param comparator_no the treatment to treat as the "comparator" to  produce the results breakdown - this should be the closest comparator on the incremental analysis


ff_QALY_table <- function(ft_QALY_breakdown, comparator_no) {
  ft_QALY_tab  <- ft_QALY_breakdown  %>%
    theme_box() |>
    set_header_labels(
      values = list(
        Health_state = "Health state",
        cabo_nivo = "QALY Cabozantinib plus nivolumab (X)",
        comp = paste0("QALY ", p$basic$lookup$ipd$mol[Number == comparator_no]$Description, " (Y)"),
        inc = "Increment",
        absinc = "Absolute increment",
        percentabs = "% absolute increment"
      )
    ) %>%
    flextable::colformat_double(j = c(6),
                                digits = 0,
                                suffix = "%") %>%
    flextable::colformat_double(j = c(2:5), digits = 3) %>%
    add_footer_lines(
      "Abbreviations: 1L, 1st line; 2L, 2nd line; 3L, 3rd line; 4L, 4th line; BSC, best supportive care; LY, life years; vs, versus"
    ) %>%
    add_footer_lines(
      "Discrepancies in sums due to rounding errors: totals shown are calculated on unrounded numbers"
    ) %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |>
    bold(bold = TRUE, part = "header") %>%
    fontsize(i = NULL,
             size = 10,
             part = c("header")) %>%
    fontsize(i = NULL,
             size = 10,
             part = c("body")) %>%
    fontsize(i = NULL,
             size = 9,
             part = c("footer")) %>% 
    style(i=11, pr_t = fp_text_default(bold = TRUE)) %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_QALY_tab)
  
}


#' function to create cost breakdown for automated reporting 
#' 
#' @param tab `disc_results` discounted results from the res file for the weighted population 
#' @param trt_no treatment number to produce the table for
#' @param pop population to produce the table for
#' @param cost_type list of the cost types to include

ff_cost_table <- function(disc_results, trt_no, pop, cost_type = c("drug" , "admin" , "mru" , "eol" , "ae_cost")) {
  
  cost_inputs <- disc_results[[pop]][L1 == trt_no,] %>% dplyr::select(contains(cost_type))
  
  subs_drug  <- sum(cost_inputs %>% dplyr::select(starts_with("drug")))    - cost_inputs$drug_L1
  subs_admin <- sum(cost_inputs %>% dplyr::select(starts_with("admin")))   - cost_inputs$admin_L1
  subs_ae    <- sum(cost_inputs %>% dplyr::select(starts_with("ae_cost"))) - cost_inputs$ae_cost_L1
  mru_1L     <- cost_inputs$mru_on_L1 + cost_inputs$mru_off_L1
  subs_mru   <- sum(cost_inputs %>%  select(starts_with("mru"))) - mru_1L
  
  return(data.table(Population = p$basic$lookup$pop_map[Overall.population.number == as.numeric(gsub("\\D", "", pop))]$Risk.population, 
                    Treatment  = p$basic$lookup$ipd$mol[Number == trt_no]$Description, 
                    L1_drug    = cost_inputs$drug_L1, 
                    L1_admin   = cost_inputs$admin_L1, 
                    L1_ae      = cost_inputs$ae_cost_L1,
                    subs_drug  = subs_drug, 
                    subs_admin = subs_admin, 
                    subs_ae    = subs_ae, 
                    mru_1L     = mru_1L, 
                    subs_mru   = subs_mru, 
                    eol_cost   = cost_inputs$eol_L5, 
                    Total      = sum(cost_inputs)))
}

#' function to create cost breakdown for automated reporting - combined results by risk status
#' 
#' @param tab `cost_breakdown_2` cost breakdown from combined output of ff_cost_table
#' @param comparator_no the treatment to treat as the "comparator" to  produce the results breakdown - this should be the closest comparator on the incremental analysis

ff_cost_byrisk_table <- function(cost_breakdown_2, comparator_no) {
  
  ft_cost2_tab  <- flextable(cost_breakdown_2)  %>%
    theme_box() |>
    set_header_labels(
      values = list(
        Type = "Item",
        Int = "Cost Cabozantinib plus nivolumab (X)",
        Comp = paste0("Cost ", p$basic$lookup$ipd$mol[Number == comparator_no]$Description, " (Y)"),
        Inc = "Increment",
        abs = "Absolute increment",
        abspercent = "% absolute increment"
      )
    ) %>%
    flextable::colformat_double(j = c(6), digits = 0, suffix = "%") %>%
    flextable::colformat_double(j = c(2:5), digits = 0, prefix = "£") %>%
    add_footer_lines("Abbreviations: 1L, 1st line; 2L, 2nd line; 2L+, 2nd line-plus; admin, administration; AE, adverse event; EOL, end of life; MRU, medical resource use") %>%
    add_footer_lines("Discrepancies in sums due to rounding errors: totals shown are calculated on unrounded numbers") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |>
    bold(bold = TRUE, part = "header") %>%
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9,  part = c("footer")) %>% 
    style(i=10, pr_t = fp_text_default(bold = TRUE)) %>%
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_cost2_tab)
  
}

#' function to create scenario analysis format table for easy copy pasting
#' 
#' @param tab `Scenario_table` table of scenario analysis results


ff_scenario_table <- function(Scenario_table) {
  
  ft_scenario_tab  <- flextable(Scenario_table)  %>%
    theme_box() |>
    set_header_labels(
      values = list(
        Scenario_name = "Scenario",
        next_best = "Next best comparator*",
        ic = "Incremental costs",
        iq = "Incremental QALYs",
        ICER = "ICER (£/QALY)"
      )
    ) %>%
    flextable::colformat_double(j = c(3,5), digits = 0, prefix = "£") %>%
    flextable::colformat_double(j = c(4), digits = 3) %>%
    add_footer_lines("Abbreviations: 1L, 1st line; 2L, 2nd line; 3L, 3rd line; 4L, 4th line; AEs, adverse events; AUC, area under the curve;  axi, axitinib; BSC, best supportive care; evero, everolimus; FP, fractional polynomial; ICER, incremental cost-effectiveness ratio;  IO, immune-oncology; IV, intravenous; KM, Kaplan-Meier; OS, overall survival; PD, progressed disease; PFS, progression free survival; PH, proportional hazards; PPS, post-progression survival; QALYs, quality adjusted life years; RDI, relative dosing intensity; RWE, real world evidence; tivo, tivozanib; TKI, tyrosine kinase inhibitor; TTD, time to discontinuation; TTP, time to progression ") %>%
    add_footer_lines("*Next best comparator defined as next most efficient non-dominated comparator.") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |>
    bold(bold = TRUE, part = "header") %>%
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9,  part = c("footer")) %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_scenario_tab)
  
}

ff_scenario_pairwise_table <- function(ft_all_pairwise, Word_width_inches) {
  
  ft_scenario_tab  <- ft_all_pairwise  %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>% 
    as_flextable() %>% 
    width(., width = (Word_width_inches/(ncol(ft_all_pairwise)))) %>% 
    theme_box() |> 
    set_header_labels(
      values = list(
        L1 = "Technologies",
        costs = "Costs (£)",
        qalys = "QALYs",
        ly = "LYG",
        ic = "Inc. Costs",
        iq = "Inc. QALYs",
        il = "Inc. LYG",
        ICER = "ICER cabo + nivo vs comparator"
      )
    ) %>%
    flextable::colformat_double(j=c(2,5), digits = 0, prefix = "£") %>%
    flextable::colformat_double(j=c(3,4,6,7), digits = 2) %>%
    add_footer_lines("Abbreviations: ICER, incremental cost-effectiveness ratio; inc. incremental; LYG, life-years gained; QALY, quality-adjusted life-year") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
    bold( bold = TRUE, part="header") %>% 
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9, part = c("footer")) %>%
    align(i = ~ !is.na(`Risk population`), align = "left") %>% 
    bold(i = ~ !is.na(`Risk population`))  %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_scenario_tab)
  
}


#' function to identify the closest comparator for a population
#' 
#' @param `res` results file containing the weighted incremental analysis
#' @param `pop` population to run the function for - needs to be input in format for overall population e.g. "pop_1"

ff_closest_comparator <- function(res,pop){
  
  if ( "Cabozantinib plus nivolumab" %in% res$weighted_incremental[[pop]]$L1) {
    
    # if cabo plus nivo is dominant then 
    # take the treatment with the highest QALYs
    QALYs <- res$weighted_model_disc[[pop]] %>% 
      select(starts_with("qaly"))
    total_QALYs <- rowSums(QALYs,na.rm = TRUE)
    names(total_QALYs) <- res$weighted_model_disc[[pop]]$L1
    n <- length(total_QALYs)
    second_QALY <- sort(total_QALYs,partial=n-1)[n-1]
    
    comparator_name <- names(total_QALYs[match(second_QALY, total_QALYs)])
    
    comparator_no <- p$basic$lookup$ipd$mol[Description == comparator_name_lookup]$Number
  
  } else if( "Cabozantinib plus nivolumab" %in% res$weighted_incremental[[pop]]$non_dominated$L1) {
    
    # if cabo plus nivo appears in the incremental results non dominated list
    # take the treatment it appears after in the list as the most relevant comparator
    
    position_cabonivo <- which("Cabozantinib plus nivolumab" == res$weighted_incremental[[pop]]$non_dominated$L1)[[1]]
    
    if (position_cabonivo == 1) { 
      
      comparator_name_lookup <- res$weighted_incremental[[pop]]$non_dominated$L1[position_cabonivo+1]
      
      } else {
      
       comparator_name_lookup <- res$weighted_incremental[[pop]]$non_dominated$L1[position_cabonivo-1]
    }
       
    comparator_no <- p$basic$lookup$ipd$mol[Description == comparator_name_lookup]$Number
    
  } else {
    
    # take the most effective treatment in the non dominated list
    
    position_comparator <- length(res$weighted_incremental[[pop]]$non_dominated$L1)
    comparator_name_lookup <- res$weighted_incremental[[pop]]$non_dominated$L1[position_comparator]
    comparator_no <- p$basic$lookup$ipd$mol[Description == comparator_name_lookup]$Number
    
  }
  
  return(comparator_no)
  
}

ff_closest_comparator_PartSA <- function(res,pop){
  
  if ( "Cabozantinib plus nivolumab" %in% res$incremental[[pop]]$L1) {
    
    # if cabo plus nivo is dominant
    # take the treatment with the highest QALYs
    QALYs <- as.data.table(res$disc_qaly[[pop]])
    total_QALYs <- rowSums(QALYs,na.rm = TRUE)
    names(total_QALYs) <- rownames(res$disc_qaly[[pop]])
    n <- length(total_QALYs)
    second_QALY <- sort(total_QALYs,partial=n-1)[n-1]
    
    comparator_name <- names(total_QALYs[match(second_QALY, total_QALYs)])
    
    comparator_no <- p$basic$lookup$ipd$mol[RCC_input_desc == comparator_name ]$Number
    
  } else if( "Cabozantinib plus nivolumab" %in% res$incremental[[pop]]$non_dominated$L1) {
    
    # if cabo plus nivo appears in the incremental results non dominated list
    # take the treatment it appears after in the list as the most relevant comparator
    
    position_cabonivo <- which("Cabozantinib plus nivolumab" == res$incremental[[pop]]$non_dominated$L1)[[1]]
    
    if (position_cabonivo == 1) {
      
      comparator_name_lookup <- res$incremental[[pop]]$non_dominated$L1[position_cabonivo+1]
      
    } else {
    
    comparator_name_lookup <- res$incremental[[pop]]$non_dominated$L1[position_cabonivo-1]
    
    }
    
    comparator_no <- p$basic$lookup$ipd$mol[Description == comparator_name_lookup]$Number
    
  } else {
    
    # take the most effective treatment in the non dominated list
    
    position_comparator <- length(res$incremental[[pop]]$non_dominated$L1)
    comparator_name_lookup <- res$incremental[[pop]]$non_dominated$L1[position_comparator]
    comparator_no <- p$basic$lookup$ipd$mol[Description == comparator_name_lookup]$Number
    
  }
  
  return(comparator_no)
  
}

f_res_cabonivo_SevMod <- function(res, oo_pop_string, pop_n, comp_numb) {
  # note that in calc_severity_modifier, the argument sex is the proportion male
  
  if (i$dd_age_sex_source == "Mean") {
    
    # So for this risk population, we need the baseline characteristics:
    bl_chars <- i$R_table_ptchar[Population == oo_pop_string & Treatment.line == 1,]
    bl_age  <- bl_chars$Starting.age..years..Mean
    bl_male <- 1-bl_chars$Starting...female.Mean
    
  } else {
    
    patient_sex_age_IPD <- as.data.table(i$R_table_patientagesex)
    patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="M","male")
    patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="F","female")
    
    bl_age <- patient_sex_age_IPD[Line ==1]$Age 
    bl_male <- patient_sex_age_IPD[Line ==1]$Gender
    
  }
  
  pna_txt <- names(res$wa_summarised)[pop_n]
  
  tab <- res$wa_summarised[[pna_txt]][L1 != 1,]
  
  met <- tab[L1 == comp_numb]
  
  q_met <- met$qalys
  comp_no_met <- met$L1
  
  out <-    calc_severity_modifier(
    age = bl_age,
    sex = bl_male,
    .patient_level = if(i$dd_age_sex_source == "Mean") {FALSE} else {TRUE},  
    qalys = q_met,
    .i = i,
    .p = p
  )
  
  out <- cbind(risk_pop = p$basic$lookup$pop_map$Risk.population[pop_n], out, SOC = comp_numb)
  
return(out)
  
}

ff_severity_table <- function(severity_table) {
  ft_severity_tab  <- flextable(severity_table)  %>%
    theme_box() |>
    set_header_labels(
      values = list(
        risk_pop = "Risk group",
        qaly_soc = "SOC QALYs",
        qaly_gpop = "Gen pop QALYs",
        abs_sf = "Abs SF",
        prop_sf= "Prop SF",
        modifier = "Modifier",
        SOC = "Treatment considered SOC"
      )
    ) %>%
    flextable::colformat_double(j = c(6), digits = 1) %>%
    flextable::colformat_double(j = c(2:5), digits = 3) %>%
    flextable::colformat_char(j = c(1:7)) %>%
    add_footer_lines(
      "Abbreviations: Abs, absolute; Fav, favourable; Gen, general; Int, intermediate; pop, population; Prop, proportional; QALYs, quality adjusted life years; SF, shortfall; SOC, standard of care"
    ) %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |>
    bold(bold = TRUE, part = "header") %>%
    fontsize(i = NULL,
             size = 10,
             part = c("header")) %>%
    fontsize(i = NULL,
             size = 10,
             part = c("body")) %>%
    fontsize(i = NULL,
             size = 9,
             part = c("footer")) %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  return(ft_severity_tab)
  
}

ff_PartSALY_table <- function(PartSA_Lys) {
  
  ft_PartSALY_tab  <- PartSA_Lys  %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>% 
    as_flextable() %>% 
    width(., width = (Word_width_inches/(ncol(PartSA_Lys)))) %>% 
    theme_box() |> 
    set_header_labels(
      values = list(
        L1 = "Technologies",
        PFS_on = "PFS on treatment",
        PFS_off = "PFS off treatment",
        PPS_on = "PPS on treatment",
        PPS_off = "PPS off treatment",
        Total = "Total"
      )
    ) %>%
    flextable::colformat_double(j=c(2:6), digits = 2) %>%
    add_footer_lines("Abbreviations: LYG, life years gained; PartSA, partitioned survival analysis; PFS, progression free survival; PPS, post progression survival") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
    bold( bold = TRUE, part="header") %>% 
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9, part = c("footer")) %>%
    align(i = ~ !is.na(`Risk population`), align = "left") %>% 
    bold(i = ~ !is.na(`Risk population`))  %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_PartSALY_tab)
  
}

ff_PartSAQALY_table <- function(PartSA_QALYs) {
  
  ft_PartSAQALY_tab  <- PartSA_QALYs  %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>% 
    as_flextable() %>% 
    width(., width = (Word_width_inches/(ncol(PartSA_QALYs)))) %>% 
    theme_box() |> 
    set_header_labels(
      values = list(
        L1 = "Technologies",
        PFS = "PFS",
        PPS = "PPS",
        AE = "1L AEs",
        AE_PPS = "AEs PPS",
        Total = "Total"
      )
    ) %>%
    flextable::colformat_double(j=c(2:6), digits = 2) %>%
    add_footer_lines("Abbreviations: AE, adverse event; PartSA, partitioned survival analysis; PFS, progression free survival; PPS, post progression survival; QALYs, quality adjusted life years") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
    bold( bold = TRUE, part="header") %>% 
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9, part = c("footer")) %>%
    align(i = ~ !is.na(`Risk population`), align = "left") %>% 
    bold(i = ~ !is.na(`Risk population`))  %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_PartSAQALY_tab)
  
}

ff_PartSAcost_table <- function(PartSA_costs) {
  
  ft_PartSAcost_tab  <- PartSA_costs  %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>% 
    as_flextable() %>% 
    width(., width = (Word_width_inches/(ncol(PartSA_costs)))) %>% 
    add_header_row(top = TRUE, values = c("","1L costs", "Subsequent treatment", "MRU","","",""), colwidths = c(1,3,3,2,1,1,1)) %>%
    theme_box() |> 
    set_header_labels(
      values = list(
        L1 = "Technologies",
        drug = "Drug cost",
        admin = "Admin cost",
        AE = "AE cost",
        substrt_drug_cost  = "Drug cost",
        substrt_admin_cost = "Admin cost",
        substrt_AE_cost = "AE cost",
        mru_preprog = "Pre-progression cost",
        mru_postprog = "Post-progression cost",
        EOL_cost = "EOL cost",
        prog_cost = "On progression cost",
        Total = "Total"
      )
    ) %>%
    flextable::colformat_double(j=c(2:12), digits = 0, prefix = "£") %>%
    add_footer_lines("Abbreviations: admin, administration; AE, adverse event; EOL, end of life; MRU, medical resource use; PartSA, partitioned survival analysis") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
    bold( bold = TRUE, part="header") %>% 
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9, part = c("footer")) %>%
    align(i = ~ !is.na(`Risk population`), align = "left") %>% 
    bold(i = ~ !is.na(`Risk population`))  %>% 
    align(i = 1, j = NULL, align = "center", part = "header")  %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_PartSAcost_tab)
  
}

ff_PartSAresults_table <- function(PartSA_results) {
  
  ft_PartSA_results_tab  <- PartSA_results  %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>% 
    as_flextable() %>% 
    width(., width = (Word_width_inches/(ncol(PartSA_results)))) %>% 
    add_header_row(top = TRUE, values = c("","Total", "Incremental"), colwidths = c(1,3,5)) %>%
    theme_box() |> 
    set_header_labels(
      values = list(
        L1 = "Technologies",
        costs = "Costs",
        lys = "LYs",
        qalys = "QALYs",
        ic = "Costs",
        il = "LYG",
        iq = "QALYs",
        Pairwise_ICER = "ICER cabo + nivo vs comparator (£/QALY)",
        ICER = "ICER incremental (£/QALY)"
        )
    ) %>%
    flextable::colformat_double(j=c(3,4,6,7), digits = 2, na_str = "") %>%
    flextable::colformat_double(j=c(2,5), digits = 0, prefix = "£", na_str = "") %>%
    add_footer_lines("Abbreviations: ext, extended; ICER, incremental cost-effectiveness ratio; LYG, life year gained; PartSA, partitioned survival analysis; QALY, quality adjusted life year; SW, south west; vs, versus") %>%
    # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
    bold( bold = TRUE, part="header") %>% 
    fontsize(i = NULL, size = 10, part = c("header")) %>%
    fontsize(i = NULL, size = 10, part = c("body")) %>%
    fontsize(i = NULL, size = 9, part = c("footer")) %>%
    align(i = ~ !is.na(`Risk population`), align = "left") %>% 
    bold(i = ~ !is.na(`Risk population`))  %>% 
    autofit() %>% 
    set_table_properties(layout = "autofit") 
  
  
  return(ft_PartSA_results_tab)
  
}

ff_scenario_output <- function(res, Scenario_name, closest_comparator, pop, structure) {
  
  if (structure=="Partitioned survival") {location = "incremental"} else {location = "weighted_incremental"}
  
  Scenario_table <- data.table()
  
  Scenario_table$Scenario <- Scenario_name
  
  Scenario_table$next_best <- p$basic$lookup$ipd$mol[Number == closest_comparator]$Description
  
  if(sum(str_detect(res[[location]][[pop]]$non_dominated$L1,"Cabozantinib plus nivolumab"))>0) {
    
    position_cabonivo <- which("Cabozantinib plus nivolumab" == res[[location]][[pop]]$non_dominated$L1)[[1]]
    
    if(nrow(res[[location]][[pop]]$non_dominated) == 1) {
      
      Scenario_table$ic <- 0
      Scenario_table$iq <- 0
      Scenario_table$ICER <- "Cabo+nivo dominant vs all"
      
    } else if (position_cabonivo == 1){ 
      
      Scenario_table$ic <- res[[location]][[pop]]$non_dominated[2]$ic
      Scenario_table$iq <- res[[location]][[pop]]$non_dominated[2]$iq
      Scenario_table$ICER <- paste0("£",round(res[[location]][[pop]]$non_dominated[2]$ICER,0)," SW quadrant comp vs cabo+nivo")
      
      } else {
      
      Scenario_table$ic <- res[[location]][[pop]]$non_dominated[L1 == "Cabozantinib plus nivolumab"]$ic
      Scenario_table$iq <- res[[location]][[pop]]$non_dominated[L1 == "Cabozantinib plus nivolumab"]$iq
      Scenario_table$ICER <- res[[location]][[pop]]$non_dominated[L1 == "Cabozantinib plus nivolumab"]$ICER
    } 
  }  else if  ("Cabozantinib plus nivolumab" %in% res[[location]][[pop]]$L1) { 
    
    Scenario_table$ic <- 0
    Scenario_table$iq <- 0
    Scenario_table$ICER <- "Cabo+nivo dominant vs all"
    
    } else {
    
    Scenario_table$ic <- 0
    Scenario_table$iq <- 0
    Scenario_table$ICER <- "Cabo+nivo dominated"
    
  }
  
  return(Scenario_table)
  
}