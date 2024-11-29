# output to word final.R
# 
# Some edits were made in a temp script and these should be the final code

#' Overarching function which produces the word document irrespective of model
#' structure. this function routes to either the ST or PS function depending on
#' `model_structure`, expecting the corresponding results set
#' 
#' @param p input set for scenario. contains much of the pertinent information.
#' @param res results object. can be loaded directly from results `.rds` file
#' @param Scenario_name Usually taken from the excel book, named range `R_Scenario_name`
#' @param Scenario_number Usually taken from the excel book, named range `R_Scenario_num`
#' @param price_options Usually taken from the excel book, named range `dd_drug_price_options`
#' @param Run_date default `date()`. Allows custom string to be used instead
#' @param word_template_location location of the word document to start from. does not control styles (see `flextable` package). default `./3_Functions/reporting/empty results doc.docx`
#' @param Word_width_inches paragraph width in word document, used for table column distribution. default `29.7*0.3937=11.6721`
#' @param auto_save whether or not to automatically output a word document with automatic naming. if not, the word document object is returned within the R session
#' @param folder string path to output folder to save word document in
#' 
#' 
f_res_ProduceWordDoc <- function(
    p,
    res, 
    Scenario_name,
    Scenario_number,
    price_options,
    Run_date = date(),
    word_template_location = "./3_Functions/reporting/empty results doc.docx",
    Word_width_inches = 29.7*0.3937,
    auto_save = FALSE,
    verbose = FALSE,
    folder = "./4_Output"
) {
  
  # preamble
  model_structure <- p$basic$structure
  lookups         <- p$basic$lookup
  
  for(pkg in c(
    "shiny","gtools","openxlsx","flexsurv","tidyverse","data.table","heemod",
    "logOfGamma","ggplot2","survminer","officer","officedown","magrittr","Hmisc",
    "future.apply","crosstable","flextable","stringr","BCEA","collapse",
    "scales","Matrix","dplyr")) {require(pkg,character.only = TRUE)}
  
  # Word document tempalte to start from:
  doc_res <- read_docx(word_template_location)
  
  # Add a 1st level header for the overall population:
  doc_res <- doc_res %>% 
    body_add_par(paste0("Results of Model Run in R Scenario name: " , Scenario_name),style = "heading 1")  %>%
    body_add_par(paste0("Date and time run: ", Run_date))
  
  # Producing report tables (state transition model) ------------------------------------------------------
  
  # Make a word document containing results tables using the object res
  
  # Produces a different format depending on model structure
  
  if(model_structure=="State transition") {
    doc_res <- f_res_ProduceWordDoc_ST(
      doc_res           = doc_res,
      res               = res,
      Scenario_name     = Scenario_name,
      Scenario_number   = Scenario_number,
      model_structure   = model_structure,
      pops_to_run       = p$basic$pops_to_run,
      ptchar            = as.data.table(i$R_table_ptchar),
      age_sex_source    = i$dd_age_sex_source,
      patientagesex     = i$R_table_patientagesex,
      lookups           = lookups,
      Run_date          = Run_date,
      Word_width_inches = Word_width_inches,
      verbose = verbose
    )
  } else {
    doc_res <-
      f_res_ProduceWordDoc_PS(
        doc_res           = doc_res,
        res               = res,
        Scenario_name     = Scenario_name,
        Scenario_number   = Scenario_number,
        model_structure   = model_structure,
        pops_to_run       = p$basic$pops_to_run,
        ptchar            = as.data.table(i$R_table_ptchar),
        age_sex_source    = i$dd_age_sex_source,
        patientagesex     = i$R_table_patientagesex,
        lookups           = lookups,
        Run_date          = Run_date,
        Word_width_inches = Word_width_inches,
        verbose = verbose
      )
  }
  
  # If the user wants to save the document with automatic naming, then do so. Otherwise
  # return the updated doc so the user can save it themselves.
  if (auto_save) {
    doc_target <- gsub(" ", "_", paste0(folder, "/Scenario ",Scenario_number,"_",price_options,"_",gsub(":","_",Run_date),".docx"))
    print(doc_res, target = doc_target)
    rm(doc_res)
    cat(paste0("Document automatically saved in location: ", doc_target,"\n"))
    return(NULL)
  } else {
    return(doc_res)
  }
}



# Structure specific functions --------------------------------------------


# ~ State transition (ST) -------------------------------------------------

#' Function to add results tables specific to cabo nivo no adjuvant population word document output
#' 
#' @param doc_res initial document with first header added (see function `f_res_ProduceWordDoc`)
#' @param res results object. can be loaded directly from results `.rds` file
#' @param Scenario_name Usually taken from the excel book, named range `R_Scenario_name`
#' @param Scenario_number Usually taken from the excel book, named range `R_Scenario_num`
#' @param model_structure taken from `p`, in location `p$basic$structure`. Make sure it's correct!
#' @param pops_to_run taken from `p`, in location `p$basic$pops_to_run`
#' @param ptchar taken from `i`, in location `i$R_table_ptchar`, as a data.table. i.e., `as.data.table(i$R_table_ptchar)`
#' @param age_sex_source taken from `i`, in location `i$dd_age_sex_source`
#' @param patientagesex taken from `i`, in location `i$R_table_patientagesex`
#' @param lookups taken from `p` in location `p$basic$lookup`
#' @param Run_date default `date()`. Allows custom string to be used instead. No default as expected to be added in uses of `f_res_ProduceWordDoc`
#' @param Word_width_inches paragraph width in word document, used for table column distribution. No default as expected to be added in uses of `f_res_ProduceWordDoc`
#' 
#' 
#' 
f_res_ProduceWordDoc_ST <- function(
    doc_res,
    res,
    Scenario_name,
    Scenario_number,
    model_structure,
    pops_to_run,
    ptchar,
    age_sex_source,
    patientagesex,
    lookups,
    Run_date,
    Word_width_inches,
    verbose = FALSE
) {
  
  landscape <- prop_section(page_size = page_size(orient = "landscape"))
  portrait  <- prop_section(page_size = page_size(orient = "landscape"))
  
  # Shortened loookups for population and molecule
  lu_pop  <- lookups$pop_map
  lu_rpop <- lookups$ipd$pop
  lu_mol  <- lookups$ipd$mol
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): "
  ))
  
  # make a table by overall population for the summary results
  ft_basic_bop <- do.call(rbind, lapply(structure(
    names(res$wa_summarised), .Names = names(res$wa_summarised)
  ), function(popu_txt) {
    popu <- res$wa_summarised[[popu_txt]]
    popu_n <- as.numeric(gsub("pop_", "", popu_txt))
    
    # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
    rsk_popu_lab <-
      lu_pop[Overall.population.number == popu_n,]$Risk.population
    
    # popu$seq_pop <- seq_popu_lab
    popu$risk_pop <- rsk_popu_lab
    popu$L1 <- lu_mol[match(popu$L1, lu_mol$Number),]$Description
    
    return(popu)
    
  }))
  
  # Now do the same thing for the incremental analysis
  ft_wa_inc <-
    do.call(rbind, lapply(structure(
      names(res$weighted_incremental),
      .Names = names(res$weighted_incremental)
    ), function(popu_txt) {
      if (is.null(res$weighted_incremental[[popu_txt]]$non_dominated)) {
        popu <- as.data.table(res$weighted_incremental[[popu_txt]])
        popu <-
          data.table(
            popu,
            ic = 0,
            iq = 0,
            il = 0,
            ICER = "Dominant"
          )
        popu$str_dom <- NULL
        
      } else {
        popu <-
          as.data.table(res$weighted_incremental[[popu_txt]]$expanded_results)
        popu$ICER[popu$extdom == FALSE] <-
          as.character(paste0("£", round(popu$ICER[popu$extdom == FALSE] , 0)))
        popu$ICER[popu$extdom == TRUE] <- "(ext dominated)"
        popu$str_dom <- NULL
        popu$extdom <- NULL
        popu$r <- NULL
        
      }
      
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        rep(lu_pop[Overall.population.number == popu_n,]$Risk.population, nrow(popu))
      
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  # Create table combining pairwise and incremental ICERs
  setDT(ft_basic_bop)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  ft_basic_bop <- ft_basic_bop[order(risk_pop, costs)] # order by increasing costs
  setDT(ft_wa_inc)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  ft_wa_inc <- ft_wa_inc[, c(1, 3, 2, 4, 5, 6, 7, 8, 9)]
  ft_pairwise <- do.call(rbind, lapply(structure(
    names(res$pairwise_vs_mol), .Names = names(res$pairwise_vs_mol)
  ), function(popu_txt) {
    popu <- res$pairwise_vs_mol[[popu_txt]]
    popu_n <- as.numeric(gsub("pop_", "", popu_txt))
    
    # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
    rsk_popu_lab <-
      lu_pop[Overall.population.number == popu_n,]$Risk.population
    
    # popu$seq_pop <- seq_popu_lab
    popu$risk_pop <- rsk_popu_lab
    
    return(popu)
    
  }))
  setDT(ft_pairwise)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  ft_pairwise$Pairwise_ICER <- ft_pairwise$icer
  ft_pairwise$Pairwise_ICER[is.na(ft_pairwise$icer) != TRUE] <- as.character(paste0("£", round(ft_pairwise$icer[is.na(ft_pairwise$icer) != TRUE] , 0)))
  ft_pairwise[ft_pairwise$icer < 0 & ft_pairwise$iq < 0]$Pairwise_ICER <- "Cabo+nivo dominated"
  ft_pairwise[ft_pairwise$icer < 0 & ft_pairwise$iq > 0]$Pairwise_ICER <- "Cabo+nivo dominant"
  ft_pairwise[ft_pairwise$icer > 0 & ft_pairwise$iq < 0]$Pairwise_ICER <- paste0("SW quadrant ", ft_pairwise[ft_pairwise$icer > 0 & ft_pairwise$iq < 0]$Pairwise_ICER)
  ft_pairwise <- ft_pairwise[, .SD, .SDcols = c("L1", "costs", "qalys", "ly", "Pairwise_ICER", "risk_pop")]
  ft_wa_inc <- merge(ft_pairwise, ft_wa_inc, all.x = TRUE)
  ft_wa_inc <- ft_wa_inc[, c(1, 2, 4, 3, 7, 9, 8, 6, 10, 5)]
  ft_wa_inc <- ft_wa_inc[order(risk_pop, costs)] # order by increasing costs
  ft_wa_inc[is.na(ICER)]$ICER <- "(dominated)"
  
  # Pull out nearest comparators
  comparator_no_allrisk <- ff_closest_comparator(res, "pop_1")
  comparator_no_favrisk <- ff_closest_comparator(res, "pop_2")
  comparator_no_IPrisk  <- ff_closest_comparator(res, "pop_3")
  
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "LY breakdown"
  ),31)
  
  # Create table for LY breakdown
  
  # All risk
  cabo_nivo_LY_allrisk <- res$weighted_model_undisc$pop_1[L1 == 1] %>% select(starts_with("ly"))
  comparator_LY_allrisk <- res$weighted_model_undisc$pop_1[L1 == comparator_no_allrisk] %>% select(starts_with("ly"))
  ft_LY_all_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_LY_allrisk,
    comparator_outcomes = comparator_LY_allrisk,
    comparator_no = comparator_no_allrisk,
    LYorQALY = "LY"
  )
  
  # Fav risk
  
  cabo_nivo_LY_favrisk <- res$weighted_model_undisc$pop_2[L1 == 1] %>% select(starts_with("ly"))
  comparator_LY_favrisk <- res$weighted_model_undisc$pop_2[L1 == comparator_no_favrisk] %>% select(starts_with("ly"))
  ft_LY_fav_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_LY_favrisk,
    comparator_outcomes = comparator_LY_favrisk,
    comparator_no = comparator_no_favrisk,
    LYorQALY = "LY"
  )
  
  # Int/poor risk
  cabo_nivo_LY_IPrisk <- res$weighted_model_undisc$pop_3[L1 == 1] %>% select(starts_with("ly"))
  comparator_LY_IPrisk <- res$weighted_model_undisc$pop_3[L1 == comparator_no_IPrisk] %>% select(starts_with("ly"))
  ft_LY_IP_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_LY_IPrisk,
    comparator_outcomes = comparator_LY_IPrisk,
    comparator_no = comparator_no_IPrisk,
    LYorQALY = "LY"
  )
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "QALY breakdown"
  ),33)
  
  # Create tables for QALY breakdown
  
  # All risk
  
  cabo_nivo_QALY_allrisk    <- res$weighted_model_disc$pop_1[L1 == 1] %>% select(starts_with("qaly"))
  cabo_nivo_AEQALY_allrisk  <- cbind(res$weighted_model_disc$pop_1[L1 == 1] %>% select(starts_with("ae_qaly")),BSC = 0)
  cabo_nivo_QALY_allrisk    <- colSums(rbind(cabo_nivo_QALY_allrisk, cabo_nivo_AEQALY_allrisk, use.names = FALSE))
  cabo_nivo_QALY_allrisk    <- cabo_nivo_QALY_allrisk[c(1:2, 9, 3:8)]
  comparator_QALY_allrisk   <- res$weighted_model_disc$pop_1[L1 == comparator_no_allrisk] %>% select(starts_with("qaly"))
  comparator_AEQALY_allrisk <- cbind(res$weighted_model_disc$pop_1[L1 == comparator_no_allrisk] %>% select(starts_with("ae_qaly")), BSC = 0)
  comparator_QALY_allrisk   <- colSums(rbind(comparator_QALY_allrisk,comparator_AEQALY_allrisk,use.names = FALSE))
  comparator_QALY_allrisk   <- comparator_QALY_allrisk[c(1:2, 9, 3:8)]
  ft_QALY_all_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_QALY_allrisk,
    comparator_outcomes = comparator_QALY_allrisk,
    comparator_no = comparator_no_allrisk,
    LYorQALY = "QALY"
  )
  
  
  # Fav risk
  cabo_nivo_QALY_favrisk    <- res$weighted_model_disc$pop_2[L1 == 1] %>% select(starts_with("qaly"))
  cabo_nivo_AEQALY_favrisk  <- cbind(res$weighted_model_disc$pop_2[L1 == 1] %>% select(starts_with("ae_qaly")),BSC = 0)
  cabo_nivo_QALY_favrisk    <- colSums(rbind(cabo_nivo_QALY_favrisk, cabo_nivo_AEQALY_favrisk,use.names = FALSE))
  cabo_nivo_QALY_favrisk    <- cabo_nivo_QALY_favrisk[c(1:2, 9, 3:8)]
  comparator_QALY_favrisk   <- res$weighted_model_disc$pop_2[L1 == comparator_no_favrisk] %>% select(starts_with("qaly"))
  comparator_AEQALY_favrisk <- cbind(res$weighted_model_disc$pop_2[L1 == comparator_no_favrisk] %>% select(starts_with("ae_qaly")),BSC = 0)
  comparator_QALY_favrisk   <- colSums(rbind(comparator_QALY_favrisk,comparator_AEQALY_favrisk,use.names = FALSE))
  comparator_QALY_favrisk   <- comparator_QALY_favrisk[c(1:2, 9, 3:8)]
  ft_QALY_fav_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_QALY_favrisk,
    comparator_outcomes = comparator_QALY_favrisk,
    comparator_no = comparator_no_favrisk,
    LYorQALY = "QALY"
  )
  
  
  # Int/poor risk
  cabo_nivo_QALY_IPrisk   <- res$weighted_model_disc$pop_3[L1 == 1] %>% select(starts_with("qaly"))
  cabo_nivo_AEQALY_IPrisk <- cbind(res$weighted_model_disc$pop_1[L1 == 1] %>% select(starts_with("ae_qaly")),BSC = 0) 
  cabo_nivo_QALY_IPrisk   <- colSums(rbind(cabo_nivo_QALY_IPrisk,  cabo_nivo_AEQALY_IPrisk, use.names =FALSE))
  cabo_nivo_QALY_IPrisk   <- cabo_nivo_QALY_IPrisk[c(1:2, 9, 3:8)]
  comparator_QALY_IPrisk  <- res$weighted_model_disc$pop_3[L1 == comparator_no_IPrisk] %>% select(starts_with("qaly"))
  comparator_AEQALY_IPrisk <- cbind(res$weighted_model_disc$pop_3[L1 == comparator_no_IPrisk] %>% select(starts_with("ae_qaly")),BSC = 0)
  comparator_QALY_IPrisk  <- colSums(rbind(comparator_QALY_IPrisk,comparator_AEQALY_IPrisk,use.names = FALSE))
  comparator_QALY_IPrisk  <- comparator_QALY_IPrisk[c(1:2, 9, 3:8)]
  ft_QALY_IP_tab <- ff_report_outcomes_breakdown(
    cabo_nivo_outcomes = cabo_nivo_QALY_IPrisk,
    comparator_outcomes = comparator_QALY_IPrisk,
    comparator_no = comparator_no_IPrisk,
    LYorQALY = "QALY"
  )
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Cost breakdown"
  ),34)
  
  # Create tables for overall cost breakdown
  cost_type <- c("drug" , "admin" , "mru" , "eol" , "ae_cost")
  populations <- names(res$weighted_model_disc)
  summary_costs_table <-
    rbindlist(lapply(populations, function(popu) {
      treatments <- res$weighted_model_disc[[popu]]$L1
      rbindlist(lapply(treatments, function(mol) {
        ff_cost_table(
          disc_results = res$weighted_model_disc,
          trt_no       = mol,
          pop          = popu
        )
      }))
    }))
  
  summary_costs_table[, risk_pop := str_replace(Population, "Int/poor", "Intermediate / poor risk")]
  summary_costs_table[, Population := NULL]
  summary_costs_table <- summary_costs_table[order(risk_pop, Total)] # order by increasing costs
  ft_cost_tab <- summary_costs_table %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>%
    as_flextable() %>%
    width(., width = (Word_width_inches / (ncol(
      summary_costs_table
    )))) %>%
    add_header_row(
      top = TRUE,
      values = c("", "1L costs", "Subsequent treatment", "MRU", "", ""),
      colwidths = c(1, 3, 3, 2, 1, 1)
    ) %>%
    theme_box() |>
    set_header_labels(
      values = list(
        Treatment = "Technologies",
        L1_drug = "Drug cost",
        L1_admin = "Admin cost",
        L1_ae = "AE cost",
        subs_drug = "Drug cost",
        subs_admin = "Admin cost",
        subs_ae = "AE cost",
        mru_1L = "1L",
        subs_mru = "Subsequent treatment",
        eol_cost = "EOL cost",
        Total = "Total cost"
      )
    ) %>%
    colformat_double(j = c(2:11),
                     digits = 0,
                     prefix = "£") %>%
    add_footer_lines(
      "Abbreviations: admin, administration; AE, adverse event; EOL, end of life;  MRU, medical resource use"
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
    align(i = ~ !is.na(`Risk population`), align = "left") %>%
    align(i = NULL,
          align = "center",
          part = c("header")) %>%
    bold(i = ~ !is.na(`Risk population`))  %>%
    autofit() %>%
    set_table_properties(layout = "autofit")
  
  
  # produce break downs by population
  intervention_name <- lu_mol[Number == 1]$Description
  
  # all risk
  comparator_name <- lu_mol[Number == comparator_no_allrisk]$Description
  cost_breakdown_2 <- rbind(
    summary_costs_table[risk_pop == "All risk" & Treatment == intervention_name], 
    summary_costs_table[risk_pop == "All risk" & Treatment == comparator_name]
  )
  
  # reshape the data:
  cb2 <- melt.data.table(cost_breakdown_2, id.vars = c("Treatment", "risk_pop"))
  cb2$risk_pop <- NULL
  cb2 <- dcast.data.table(cb2, variable ~ Treatment)
  colnames(cb2) <- c("Type", "Int", "Comp")
  cb2$Inc <- cb2[, Int] - cb2[, Comp]
  cb2$abs <- abs(cb2$Inc)
  cb2$abs[10] <- sum(cb2$abs[1:9])
  cb2$abspercent <-  cb2$abs / cb2$abs[10] * 100
  cb2[, 1] <-
    c(
      "Drug acquisition cost (1L)",
      "Admin cost (1L)",
      "AE cost (1L)",
      "Drug acquisition cost (2L+)",
      "Admin cost (2L+)",
      "AE cost (2L+)",
      "MRU 1L",
      "MRU 2L+",
      "EOL",
      "Total"
    )
  cost_table_2_allrisk <- ff_cost_byrisk_table(cb2, comparator_no_allrisk)
  
  # favourable risk
  comparator_name <- lu_mol[Number == comparator_no_favrisk]$Description
  cost_breakdown_2 <- rbind(
    summary_costs_table[risk_pop == "Favourable risk" & Treatment == intervention_name], 
    summary_costs_table[risk_pop == "Favourable risk" & Treatment == comparator_name]
  )
  cost_breakdown_2 <- as.data.table(x = t(cost_breakdown_2),stringsAsFactors = FALSE)
  cost_breakdown_2 <- cbind(colnames(summary_costs_table), cost_breakdown_2)
  cost_breakdown_2 <- cost_breakdown_2[2:11,]
  cost_breakdown_2[, 2:3] <- lapply(cost_breakdown_2[, 2:3], as.numeric)
  colnames(cost_breakdown_2) <- c("Type", "Int", "Comp")
  cost_breakdown_2$Inc <- cost_breakdown_2[, Int] - cost_breakdown_2[, Comp]
  cost_breakdown_2$abs <- abs(cost_breakdown_2$Inc)
  cost_breakdown_2$abs[10] <- sum(cost_breakdown_2$abs[1:9])
  cost_breakdown_2$abspercent <- cost_breakdown_2$abs / cost_breakdown_2$abs[10] * 100
  cost_breakdown_2[, 1] <-
    c(
      "Drug acquisition cost (1L)",
      "Admin cost (1L)",
      "AE cost (1L)",
      "Drug acquisition cost (2L+)",
      "Admin cost (2L+)",
      "AE cost (2L+)",
      "MRU 1L",
      "MRU 2L+",
      "EOL",
      "Total"
    )
  cost_table_2_favrisk <- ff_cost_byrisk_table(cost_breakdown_2, comparator_no_favrisk)
  
  # int / poor risk
  comparator_name <- lu_mol[Number == comparator_no_IPrisk]$Description
  cost_breakdown_2 <- rbind(
    summary_costs_table[risk_pop == "Intermediate / poor risk" & Treatment == intervention_name], 
    summary_costs_table[risk_pop == "Intermediate / poor risk" & Treatment == comparator_name]
  )
  cost_breakdown_2           <- as.data.table(x = t(cost_breakdown_2), stringsAsFactors = FALSE)
  cost_breakdown_2           <- cbind(colnames(summary_costs_table), cost_breakdown_2)
  cost_breakdown_2           <- cost_breakdown_2[2:11,]
  cost_breakdown_2[, 2:3]    <- lapply(cost_breakdown_2[, 2:3], as.numeric)
  colnames(cost_breakdown_2) <- c("Type", "Int", "Comp")
  cost_breakdown_2$Inc       <- cost_breakdown_2[, Int] - cost_breakdown_2[, Comp]
  cost_breakdown_2$abs       <- abs(cost_breakdown_2$Inc)
  cost_breakdown_2$abs[10]   <- sum(cost_breakdown_2$abs[1:9])
  cost_breakdown_2$abspercent <- cost_breakdown_2$abs / cost_breakdown_2$abs[10] * 100
  cost_breakdown_2[, 1] <-
    c(
      "Drug acquisition cost (1L)",
      "Admin cost (1L)",
      "AE cost (1L)",
      "Drug acquisition cost (2L+)",
      "Admin cost (2L+)",
      "AE cost (2L+)",
      "MRU 1L",
      "MRU 2L+",
      "EOL",
      "Total"
    )
  cost_table_2_IPrisk <- ff_cost_byrisk_table(cost_breakdown_2, comparator_no_IPrisk)
  
  #### Scenario analysis tables
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Scenario tables"
  ),35)
  
  # all risk
  Scenario_table <- ff_scenario_output(res,Scenario_name,comparator_no_allrisk,"pop_1",model_structure)
  Scenario_table_allrisk <- ff_scenario_table(Scenario_table)
  
  # favourable risk
  Scenario_table <- ff_scenario_output(res,Scenario_name,comparator_no_favrisk,"pop_2",model_structure)
  Scenario_table_favrisk <- ff_scenario_table(Scenario_table)
  
  # int/poor risk
  Scenario_table <- ff_scenario_output(res,Scenario_name,comparator_no_IPrisk,"pop_3",model_structure)
  Scenario_table_IPrisk <- ff_scenario_table(Scenario_table)
  
  # base case table
  ft_basecase <- ft_wa_inc %>%
    rename(`Risk population` = risk_pop) %>%
    as_grouped_data(groups = "Risk population") %>%
    as_flextable() %>%
    width(., width = (Word_width_inches / (ncol(ft_wa_inc)))) %>%
    theme_box() |>
    set_header_labels(
      values = list(
        L1 = "Technologies",
        costs = "Costs (£)",
        ly = "LYG",
        qalys = "QALYs",
        ic = "Inc. Costs",
        il = "Inc. LYG",
        iq = "Inc. QALYs",
        Pairwise_ICER = "ICER cabo + nivo vs comparator",
        ICER = "ICER incremental"
      )
    ) %>%
    flextable::colformat_double(j = c(2, 5, 8, 9),
                                digits = 0,
                                prefix = "£") %>%
    flextable::colformat_double(j = c(3, 4, 6, 7), digits = 2) %>%
    add_footer_lines(
      "Abbreviations: ICER, incremental cost-effectiveness ratio; LYG, life-years gained; QALY, quality-adjusted life-year"
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
    align(i = ~ !is.na(`Risk population`), align = "left") %>%
    bold(i = ~ !is.na(`Risk population`))
  
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Severity modifier"
  ),36)
  
  # Severity modifier
  if(p$basic$decision_problem == "cabo+nivo") {pops_to_run <- 1:3} else {pops_to_run <- 1:3}
  
  population_numbers <- if(sum(pops_to_run == 1:3)>0){1:3} else{1:6}
  res$mk$qaly_shortfall_1_to_3 <- lapply(population_numbers, function(npa_pop) {
    
    # npa_pop is overall population, we need to look up risk population from it:
    
    risk_pop_n <- lu_pop[match(npa_pop,lu_pop$Overall.population.number),]$Risk.population.number
    risk_pop <- lu_rpop[match(risk_pop_n,lu_rpop$Number),]$Description  
    
    if (age_sex_source == "Mean") {
      
      # So for this risk population, we need the baseline characteristics:
      bl_chars <- ptchar[Population == risk_pop & Treatment.line == 1,]
      bl_age  <- bl_chars$Starting.age..years..Mean
      bl_male <- 1-bl_chars$Starting...female.Mean
      
    } else {
      
      patient_sex_age_IPD <- as.data.table(patientagesex)
      patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="M","male")
      patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="F","female")
      
      bl_age <- patient_sex_age_IPD[Line ==1]$Age 
      bl_male <- patient_sex_age_IPD[Line ==1]$Gender
      
    }
    
    pna_txt <- names(res$wa_summarised)[npa_pop]
    
    tab <- res$wa_summarised[[pna_txt]][L1 != 1,]
    
    met <- tab[which.max(qalys),]
    
    q_met <- met$qalys
    comp_no_met <- met$L1
    
    out <-    calc_severity_modifier(
      age = bl_age,
      sex = bl_male,
      .patient_level = if(age_sex_source == "Mean") {FALSE} else {TRUE},  
      qalys = q_met,
      .i = i,
      .p = p
    )
    
    out <- cbind(out, SOC = comp_no_met)
    
    return(out)
    
  })
  
  
  severity_table <- data.table(do.call(rbind, res$mk$qaly_shortfall_1_to_3))
  severity_table <-
    cbind(risk_pop = lu_pop$Risk.population[1:3], severity_table)
  severity_table <- rbind(
    severity_table,
    f_res_cabonivo_SevMod(
      res = res,
      oo_pop_string = "Poor / intermediate risk",
      pop_n = 3,
      comp_numb = 5
    )
  )
  severity_table <- rbind(
    severity_table,
    f_res_cabonivo_SevMod(
      res = res,
      oo_pop_string = "Poor / intermediate risk",
      pop_n = 3,
      comp_numb = 8
    )
  )
  setDT(severity_table)[, risk_pop := str_replace(risk_pop, "Favourable risk", "Fav")]
  setDT(severity_table)[, risk_pop := str_replace(risk_pop, "All risk", "All")]
  severity_table$SOC <-
    unlist(lapply(1:nrow(severity_table), function(mol) {
      lu_mol[Number == severity_table$SOC[mol]]$Description
    }))
  ft_severity_mod <- ff_severity_table(severity_table)
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Pairwise results"
  ),37)
  
  # Scenario analysis pairwise results
  ft_all_pairwise <-
    do.call(rbind, lapply(structure(
      names(res$pairwise_vs_mol), .Names = names(res$pairwise_vs_mol)
    ), function(popu_txt) {
      popu <- res$pairwise_vs_mol[[popu_txt]]
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        lu_pop[Overall.population.number == popu_n,]$Risk.population
      
      # popu$seq_pop <- seq_popu_lab
      popu$risk_pop <- rsk_popu_lab
      
      return(popu)
      
    }))
  
  ft_all_pairwise$ICER[is.na(ft_all_pairwise$icer) != TRUE] <-
    as.character(paste0("£", round(ft_all_pairwise$icer[is.na(ft_all_pairwise$icer) != TRUE] , 0)))
  ft_all_pairwise[ft_all_pairwise$icer < 0 &
                    ft_all_pairwise$iq < 0]$ICER <-
    "Cabo+nivo dominated"
  ft_all_pairwise[ft_all_pairwise$icer < 0 &
                    ft_all_pairwise$iq > 0]$ICER <-
    "Cabo+nivo dominant"
  ft_all_pairwise[ft_all_pairwise$icer > 0 &
                    ft_all_pairwise$iq < 0]$ICER <-
    paste0("SW quadrant ", ft_all_pairwise[ft_all_pairwise$icer > 0 &
                                             ft_all_pairwise$iq < 0]$ICER)
  ft_all_pairwise[, icer := NULL]
  setDT(ft_all_pairwise)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  ft_all_pairwise_tab <- ff_scenario_pairwise_table(ft_all_pairwise, Word_width_inches)
  
  
  # Outputting report (state transition) ------------------------------------------------------
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Generating word document from results..."
  ),40)
  
  # Add base case results.
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("Base-case results (ordered in increasing cost)"),
                          bookmark = "tab1") %>%
    body_add_flextable(ft_basecase,
                       align = "left",
                       topcaption = TRUE,
                       split = TRUE) %>%
    
    body_add_break()
  
  
  doc_res <- body_end_section_landscape(doc_res)
  doc_res <- doc_res %>%
    body_add_par("Qualification for the severity modifier", style = "heading 2")  %>%
    body_add_table_legend(paste0("Application of the severity modifier to the base case"),
                          bookmark = "tab2") %>%
    body_add_flextable(
      ft_severity_mod,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_par("Breakdowns by health state and cost category", style = "heading 2")  %>%
    body_add_table_legend(
      paste0(
        "Summary of LY gain by health state (all risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_allrisk]$Description,
        ")"
      ),
      bookmark = "tab3"
    ) %>%
    body_add_flextable(
      ft_LY_all_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of LY gain by health state (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_favrisk]$Description,
        ")"
      ),
      bookmark = "tab4"
    ) %>%
    body_add_flextable(
      ft_LY_fav_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of LY gain by health state (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_IPrisk]$Description,
        ")"
      ),
      bookmark = "tab5"
    ) %>%
    body_add_flextable(ft_LY_IP_tab,
                       align = "left",
                       topcaption = TRUE,
                       split = TRUE) %>%
    body_add_break()
  
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of QALY gain by health state (all risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_allrisk]$Description,
        ")"
      ),
      bookmark = "tab1"
    ) %>%
    body_add_flextable(
      ft_QALY_all_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of QALY gain by health state (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_favrisk]$Description,
        ")"
      ),
      bookmark = "tab6"
    ) %>%
    body_add_flextable(
      ft_QALY_fav_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of QALY gain by health state (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_IPrisk]$Description,
        ")"
      ),
      bookmark = "tab7"
    ) %>%
    body_add_flextable(
      ft_QALY_IP_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  
  doc_res <- body_end_section_portrait(doc_res)
  
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("Summary of costs by health state"),
                          bookmark = "tab8") %>%
    body_add_flextable(ft_cost_tab,
                       align = "left",
                       topcaption = TRUE,
                       split = TRUE) %>%
    body_add_break()
  
  doc_res <- body_end_section_landscape(doc_res)
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of predicted resource use by category of cost (all risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_allrisk]$Description,
        ")"
      ),
      bookmark = "tab1"
    ) %>%
    body_add_flextable(
      cost_table_2_allrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of predicted resource use by category of cost (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_favrisk]$Description,
        ")"
      ),
      bookmark = "tab9"
    ) %>%
    body_add_flextable(
      cost_table_2_favrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      paste0(
        "Summary of predicted resource use by category of cost (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,
        lu_mol[Number == comparator_no_IPrisk]$Description,
        ")"
      ),
      bookmark = "tab10"
    ) %>%
    body_add_flextable(
      cost_table_2_IPrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  
  population_table <-
    as.data.table(cbind(
      risk_pop = lu_pop$Risk.population[1:3],
      pop_labels = c("pop_1", "pop_2", "pop_3")
    ))
  setDT(population_table)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  
  
  
  for (popu in population_table$pop_labels) {
    for (mol in names(res$weighted_trace_plots[[popu]]$plots)) {
      doc_res <- doc_res %>% body_add_figure_legend(
        legend = paste0(
          "Markov trace: ",
          population_table[pop_labels == popu]$risk_pop ,
          ", ",
          lu_mol[Number == str_sub(mol, -1, -1)]$Description
        ),
        bookmark = "fig1"
      ) %>%
        body_add_plot(print(res$weighted_trace_plots[[popu]]$plots[mol]), width = 6) %>%
        body_add_par(
          paste0(
            "Abbreviations: L1, 1st line; L2, 2nd line; L3, 3rd line; L4, 4th line; L5, 5th line"
          ),
          style = "Table footnote"
        ) %>% body_add_break()
      
    }
  }
  
  doc_res <- doc_res  %>% body_add_break() %>%
    body_add_par("Cost-effectiveness acceptability frontiers", style = "heading 2")  %>%
    body_add_par(
      paste0(
        "Cost-effectiveness acceptability frontiers are presented for all non-dominated treatments for each of the risk groups"
      )
    )  %>%
    body_add_break() %>%
    body_add_figure_legend(
      legend = paste0("Cost-effectiveness acceptability frontier – all risk"),
      bookmark = "fig2"
    ) %>%
    body_add_plot(print(res$weighted_incremental$pop_1$p), height = 4) %>%
    body_add_par(paste0("Abbreviations: QALYs, quality-adjusted life-years"),
                 style = "Table footnote")
  
  doc_res <- doc_res %>%
    body_add_figure_legend(
      legend = paste0("Cost-effectiveness acceptability frontier – favourable risk"),
      bookmark = "fig3"
    ) %>%
    body_add_plot(print(res$weighted_incremental$pop_2$p), height = 4) %>%
    body_add_par(paste0("Abbreviations: QALYs, quality-adjusted life-years"),
                 style = "Table footnote")
  
  doc_res <- doc_res %>%
    body_add_figure_legend(
      legend = paste0(
        "Cost-effectiveness acceptability frontier – intermediate / poor risk"
      ),
      bookmark = "fig4"
    ) %>%
    body_add_plot(print(res$weighted_incremental$pop_3$p), height = 4) %>%
    body_add_par(paste0("Abbreviations: QALYs, quality-adjusted life-years"),
                 style = "Table footnote")   %>%
    
    body_add_break()
  
  doc_res <- body_end_section_portrait(doc_res)
  
  
  doc_res <- doc_res %>%
    body_add_par("Scenario analysis style tables", style = "heading 1")  %>%
    body_add_table_legend(legend = paste0("Scenario analysis - all risk"),
                          bookmark = "tab11") %>%
    body_add_flextable(
      Scenario_table_allrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(legend = paste0("Scenario analysis - favourable risk"),
                          bookmark = "tab12") %>%
    body_add_flextable(
      Scenario_table_favrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      legend = paste0("Scenario analysis - intermediate / poor risk"),
      bookmark = "tab13"
    ) %>%
    body_add_flextable(
      Scenario_table_IPrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      legend = paste0("Scenario analysis pairwise comparison table"),
      bookmark = "tab14"
    ) %>%
    body_add_flextable(
      ft_all_pairwise_tab ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    body_add_break()
  
  
  # return the updated document:
  return(doc_res)
}



# ~ Partitioned survival (PS) ---------------------------------------------

#' Function to add results tables specific to cabo nivo no adjuvant population word document output
#' 
#' @param doc_res initial document with first header added (see function `f_res_ProduceWordDoc`)
#' @param res results object. can be loaded directly from results `.rds` file
#' @param Scenario_name Usually taken from the excel book, named range `R_Scenario_name`
#' @param Scenario_number Usually taken from the excel book, named range `R_Scenario_num`
#' @param model_structure taken from `p`, in location `p$basic$structure`. Make sure it's correct!
#' @param pops_to_run taken from `p`, in location `p$basic$pops_to_run`
#' @param ptchar taken from `i`, in location `i$R_table_ptchar`, as a data.table. i.e., `as.data.table(i$R_table_ptchar)`
#' @param age_sex_source taken from `i`, in location `i$dd_age_sex_source`
#' @param patientagesex taken from `i`, in location `i$R_table_patientagesex`
#' @param lookups taken from `p` in location `p$basic$lookup`
#' @param Run_date default `date()`. Allows custom string to be used instead. No default as expected to be added in uses of `f_res_ProduceWordDoc`
#' @param Word_width_inches paragraph width in word document, used for table column distribution. No default as expected to be added in uses of `f_res_ProduceWordDoc`
#' 
#' 
f_res_ProduceWordDoc_PS <- function(
    doc_res,
    res,
    Scenario_name,
    Scenario_number,
    model_structure,
    pops_to_run,
    ptchar,
    age_sex_source,
    patientagesex,
    lookups,
    Run_date,
    Word_width_inches,
    verbose = FALSE
) {
  
  
  landscape <- prop_section(page_size = page_size(orient = "landscape"))
  portrait  <- prop_section(page_size = page_size(orient = "landscape"))
  
  # Shortened loookups for population and molecule
  lu_pop <- lookups$pop_map
  lu_mol <- lookups$ipd$mol
  
  # Producing report tables (PartSA) ------------------------------------------------------
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): "
  ))
  
  # Make LY table
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "LYs"
  ),col_num = 31)
  
  PartSA_Lys <-
    do.call(rbind, lapply(structure(names(res$ly), .Names = names(res$ly)), function(popu_txt) {
      popu <- as.data.table(res$ly[[popu_txt]])
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      popu <-  data.table(L1 = rownames(res$ly[[popu_txt]]), popu)
      popu$L1 <-
        lu_mol[match(popu$L1, lu_mol$RCC_input_desc), ]$Description
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  setDT(PartSA_Lys)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  PartSA_Lys$Total <- rowSums(PartSA_Lys[, 2:5])
  
  PartSA_LYs_table <- ff_PartSALY_table(PartSA_Lys)
  
  # Make QALYs table
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "QALYs"
  ),col_num = 33)
  
  PartSA_QALYs <-
    do.call(rbind, lapply(structure(
      names(res$disc_qaly), .Names = names(res$disc_qaly)
    ), function(popu_txt) {
      popu <- as.data.table(res$disc_qaly[[popu_txt]])
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      popu <-
        data.table(L1 = rownames(res$disc_qaly[[popu_txt]]), popu)
      popu$L1 <-
        lu_mol[match(popu$L1, lu_mol$RCC_input_desc), ]$Description
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  setDT(PartSA_QALYs)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  PartSA_QALYs$Total <- rowSums(PartSA_QALYs[, 2:5])
  
  PartSA_QALYs_table <- ff_PartSAQALY_table(PartSA_QALYs)
  
  # Make costs table
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Costs"
  ),col_num = 34)
  
  PartSA_costs <-
    do.call(rbind, lapply(structure(
      names(res$disc_cost), .Names = names(res$disc_cost)
    ), function(popu_txt) {
      popu <- as.data.table(res$disc_cost[[popu_txt]])
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      popu <-
        data.table(L1 = rownames(res$disc_cost[[popu_txt]]), popu)
      popu$L1 <-
        lu_mol[match(popu$L1, lu_mol$RCC_input_desc), ]$Description
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  setDT(PartSA_costs)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  PartSA_costs$Total <- rowSums(PartSA_costs[, 2:11])
  
  PartSA_costs_table <- ff_PartSAcost_table (PartSA_costs)
  
  # Make results table
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Results tables"
  ),col_num = 35)
  
  PartSA_wa <-
    do.call(rbind, lapply(structure(
      names(res$incremental), .Names = names(res$incremental)
    ), function(popu_txt) {
      if (is.null(res$incremental[[popu_txt]]$non_dominated)) {
        popu <- as.data.table(res$incremental[[popu_txt]])
        popu <-
          data.table(
            popu,
            ic = 0,
            iq = 0,
            il = 0,
            ICER = "Dominant"
          )
        popu$str_dom <- NULL
        
      } else {
        popu <- as.data.table(res$incremental[[popu_txt]]$expanded_results)
        popu$ICER[popu$extdom == FALSE] <-
          as.character(paste0("£", round(popu$ICER[popu$extdom == FALSE] , 0)))
        popu$ICER[popu$extdom == TRUE] <- "(ext dominated)"
        popu$str_dom <- NULL
        popu$extdom <- NULL
        popu$r <- NULL
        
      }
      
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  PartSA_wa <- PartSA_wa[, c(2, 3, 1, 4, 5, 7, 6, 8, 9)]
  
  
  
  PartSA_totals <-
    do.call(rbind, lapply(structure(
      names(res$tables$top_line), .Names = names(res$tables$top_line)
    ), function(popu_txt) {
      popu <- as.data.table(res$tables$top_line[[popu_txt]])
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <-
        rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      popu$L1 <- lu_mol[match(popu$L1, lu_mol$Number), ]$Description
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  PartSA_totals <- PartSA_totals[order(risk_pop, costs)] # order by increasing costs
  
  PartSA_totals$L1_risk  <- paste(PartSA_totals$L1, PartSA_totals$risk_pop)
  PartSA_wa$L1_risk  <- paste(PartSA_wa$L1, PartSA_wa$risk_pop)
  
  PartSA_results <- merge(PartSA_totals, PartSA_wa, all.x = TRUE)
  PartSA_results[is.na(ICER)]$ICER <- "(dominated)"
  
  PartSA_results <-  PartSA_results[, c(4, 1, 3, 2, 7, 8, 9, 10, 5)]
  
  PartSA_results <-
    PartSA_results[order(risk_pop, costs)] # order by increasing costs
  
  
  PartSA_Pairwise <-
    do.call(rbind, lapply(structure(
      names(res$tables$top_line), .Names = names(res$tables$top_line)
    ), function(popu_txt) {
      popu   <- f_res_ICER_pairwiseVsoneTrt(res$tables$top_line[[popu_txt]], 1, lu_mol)
      popu_n <- as.numeric(gsub("pop_", "", popu_txt))
      
      # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
      rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n, ]$Risk.population, nrow(popu))
      
      # popu$seq_pop <- seq_popu_lab
      popu <- cbind(popu, risk_pop = rsk_popu_lab)
      
      return(popu)
      
    }))
  
  
  PartSA_Pairwise$Pairwise_ICER[is.na(PartSA_Pairwise$icer) != TRUE] <- as.character(paste0(
    "£",
    round(PartSA_Pairwise$icer[is.na(PartSA_Pairwise$icer) != TRUE] , 0))
  )
  
  
  PartSA_Pairwise[PartSA_Pairwise$icer < 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER <- "Cabo+nivo dominated"
  PartSA_Pairwise[PartSA_Pairwise$icer < 0 & PartSA_Pairwise$iq > 0]$Pairwise_ICER <- "Cabo+nivo dominant"
  PartSA_Pairwise[PartSA_Pairwise$icer > 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER <- paste0(
    "SW quadrant ", 
    PartSA_Pairwise[PartSA_Pairwise$icer > 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER
  )
  
  PartSA_Pairwise_Scen <- PartSA_Pairwise
  PartSA_Pairwise <-  PartSA_Pairwise[, c(4, 9, 10)]
  
  
  PartSA_results <- merge(PartSA_results, PartSA_Pairwise, all.x = TRUE)
  
  PartSA_results <- PartSA_results[, c(1:8, 10, 9)]
  
  PartSA_results[ICER == 0]$ICER <- "(dominated)"
  PartSA_results[, 6] <- as.numeric(unlist(PartSA_results[, 6][[1]]))
  PartSA_results[, 7] <- as.numeric(unlist(PartSA_results[, 7][[1]]))
  PartSA_results[, 8] <- as.numeric(unlist(PartSA_results[, 8][[1]]))
  
  PartSA_results <- PartSA_results[order(risk_pop, costs)] # order by increasing costs
  setDT(PartSA_results)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  
  PartSA_results_tab <- ff_PartSAresults_table(PartSA_results)
  
  #### Scenario analysis tables
  comparator_no_allrisk <- ff_closest_comparator_PartSA(res, "pop_1")
  comparator_no_favrisk <- ff_closest_comparator_PartSA(res, "pop_2")
  comparator_no_IPrisk  <- ff_closest_comparator_PartSA(res, "pop_3")
  
  # all risk
  Scenario_table <- ff_scenario_output(res,Scenario_name,comparator_no_allrisk,"pop_1",model_structure)
  Scenario_table_allrisk <- ff_scenario_table(Scenario_table)
  
  # favourable risk
  Scenario_table <-ff_scenario_output(res,Scenario_name,comparator_no_favrisk,"pop_2",model_structure)
  Scenario_table_favrisk <- ff_scenario_table(Scenario_table)
  
  # int/poor risk
  
  Scenario_table <- ff_scenario_output(res,Scenario_name,comparator_no_IPrisk,"pop_3",model_structure)
  Scenario_table_IPrisk <- ff_scenario_table(Scenario_table)
  
  # Scenario analysis pairwise results
  
  setDT(PartSA_Pairwise_Scen)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
  
  PartSA_Pairwise_Scen <- PartSA_Pairwise_Scen[, c(4, 1, 2, 3, 5, 6, 7, 10, 9)]
  PartSA_Pairwise_Scen$ICER <- PartSA_Pairwise_Scen$Pairwise_ICER
  PartSA_Pairwise_Scen$Pairwise_ICER <- NULL
  
  ft_all_pairwise_tab <- ff_scenario_pairwise_table(PartSA_Pairwise_Scen, Word_width_inches)
  
  
  # Outputting report (PartSA) ------------------------------------------------------
  
  if (verbose) f_misc_colcat(paste0(
    "Word output for ",
    Scenario_name,
    ". (scen #", Scenario_number,"): ",
    "Generating word document..."
  ),col_num = 36)
  
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("PartSA life years"),
                          bookmark = "tab1") %>%
    body_add_flextable(
      PartSA_LYs_table,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("PartSA QALYs"),
                          bookmark = "tab2") %>%
    body_add_flextable(
      PartSA_QALYs_table,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("PartSA costs"),
                          bookmark = "tab3") %>%
    body_add_flextable(
      PartSA_costs_table,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(paste0("PartSA results (ordered in increasing cost)"),
                          bookmark = "tab4") %>%
    body_add_flextable(
      PartSA_results_tab,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_par("Scenario analysis style table", style = "heading 1")  %>%
    body_add_table_legend(legend = paste0("Scenario analysis style table - all risk"),
                          bookmark = "tab5") %>%
    body_add_flextable(
      Scenario_table_allrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      legend = paste0("Scenario analysis style table - favourable risk"),
      bookmark = "tab6"
    ) %>%
    body_add_flextable(
      Scenario_table_favrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      legend = paste0("Scenario analysis style table - intermediate / poor risk"),
      bookmark = "tab7"
    ) %>%
    body_add_flextable(
      Scenario_table_IPrisk ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    )  %>%
    body_add_break()
  
  doc_res <- doc_res %>%
    body_add_table_legend(
      legend = paste0("Scenario analysis pairwise comparison table"),
      bookmark = "tab8"
    ) %>%
    body_add_flextable(
      ft_all_pairwise_tab ,
      align = "left",
      topcaption = TRUE,
      split = TRUE
    ) %>%
    body_add_break()
  
  return(doc_res)
  
}



# TESTING AREA ------------------------------------------------------------


# Testing area, set to not run but can be used to assign variables and run through
# the function definition line by line:

if (FALSE) {
  # Load all libraries
  library(shiny, quiet = TRUE)   
  library(gtools, quiet = TRUE)
  library(openxlsx, quiet = TRUE)
  library(flexsurv, quiet = TRUE)
  library(tidyverse, quiet = TRUE)
  library(data.table, quiet = TRUE)
  library(heemod, quiet = TRUE)
  library(logOfGamma, quiet = TRUE)
  library(ggplot2, quiet = TRUE)
  library(survminer, quiet = TRUE)
  library(officer, quiet = TRUE)
  library(officedown, quiet = TRUE)
  library(magrittr, quiet = TRUE)
  library(Hmisc, quiet = TRUE)
  library(future.apply, quiet = TRUE)
  library(crosstable, quiet = TRUE)
  library(flextable, quiet = TRUE)
  library(stringr, quiet = TRUE)
  library(BCEA, quiet = TRUE)
  library(collapse, quiet = TRUE)
  library(scales, quiet = TRUE)
  library(Matrix, quiet = TRUE)
  library(dplyr, quiet = TRUE)
  
  # Load all functions:
  source("./3_Functions/excel/extract.R")
  source("./3_Functions/sequencing/sequences.R")
  source("./3_Functions/survival/Survival_functions.R")
  source("./3_Functions/survival/other_cause_mortality.R")
  source("./3_Functions/survival/treatment_effect_waning.R")
  source("./3_Functions/misc/other.R")
  source("./3_Functions/misc/shift_and_pad.R")
  source("./3_Functions/misc/cleaning.R")
  source("./3_Functions/misc/nesting.R")
  source("./3_Functions/misc/discounting.R")
  source("./3_Functions/misc/qdirichlet.R")
  source("./3_Functions/misc/plotting.R")
  source("./3_Functions/misc/structure.R")
  source("./3_Functions/utility/age_related.R")
  source("./3_Functions/costs_and_QALYs/utility_processing.R")
  source("./3_Functions/adverse_events/AE_steps.R")
  source("./3_Functions/costs_and_QALYs/cost_processing.R")
  source("./3_Functions/markov/markov.R")
  source("./3_Functions/patient_flow/overarching.R")
  source("./3_Functions/patient_flow/partitioned_survival.R")
  source("./3_Functions/patient_flow/markov.R")
  source("./3_Functions/patient_flow/drug_costs.R")
  source("./3_Functions/patient_flow/hcru_costs.R")
  source("./3_Functions/patient_flow/qalys.R")
  source("./3_Functions/patient_flow/ae.R")
  source("./3_Functions/results/incremental_analysis.R")
  source("./3_Functions/results/model_averaging.R")
  source("./3_Functions/results/partitioned_survival.R")
  source("./3_Functions/misc/severity_modifier.R")
  source("./3_Functions/results/results_tables.R")
  
  # Load required variables from file:
  i   <- readRDS("./2_Scripts/standalone scripts/QC/i.rds")
  p   <- readRDS("./2_Scripts/standalone scripts/QC/p.rds")
  res <- readRDS("~/Downloads/PATT-Pathways-RCC/4_Output/results_scenario_RevisedBC_2Aug2023.rds")
  
  # Test the function
  f_res_ProduceWordDoc(
    p                      = p,
    res                    = res,
    Scenario_name          = "Test scenario (base-case)",
    Scenario_number        = i$R_Scenario_num,
    price_options          = i$dd_drug_price_options,
    Run_date               = date(),
    word_template_location = "./3_Functions/reporting/empty results doc.docx",
    Word_width_inches      = 29.7*0.3937,
    auto_save              = TRUE
  )
}


