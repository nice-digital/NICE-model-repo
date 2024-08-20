#' Function to compute "patient flow" i.e. either a Markov trace or state transition 
#' This top level function simply routes to the appropriate underlying pf generating
#' function. Either for a state transition or ps model.
#' 
#' @param p object p
#' @param struct the model structure. must be one of the dropdown options for named range i$dd_model_struct
#' @param verbose whether to print more to the console or not
#' @param plots whether to generate plots in the output
#' @param just_pop just run a subset of the risk populations. default is to run them ALL (a lot of iteration!)
#' @param just_nlines just run a subset of the pathways according to how many ACTIVE treatment lines. Markov only.
#' 
f_pf_computePF <-
  function(p,
           struct = "State transition",
           verbose = FALSE,
           plots = FALSE,
           just_pop = NULL,
           just_nlines = NULL,
           just_seq = NULL) {
    
  
  # list the populations to simulate:
  if (!is.null(just_pop)) {
    if(0%in% just_pop) stop ("this is overall population, not risk population, it can't be 0.")
    overall_pops <- structure(paste0("pop_",just_pop),.Names=paste0("pop_",just_pop))
  } else {
    overall_pops <- structure(
      paste0("pop_",p$basic$lookup$pop_map$Overall.population.number),
      .Names=paste0("pop_",p$basic$lookup$pop_map$Overall.population.number)
    )
  }
  
  # map pops for risk for the demographics:
  rpop <- paste0("pop_",p$basic$lookup$pop_map[match(as.numeric(gsub("pop_","",overall_pops)),p$basic$lookup$pop_map$Overall.population.number),]$Risk.population.number)
  
  
  
  # If the dropdown is broken, then return an error
  if(!struct %in% c("State transition", "Partitioned survival")) {
    stop ("structure must be at least one of the dropdown options in Excel (named range List_model_structures)")
  }
  
  # Otherwise it has to be one of the below:
  if (struct == "State transition") {
    out <- f_pf_computePF_mk(
      pops          = overall_pops,
      basic         = p$basic,
      demo          = p$demo$live[which(names(p$demo$live) %in% rpop)],
      sequences     = p$seq,
      survival      = p$surv,
      costs         = list(per_cycle = p$costs$mk, one_off = p$costs$oneoff_mk),
      util          = list(hsuv  = p$util$mk, gpop = p$util$gpop),
      ae            = list(one_off = p$ae$duration, per_cycle = p$ae$mk$per_cycle, approach = p$ae$approach),
      eff_table     = p$releff$excel_table,
      verbose       = verbose,
      include_plots = plots,
      just_nlines   = just_nlines,
      just_seq      = just_seq
    )
  } else {
    out <- f_pf_computePF_ps(
      pops          = overall_pops,
      basic         = p$basic,
      p             = p,
      cyclecosts    = p$costs$mk,
      oneoffcosts   = p$costs$oneoff,
      util          = list(hsuv  = p$util$mk, gpop = p$util$gpop),
      ae            = list(one_off = p$ae$duration, per_cycle = p$ae$mk$per_cycle, approach = p$ae$approach),
      substrt       = p$substrt
    )
  }
  return(out)
}