#' Function which creates an empty version of the parameters list p and populates
#' it using the inputs from i
#' 
#' @param i raw extracted inputs from Excel
#' 
f_misc_param_generate_p <- function(i) {
  
  p <- list(
    basic = list(
      th   = ceiling(i$ui_time_horizon * 365.25 / 7), 
      th_y = i$ui_time_horizon, 
      cl_d = i$i_cycle_length_weeks*7, 
      cl_w = i$i_cycle_length_weeks, 
      cl_y = i$i_cycle_length_weeks*7/365.25, 
      discQ = i$ui_disc_qaly,
      discC = i$ui_disc_cost,
      structure = str_trim(i$dd_model_struct)
    ),
    demo = list(
      table = data.table(i$R_table_ptchar)
    ), # patient demographic data
    seq   = list(), # treatment sequences. same for both model structures.
    surv  = list(), # Survival extrapolations (if needed for iteration, i.e. in part surv model). same for both model structures
    releff = list(), # Relative efficacy network, for use populating the disease model. same for both model structures
    costs  = list(
      mk = list(),
      ps = list(),
      settings = list(
        subsTx = data.table(i$R_table_sub_txts_prop_n_costs)
      )
    ), # All drug inputs after they've been processed and tidied up
    util  = list(
      mk = list(),
      ps = list()
    ), # All inputs to apply to the disease model
    ae    = list(
      mk = list(),
      ps = list()
    ), # Inputs to generate AE matrices AC and AQ
    misc  = list(
      mk = list(),
      ps = list(),
      plot = list(
        xlim_survplots_yr = 20
      )
    )
  )
  
  p$basic$n_cyc <- ceiling(p$basic$th) # rounding up to a whole number
  p$basic$t_cyc <- rep(0:p$basic$n_cyc)
  p$basic$t_yr  <- p$basic$t_cyc*7/365.25
  
  # discount factors - base-case edition. These will get replaced for scenarios affecting
  # discount rates or the time horizon
  p$basic$discFacQ <- f_discFac(p$basic$discQ,p$basic$t_yr)
  p$basic$discFacC <- f_discFac(p$basic$discC,p$basic$t_yr)
  
  return(p)
  
}


#' function to pull out the HRs being used in the network. Only for time-invariant HRs of course
f_releff_extract_all_HRs <- function(network) {
  lapply(network,function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        do.call(
          rbind,
          lapply(mol, function(tr) {
            unlist(lapply(tr, function(endp) {
              endp$hr
            }))
          })
        )
      })
    })
  })
}



