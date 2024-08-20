#' Produce a summary of the PS model
#' 
#' @param pf_ps patient flow output for the PS mdoel
#' @param lookups the lookup tables i.e., `p$basic$lookup`
#' @param top_line only top line results or full results. should be TRUE in the PSA.
#' 
f_res_summary_ps <- function(pf_ps,lookups, top_line = FALSE) {
  
  lu_mol <- lookups$ipd$mol
  
  # Produce a breakdown table for discounted everything except for life years(undiscounted)
  bdt <- lapply(pf, function(popu) {
    # popu <- pf$pop_1
    rbindlist(
      lapply(popu, function(L1_treatment) {
        # L1_treatment <- popu[[1]]
        
        lys <- as.data.table(L1_treatment$lys)
        setnames(lys,c("PFS_on","PFS_off","PPS_on","PPS_off"),c("ly_PFS_on","ly_PFS_off","ly_PPS_on","ly_PPS_off"))
        qalys <- as.data.table(L1_treatment$qalys$disc)
        setnames(qalys,c("PFS","PPS","AE","AE_PPS"),c("qaly_PFS","qaly_PPS","qaly_AE","qaly_AE_PPS"))
        costs <- as.data.table(L1_treatment$costs$disc)
        setnames(costs,c("drug",
                         "admin",
                         "AE",
                         "substrt_drug_cost",
                         "substrt_admin_cost",
                         "substrt_AE_cost",
                         "mru_preprog",
                         "mru_postprog",
                         "EOL_cost",
                         "prog_cost"),c("cost_drug",
                                        "cost_admin",
                                        "cost_AE",
                                        "cost_substrt_drug",
                                        "cost_substrt_admin",
                                        "cost_substrt_AE",
                                        "cost_mru_preprog",
                                        "cost_mru_postprog",
                                        "cost_EOL",
                                        "cost_prog"))
        
        return(data.table(t(colSums(data.table(lys,qalys,costs)))))
      })
    )
  })
  
  # summarise the breakdown table into final results (tl = top-line)
  tl <- lapply(bdt, function(popu) {
    
    # add up some columns:
    popu[, `:=`(
      costs = cost_drug +
        cost_admin +
        cost_AE +
        cost_substrt_drug +
        cost_substrt_admin +
        cost_substrt_AE +
        cost_mru_preprog +
        cost_mru_postprog +
        cost_EOL +
        cost_prog,
      qalys = qaly_PFS +
        qaly_PPS +
        qaly_AE +
        qaly_AE_PPS,
      lys  = ly_PFS_on +
        ly_PFS_off +
        ly_PPS_on +
        ly_PPS_off
    )]
    
    # get just summary columns:
    return(popu[,list(costs,qalys,lys)])
    
  })
  
  # Get the 1L therapies for each
  trts <- lapply(pf, names)
  trt_n <- lapply(trts, function(pop_oo) {lu_mol[match(pop_oo,RCC_input_desc),]$Number})
  
  bdt <- lapply(structure(1:length(bdt),.Names = names(bdt)), function(pop_n){
    dat <- bdt[[pop_n]]
    dat$L1 <- trt_n[[pop_n]]
    return(dat)
  })
  tl <- lapply(structure(1:length(tl),.Names = names(tl)), function(pop_n){
    dat <- tl[[pop_n]]
    dat$L1 <- trt_n[[pop_n]]
    return(dat)
  })
  
  if(top_line) {
    return(tl)
  } else {
    return(list(
      breakdown = bdt,
      top_line  = tl
    ))
  }
  
}