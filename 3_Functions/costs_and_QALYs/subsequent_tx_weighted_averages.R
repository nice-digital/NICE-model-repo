# code to calculate weighted averages of subsequent therapies

#' an early attempt at computing model averaging for an overall population (see `p$basic$lookup$pop_map`).
#' Not used as averaging was incorporated into results processing in `Model_Structure.R`
#' 
f_calculate_weighted_average_of_subsequent_treatments <- function(proportions = i$R_table_sub_txts_prop_n_costs,
                                                          mol_lookup = i$r_pld_lookup_mol,
                                                          popn_lookup = i$r_overall_lookup_pop) {

  #temp line at present to read in a res object to process
  res <- readRDS("D:/Downloads/markov_model_results_bc_PHNMA.rds")
  
  proportions <- proportions[,c(1:5,11)]
  mol_lookup  <- rbind(mol_lookup, c("BSC","BSC",999))
  results     <- list()
  
  for (i in 1:6) {
    seq_popn <- popn_lookup$Sequencing.population.number[i]
  
    results[[i]] <- list(population = i,
                         seq_popn = seq_popn,
                         proportions = proportions[proportions$Population == paste0("pop",seq_popn),])
    
    results[[i]]$proportions$trt_n <- paste0(mol_lookup[match(results[[i]]$proportions$Line.1, mol_lookup[,2]),3], "→",
                                mol_lookup[match(results[[i]]$proportions$Line.2, mol_lookup[,2]),3], "→",
                                mol_lookup[match(results[[i]]$proportions$Line.3, mol_lookup[,2]),3], "→",
                                mol_lookup[match(results[[i]]$proportions$Line.4, mol_lookup[,2]),3])
    results[[i]]$proportions$trt_n <- gsub("→NA","",results[[i]]$proportions$trt_n)
  
    results[[i]]$proportions$costs_undisc   <- res$undisc[[i]]$res$costs[match(results[[i]]$proportions$trt_n, res$undisc[[i]]$res$trt_n)]
    results[[i]]$proportions$qalys_undisc   <- res$undisc[[i]]$res$qalys[match(results[[i]]$proportions$trt_n, res$undisc[[i]]$res$trt_n)]
    results[[i]]$proportions$ly_undisc      <- res$undisc[[i]]$res$ly[match(results[[i]]$proportions$trt_n, res$undisc[[i]]$res$trt_n)]
    
    results[[i]]$proportions$costs_disc     <- res$disc[[i]]$res$costs[match(results[[i]]$proportions$trt_n, res$disc[[i]]$res$trt_n)]
    results[[i]]$proportions$qalys_disc     <- res$disc[[i]]$res$qalys[match(results[[i]]$proportions$trt_n, res$disc[[i]]$res$trt_n)]
    
    results[[i]]$proportions$p.costs_undisc <- results[[i]]$proportions$Adj.proportion.given.line.1 * results[[i]]$proportions$costs_undisc
    results[[i]]$proportions$p.qalys_undisc <- results[[i]]$proportions$Adj.proportion.given.line.1 * results[[i]]$proportions$qalys_undisc
    results[[i]]$proportions$p.ly_undisc    <- results[[i]]$proportions$Adj.proportion.given.line.1 * results[[i]]$proportions$ly_undisc
    
    results[[i]]$proportions$p.costs_disc   <- results[[i]]$proportions$Adj.proportion.given.line.1 * results[[i]]$proportions$costs_disc
    results[[i]]$proportions$p.qalys_disc   <- results[[i]]$proportions$Adj.proportion.given.line.1 * results[[i]]$proportions$qalys_disc
    
    results[[i]]$proportions <- as.data.table(results[[i]]$proportions)
    
    results[[i]]$summaries <- results[[i]]$proportions[,list(costs_undisc = sum(p.costs_undisc),
                                                             qalys_undisc = sum(p.qalys_undisc),
                                                             ly_undisc    = sum(p.ly_undisc),
                                                             costs_disc   = sum(p.costs_disc),
                                                             qalys_disc   = sum(p.qalys_disc)), 
                                                       by = results[[i]]$proportions$Line.1]
  }
  
  return(results)
}          

