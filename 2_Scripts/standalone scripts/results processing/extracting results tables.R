#Ed's code for extracting results tables



# res <-  readRDS(rstudioapi::selectFile(
#   path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files",
#   existing = TRUE))

#to get total cost, QALYS and LY:

dir.create("4_Output/tables_for_report_PAS_TO BE DELETED/basecase")
for (n in 1:3) {
  summary <- f_res_sum_weighted_model(
    rd  = res$weighted_model_disc[[n]],
    rud = res$weighted_model_undisc[[n]]
  )
  summary$drug <- i$r_pld_lookup_mol$Description[match(summary$L1, i$r_pld_lookup_mol$Number)]
  summary <- summary[,c(1,5,2,4,3)]
  summary <- summary[order(summary$costs),]
  summary[,c("ly","qalys")] <- round(summary[,c("ly","qalys")],3)
  write.csv(summary, file = paste0("4_Output/tables_for_report_PAS_TO BE DELETED/basecase/summary",n,".csv"), row.names = F)
  
  #to get incremental analyses (note excludes dominated):
  nondom <- res$weighted_incremental[[n]]$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(1,9,3,2,4,5,7,6,8)]
  nondom[,c(4:5,7:8)] <- round(nondom[,c(4:5,7:8)],3)
  nondom[,c(3,6,9)] <- round(nondom[,c(3,6,9)],2)
  write.csv(nondom, file = paste0("4_Output/tables_for_report_PAS_TO BE DELETED/basecase/nondom",n,".csv"), row.names = F)
  
  LY <- res$weighted_model_undisc[[n]][,c(1,45:53)]
  LY$drug <-  i$r_pld_lookup_mol$Description[match(LY$L1, i$r_pld_lookup_mol$Number)]
  LY <- LY[,c(1,11,2:3, 5:10, 4)]
  LY <- t(LY)
  write.csv(LY, file = paste0("4_Output/tables_for_report_PAS_TO BE DELETED/basecase/LY",n,".csv"), row.names = T)
  
  QALY <- res$weighted_model_disc[[n]][,c(1,28:44)]
  QALY$drug <-  i$r_pld_lookup_mol$Description[match(QALY$L1, i$r_pld_lookup_mol$Number)]
  QALY <- as.data.frame(QALY)
  for (m in seq(2,14,4)) {
    QALY[,m] <- QALY[,m] + QALY[,(m+1)]
  }
  QALY <- t(QALY[,c(1,19,seq(2,17,2),18)])
  write.csv(QALY, file = paste0("4_Output/tables_for_report_PAS_TO BE DELETED/basecase/QALY",n,".csv"), row.names = T)
  
  cost <- res$weighted_model_disc[[n]][,c(1:27)]
  cost$drug <-  i$r_pld_lookup_mol$Description[match(cost$L1, i$r_pld_lookup_mol$Number)]
  cost <- cost[,c(1,28,2:27)]
  cost$substx_drug <- apply(cost[,c("drug_L2","drug_L3","drug_L4","drug_L5")],1,sum)
  cost <- cost[,-c("drug_L2","drug_L3","drug_L4","drug_L5")]
  cost$substx_admin <- apply(cost[,c("admin_L2","admin_L3","admin_L4","admin_L5")],1,sum)
  cost <- cost[,-c("admin_L2","admin_L3","admin_L4","admin_L5")]
  cost$substx_ae <- apply(cost[,c("ae_cost_L2","ae_cost_L3","ae_cost_L4","ae_cost_L5")],1,sum)
  cost <- cost[,-c("ae_cost_L2","ae_cost_L3","ae_cost_L4","ae_cost_L5")]
  cost$mru_L1 <- apply(cost[,c("mru_on_L1","mru_off_L1")],1,sum)
  cost <- cost[,-c("mru_on_L1","mru_off_L1")]
  cost$mru_substx <- apply(cost[,c("mru_on_L2","mru_off_L2",
                                   "mru_on_L3","mru_off_L3",
                                   "mru_on_L4","mru_off_L4",
                                   "mru_on_L5","mru_off_L5")],1,sum)
  cost <- cost[,-c("mru_on_L2","mru_off_L2",
                                   "mru_on_L3","mru_off_L3",
                                   "mru_on_L4","mru_off_L4",
                                   "mru_on_L5","mru_off_L5")]
  cost <- cost[,c(1:5,7:11,6)]
  cost$total <- apply(cost[,3:11],1,sum)
  cost <- cost[order(cost$total),]
  write.csv(cost, file = paste0("4_Output/tables_for_report_PAS_TO BE DELETED/basecase/cost",n,".csv"), row.names = F)
  
}

dir.create("4_Output/tables_for_report_PAS_TO BE DELETED/scenarios")

rds_files <- list.files("E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files",
                        full.names = TRUE)

rds_files <- rds_files[grepl("\\.rds", rds_files)]
rds_files <- rds_files[grepl("Scenario", rds_files)]




# for (n in 1:length(rds_files)) {
#   res <-  readRDS(rds_files[n])
#   
#   summary <- res$mk$wa_summarised[1:3]
#   
#   summary <- lapply(summary, function(x) {
#     x$drug <- i$r_pld_lookup_mol$Description[match(x$L1, i$r_pld_lookup_mol$Number)]
#     x <- x[,c(1,5,2:4)]
#     x <- x[order(x$costs)]
#     return(x)})
# }
scenarios <- list(pop_1 = data.table(),
                  pop_2 = data.table(),
                  pop_3 = data.table())

for (n in 1:length(rds_files)) {
  res <-  readRDS(rds_files[n])

  scenario <- substring(rds_files[n], gregexpr("Scenario", rds_files[n])[[1]][1], (gregexpr("\\.rds", rds_files[n])[[1]][1]-1))
  
  non_dom <- lapply(res$mk$weighted_incremental[1:3], function(x) x$non_dominated)
  
  if (length(non_dom) == 0) {
    warning("No data recorded for ", scenario)
    next
  }
  
  next_best <- 
    lapply(non_dom, function(x) 
      x$L1[which(x$L1 == "Cabozantinib plus nivolumab")-1]
      )
  
  incrementals <- lapply(non_dom, function(x)
    x[which(x$L1 == "Cabozantinib plus nivolumab"),c("ic","iq","ICER")])

  for (j in 1:3) {
    scenarios[[j]] <- rbind(scenarios[[j]], cbind(scenario, next_best = next_best[j], incrementals[[j]]))
  }
}

write.csv(as.matrix(scenarios[[1]]), file = "4_Output/tables_for_report_PAS_TO BE DELETED/scenarios/pop1.csv", row.names = F)
write.csv(as.matrix(scenarios[[2]]), file = "4_Output/tables_for_report_PAS_TO BE DELETED/scenarios/pop2.csv", row.names = F)
write.csv(as.matrix(scenarios[[3]]), file = "4_Output/tables_for_report_PAS_TO BE DELETED/scenarios/pop3.csv", row.names = F)


#population 3 analyses fail in most cases - correcting manually.

#Manual extractions
#Scenario 1
  #pop_1
  nondom <- Scenario_1_PAS$ps$incremental$pop_1$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]

  ce <- Scenario_1_PAS$ps$tables$top_line$pop_1
  ce$drug <- i$r_pld_lookup_mol$Description[match(ce$L1, i$r_pld_lookup_mol$Number)]
  ce <- ce[,c(4,5,1,2,3)]
  ce <- ce[order(ce$costs),]
  
  round(ce[ce$drug == "Cabozantinib plus nivolumab",costs] - ce[ce$drug == "Pazopanib",costs],2)
  round(ce[ce$drug == "Cabozantinib plus nivolumab",qalys] - ce[ce$drug == "Pazopanib",qalys],3)
  #pop_2

  nondom <- Scenario_1_PAS$ps$incremental$pop_2$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]

  #pop_3
  
  nondom <- Scenario_1_PAS$ps$incremental$pop_3$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]
  
  nondom
  
#Scenario 2
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  nondom <- pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  nondom <- pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  nondom <- pas$mk$weighted_incremental$pop_3$non_dominated
  
#Scenario 3
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated
  
#Scenario 4
  #use base_Case for this (undisc)
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  #pop_1
  pop1 <- pas$weighted_model_undisc$pop_1
  pop1$costs <- apply(pop1[,2:27],1,sum)
  pop1$qalys <- apply(pop1[,28:44],1,sum)
  pop1 <- pop1[,c(1,54,55)]
  pop1$drug <- i$r_pld_lookup_mol$Description[match(pop1$L1, i$r_pld_lookup_mol$Number)]
  pop1 <- pop1[,c(1,4,2,3)]
  pop1 <- pop1[order(pop1$cost),]
  
  pop1
  #comparator = sunitinib
  round(pop1[pop1$drug == "Cabozantinib plus nivolumab",costs] - pop1[pop1$drug == "Sunitinib",costs])
  round(pop1[pop1$drug == "Cabozantinib plus nivolumab",qalys] - pop1[pop1$drug == "Sunitinib",qalys],3)
  
  round((pop1[pop1$drug == "Cabozantinib plus nivolumab",costs] - pop1[pop1$drug == "Sunitinib",costs]) /
  (pop1[pop1$drug == "Cabozantinib plus nivolumab",qalys] - pop1[pop1$drug == "Sunitinib",qalys]))
  
  #pop_2
  pop <- pas$weighted_model_undisc$pop_2
  pop$costs <- apply(pop[,2:27],1,sum)
  pop$qalys <- apply(pop[,28:44],1,sum)
  pop <- pop[,c(1,54,55)]
  pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  pop <- pop[,c(1,4,2,3)]
  pop <- pop[order(pop$cost),]
  
  pop
  #comparator = sunitinib
  round(pop[pop$drug == "Cabozantinib plus nivolumab",costs] - pop[pop$drug == "Sunitinib",costs])
  round(pop[pop$drug == "Cabozantinib plus nivolumab",qalys] - pop[pop$drug == "Sunitinib",qalys],3)
  
  round((pop[pop$drug == "Cabozantinib plus nivolumab",costs] - pop[pop$drug == "Sunitinib",costs]) /
          (pop[pop$drug == "Cabozantinib plus nivolumab",qalys] - pop[pop$drug == "Sunitinib",qalys]))
  
  #pop_3
  pop <- pas$weighted_model_undisc$pop_3
  pop$costs <- apply(pop[,2:27],1,sum)
  pop$qalys <- apply(pop[,28:44],1,sum)
  pop <- pop[,c(1,54,55)]
  pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  pop <- pop[,c(1,4,2,3)]
  pop <- pop[order(pop$cost),]
  
  pop
  #cabo+nivo dominated

#Scenario 5
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 6
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 7
  pas <- readRDS(rstudioapi::selectFile(
               path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))

  #pop_1
  nondom <- pas$ps$incremental$pop_1$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]
  nondom
  #pop_2
  nondom <- pas$ps$incremental$pop_2$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]
  nondom
  #pop_3
  nondom <- pas$ps$incremental$pop_3$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(2,9,3,4,1,5:8)]
  nondom
  
#Scenario 8
  # n/a'd in table?
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated  
  
  
#Scenario 9
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 10
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  nondom <- pas$mk$weighted_incremental$pop_1$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(1,9,2:8)]
  nondom
  #pop_2
  nondom <- pas$mk$weighted_incremental$pop_2$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(1,9,2:8)]
  nondom  
  #pop_3
  nondom <- pas$mk$weighted_incremental$pop_3$non_dominated
  nondom$drug <- i$r_pld_lookup_mol$Description[match(nondom$L1, i$r_pld_lookup_mol$Number)]
  nondom <- nondom[,c(1,9,2:8)]
  nondom  

#Scenario 11
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated    
  
#Scenario 12
  #TO FOLLOW?

#Scenario 13
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated    

#Scenario 14
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))

  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated
  
#Scenario 15
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 16
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 17
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 18
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 19
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 20
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated  
  
  
#Scenario 21
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 22
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  
  
#Scenario 23
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 24
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 25
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated  

#Scenario 26
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pop <- pas$tables$top_line$pop_1
  pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  pop <- pop[,c(4:5,1:3)]
  pop
  #pop_2
  pas$incremental$pop_2$non_dominated
  #pop_3
  pas$incremental$pop_3$non_dominated
  
#Scenario 27
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated
  
# Scenarios 28-37 
  #TO FOLLOW

#Scenario 38
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated
  
#Scenario 39
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated
  
#Scenario 40
  #TO FOLLOW

#Scenario 41
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated
  
#Scenario 42
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated

#Scenarios 43-51
  #TO FOLLOW

#Scenario 52
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files"))
  
  #pop_1
  pas$mk$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$mk$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$mk$weighted_incremental$pop_3$non_dominated
  
  
  
# LIST PRICES
#Base Case
  #RERUNNING
  # pas <- readRDS(rstudioapi::selectFile(
  #   path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  # #pop_1
  # pop <- pas$weighted_model_disc$pop_1
  # colnames(pop)[2:27]
  # pop$costs <- apply(pop[,2:27],1,sum)
  # pop$qalys <- apply(pop[,28:44],1,sum)
  # pop <- pop[,c(1,45,46)]
  # pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  # pop <- pop[,c(1,4,2,3)]
  # pop <- pop[order(pop$cost),]
  # pop
  # 
  # pas$weighted_incremental$pop_3$non_dominated
  # 
  # 
  # #pop_2
  # pop <- pas$weighted_model_undisc$pop_2
  # pop$costs <- apply(pop[,2:27],1,sum)
  # pop$qalys <- apply(pop[,28:44],1,sum)
  # pop <- pop[,c(1,54,55)]
  # pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  # pop <- pop[,c(1,4,2,3)]
  # pop <- pop[order(pop$cost),]
  # 
  # pop
  # #comparator = sunitinib
  # round(pop[pop$drug == "Cabozantinib plus nivolumab",costs] - pop[pop$drug == "Sunitinib",costs])
  # round(pop[pop$drug == "Cabozantinib plus nivolumab",qalys] - pop[pop$drug == "Sunitinib",qalys],3)
  # 
  # round((pop[pop$drug == "Cabozantinib plus nivolumab",costs] - pop[pop$drug == "Sunitinib",costs]) /
  #         (pop[pop$drug == "Cabozantinib plus nivolumab",qalys] - pop[pop$drug == "Sunitinib",qalys]))
  # 
  # #pop_3
  # pop <- pas$weighted_model_undisc$pop_3
  # pop$costs <- apply(pop[,2:27],1,sum)
  # pop$qalys <- apply(pop[,28:44],1,sum)
  # pop <- pop[,c(1,54,55)]
  # pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  # pop <- pop[,c(1,4,2,3)]
  # pop <- pop[order(pop$cost),]
  # 
  # pop
  # #cabo+nivo dominated
  

#Scenario 1
  # NOTE: use file ps_model_Scenario_1.rds for this scenario
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))

  #pop_1
  pop <- pas$tables$top_line$pop_1
  pop$drug <- i$r_pld_lookup_mol$Description[match(pop$L1, i$r_pld_lookup_mol$Number)]
  pop <- pop[,c(4:5,1:3)]
  pop <- pop[order(pop$cost),]
  pop
  
  pas$incremental$pop_1
  
  #pop_2
  pas$incremental$pop_2$non_dominated
  
  #pop_3
  pas$incremental$pop_3$non_dominated
  
#Scenario 2
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated

#Scenario 3
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated

#Scenario 4 (Use Base Case file)
  #RERUNNING BASE CASE
  # pas <- readRDS(rstudioapi::selectFile(
  #   path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  # 
  # #pop_1
  # pas$weighted_incremental$pop_1$non_dominated
  # #pop_2
  # pas$weighted_incremental$pop_2$non_dominated
  # #pop_3
  # pas$weighted_incremental$pop_3$non_dominated
  
#Scenario 5
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated

#Scenario 6
  pas <- readRDS(rstudioapi::selectFile(
    path = "E:/University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files/list prices"))
  
  #pop_1
  pas$weighted_incremental$pop_1$non_dominated
  #pop_2
  pas$weighted_incremental$pop_2$non_dominated
  #pop_3
  pas$weighted_incremental$pop_3$non_dominated
  