#install.packages("devtools")
#devtools::install_github('Sheffield-Accelerated-VoI/SAVI-package')


library(BCEA)
library(SAVI)
library(dplyr)
library(ggplot2)
library(purrr)

#run Model_structure.R 1-1801 to generate p
lu_mol = p$basic$lookup$ipd$mol

source("./3_Functions/results/incremental_analysis.R")
wa_model <- readRDS("./4_Output/PSA_output_weighted_summaryThu Aug 24 10_59_57 2023.rds")

results <- wa_model$weighted_average
results$cost <- apply(results[, c("mol_0", "mol_1", "mol_2", "mol_3","mol_4", "mol_5", 
                                  "mol_6", "mol_7", "mol_8", "mol_9", "mol_10",
                                  "mol_11", "mol_12", "mol_999", "other_costs")],1, sum)

results$trt <- NA

results$trt[results$L1 == 1] <- "Cabozantinib plus nivolumab"
results$trt[results$L1 == 2] <- "Nivolumab plus ipilimumab"
results$trt[results$L1 == 3] <- "Lenvatinib plus pembrolizumab"
results$trt[results$L1 == 5] <- "Pazopanib"
results$trt[results$L1 == 6] <- "Tivozanib"
results$trt[results$L1 == 7] <- "Sunitinib"
results$trt[results$L1 == 8] <- "Cabozantinib"


results <- results[,c("dd_drug_price_options", "oo_pop", "L1", "trt", "iteration", "cost", "LY", "qaly")]

#set up tables list
report_tables <- list()

# read in formatted data
for (pricing in unique(results$dd_drug_price_options)) {
  for (pop in unique(results$oo_pop)) {
    results2 <- results[(results$dd_drug_price_option == pricing & 
                         results$oo_pop == pop),]
    
    cost <- results2[,c(3,5,6)]
    cost <- spread(cost, key = "L1", value = "cost")
    cost <- cost[,-1]
    
    LY <- results2[,c(3,5,7)]
    LY <- spread(LY, key = "L1", value = "LY")
    LY <- LY[,-1]

    qaly <- results2[,c(3,5,8)]
    qaly <- spread(qaly, key = "L1", value = "qaly")
    qaly <- qaly[,-1]
    
    data <- list(raw = results2,
                 cost = cost,
                 LY = LY,
                 qaly = qaly)
    rm(cost, LY, qaly)
    
    report_tables[[pricing]][[pop]] <- list(data = data)
    names(report_tables[[pricing]])[pop] <- paste0("pop",pop)
  }
}

# create BCEA objects
for (pricing in unique(names(report_tables))) {
  #pricing <- unique(names(report_tables))[1]
  for (pop in unique(names(report_tables[[pricing]]))) {
    #pop <- unique(names(report_tables[[pricing]]))[1]
    interventions <- lu_mol[match(colnames(report_tables[[pricing]][[pop]][["data"]][["cost"]]), lu_mol$Number)]$Description 
    report_tables[[pricing]][[pop]][["bcea"]][["LY"]] <- 
      bcea(eff = as.matrix(report_tables[[pricing]][[pop]][["data"]][["LY"]]),
      cost = as.matrix(report_tables[[pricing]][[pop]][["data"]][["cost"]]),
      ref = 1, 
      interventions = interventions,
      Kmax = 60000) 
    
    report_tables[[pricing]][[pop]][["bcea"]][["qaly"]] <- 
      bcea(eff = as.matrix(report_tables[[pricing]][[pop]][["data"]][["qaly"]]),
           cost = as.matrix(report_tables[[pricing]][[pop]][["data"]][["cost"]]),
           ref = 1, 
           interventions = interventions,
           Kmax = 60000) 
  }
}

#create tables
for (pricing in unique(names(report_tables))) {
  #pricing <- unique(names(report_tables))[1]
  for (pop in unique(names(report_tables[[pricing]]))) {
    #pop <- unique(names(report_tables[[pricing]]))[1]
    
    data <- report_tables[[pricing]][[pop]][["data"]]
    totals <- cbind(t(rbind(apply(data$cost,2,mean),
                    apply(data$cost,2,quantile, probs = c(0.025, 0.975)))),
                   t(rbind(apply(data$LY,2,mean),
                           apply(data$LY,2,quantile, probs = c(0.025, 0.975)))),
                   t(rbind(apply(data$qaly,2,mean),
                           apply(data$qaly,2,quantile, probs = c(0.025, 0.975)))))
    totals <- cbind(as.numeric(rownames(totals)), totals)
    
    colnames(totals) <- c("L1","Cost_mean","Cost_2.5","Cost_97.5",
                         "LY_mean","LY_2.5","LY_97.5",
                         "qaly_mean","qaly_2.5","qaly_97.5")
    table_for_increments <- totals[,c("L1","Cost_mean","qaly_mean","LY_mean")]
    colnames(table_for_increments) <- c("L1", "costs", "qalys", "ly")
    table_for_increments <- as.data.table(table_for_increments)

    incremental <- f_res_mk_incremental(res_d = table_for_increments, 
                         res_ud = table_for_increments,
                         lu_mol = lu_mol, 
                         produce_plot = FALSE, 
                         no_active_lines = 4, 
                         output_weighted = "Yes")$non_dominated
    
    colnames(incremental)[colnames(incremental) == "L1"] <- "trt"
    incremental$L1 <- lu_mol[match(incremental$trt, lu_mol$Description)]$Number
    table_for_increments$trt <- lu_mol[match(table_for_increments$L1, lu_mol$Number)]$Description
    table_for_increments <- merge(table_for_increments, incremental, all.x = TRUE)
    table_for_increments <- table_for_increments[order(table_for_increments$costs),]
    table_for_increments <- table_for_increments[,c("L1","trt", "costs", "qalys", "ly", "ic", "iq", "il", "ICER")]
    
    #for each increment calculate 95%CrI
    #get comparisons
      table_for_increments$il_ub <- table_for_increments$il_lb <- 
        table_for_increments$iq_ub <- table_for_increments$iq_lb <- 
        table_for_increments$ic_ub <- table_for_increments$ic_lb <- NA 
          
    comparisons <- table_for_increments$L1[!is.na(table_for_increments$ic)]
    
    if (length(comparisons)<2) stop("Less than one valid incremental pair")
    for (n in 2:length(comparisons)) {
      inc <- data.frame(cost = data$cost[[which(comparisons[n] == names(data$cost))]] - 
                    data$cost[[which(comparisons[n-1] == names(data$cost))]],
                  
                  ly = data$LY[[which(comparisons[n] == names(data$LY))]] - 
                    data$LY[[which(comparisons[n-1] == names(data$LY))]],
                  
                  qaly = data$qaly[[which(comparisons[n] == names(data$qaly))]] - 
                    data$qaly[[which(comparisons[n-1] == names(data$qaly))]])
      
      if (round(mean(inc$cost),0) != round(table_for_increments$ic[table_for_increments$L1 == comparisons[n]],0)) stop("Error in mean increments")
      
      inc <- apply(inc, 2, quantile,c(0.025, 0.975))
      table_for_increments$ic_lb[table_for_increments$L1 == comparisons[n]] <- inc["2.5%","cost"]
      table_for_increments$ic_ub[table_for_increments$L1 == comparisons[n]] <- inc["97.5%","cost"]
      table_for_increments$il_lb[table_for_increments$L1 == comparisons[n]] <- inc["2.5%","ly"]
      table_for_increments$il_ub[table_for_increments$L1 == comparisons[n]] <- inc["97.5%","ly"]
      table_for_increments$iq_lb[table_for_increments$L1 == comparisons[n]] <- inc["2.5%","qaly"]
      table_for_increments$iq_ub[table_for_increments$L1 == comparisons[n]] <- inc["97.5%","qaly"]
    }
    totals <- merge(totals,table_for_increments)

    #some final tests for consistency
    if (sum(totals$Cost_mean != totals$costs) + sum(totals$LY_mean != totals$ly) +
        sum(totals$qaly_mean != totals$qalys)) stop("Inconsistency in incremental results") 
    
    totals <- totals[,c("L1", "trt", "Cost_mean", "Cost_2.5", "Cost_97.5",
                        "LY_mean", "LY_2.5", "LY_97.5",
                        "qaly_mean", "qaly_2.5", "qaly_97.5",
                        "ic","ic_lb","ic_ub",
                        "il", "il_lb", "il_ub",
                        "iq", "iq_lb", "iq_ub", "ICER")]
    totals <- totals[order(totals$Cost_mean),]
    
    compact <- totals
    compact[,c("Cost_mean", "Cost_2.5", "Cost_97.5",
               "ic","ic_lb","ic_ub",
               "ICER")] <- round(compact[,c("Cost_mean", "Cost_2.5", "Cost_97.5",
                                            "ic","ic_lb","ic_ub",
                                            "ICER")],0)
    compact[,c("LY_mean", "LY_2.5", "LY_97.5",
               "qaly_mean", "qaly_2.5", "qaly_97.5",
               "il", "il_lb", "il_ub",
               "iq", "iq_lb", "iq_ub")] <- round(compact[,c("LY_mean", "LY_2.5", "LY_97.5",
                                                            "qaly_mean", "qaly_2.5", "qaly_97.5",
                                                            "il", "il_lb", "il_ub",
                                                            "iq", "iq_lb", "iq_ub")],3)
    
    compact$cost <- paste0(compact$Cost_mean," (", compact$Cost_2.5, ", ", compact$Cost_97.5, ")")
    compact$LY <- paste0(compact$LY_mean," (", compact$LY_2.5, ", ", compact$LY_97.5, ")")
    compact$qaly <- paste0(compact$qaly_mean," (", compact$qaly_2.5, ", ", compact$qaly_97.5, ")")
    compact$ic <- paste0(compact$ic," (", compact$ic_lb, ", ", compact$ic_ub, ")")
    compact$il <- paste0(compact$il," (", compact$il_lb, ", ", compact$il_ub, ")")
    compact$iq <- paste0(compact$iq," (", compact$iq_lb, ", ", compact$iq_ub, ")")
    
    compact <- compact[,c(c("L1", "trt", "cost", "LY", "qaly",
                            "ic","il","iq", "ICER"))]
    report_tables[[pricing]][[pop]][["tables"]] <- list(totals = totals,
                                                        compact = compact)
      
  }
}

write.csv(report_tables$`List price`$pop1$tables$compact, "./4_output/table_list_pop1_.csv")
write.csv(report_tables$`List price`$pop2$tables$compact, "./4_output/table_list_pop2_.csv")
write.csv(report_tables$`List price`$pop3$tables$compact, "./4_output/table_list_pop3_.csv")
write.csv(report_tables$`PAS price`$pop1$tables$compact, "./4_output/table_PAS_pop1_.csv")
write.csv(report_tables$`PAS price`$pop2$tables$compact, "./4_output/table_PAS_pop2_.csv")
write.csv(report_tables$`PAS price`$pop3$tables$compact, "./4_output/table_PAS_pop3_.csv")

x <- multi.ce(report_tables$`List price`$pop1$bcea$qaly)
ceac.plot(x)
x <- multi.ce(report_tables$`List price`$pop2$bcea$qaly)
ceac.plot(x)
x <- multi.ce(report_tables$`List price`$pop3$bcea$qaly)
ceac.plot(x)
x <- multi.ce(report_tables$`PAS price`$pop1$bcea$qaly)
ceac.plot(x)
x <- multi.ce(report_tables$`PAS price`$pop2$bcea$qaly)
ceac.plot(x)
x <- multi.ce(report_tables$`PAS price`$pop3$bcea$qaly)
ceac.plot(x)

# plot(bcea_oo_pop1_LY)
# ceplane.plot(bcea_oo_pop1_LY, wtp = 20000)
# eib.plot(bcea_oo_pop1_LY)
# contour(bcea_oo_pop1_LY)
# 
# bcea_oo_pop1_LY <- multi.ce(bcea_oo_pop1_LY)
# ceac.plot(bcea_oo_pop1_LY)
# eib.plot(bcea_oo_pop1_LY)