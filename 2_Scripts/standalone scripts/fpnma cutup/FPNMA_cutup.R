# code to aggregate FPNMA means
library(data.table)

folder_to_FPNMA_codas <- rstudioapi::selectDirectory()
filelist <- list.files(folder_to_FPNMA_codas)
filelist <- filelist[grepl(".rds",filelist)]


# restructure the coda into complete sets of 1000 runs for all PLMTE (so 16 files for 16,000 runs)
FPNMA <- readRDS(paste0(folder_to_FPNMA_codas,"/",filelist[1]))
runs <- max(FPNMA$run)

for (i in seq(1, runs, by = 1000)) {
  cat(i, "of", runs,"\n")
  run <- FPNMA[FPNMA$run >= i & FPNMA$run < (i+1000),]
  saveRDS(run, file = paste0(folder_to_FPNMA_codas,"/by_run/FPNMA_",i,"_",i+1000,".rds"))
}

rm(FPNMA)
rm(run)

for (j in 2:length(filelist)) {
  cat("Processing file", filelist[j])
  FPNMA <- readRDS(paste0(folder_to_FPNMA_codas,"/",filelist[j]))
  runs <- max(FPNMA$run)
  
  for (i in seq(1, runs, by = 1000)) {
    cat(i, "of", runs,"\n")
    run <- readRDS(paste0(folder_to_FPNMA_codas,"/by_run/FPNMA_",i,"_",i+1000,".rds"))
    run2 <- FPNMA[FPNMA$run >= i & FPNMA$run < (i+1000),]
    run <- rbind(run, run2)
    rm(run2)
    saveRDS(run, file = paste0(folder_to_FPNMA_codas,"/by_run/FPNMA_",i,"_",i+1000,".rds"))
  }
  
  rm(FPNMA)
  rm(run)
}

#calculate means from 16,000 runs
filelist <- list.files(paste0(folder_to_FPNMA_codas,"/by_run"))
filelist <- filelist[grepl(".rds",filelist)]

means <- list()
for (i in 1:length(filelist)) {
  cat(filelist[i],"\n")
  FPNMA <- readRDS(paste0(folder_to_FPNMA_codas,"/by_run/",filelist[i]))
  FPNMA <- as.data.table(FPNMA)
  means[[i]] <- FPNMA[,mean(HR), by = .(time, intervention_code, reference_treatment_code, line, endpoint, population, ref_trial_code)]
}

means2 <- NULL
for (i in 1:length(means)) {
  means2 <- rbind(means2,means[[i]])
}
colnames(means2)[8] <- "HR"

means2 <- means2[,mean(HR), by = .(time, intervention_code, reference_treatment_code, line, endpoint, population, ref_trial_code)]

saveRDS(means2, file = paste0(folder_to_FPNMA_codas,"/FPNMA_means.rds"))
