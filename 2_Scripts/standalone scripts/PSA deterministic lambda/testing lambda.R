# this script uses a saved input set p and i to perform a controlled test
# of model comparability to the deterministic model extrapolated with "full tunnels"
# for all states.

# preamble ----------------------------------------------------------------

# ~  Packages and functions --------------------------------------------------

#library(shiny, quiet = TRUE)   
library(gtools, quiet = TRUE)
library(openxlsx, quiet = TRUE)
#library(flexsurv, quiet = TRUE)
library(tidyverse, quiet = TRUE)
library(data.table, quiet = TRUE)
#library(heemod, quiet = TRUE)
#library(logOfGamma, quiet = TRUE)
library(ggplot2, quiet = TRUE)
library(survminer, quiet = TRUE)
library(officer, quiet = TRUE)
library(officedown, quiet = TRUE)
library(magrittr, quiet = TRUE)
library(Hmisc, quiet = TRUE)
library(future.apply, quiet = TRUE)
#library(crosstable, quiet = TRUE)
#library(flextable, quiet = TRUE)
library(stringr, quiet = TRUE)
#library(BCEA, quiet = TRUE)
library(collapse, quiet = TRUE)
library(scales, quiet = TRUE)
library(Matrix, quiet = TRUE)
library(progressr)
library(pracma)


# Multi-core processing:
# 
# Instructions.
# 
# This model is highly RAM intensive. You need a lot of RAM on your computer
# to run this model due to the large amount of very large matrix multiplications
# (up to approximately 15,000 discrete health states in the model). Therefore,
# in order to efficiently run the model, it is a balancing act between RAM
# usage and CPU usage. 
# 
# Some rough guidance is:
# 
# - If you have 8GB of RAM on your computer, you can run this model with 2 cores only
#   but it may even be faster to run in series if you have other things open on your
#   computer at the same time. Therefore, please set keep_free_cores to NA and run
#   the model in series. This is because when the RAM on your computer runs out
#   your computer will use the hard-disk instead which is extremely slow.
# - If you have 16GB of RAM on your computer, parallel should be a lot faster.
#   On my laptop (I7 8th gen, 16GB RAM, running Linux for low RAM usage) I can
#   run with 5 cores whilst using about 12GB of RAM running this model. 
# - if you have 24GB or 32GB of RAM, you should be able to run the model with 8
#   and up to around 14 cores before running out of RAM whilst running the model.
# - if you are using a HPC, you should be able to run this model with many cores
#   due to the typically large amount of RAM available per core in a HPC
# 
# 
keep_free_cores <- 3
if (is.na(keep_free_cores)) {
  plan(sequential)
} else {
  plan(multisession(workers = max(availableCores()-keep_free_cores,1)))
}

# Other generic settings for the progress bar and units for table widths
handlers("progress")
options(crosstable_units="cm")

#### 2. Loading functions ###########


# This variable is used throughout the model to define whether to provide additional outputs useful for QC or not
# The model will take longer to run when this is set to TRUE
qc_mode <- FALSE



# 2.1. Excel data extraction functions -----------------------------------------

#### These functions are used to extract parameters from the Excel input workbook for use in R
#### During Phase 2 a Shiny front-end will be added to the model which will allow an alternative mechanism to upload these types of inputs

source(file.path("./3_Functions/excel/extract.R"))

# 2.2. Treatment sequencing functions ----------------------------------------

#### Function: filter to active treatments and lines
##### Takes as an input the defined sequences, evaluation type and line to start the evaluation from 
##### Other input is % receiving each subs therapy at each line dependent on previous treatments received 
##### Reweights so that the % receiving each treatment sums to 100% within each arm / line being studied
##### Outputs a matrix that has the % receiving each possible combination

source("./3_Functions/sequencing/sequences.R")

# 2.3. Survival analysis functions ---------------------------------------------

# Function: conduct survival analysis
##### by treatment, line, population and outcome fitted survival curves using Flexsurvreg (exp, Weibull, lognormal, loglog, Gompertz, gen gamma)
##### calculation of and adjustment for general population
##### adjustment for treatment effect waning

source("./3_Functions/survival/Survival_functions.R")
source("./3_Functions/survival/other_cause_mortality.R")
source("./3_Functions/survival/treatment_effect_waning.R")
source("./3_Functions/misc/fpnma_fns.R")

# 2.4 Misc functions ----------------------------------------------------------

### these functions enable smoother data cleaning and manipulation

source("./3_Functions/misc/other.R")
source("./3_Functions/misc/shift_and_pad.R")
source("./3_Functions/misc/cleaning.R")

# 2.4.1 Functions imposing list structures -----------------------------------

source("./3_Functions/misc/nesting.R")
source("./3_Functions/misc/discounting.R")
source("./3_Functions/misc/qdirichlet.R")
source("./3_Functions/misc/plotting.R")
source("./3_Functions/misc/structure.R")


# 2.5 Utility functions -------------------------------------------------------

source("./3_Functions/utility/age_related.R")
source("./3_Functions/costs_and_QALYs/utility_processing.R")

# 2.6 AE functions --------------------------------------------------------

source("./3_Functions/adverse_events/AE_steps.R")

# 2.7 Cost calculation functions --------------------------------------------

source("./3_Functions/costs_and_QALYs/cost_processing.R")


# 2.8 State transition modelling functions --------------------------------

source("./3_Functions/markov/markov.R")

# 2.9 Patient flow functions ----------------------------------------------

source("./3_Functions/patient_flow/overarching.R")
source("./3_Functions/patient_flow/partitioned_survival.R")
source("./3_Functions/patient_flow/markov.R")
source("./3_Functions/patient_flow/drug_costs.R")
source("./3_Functions/patient_flow/hcru_costs.R")
source("./3_Functions/patient_flow/qalys.R")
source("./3_Functions/patient_flow/ae.R")



# 2.10 Results processing functions ---------------------------------------

source("./3_Functions/results/incremental_analysis.R")
source("./3_Functions/results/model_averaging.R")
source("./3_Functions/results/partitioned_survival.R")
source("./3_Functions/misc/severity_modifier.R")
source("./3_Functions/results/results_tables.R")

# PSA RELATED FUNCTIONS ---------------------------------------------------

source("./3_Functions/psa/psa functions.R")



# Load in i and p ---------------------------------------------------------

# Load a pre-saved i and p set so we know we're starting from the right place
i <- readRDS("./2_Scripts/standalone scripts/QC/i.rds")
p <- readRDS("./2_Scripts/standalone scripts/QC/p.rds")

# Now we have all the inputs required, EXCEPT for lambdas. we can use the functions
# from 



# Replicate prob model but det --------------------------------------------

i$psa_psm <- f_PSA_drawFSParams(
  surv_regs = i$surv$reg,
  n_psa = 1,
  return_rands = FALSE,
  lookups = p$basic$lookup$ipd,
  verbose = FALSE
)

# Now that the structure is laid down, you can see that the entries containing
# something already have "id" in them, we can use this:
# get_elem(i$psa_psm,"id")

i$psa_psm <- lapply(i$psa_psm, function(risk_pop) {
  lapply(risk_pop, function(tr_line) {
    lapply(tr_line, function(mol) {
      lapply(mol, function(trial) {
        lapply(trial, function(endpoint) {
          if (is.null(endpoint)) {
            return(NULL)
          } else {
            id <- endpoint$id
            fits <- f_misc_get_plmte(i$surv$reg,id)$fs_fits
            draws <- lapply(fits, function(x) x$coefs)
            namdist <- names(fits)
            names(namdist) <- namdist
            draws <- lapply(namdist, function(dis) {
              if(dis=="exp") {
                as.matrix(draws[[dis]])
              } else {
                t(as.matrix(draws[[dis]]))
              }
            })
            return(list(
              id = id,
              draws = draws
            ))
          }
        })
      })
    })
  })
})

i$psa_psm_filtered <- f_psa_surv_params_filter(
  psa_psm_param = i$psa_psm,
  excel_eff_table = data.table(i$R_table_eff_data_settings),
  lookups = p$basic$lookup
)

i$psa_psm <- NULL

# Generate lambdas for all
i$PSA_est_Lambdas <- f_psa_approx_lambda(
  psa_params = i$psa_psm_filtered,
  method = "sum",
  th = p$basic$th
)


i_psa <- list(
  ps = 1,
  releff = list(
    coda = list(
      ph = i$PHNMA$data,
      fp = i$FPNMA$data
    ),
    table = p$releff$excel_table
  ),
  surv = list(
    lambda = i$PSA_est_Lambdas,
    ref_curves = f_psa_lambda2St(i$PSA_est_Lambdas,0:p$basic$th),
    rc_params = i$psa_psm_filtered
  ),
  cost = list(),
  hrql = list(),
  ae   = list()
)

# Make blank network - we dont' need to mess about with PSA style HR generation
# in this table as we're using the determinstic values:
i_psa$releff$blank_network <- f_NMA_generateNetwork(p$basic$id$ipd, p$basic$lookup$ipd)


p_psa <- list(
  demo = p$demo$live,
  util =  f_process_utilities(
    raw_utilities = i$R_table_util,
    PSA = FALSE
  ),
  util_gpop_coefs = {
    .p <- add_population_utility_params(list(), psa = FALSE, .i = i)
    .p$util$pop_norms
  },
  costs = f_psa_lambda_cost(f_process_cost_data(
    drug_and_admin  = i$R_table_drug_admin_costs,
    per_cycle_costs = i$R_table_MRU,
    time_horizon    = p$basic$th,
    max_trt_lines   = p$basic$R_maxlines,
    RDI_source      = i$dd_sc_RDI,
    verbose         = FALSE,
    samples         = p$basic$npsa,
    PSA             = FALSE)),
  releff = list()
)

base_utility <- data.frame(cycle = 0:p$basic$th, utility = 1)
if (i$dd_ageadjuutilities == "Yes") {
  if (i$dd_age_sex_source == "Mean") {
    # We find the row corresponding to line 1 for each relevant population
    
    # Do a lot of wrangling to get in the format we want...
    ptc_L1 <- i$R_table_ptchar[Treatment.line == 1, c(1, 3, 4)]
    colnames(ptc_L1) <- c("Population", "age", "sex")
    ptc_L1$sex <- 1 - ptc_L1$sex
    ptc_L1 <- merge(ptc_L1, i$lookup$ipd$pop, by.x = "Population", by.y = "Description")
    ptc_L1 <- ptc_L1[order(ptc_L1$Number), c("age", "sex", "Number")]
    ptc_L1 <- split(ptc_L1[, c("age", "sex")], paste0("pop_", ptc_L1$Number))
    
    p_psa$util_gpop <- lapply(
      X = 1,
      FUN = function(nested_psa_iteration) lapply(ptc_L1, function(pop) adjust_utility(
        age            = pop$age,
        sex            = pop$sex,
        utilities      = base_utility,
        .patient_level = FALSE,
        .p             =
          list(
            basic = list(cl_y = p$basic$cl_y),
            util = list(pop_norms = p_psa$util_gpop_coefs)
          )
      ))
    )
    
  } else {
    # We will only include IPD from line 1, since the population
    # norm is applied according to absolute model time rather than
    # than time in state. We don't know which population they are
    # in, so we will replicate for pop_0, pop_1 and pop_2.
    ipd_L1 <- i$R_table_patientagesex$Line == 1
    p_psa$util_gpop <- lapply(
      X   = 1,
      FUN = function(nested_psa_iteration) {
        pop_0 <- adjust_utility(
          age            = i$R_table_patientagesex$Age[ipd_L1],
          sex            = if_else(i$R_table_patientagesex$Gender[ipd_L1] == "M", "male", "female"),
          utilities      = base_utility,
          .patient_level = TRUE,
          .p             =
            list(
              basic = list(cl_y = p$basic$cl_y),
              util = list(pop_norms = p_psa$util_gpop_coefs)
            )
        )
        list(
          pop_0 = pop_0,
          pop_1 = pop_0,
          pop_2 = pop_0
        )
      }
    )
  }
} else {
  p_psa$util_gpop <- lapply(1:p$basic$npsa, function(nested_psa_iteration) list(pop_0 = 1, pop_1 = 1, pop_2 = 1))
} 







# Run model deterministic -------------------------------------------------


# PHNMA use means
phnma_coda_run       <- p$releff$CODA$PH

# FPNMA use means
p_psa$releff$CODA$FP <- p$releff$CODA$FP

# link releff
network <- f_NMA_linkPHNMA(
  network       = i_psa$releff$blank_network,
  hr_table      = phnma_coda_run
)
network <- f_NMA_linkFPNMA(
  network       = network,
  destinations  = p$releff$fp_dest,
  hr_table      = p$releff$CODA$FP,
  time_horizon  = p$basic$th
)

network <- f_NMA_AddAssumptionsToNetwork(
  network            = network,
  phnma_table        = phnma_coda_run,
  fpnma_table        = p$releff$CODA$FP,
  fpnma_destinations = p$releff$fp_dest,
  excel_table        = i_psa$releff$table,
  trial_flag         = i$List_eff_datasources[1],
  fpnma_flag         = i$List_eff_datasources[3],
  phnma_flag         = i$List_eff_datasources[2],
  et_flag            = i$List_eff_datasources[4],
  ahr_flag           = i$List_eff_datasources[5],
  verbose            = FALSE,
  psa_flag           = TRUE
)

st <- f_releff_PropNetwork(
  network         = network,
  extraps         = i_psa$surv$lambda,
  dos             = 10,
  verbose         = FALSE,
  dist_lookups    = p$basic$lookup$dist,
  excel_table     = i_psa$releff$table,
  psa_lambda_flag = TRUE,
  psa_iteration   = 1,
  psa_params      = i_psa$surv$rc_params,
  th              = p$basic$th
)


# Adjust curves:
if (sum(i$R_table_TE_waning_settings$apply.waning == "Yes") > 0) {
  st <- f_surv_twaning_apply(
    st_list     = st,
    tab_waning  = data.table(i$R_table_TE_waning_settings),
    tab_eff_set = data.table(i$R_table_eff_data_settings),
    verbose     = FALSE
  )
}
st <- f_surv_gpopadjust(
  st      = st,
  gpop    = p$surv$gpop,
  method  = "hazardmax",
  verbose = FALSE
)
st <- f_surv_PFSxOS(st = st, method = if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})
st <- f_surv_TTDxOS(st, if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})
st <- f_surv_PFSxTTP(st = st,method =  "abs")
st <- lapply(st, function(popu) {
  popu$line_5 <- popu$line_4
  return(popu)
})


# Make transition probs for lambda approach:
tp <- f_psa_collapse_st_lambda2lplus(st = st, th = p$basic$th, disc = FALSE)

# top-line inputs:
struct      = p$basic$structure
verbose     = FALSE
plots       = FALSE
just_pop    = p$basic$pops_to_run
just_nlines = NULL
just_seq    = NULL

if (!is.null(just_pop)) {
  if(0 %in% just_pop) stop ("this is overall population, not risk population, it can't be 0.")
  overall_pops <- structure(paste0("pop_",just_pop),.Names=paste0("pop_",just_pop))
} else {
  overall_pops <- structure(
    paste0("pop_",p$basic$lookup$pop_map$Overall.population.number),
    .Names=paste0("pop_",p$basic$lookup$pop_map$Overall.population.number)
  )
}
rpop <- paste0("pop_",p$basic$lookup$pop_map[match(as.numeric(gsub("pop_","",overall_pops)),p$basic$lookup$pop_map$Overall.population.number),]$Risk.population.number)


PATIENT_FLOW <- f_psa_pf_computePF_mkLambda(
  pops          = overall_pops,
  basic         = p$basic,
  demo          = p_psa$demo,
  sequences     = p$seq,
  survival      = list(gpop = p$surv$gpop, tp = tp),
  costs         = list(per_cycle = p_psa$costs, one_off = p$costs$oneoff_mk),
  util          = list(hsuv  = p$util$mk, gpop = p_psa$util_gpop[[1]]),
  ae            = list(one_off = p$ae$duration, per_cycle = p$ae$mk$per_cycle, approach = p$ae$approach),
  eff_table     = p$releff$excel_table,
  verbose       = TRUE,
  include_plots = FALSE,
  just_nlines   = NULL,
  just_seq      = NULL
)

res_undisc <- f_pf_mk_summary(
  pf_list = PATIENT_FLOW,
  disc_undisc = "undisc",
  lookups = p$basic$lookup,
  full_breakdown = TRUE,
  breakdown = TRUE,
  ypc = p$basic$cl_y
)
res_disc <- f_pf_mk_summary(
  pf_list = PATIENT_FLOW,
  disc_undisc = "disc",
  lookups = p$basic$lookup,
  full_breakdown = TRUE,
  breakdown = TRUE,
  ypc = p$basic$cl_y
)

# Cost per drug per sequence per population

# Target - do this for all pathways!
empty_cost_mol_list <- paste0("mol_",p$basic$lookup$trt)
empty_cost_mol_list <- structure(numeric(length(empty_cost_mol_list)),.Names=empty_cost_mol_list)

# undiscounted and discounted results by treatment pathway:
results_undsic <- lapply(res_undisc, function(popu) {
  fbd <- popu$full_breakdowns
  
  # Get the drug costs per molecule per sequence:
  trt_list_L1  <- unlist(lapply(fbd, function(trt_sq) trt_sq$numb[1]))
  trt_list  <- lapply(fbd, function(trt_sq) paste0("mol_",trt_sq$numb))
  cost_list <- lapply(fbd, function(trt_sq) trt_sq$cost[,"drug"])
  drug_cost <- data.table(do.call(
    rbind,
    lapply(structure(1:length(cost_list),.Names=names(trt_list)), function(trt_sq) {
      cl <- cost_list[[trt_sq]]
      names(cl) <- trt_list[[trt_sq]]
      empty_cost_mol_list[names(cl)] <- cl
      return(empty_cost_mol_list)
    })
  ))
  
  # now we need other costs, which are nicely in breakdowns
  bdt <- cbind(popu$breakdowns,drug_cost)
  bdt$drug <- NULL
  bdt$L1 <- trt_list_L1
  
  bdt[,`:=`(other_costs = admin + mru_on + mru_off + ae_cost + eol, qaly = qaly + ae_qaly)]
  bdt[,`:=`(admin = NULL, mru_on = NULL, mru_off = NULL, ae_cost = NULL, eol = NULL, ae_qaly = NULL)]
  
  # Now we just need life years and weightings:
  lybd <- popu$ly$breakdown[,.(LY = sum(L1_on,L1_off,BSC,L2_on,L2_off,L3_on,L3_off,L4_on,L4_off,na.rm = TRUE)),by="trt"]
  
  model_breakdown <- merge.data.table(bdt,lybd)
  rm(bdt)
  rm(lybd)
  
  # Now we have all the results we need for this PSA iteration consolidated
  # together in one table. However, it's too granular and we haven't
  # merged in the weightings yet.
  model_breakdown
  
})
results_disc <- lapply(res_disc, function(popu) {
  fbd <- popu$full_breakdowns
  
  # Get the drug costs per molecule per sequence:
  trt_list_L1  <- unlist(lapply(fbd, function(trt_sq) trt_sq$numb[1]))
  trt_list  <- lapply(fbd, function(trt_sq) paste0("mol_",trt_sq$numb))
  cost_list <- lapply(fbd, function(trt_sq) trt_sq$cost[,"drug"])
  drug_cost <- data.table(do.call(
    rbind,
    lapply(structure(1:length(cost_list),.Names=names(trt_list)), function(trt_sq) {
      cl <- cost_list[[trt_sq]]
      names(cl) <- trt_list[[trt_sq]]
      empty_cost_mol_list[names(cl)] <- cl
      return(empty_cost_mol_list)
    })
  ))
  
  # now we need other costs, which are nicely in breakdowns
  bdt <- cbind(popu$breakdowns,drug_cost)
  bdt$drug <- NULL
  bdt$L1 <- trt_list_L1
  
  bdt[,`:=`(other_costs = admin + mru_on + mru_off + ae_cost + eol, qaly = qaly + ae_qaly)]
  bdt[,`:=`(admin = NULL, mru_on = NULL, mru_off = NULL, ae_cost = NULL, eol = NULL, ae_qaly = NULL)]
  
  # Now we just need life years and weightings:
  
  model_breakdown <- bdt
  rm(bdt)
  
  # Now we have all the results we need for this PSA iteration consolidated
  # together in one table. However, it's too granular and we haven't
  # merged in the weightings yet.
  model_breakdown
  
})

# Return a list for this PSA iteration containing lifetime outcome table
# which collapses to the following columns
# 
# Run | L1 trt | overall population | dcost mol 1, 2, 3,... | other costs | QALYs | LYs
# 
# 
# 

pop_lab <- names(results_disc)
names(pop_lab) <- pop_lab
# 
# Return
deterministic_lambda_results <- rbindlist(lapply(pop_lab, function(pop_name) {
  tab_with_ly <- merge.data.table(results_disc[[pop_name]],results_undsic[[pop_name]][,list(trt_n,LY)])
  tab_with_ly$oo_pop <- as.numeric(gsub("pop_","",pop_name))
  tab_with_ly
}))

deterministic_lambda_results$dd_drug_price_options <- rep(i$dd_drug_price_options, nrow(deterministic_lambda_results))
deterministic_lambda_results$iteration <- rep(1, nrow(deterministic_lambda_results))

saveRDS(deterministic_lambda_results, paste0("./4_Output/lambda_det_output_",gsub(":","_",date()),".rds"))

# Weight results -------------------------------------------------

wa_model <- f_psa_computeWAModelRes(
  R_table_sub_txts_prop_n_costs = i$R_table_sub_txts_prop_n_costs,
  sims = 1,
  lookups = p$basic$lookup,
  psa_results = deterministic_lambda_results, 
  PSA = FALSE
)

# Report results -------------------------------------------------

lu_mol = p$basic$lookup$ipd$mol

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

#create tables
for (pricing in unique(names(report_tables))) {
  #pricing <- unique(names(report_tables))[1]
  for (pop in unique(names(report_tables[[pricing]]))) {
    #pop <- unique(names(report_tables[[pricing]]))[1]
    
    data <- report_tables[[pricing]][[pop]][["data"]]
    totals <- t(rbind(t(apply(data$cost,2,mean)),
                    t(apply(data$LY,2,mean)),
                    t(apply(data$qaly,2,mean))))
    totals <- cbind(as.numeric(rownames(totals)), totals)
    
    colnames(totals) <- c("L1","costs",
                          "ly",
                          "qalys")
    
    table_for_increments <- as.data.table(totals)
    
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
      
    
    }
    totals <- merge(totals,table_for_increments)
    
    totals <- totals[,c("trt", "costs", 
                        "ly", 
                        "qalys", 
                        "ic",
                        "il", 
                        "iq", "ICER")]
    totals <- totals[order(totals$costs),]
    
    report_tables[[pricing]][[pop]][["tables"]] <- list(totals = totals)
    
  }
}


pricing <- i$dd_drug_price_options

lu_pop <- p$basic$lookup$pop_map

ft_basic_bop <- do.call(rbind, lapply(structure(
  names(report_tables[[pricing]]), .Names = names(report_tables[[pricing]])
), function(popu_txt) {
  popu <- report_tables[[pricing]][[popu_txt]]$tables$totals
  popu_n <- as.numeric(gsub("pop", "", popu_txt))
  
  # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
  rsk_popu_lab <-
    lu_pop[Overall.population.number == popu_n,]$Risk.population
  
  # popu$seq_pop <- seq_popu_lab
  popu$risk_pop <- rsk_popu_lab
  
  
  return(popu)
  
}))

setDT(ft_basic_bop)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]

Word_width_inches = 29.7*0.3937

ft_det_lamda_output <- ft_basic_bop %>%
  rename(`Risk population` = risk_pop) %>%
  as_grouped_data(groups = "Risk population") %>%
  as_flextable() %>%
  width(., width = (Word_width_inches / (ncol(ft_basic_bop)))) %>%
  theme_box() |>
  set_header_labels(
    values = list(
      trt = "Technologies",
      costs = "Costs (£)",
      ly = "LYG",
      qalys = "QALYs",
      ic = "Inc. Costs",
      il = "Inc. LYG",
      iq = "Inc. QALYs",
      ICER = "ICER incremental"
    )
  ) %>%
  flextable::colformat_double(j = c(2, 5, 8),
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

ft_det_lamda_output
