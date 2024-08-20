
# The primary difference with the probabilistic model is that the survival
# analysis (along with several of the other input parameters) are drawn
# probabilistically repeatedly and those randomly drawn numbers
# are fed through the SAME functions that have been used in Model_structure.R
# to ensure consistency.
# 
# This is analogous to generating values in an excel model and passing it through
# the same calculation chain as the deterministic model via the patient flow
# sheet.
# 
# To achieve this, we generate the list "p" repeatedly, with each "p" being
# a "copy" of the p used in the determinstic model, following exactly the same
# structure, but with random numbers instead of point estimate numbers each time.
# 
# By running the same structure of p through the same function to run the model
# we can generate a probabilistic version of the model results, which can then
# be used to capture the uncertainty surrounding the ICER estimates.
# 
# !!!PLEASE READ!!! THIS IS IMPORTANT:
# 
# However, the ST model is highly computationally intensive. On a typical laptop,
# one run of the determinstic model across 3 different populations for ~150 pathways
# for 1-4 active treatment line on and off treatment tunnel states for 2000+ time cycles
# results in a lot of computations to get the most precise possible estimate
# of the time that patients will spend in each treatment line within a given
# pathway given a risk/sequencing population pairing. This is because it takes
# into account both absolute time and time within treatment line given treatment
# status, resulting in up to 18,000 discrete health states per treatment pathway.
# 
# The 4 active treatment line plus BSC pathways have f_markov_calcN(4,p$basic$th)
# or 14619 discrete health states. This is 14619 * p$basic$th or 30,524,472 
# state residency calculations per treatment pathway.
# 
# sequencing population 0 has p$seq$n$pop_0[!is.na(line_5),] (54) 4 active treatment
# line pathways, all plus BSC. 14619 * p$basic$th * 54 = 1,648,321,488 state
# residency values, just for overall population 1, and just for those pathways
# with 4 active treatment lines. The total number of values is much larger.
# 
# As can be seen from this example, the reason that the model is so computationally
# intensive is the inclusion of the tunnel states for all 2L+ living health states.
# 
# The computational time for each 4 ATL treatment pathway is around 2-3 seconds.
# Given a PSA of 1000 iterations, this would likely be something close to:
# 
# 2.5 seconds * 54 pathways * 1000
# 
# which would be (2.5*54*1000)/60/60 = 37.5 hours, only for one population 
# and only for the 4atl pathways in that population, and only to compute the Markov
# trace...
#  
# On top of that, due to the large matrix required to incorporate full memory
# into a sequencing Markov model, the RAM requirement for such a process would
# be extreme. Currently, to run the model using 8 cores requires something in the
# region of 24GB of RAM, but of course considerably improves runtime.
#  
# Consequently at this point there are three options which are available for
# generating probabilistic results from the state-transition model. These are:
#  
#  - Use a supercomputer. (or high-performance computing cluster, HPC) to provide 
#    hundreds of cores and TBs of RAM. This would reduce runtime to a few hours
#    but requires outsourcing of the computations. There are many options available
#    including commercial (Microsoft Azure, Amazon AWS), academic (Exeter university,
#    Sheffield university, UEA, others) to satisfy this. This would produce
#    the most accurate results and would ensure the use of identical code for
#    determinstic and iterative analyses.
#  - Approximate time in state using exponential approximation for the CURRENT
#    decision problem. In this specific case, there are no stopping rules or
#    complicated cost dynamics in any 2L+ state. Consequently the tunnels
#    are not actually required in this specific simple case. We therefore
#    define a non-tunnel trace calculator which operates with vastly reduced
#    computation, resulting in a PSA taking a few seconds rather than multiple
#    days.
#  - Forego a probabilsitic analysis
# 
# We would like to avoid not doing a probabilistic analysis, so options 1 and
# 2 are preferable. However, we cannot ensure that a HPC is available for running
# of this application. Therefore, we provide the alternative within this script
# for approximation.
# 


# Packages and functions --------------------------------------------------

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
library(optparse)

parser_results <- OptionParser() |>
  add_option(c("-k", "--job-id"), default = 1L, type = "integer", help = "The job ID in an array") |>
  add_option(c("-N", "--iterations-per-job"), default = 2L, type = "integer", help = "The number of PSA iterations for this job") |>
  parse_args()

job_id <- parser_results$`job-id`
n_psa <- parser_results$`iterations-per-job`

set.seed(job_id)

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
keep_free_cores <- NA
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




# i and p ---------------------------------------------------------------

i <- readRDS("./2_Scripts/standalone scripts/QC/i.rds")
p <- readRDS("./2_Scripts/standalone scripts/QC/p.rds")


# Replace the components of p with their probabilistic versions, one section
# at a time, starting with the disease model:

p$basic$jobid <- job_id
p$basic$npsa <- n_psa

# Component generation ----------------------------------------------------
# ~ Disease model ---------------------------------------------------------

# To test - generate all of the parameters for all distributions for all reference
# curve PLMTEs for npsa PSA iterations
i$psa_psm <- f_PSA_drawFSParams(
  surv_regs = i$surv$reg,
  n_psa = p$basic$npsa,
  return_rands = FALSE,
  lookups = p$basic$lookup$ipd,
  verbose = FALSE
)

# Great, that works. Now let's introduce a version which filters down to only
# reference curves which are included in the analysis for THIS scenario (i.e.
# according to THIS excel file)


# filter down the parameters to only the distributions which have been
# selected in Excel, and remove all reference curves that aren't used in the 
# model.
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

i$PSA_est_Lambdas_disc <- f_psa_approx_lambda(
  psa_params = i$psa_psm_filtered,
  method = "sum",
  th = p$basic$th,
  disc = p$basic$discFacQ
)

# So, when using a lambda computed to match the area under the curve of the 
# discounted time in state plot, one comes to a similar time in state when
# just using the lambda approximation method, which is great

# Using lambdas in the model ----------------------------------------------

# First, some simple maths:
# 
# Lambda can be translated to transition probability TP like so:
# 
#   TP_t = 1-s(t)/s(t-1)
#   
# Given that the rate is constant with exponential, all TP for a given endpoint
# are equal to TP_t, that is TP = 1-s(t)/s(t-1). 
# 
# At cycle 1 (given cycle 0 is model start), TP_t = 1-(s(t) / s(t-1)), but s(t-1)
# is known to be 1. Therefore TP = 1-s(1).
# 
# Thus, the entire set of lambdas can be converted to TPs by cycling through them
# and applying 1-f_psa_exp(1,lambda)

# Exponential lines -------------------------------------------------------

# The reference curves will all get exponential fits for undiscounted and 
# discounted. This will be done using lambda.
# 
# then, the relative efficacy network will be applied to them all to get
# the full st set
# 
# Then, the "full set" of s(t)  will be translated into TPs as normal
# 
# Then, instead of M there will be TH small M's which will only have 4, 6, 8 or 10
# as dim.
# 
# These will be used to compute a trace for each PSA run, for each 
# 
# 

# All of these steps happen per PSA run, and the results that are returned are 
# Markov traces. this way the individual traces can be calculated in parallel

# PSA input set (temporarily just for disease model):
i_psa <- list(
  ps = 1:p$basic$npsa,
  releff = list(
    coda = list(
      ph = i$PHNMA$data,
      fp = i$FPNMA$data
    ),
    table = p$releff$excel_table
  ),
  surv = list(
    lambda = i$PSA_est_Lambdas,
    lambda_dsic = i$PSA_est_Lambdas_disc,
    ref_curves = f_psa_lambda2St(i$PSA_est_Lambdas,0:p$basic$th),
    rc_params = i$psa_psm_filtered
  ),
  cost = list(),
  hrql = list(),
  ae   = list()
)

# At this point we have all the probabilistic reference curves as well as the CODA
# samples and the relative efficacy table.

# excel_efficacy_table <- i$R_table_eff_data_settings

i_psa$releff$table_hr <- f_psa_assumptionsTab_genHRs(excel_efficacy_table = i_psa$releff$table,p$basic$npsa)
i_psa$releff$table_noHR <- i_psa$releff$table[Include.in.this.analysis. == "Yes" & Effectiveness.data.source != "Apply HR to",]

# Generate the network p$basic$npsa times (this is a big object, and uses a lot of RAM)

i_psa$releff$blank_network <- f_NMA_generateNetwork(p$basic$id$ipd, p$basic$lookup$ipd)

# P_PSA -------------------------------------------------------------------

# Generate a version of p which contains probabilistic versions of the inputs
# per iteration. Where things have no probabilistic uncertainty they can
# come from the base p as why generate that data thousands of times
# 
# Everything in basic is base
# Demo gets generated
# 
# 
# Notes - 
#  - The cost object is large (around 3.4GB for 1000 iterations at 40 year TH)
#    HOWEVER, we are doing lambda approximation so only require one number for
#    all 2L+. This DRASTICALLY reduces the size of this object (to about 850MB)
# 
p_psa <- list(
  demo = lapply(p$demo$agg, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(category) {
        if (is.null(category$mean)) {
          return(category)
        } else if (is.null(category$se)) {
          rep(category$mean,p$basic$npsa)
        } else {
          rep(category$mean,p$basic$npsa)   #### note we do not include uncertainty around the cohort makeup in PSA
        }
      })
    })
  }),
  util =  f_process_utilities(
    raw_utilities = i$R_table_util,
    PSA = TRUE,
    samples = p$basic$npsa
  ),
  util_gpop_coefs = lapply(1:p$basic$npsa, function(nested_psa_iteration) {
    .p <- add_population_utility_params(list(), psa = TRUE, .i = i)
    .p$util$pop_norms
  }),
  costs = lapply(f_process_cost_data(
    drug_and_admin  = i$R_table_drug_admin_costs,
    per_cycle_costs = i$R_table_MRU,
    time_horizon    = p$basic$th,
    max_trt_lines   = p$basic$R_maxlines,
    RDI_source      = i$dd_sc_RDI,
    verbose         = FALSE,
    samples         = p$basic$npsa,
    PSA             = TRUE)
    ,f_psa_lambda_cost),
  releff = list()
)

# Pre-calculate the population utility norms since they will be the same across
# all sequences (though may vary across populations), and store in p

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
      X = 1:p$basic$npsa,
      FUN = function(nested_psa_iteration) lapply(ptc_L1, function(pop) adjust_utility(
        age            = pop$age,
        sex            = pop$sex,
        utilities      = base_utility,
        .patient_level = FALSE,
        .p             =
          list(
            basic = list(cl_y = p$basic$cl_y),
            util = list(pop_norms = p_psa$util_gpop_coefs[[nested_psa_iteration]])
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
      X   = 1:p$basic$npsa,
      FUN = function(nested_psa_iteration) {
        pop_0 <- adjust_utility(
          age            = i$R_table_patientagesex$Age[ipd_L1],
          sex            = if_else(i$R_table_patientagesex$Gender[ipd_L1] == "M", "male", "female"),
          utilities      = base_utility,
          .patient_level = TRUE,
          .p             =
            list(
              basic = list(cl_y = p$basic$cl_y),
              util = list(pop_norms = p_psa$util_gpop_coefs[[nested_psa_iteration]])
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


# QC stuff - un-comment to investigate:

# util[, .(
#   mean_on = mean(OnTxt), 
#   lb_on = quantile(OnTxt, 0.025),
#   ub_on = quantile(OnTxt, 0.975),
#   mean_off = mean(OffTxt), 
#   lb_off = quantile(OffTxt, 0.025),
#   ub_off = quantile(OffTxt, 0.975),
#   mean_PFS = mean(PFS), 
#   lb_PFS = quantile(PFS, 0.025),
#   ub_PFS = quantile(PFS, 0.975),
#   mean_PD = mean(PD), 
#   lb_PD = quantile(PD, 0.025),
#   ub_PD = quantile(PD, 0.975)
#   ), 
#   by = list(Population, Treatment.line, Molecule)]

# ~ FPNMA CODA sample -----------------------------------------------------

#read in FPNMA PSA parameters
p_psa$releff$PSAcoefficients <- read.csv("./1_Data/osipsa.csv")
p_psa$releff$PSAcoefficients <- rbind(p_psa$releff$PSAcoefficients , read.csv("./1_Data/osopsa.csv"))
p_psa$releff$PSAcoefficients <- rbind(p_psa$releff$PSAcoefficients , read.csv("./1_Data/pfipsa.csv"))
p_psa$releff$PSAcoefficients <- rbind(p_psa$releff$PSAcoefficients , read.csv("./1_Data/pfopsa.csv"))

#tidy and add exponents
p_psa$releff$PSAcoefficients <- f_FPNMA_tidy_and_add_exponents(
  PSAcoefficients = p_psa$releff$PSAcoefficients,
  exponents = i$R_table_FPNMA_coefficients
)

# The coefficients are used to generate HRs for each PSA iteration one at a time
# to avoid generating a huge amount of data, which will increase the RAM requirements
# of the PSA even more, slowing down the model considerably.


# ~ one-off costs ---------------------------------------------------------

# Calculating PSA samples for one off costs

Oneoffcost_dat <- data.table(
  p$costs$oneoff,
  SE_cost = as.numeric(data.table(i$R_table_MRU_oneoff)[Type.of.cost != "Treatment initiation\r\n",]$X5)
)
 
p_psa$costs$oneoff <- rbindlist(lapply(1:nrow(Oneoffcost_dat),  function(param_row) {
  id     <- Oneoffcost_dat[param_row,]
  ooc_m  <- id$cost
  ooc_se <- id$SE_cost
  id$cost <- NULL
  id$SE_cost <- NULL
  id     <- as.data.table(lapply(id, function(x) rep(x,p$basic$npsa)))
  id$iteration <- 1:p$basic$npsa
  id$cost      <- rnorm(p$basic$npsa,ooc_m,ooc_se)
  return(id)
}))

p_psa$costs$oneoff_mk <- p_psa$costs$oneoff[,.(cost = sum(cost)),by=list(Apply.to, iteration)]
Apply_to <- p$costs$oneoff[,.(cost = sum(cost)),by=list(Apply.to)]

p_psa$costs$oneoff_mk <- lapply(structure(
  Apply_to$Apply.to,
  .Names = Apply_to$Apply.to
), function(health_state) {
  p_psa$costs$oneoff_mk[Apply.to == health_state, ]$cost
})

# ~ AEs ---------------------------------------------------------------------

# Calculating PSA samples for AEs
AE_dat <- data.table(
  p$ae$mk$per_cycle,
  SE_cost = p$ae$mk$per_cycle$cost * i$i_SE_costs,
  SE_util = abs(p$ae$mk$per_cycle$QALYs * i$i_SE_util)
)

p_psa$ae <- rbindlist(lapply(1:nrow(AE_dat),  function(param_row) {
  id     <- AE_dat[param_row,]
  aec_m  <- id$cost
  aec_se <- id$SE_cost
  aeq_m  <- id$QALYs
  aeq_se <- id$SE_util
  id     <- id[,list(trt,line,molecule)]
  id     <- as.data.table(lapply(id, function(x) rep(x,p$basic$npsa)))
  id$iteration <- 1:p$basic$npsa
  id$cost      <- rnorm(p$basic$npsa,aec_m,aec_se)
  id$QALYs     <- rnorm(p$basic$npsa,aeq_m,aeq_se)
  return(id)
}))


# Dropping some temp stuff:
rm(AE_dat, Oneoffcost_dat, Apply_to)


# Minor prep functions for the patient flow calcs -------------------------




# PSA trace process -------------------------------------------------------


# The pre-generation method is extremely RAM intensive, with the network
# taking up 11GB of memory before even extrapolating the curves.
# 
# instead, it may be more optimal to produce st on the PSA level, so we
# here attempt to go all the way to that point:
tick <- Sys.time()

psa_results <- with_progress({
  pr <- progressr::progressor(along = i_psa$ps)
  rbindlist(future_lapply(i_psa$ps, future.seed = TRUE, future.chunk.size = 1, function(nested_psa_iteration) {
    
    psa_iteration <- nested_psa_iteration + (job_id - 1) * n_psa
    
    # PHNMA CODA sample:
    phnma_coda_run <- data.table(data.frame(i_psa$releff$coda$ph)[i_psa$releff$coda$ph$Run == psa_iteration,])
    
    # FPNMA CODA sample and HR extrapolation:
    p_psa$releff$CODA$FP <- f_generate_FPNMA_coda(
      coeffs = p_psa$releff$PSAcoefficients[run == psa_iteration, ],
      TH = p$basic$th,
      wks_per_month = i$i_mon_to_weeks
    )
    
    #add in deterministic 2L coda to PSA (PSA only conducted on 1L due to unstable 2L results in FPNMA)
    i$FPNMA$means_for_PSA_L2 <- i$FPNMA$means[Line==2]
    p_psa$releff$CODA$FP <- rbind(p_psa$releff$CODA$FP, i$FPNMA$means_for_PSA_L2)
    
    #rebase to add in cabo as 2L
    p_psa$releff$CODA$FP <- f_3L_rel_effect_same_as_2L(FPNMAdata = p_psa$releff$CODA$FP)
    p_psa$releff$fp_dest <- f_gen_destinations(fp_data = p_psa$releff$CODA$FP)
    
    # add in reference.trial 2
    p_psa$releff$CODA$FP <- f_add_reference_trial_2(fp_data = p_psa$releff$CODA$FP)
    
    #eliminate NAs in molecule
    p_psa$releff$fp_dest <- i$FPNMA$destinations[!is.na(Molecule), ]
    
    pr(paste0("PSA iteration #",psa_iteration))
    # Make empty network
    network <- f_NMA_linkPHNMA(
      network       = i_psa$releff$blank_network,
      hr_table      = phnma_coda_run
    )
    
    # TEMPORARY: put in the FPNMA sample
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
    
    # We are not going to keep the network. We are going to use it to propagate
    # and then return only the st object for efficiency
    st <- f_releff_PropNetwork(
      network = network,
      extraps = i_psa$surv$lambda,
      dos     = 10,
      verbose = FALSE,
      dist_lookups = p$basic$lookup$dist,
      excel_table = i_psa$releff$table,
      psa_lambda_flag = TRUE,
      psa_iteration = nested_psa_iteration,
      psa_params    =  i_psa$surv$rc_params,
      th = p$basic$th
    )
    
    # Now that S(t) is calculated for this PSA iteration for all PLMTEs, we can now adjust
    # all the curves ready for computation of the Markov traces
    
    # Commented out as adjuvant adjustment not used in base case / PSA
    
    # if(i$dd_adjforprioradjuvant == "Yes") {
    #   p$surv$st <- f_surv_adjuvant_HR(
    #     st              = st,
    #     adjuvant_impact = i$R_table_prior_IO_impact_eff,
    #     demo_table      = p$demo$table,
    #     lookup          = p$basic$lookup,
    #     verbose         = TRUE)
    # }
    
    
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
    
    
    # Now, we want to reduce the amount of data being stored as much as possible as 
    # the st object is huge. Each st is about 16MB, and 1000 PSA iterations is therefore
    # 16GB!
    # 
    # If we can re-compute lambdas for all 2L+ PLMTEs, we can collapse much of
    # this data back to being one value from 2000+ This should scale up tremendously.
    # 
    # To do this we will compute 2l+ lambdas and 1L  (TP=1-(s(t)/s(t-1))) to reduce
    # data need as much as possible!
    
    tp      <- f_psa_collapse_st_lambda2lplus(st = st, th = p$basic$th, disc = FALSE)
    
    # So, with these tp objects we can now cycle through all the treatment pathways
    # one at a time generating the TP matrices. these are time-varying because
    # of first-line, but do not include tunnels so are drastically simplified compared
    # to the full tunnels.
    
    # We now do some of the preamble for the deterministic model running in ST mode:
    
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
    
    # We essentially compute the PF object in the same way as we did for the
    # deterministic ST model, but using lambda approximation to simplify the computation
    # of the Markov trace. this removes all tunnels from the model, reducing
    # potentially 14,000+ health states to just 10. No sparse matrix multiplication
    # is required either. 
    
    PATIENT_FLOW <- f_psa_pf_computePF_mkLambda(
      pops          = overall_pops,
      basic         = p$basic,
      demo          = f_psa_get_it_demo(p_psa$demo, nested_psa_iteration),
      sequences     = p$seq,
      survival      = list(gpop = p$surv$gpop, tp = tp),
      costs         = list(per_cycle = p_psa$costs[[nested_psa_iteration]], one_off = lapply(p_psa$costs$oneoff_mk, "[[", nested_psa_iteration)),
      util          = list(hsuv  = p_psa$util[iteration == nested_psa_iteration,], gpop = p_psa$util_gpop[[nested_psa_iteration]]),
      ae            = list(one_off = p$ae$duration, per_cycle = p_psa$ae[iteration == nested_psa_iteration], approach = p$ae$approach),
      eff_table     = p$releff$excel_table,
      verbose       = TRUE,
      include_plots = FALSE,
      just_nlines   = NULL,
      just_seq      = NULL
    )
    
    # We summarise these individual patient flow objects:
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
    results_disc   <- lapply(res_disc, function(popu) {
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
    
    # These two tables are by drug costs, other costs, QALYs, LYs and some ids.
    # We can cycle through the different overall populations and merge it all together
    # into one big table for this PSA iteration. This table then gets merged with
    # the tables for all the other iterations into one table containing all the 
    # necessary PSA results for post evaluation.
    pop_lab <- names(results_disc)
    names(pop_lab) <- pop_lab
    return(rbindlist(lapply(pop_lab, function(pop_name) {
      tab_with_ly <- merge.data.table(results_disc[[pop_name]],results_undsic[[pop_name]][,list(trt_n,LY)])
      tab_with_ly$oo_pop <- as.numeric(gsub("pop_","",pop_name))
      tab_with_ly$iteration <- psa_iteration
      tab_with_ly
    })))
  }))
})
tock <- Sys.time() - tick
print(tock)

# Save the full table which contains results for ALL individual treatment pathways:
saveRDS(psa_results, paste0("./4_Output/PSA_output_",job_id,gsub(":","_",date()),".rds"))



# Model averaging ---------------------------------------------------------

# Now that the PSA has finished for all individual pathways for all PSA iterations,
# we can use that along with information on the market shares of different 
# subsequent treatments, AND the uncertainty around those to generate
# a probabilistic weighted average model results representing ALL treatment pathways
# that is consolidated by first-line therapy relevant to this decision problem.


# Note - the Exeter HPC cannot deal with the character "→" and encodes it as
# correct encoding error in output from HPC. It's fine to leave this uncommented
# as if the command doesn't find the character it won't correct it.
psa_results$trt_n <- gsub("b\006\022","\u2192",psa_results$trt_n)
psa_results$trt <- gsub("b\006\022","\u2192",psa_results$trt)

# Compute weighted average summary results, including weightings, weighted, weighted average
# and mean weighted average results:
wa_model <- f_psa_computeWAModelRes(
  R_table_sub_txts_prop_n_costs = i$R_table_sub_txts_prop_n_costs,
  sims = n_psa,
  lookups = p$basic$lookup,
  psa_results = psa_results, 
  PSA = TRUE
)

saveRDS(wa_model, paste0("./4_Output/PSA_output_weighted_",job_id,gsub(":","_",date()),".rds"))



