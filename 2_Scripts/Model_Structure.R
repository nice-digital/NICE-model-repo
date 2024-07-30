#### 1. Installation ###########
#### This code has been created using R version 4.3.1
#### All packages used by this model are provided here

#### Comment out the below section which installs the relevant packages after the first run of the model
# install.packages("shiny", quiet = TRUE)   ### the quiet argument is used to avoid warnings appearing in the console (useful for later conversion to web app)
# install.packages("gtools", quiet = TRUE)
# install.packages("openxlsx", quiet = TRUE)
# install.packages("flexsurv", quiet = TRUE)
# install.packages("tidyverse", quiet = TRUE)
# install.packages("data.table", quiet = TRUE)
# install.packages("heemod", quiet = TRUE)
# install.packages("logOfGamma", quiet = TRUE)
# install.packages("ggplot2", quiet = TRUE)
# install.packages("survminer", quiet = TRUE)
# install.packages("officer", quiet = TRUE)
# install.packages("officedown", quiet = TRUE)
# install.packages("magrittr", quiet = TRUE)
# install.packages("Hmisc", quiet = TRUE)
# install.packages("future.apply", quiet = TRUE)
# install.packages("crosstable", quiet = TRUE)
# install.packages("flextable", quiet = TRUE)
# install.packages("stringr", quiet = TRUE)
# install.packages("BCEA", quiet = TRUE)
# install.packages("collapse", quiet = TRUE)
# install.packages("scales", quiet = TRUE)
# install.packages("Matrix", quiet = TRUE)
# install.packages("dplyr", quiet = TRUE)
# install.packages("progressr", quiet = TRUE)
# install.packages("microbenchmark", quiet = TRUE)

### Loading libraries 

#### This section needs to be run every time and calls each package from the library 
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
library(progressr, quiet = TRUE)
library(microbenchmark, quiet = TRUE)


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
# IF YOU DO NOT WANT MULTICORE SET keep_free_cores TO NA
# 
# 
keep_free_cores <- 4
if (any(is.na(keep_free_cores), keep_free_cores<0)) {
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

# 2.4.2 Functions calculating HRs from FPNMA coefficients and other FPNMA manipulation ------

source("./3_Functions/misc/fpnma_fns.R")


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
source("./3_Functions/psa/psa functions.R")



# 2.11 Office software outputs --------------------------------------------

source("./3_Functions/reporting/word_document_output.R")



# 3. Model inputs structure --------------------------------------------------

# Model inputs should be in a list called i. This list then contains all of the
# inputs for the model, NOT the parameters used to calculate the model. In effect,
# this is a place to store all model information BEFORE it gets boiled down to
# what's needed to run 1 model.
# 
# using i allows subsetting by categorisation, which makes things a lot easier
# to find and avoids all long variable names
# 
# the structure of i should be by category. There are the following 
# categories:
# 
# dd - dropdown inputs taken from Excel
# i - parameter inputs taken from Excel
# r_ tables taken from Excel
# List, id and lookup - lists defined and used within the code
# basic - basic inputs (time horizon, cycle length, discount rates, so on so forth)
# surv  - survival analysis inputs including raw data
# sequences and seq - inputs and outputs related to the possible sequences of treatments
# cost  - drug and hcru costs. All costs are here to keep things together (dosing is not cost)
# util and QALYs  - utility and QALY inputs
# misc  - misc inputs e.g. graph labelling
# 

#### 3.1 Loading input parameters ###########

# This model allows two possible structures to be analysed: state transition with a user definable number of lines
# with health states based on time to discontinuation (drug costs) and progression status (quality of life and movement 
# between lines) and PartSA with 3 health states (pre-progression, post-progression and death)

# During Phase 1 of this pilot we use the model to evaluate the decision problem for a single therapy 
# (cabo+nivo, defined as molecule 1) starting at 1st line
# During Phase 2 we will adapt this code to evaluate the cost-effectiveness of sequences starting at a user-defined line

# Inputs to this model need to be downloaded from NICEdocs 

User_types <- c("Submitting company", "NICE", "EAG", "Committee", "NHSE", "Clinical expert", "Patient expert", "Non-intervention stakeholder", "Public")

# The submitting company are able to see their own CIC and AIC data (marked up blue / yellow in reporting but not anything else: green marking
# green marked data has been either be replaced with 0 [PAS discounts, RWE IPD] or dummy data)
# NICE users will be able to see everything
# Other users will not be able to see any marked data, this is replaced with dummy data

# The way raw data is fed into the model currently works as follows
# Define the path to where the data file lives using the select file functionality

# The model then processes the file the user selected

# There are a number of files which contain raw or intermediate inputs:
# 1. The Excel user interface - this contains information from company data and the UK RWE
# 2. The proportional hazards NMA CODA RDS file - this contains information from company data
# 3. The fractional polynomials NMA RDS file - this contains information from company data 
# 4. Either the raw data file containing the pseudo-IPD for all trials for survival analysis (RWE and company data included); or
# 5. The RDS output from the survival analysis using both RWE and company data

# You will need to manually select the inputs file relevant to your user type, this is not stored on Github as access to CIC information differs by user type

# The first part of this code pulls all of the named ranges from the excel workbook, expand the parameters table

#Option to define Excel path on local machine - comment in this and comment out the code below to select file
excel_path <- "./1_Data/ID6186_RCC_model inputs FAD version [ACIC cPAS redacted and CIC redacted].xlsm"
#i <- f_excel_extract(excel_path, verbose = TRUE)

if (file.exists(excel_path)) {
  i <- f_excel_extract(excel_path, verbose = TRUE)
} else {
  i <- f_excel_extract(rstudioapi::selectFile(
    caption = "Select the Excel inputs file (ID6184_RCC_model inputs....xlsm)",
    label = "ID6184_RCC_model inputs....xlsm",
    path = "./1_Data/",
    filter = "Excel Files (*.xlsm)",
    existing = TRUE
  ), verbose = TRUE)
}

i <- c(i,f_excel_cleanParams(i$R_table_param))


# Set which decision problem to look at, initially functionality has been geared towards the decision problem for cabozantinib plus nivolumab
i$decision_problem <- "cabo+nivo"

# We then create a place for identifiers. Adding in an object to i full of lookup tables makes automated translation
# possible even when one doesn't know the number of items ex ante, or how they combine.
# 
# If the lookup table is correct one can translate id numbers to text strings which are
# consistent throughout the entire model. This is extremely useful as the model can
# be expanded to any number of treatments and potentially even any number of lines 
# (up to a reasonable maximum)

i$id     <- list(ipd = list())
i$lookup <- list(ipd = list())

# Add distribution names to i
# This model only includes standard parametric distributions as more complex distributions were not deemed to be required for the included treatments

i$distnames <- 
  c(
    gengamma      = "gengamma",
    exp           = "exp",
    weibull       = "weibull",
    lnorm         = "lnorm",
    gamma         = "gamma",
    gompertz      = "gompertz",
    llogis        = "llogis"
  )



# The next step is to then "tidy up" i into another object, p. p doesn't necessarily
# have to house everything, only things that will change in PSA

p <- f_misc_param_generate_p(i)

# Set seed for PSA - note this is done in the script to run the PSA, not here!
# set.seed(1475)

# Max lines within the R model
p$basic$R_maxlines <- 4

# Pass this into p so that p can be used to exclusively compute the model:
p$basic$decision_problem <- i$decision_problem

#### 3.2 Define sequences  ###########

#### This code produces a list of possible sequences per population based upon the rules defined for RCC
#### and the user input number of lines


# Add drug names to comparators vector extracted from inputs list.

i$sequences <- f_generate_sequences(
  comparators = i$List_comparators, 
  maxlines    = p$basic$R_maxlines
)

# restrict the pathways to those that are possible and permitted.
i$sequences <- as.data.frame(i$sequences)

populations <- i$i_nr_populations

seqs <- NULL
for (population in 1:populations) {
  cat("Applying sequence restrictions to population", population,"\n")
  
  s <- f_path_tx_restrict(
    sequences                = i$sequences,
    allowed                  = f_get_allowed_lists(i, population), #overall list of allowed drugs in this popn
    L1                       = f_get_L1_lists(i, population), # 1L drugs allowed in this popn
    L2                       = f_get_L2_lists(i, population), # 2L drugs allowed in this popn
    L3                       = f_get_L3_lists(i, population), # 3L drugs allowed in this popn
    L4                       = f_get_L4_lists(i, population), # 4L drugs allowed in this popn
    only_after               = f_get_only_after_lists(i, population), #list of restrictions where tx can be only after the listed txs
    not_immediate_after      = f_get_not_immediate_after_lists(i, population), #list of restrictions where tx can be only immediately before the listed txs
    one_in_list              = f_get_one_in_list_lists(i, population), #list of restrictions where only one of the tx in each list is allowed 
    only_after_one           = f_get_only_after_one_lists(i, population), #list of restrictions where only one of the listed treatments is allowed prior to current therapy 
    L2_only_after            = f_get_2L_only_after_lists(i, population), #list of 2L+ restrictions: if drug is used 2L, 3L or 4L, can only be after drug x
    L2_only_immediate_after  = f_get_2L_only_immediate_after_lists(i, population), #list of 2L+ restrictions: if drug is used 2L, 3L or 4L, can only be immediately after drug x
    L2_only_one              = f_get_2L_only_one_lists(i, population) #list of 2L+ drugs where only one of them allowed in a given sequence
  )
  s <- cbind(rep(paste0("pop", population),nrow(s)), s)
  colnames(s) <- paste0('V', seq_len(ncol(s))) # rbind no longer likes un-named columns so added this
  seqs <- rbind(seqs, s)
}
rownames(seqs) <- NULL

i$sequences <- seqs

#### Uncomment this code to view the sequences and write the sequences defined to csv

# i$sequences
# write.csv(seqs, "4_Output/sequences.csv", row.names = F)
rm(s, seqs, populations)

# define number of cycles and a vector of the cycles 


# 3.3. Survival analysis -------------------------------------------------------

# All objects here go in i$surv initially, and are then streamlined down to 
# what's needed to run models in the transition from i to p.
# 
# Some values of p are used during the below (primarily p$surv$distNames, which
# controls which distributions are included in the flexsurv runs)


# 3.3.1 Survival input structure ------------------------------------------

i$surv <- list()

#### Read in survival data from Excel workbook 

# Pull out the raw data from the IPD excel book - one named range per treatment at each line
# Each reference curve is defined in Excel as time (weeks), event/censor (event coded as 1, censor as 0), patient group, line, molecule, trial and endpoint
# Pull all of the named ranges from the excel workbook, expand the parameters table

excel_path2 <- "./1_Data/IPD_R_input_noACIC.xlsx"
if (file.exists(excel_path2)) {
  wb <- f_excel_extract(excel_path2, verbose = TRUE)
  i$surv$pld <- as.data.table(wb$`_xlnm._FilterDatabase`)
  rm(wb)
} else {
  wb <- f_excel_extract(rstudioapi::selectFile(
    caption = "Select the IPD file (IPD_R_input_noACIC.xlsx)",
    label = "IPD_R_input_noACIC.xlsx",
    path = "./1_Data/",
    filter = "Excel Files (*.xlsx)",
    existing = TRUE
  ), verbose = TRUE)
  i$surv$pld <- as.data.table(wb$`_xlnm._FilterDatabase`)
  
}


# Some small cleaning of the PLD.
i$surv$pld <- i$surv$pld[,list(population,line,molecule,trial,endpoint,timew,event_censor)]

# Do not allow zero survival times, they have to be at least 1 day. the TUotA is
# weeks, so 1 day is 1/7 weeks:
i$surv$pld[timew ==0,"timew"] <- 1/7

# The named range r_pld has numeric identifiers for:
# 
# - pop
# - line
# - mol (i.e., regimen - combination therapies are under the same number)
# - trial (trial id WITHIN population line and molecule to set them apart from each other - usually just 1!)
# - endpoint

# These numeric identifiers are then used to create a nested list of survival regression models and
# extrapolations. The extrapolations are filtered down to the extrapolations that are selected
# within the excel input sheet, but the rest are kept here in i in case of scenario analysis.
# 
# Note that the lookup tables in the next section are used to translate these numbers
# into human-readable identifiers.

# 3.3.2 Data identification ------------------------------------------

# There is a lot of nesting involved in this part of the analysis, with population line, regimen trial and endpoint
# making a total of 5 layers of nesting to automatically go through each endpoint for each trial for
# each regimen for each line for each population, perform all regression analyses, produce parameters
# and have an easily identifiable (and therefore programmable) spaces for the results of each analysis
# which can then be spat out into reporting.

# The first step is to break up r_pld into separate datasets depending on the identifiers. A function
# is used to do this which returns nothing if such data for one id set doesn't exist. 
# 
# Note that at this stage it is just those contexts which HAVE got PLD which are to be organised.
# For those endpoints and so on that do not have data, a separate step after this one to populate 
# every endpoint for every treatment line for every treatment sequence is performed.

i$id$ipd <- list(
  pop      = i$r_pld_lookup_pop$Number[!is.na(i$r_pld_lookup_pop$Number)],
  line     = i$r_pld_lookup_line$Number[!is.na(i$r_pld_lookup_line$Number)],
  mol      = i$r_pld_lookup_mol$Number[!is.na(i$r_pld_lookup_mol$Number)],
  trial    = i$r_pld_lookup_trial$Number[!is.na(i$r_pld_lookup_trial$Number)],
  endpoint = i$r_pld_lookup_endpoint$Number[!is.na(i$r_pld_lookup_endpoint$Number)]
)

names(i$id$ipd$pop)      <- paste0("pop_"     , i$id$ipd$pop)
names(i$id$ipd$line)     <- paste0("line_"    , i$id$ipd$line)
names(i$id$ipd$mol)      <- paste0("mol_"     , i$id$ipd$mol)
names(i$id$ipd$trial)    <- paste0("trial_"   , i$id$ipd$trial)
names(i$id$ipd$endpoint) <- paste0("endpoint_", i$id$ipd$endpoint)


# to see this, we have:
#i$id$ipd

# Generating the same structure but with the translation table from number to
# text:

i$lookup$ipd <- list(
  pop      = data.table(i$r_pld_lookup_pop)[Description != 0],
  line     = data.table(i$r_pld_lookup_line)[Description != 0],
  mol      = data.table(i$r_pld_lookup_mol)[Description != 0],
  trial    = data.table(i$r_pld_lookup_trial)[Description != 0],
  endpoint = data.table(i$r_pld_lookup_endpoint)[Description != 0]
)

# For treatment line, add a translator for the column in the sequences output:

i$lookup$ipd$line$seq_col <- paste0("V",2:(nrow(i$lookup$ipd$line)+1))
i$lookup$ipd$line$R_id    <- paste0("line_",1:nrow(i$lookup$ipd$line))

i$lookup$dist <- i$r_pld_lookup_dist


# This means that you can easily look up things like so:

# i$lookup$ipd$mol[Number == 1,list(Description,RCC_input_desc)]
# i$lookup$ipd$mol[Number == 2,list(Description,RCC_input_desc)]
# i$lookup$ipd$line[Number == 1,list(Description,RCC_input_desc)]
# i$lookup$ipd$pop[Number == 0,list(Description,RCC_input_desc)]

# One can also do the opposite, translating input file descriptions into numbers:

# i$lookup$ipd$mol[RCC_input_desc == "ipi_nivo",list(Description,Number)]

i$lookup$trt <- i$lookup$ipd$mol$Number
names(i$lookup$trt) <- i$lookup$ipd$mol$RCC_input_desc
names(i$lookup$trt)[length(i$lookup$trt)] <- "BSC"

# pass to p whenever i$lookup has been populated/updated.
p$basic$lookup <- i$lookup
p$basic$id <- i$id

# one can then simply i$lookup$trt["nivolumab"] or i$lookup$trt["sorafenib"] to 
# get the id numbers.

# This then means that one can translate the treatment sequence data generated earlier
# into numerical versions in one go:

# Start by making the id for population fit with the rest of the model (pop_ with pop
# starting from 0). NOTE that there is 1 more population in treatment sequences than
# in the rest of the model...

i$seq_clean <- data.table(i$sequences)

i$seq_clean$V1 <- paste0("pop_",as.numeric(substr(i$seq_clean$V1,4,4)) - 1)

i$seq_pops <- unique(i$seq_clean$V1)
names(i$seq_pops) <- i$seq_pops

# The "clean" version of sequences - first with words, then with numbers, then references

i$seq_clean <- lapply(i$seq_pops, function(popu) {
  tmp <- i$seq_clean[V1 == popu,-1]
  colnames(tmp) <- i$lookup$ipd$line$R_id[1:(p$basic$R_maxlines + 1)]
  tmp
})

# It's pretty nested this but simplifies upon explanation: lapply on a data.frame
# or data.table goes column-wise, so going across columns substitute the values
# for the values in i$lookup$trt which have corresponding names, returning the numbers
# which are consistent throughout the model. The way of looking inside e.g. network
# is e.g. pop_2$line_5$mol_2$endpoint_1, so now we can use the tables produced below
# to "order" the inputs for a treatment pathway 
i$seq_n <- lapply(i$seq_clean, function(popu) {
  as.data.table(lapply(popu, function(co) i$lookup$trt[co]))
})
i$seq_ref <- lapply(i$seq_clean, function(popu) {
  tmp <- as.data.table(lapply(popu, function(co) {
    vals <- paste0("mol_",i$lookup$trt[co])
    ifelse(vals == "mol_NA",NA,vals)
  }))
})


# Now that we have the final sequence list, we can add them to p:

p$seq$n   <- i$seq_n
p$seq$ref <- i$seq_ref
p$seq$qc <- i$seq_clean

# NOTE: QC check here is for NAs that are not beyond a 999 (i.e. past BSC)

# We now have all the treatment sequences in the form of the molecule
# number and the consistent reference linking right back to the named range
# r_pld_lookup_mol in the excel front end. This ensures that the R model is
# consistent with the R model in terms of which drugs are feeding through
# to different places, as manually checking that is a very difficult and time 
# consuming task.
# 
# Long story short:
# 
#  - i$seq_clean: names of treatments per excel front end in order for all populations. use i$lookup$ipd$mol as reference table.
#  - i$seq_n: corresponding treatment numbers per named range r_pld_lookup_mol in excel
#  - i$seq_ref: reference name for pulling things out of R lists (e.g. p$drug[unlist(i$seq_ref$pop_0[1,])]) pulls pop 0 first sequence drug info IN ORDER :)
#
# This is automatically in line with the reference tables in the excel front end
# loaded at the time. If the ordering is changed there it needs updating in the IPD
# and in the lookup tables in the lists sheet of excel (and throughout excel!)
# 
# 
# 
# i.e, if Excel lookup tables are wrong, this will be wrong!!!
# 
# 


# 3.3.3 TSD14 survival analysis ------------------------------------------

# Now that treatment sequences are brought in and cleaned up ready for use, we
# can perform the survival analysis.
# 
# Use the function in Survival_functions.R to perform "simple" extrapolations
# on all pop line mol trial endpoint combinations with available data and return
# NULL for the rest

# Let's  perform some labelling like we did for treatment sequences for convenience/QC

i$surv$lab_pld <- list()

i$surv$lab_pld$population <- i$lookup$ipd$pop$Number
names(i$surv$lab_pld$population) <- i$lookup$ipd$pop$Description

i$surv$lab_pld$line <- i$lookup$ipd$line$Number
names(i$surv$lab_pld$line) <- i$lookup$ipd$line$Description

i$surv$lab_pld$molecule <- i$lookup$ipd$mol$Number
names(i$surv$lab_pld$molecule) <- i$lookup$ipd$mol$Description

i$surv$lab_pld$trial <- i$lookup$ipd$trial$Number
names(i$surv$lab_pld$trial) <- i$lookup$ipd$trial$Description

i$surv$lab_pld$endpoint <- i$lookup$ipd$endpoint$Number
names(i$surv$lab_pld$endpoint) <- i$lookup$ipd$endpoint$Description


# Now, put the data in a space and replace numbers with labels:

i$surv$lab_pld$dat <- i$surv$pld
i$surv$lab_pld$dat$population <- names(i$surv$lab_pld$population)[match(i$surv$lab_pld$dat$population,i$surv$lab_pld$population)]
i$surv$lab_pld$dat$line       <- names(i$surv$lab_pld$line)[match(i$surv$lab_pld$dat$line,i$surv$lab_pld$line)]
i$surv$lab_pld$dat$molecule   <- names(i$surv$lab_pld$molecule)[match(i$surv$lab_pld$dat$molecule,i$surv$lab_pld$molecule)]
i$surv$lab_pld$dat$trial      <- names(i$surv$lab_pld$trial)[match(i$surv$lab_pld$dat$trial,i$surv$lab_pld$trial)]
i$surv$lab_pld$dat$endpoint   <- names(i$surv$lab_pld$endpoint)[match(i$surv$lab_pld$dat$endpoint,i$surv$lab_pld$endpoint)]

# Now we have a labelled version which is a bit easier to QC.

# Note to debug it is very helpful to set verbose to TRUE below so that you can identify
# the datasets which are problematic (e.g. not converging, 0 time values)

i$surv$n_by_plmte <- i$surv$pld[, .N, by = list(population, line,molecule,trial,endpoint)] %>%
  arrange(population,line, molecule,trial,endpoint)

i$surv$n_by_plmte$population <- i$lookup$ipd$pop[match(i$surv$n_by_plmte$population       ,i$lookup$ipd$pop$Number),Description]
i$surv$n_by_plmte$line       <- i$lookup$ipd$line[match(i$surv$n_by_plmte$line       ,i$lookup$ipd$line$Number),Description]
i$surv$n_by_plmte$molecule   <- i$lookup$ipd$mol[match(i$surv$n_by_plmte$molecule       ,i$lookup$ipd$mol$Number),Description]
i$surv$n_by_plmte$molecule[which(is.na(i$surv$n_by_plmte$molecule))]   <- "Non-UK treatments (pooled)"
i$surv$n_by_plmte$trial      <- i$lookup$ipd$trial[match(i$surv$n_by_plmte$trial       ,i$lookup$ipd$trial$Number),Description]
i$surv$n_by_plmte$endpoint   <- i$lookup$ipd$endpoint[match(i$surv$n_by_plmte$endpoint       ,i$lookup$ipd$endpoint$Number),Description]

# The number of rows in this table is the number of SETS of regression analyses
# that are going to be run (each is 7 regressions)


# The below code runs the survival analysis and saves as an RDS file for upload, this will only run if you set 
# i$dd_run_surv_reg to "Yes" either in the Excel input or here

if (i$dd_run_surv_reg == "Yes") {
  
  i$surv$reg <- f_surv_runAllTSD14(
    r_pld             = i$surv$pld,
    id                = i$id$ipd,
    lookups           = i$lookup$ipd,
    draw_plots        = FALSE,
    distnames         = i$distnames,
    cl_y              = p$basic$cl_y,
    t_cyc             = p$basic$t_cyc,
    xlim_survplots_yr = p$misc$plot$xlim_survplots_yr,
    t_yr              = p$basic$t_yr,
    verbose           = qc_mode,
    min_obs           = 28
  )
  
  # now, there is very little information available on BSC overall survival,
  # for those people that decide they do not want further treatment
  #
  # The best data available currently is pooled PPS data on 4L patients, these
  # are then 5L+ patients and given that there are currently 4 lines of therapy
  # the proportion that receive something active subsequently is likely to be
  # small. Consequently, this is likely a pooled analysis which can inform
  # early (and 5L) BSC OVERALL SURVIVAL.
  #
  # Therefore the molecule 999 4th line PPS should be informed by a pooled analysis
  # of all molecules' PPS at 4th line. That is, i$surv$pld[line == 4 & trial == 2 & endpoint == 4,]
  # is the data that should inform endpoint 0 for all BSC.
  
  # MANUALLY RUN SURVIVAL FOR BSC PPS AS POOLED!!!
  
  i$surv$reg$pop_0$line_4$mol_999$trial_2$endpoint_4 <- lapply(1:1,function(x) {
    
    # Filter down to the parameters above associated with this combination:
    ipd <- i$surv$pld[line==4 & endpoint==4,list(timew,event_censor)]
    
    names(ipd) <- c("t","e")
    
    cat(paste0(
      "Survival analysis - population: ", i$lookup$ipd$pop[Number      == 0, Description],
      "\t line: "                       , i$lookup$ipd$line[Number     == 4, Description],
      "\t molecule: "                   , i$lookup$ipd$mol[Number     == 999, Description],
      "\t trial: "                      , i$lookup$ipd$trial[Number     == 2, Description],
      "\t endpoint: "                   , i$lookup$ipd$endpoint[Number == 4, Description], "\n"
    ))
    
    fs_fits <- lapply(i$distnames, function(dist) {  # applying all parametric survival curves in the list of distNames
      fs_fit <- flexsurvreg(
        formula = Surv(t, e) ~ 1,
        data = ipd,
        dist = dist
      )
      return(list(
        coefs = coefficients(fs_fit),                                         # coefficients for the fitted model
        vcov  = vcov(fs_fit),                                                 # variance covariance matrix for the fitted model
        fit   = c(AIC= AIC(fs_fit), BIC=BIC(fs_fit), logLik = logLik(fs_fit)) # goodness of fit statistics for the fitted model
      ))
    })
    
    gof <- do.call(rbind, lapply(i$distnames, function(dist) fs_fits[[dist]]$fit))
    
    st <- matrix(
      unlist(lapply(i$distnames, function(dist) {
        f_extrapolate(p$basic$t_cyc, fs_fits[[dist]]$coefs, dist)
      })),
      ncol = length(i$distnames),
      dimnames = list(NULL, i$distnames),
      byrow = FALSE
    )
    
    
    # curly braces on their own mean do this stuff and only return the last thing
    # or what's in a return call
    plot <- {
      # First the IPD is produced in a format that survminer will accept. Data must all be
      # the same format with the same column names.
      # this assumes no covariate adjustment
      
      sm_ipd <- f_ce_km_MakeDatSurvFriendly(
        Data_required = ipd,
        time_column   = "t",                 # note that this is taking IPD in weeks
        event_column  = "e",
        t_multiplier  = p$basic$cl_y             # data in weeks, cycle length in plot years
      )
      
      # get the survival analysis in the form we need for survminer
      # and make the extrapolations we need for survminer
      
      form          <- Surv(t, ec) ~ 1
      sm_surv_est   <- surv_fit(formula = form, data = sm_ipd)
      
      # make the plot with the input data:
      survival_plot <- suppressMessages(f_extrap_plot(
        SurvEstimate   = sm_surv_est,
        Data_required  = sm_ipd,
        curvefits_data = st,
        time_vector    = p$basic$t_yr,
        xlim           = p$misc$plot$xlim_survplots_yr,   #### this will need replacing dependent on how many years we decide to show per time horizon
        break_by       = round(20/8,0) #### this will need replacing dependent on how many years we decide to show per time horizon
      ))
      list(
        ipd     = sm_ipd,
        formula = form,
        plot    = survival_plot
      )
    }
    
    
    # Now that we've done everything for this dataset, return a list of the stuff
    # we need for it:
    return(list(
      pop      = i$lookup$ipd$pop[     Number == 0,Description],
      line     = i$lookup$ipd$line[    Number == 4,Description],
      mol      = i$lookup$ipd$mol[     Number == 999,Description],
      tr       = i$lookup$ipd$trial[     Number == 2,Description],
      endpoint = i$lookup$ipd$endpoint[Number == 4,Description],
      ipd      = ipd,
      fs_fits  = fs_fits,
      gof      = gof,
      st       = st,
      plot     = plot
    ))
  })[[1]]
  
  saveRDS(i$surv$reg, file = "./1_Data/Survival_analysis.rds")
  
}

# to load in pre-run survival analysis select the RDS file here

# option to load from pre-specified file path on local machine, uncomment this and comment out the line below to use

RDS_path <- "./1_Data/Survival_analysis_noTTDorTTPorPPS[NoACIC].rds"
if (file.exists(RDS_path)) {
  i$surv$reg <- readRDS(RDS_path)
} else {
  i$surv$reg <- readRDS(rstudioapi::selectFile(
    caption = "Please select 'Survival_analysis_noTTDorTTPorPPS[NoACIC].rds'",
    label = "Survival_analysis_noTTDorTTPorPPS[NoACIC].rds",
    path = "./1_Data/",
    filter = "R Files (*.rds)",
    existing = TRUE
  ))
}


# Limit to model time horizon

TH <- p$basic$th + 1

i$surv$reg <-lapply(i$surv$reg, function(popu) {
  lapply(popu, function(li) {
    lapply(li, function(mol) {
      lapply(mol, function(tr) {
        lapply(tr, function(endp) {
          if (is.null(endp$st)) {
            return(endp)
          } else {
            endp$st <- endp$st[1:TH,]
            return(endp)
          }
        })
      })
    })
  })
})


# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!
# Note: i$surv$reg$pop_0$line_4$mol_999$trial_2$endpoint_4 is used
# to inform ALL BSC OS. This will be decided in the EXCEL FILE, which 
# dropdowns for BSC OS should link to 4L PPS for mol 999
# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!



# Note that draw_plots will be a switch in the shiny application.
# In this case we draw plots because we need those plots later (for word output
# assisting with model selection)


# So that's all of the TSD14 survival analysis done. The next step is to programmatically 
# proliferate comparative efficacy 



# On a tablet with very little computational power this takes a couple of minutes to run. on a new
# laptop its not long at all


# So, to pull out the visual fit of the analysis of TTD for a population

# i$surv$reg$pop_0$line_1$mol_7$trial_0$endpoint_0$plot$plot
# i$surv$reg$pop_0$line_1$mol_1$trial_0$endpoint_3$plot$plot
# i$surv$reg$pop_0$line_1$mol_7$trial_0$endpoint_1$plot$plot

# More importantly, to pull a particular extrapolation:
# i$surv$reg$pop_0$line_1$mol_1$trial_0$endpoint_1$st[,"weibull"]

# Where the "weibull" part would come from a dropdown list in the Excel front-end of the model specific
# to that endpoint for that treatment for that line for that population (i.e. a lot of selections need to be made!)


# The stuff inside of i$surv$reg can be used to automatically populate a report presenting the full plot, regression summaries,
# goodness-of-fit results and the fit of the selected (via excel) distribution. The output can then be manually appended to include
# written justification for the selection(s) to drastically reduce the overhead associated with reporting survival analysis
# results and decisions made.


# 3.3.4 Survival analysis reporting ---------------------------------------

# the next step is to go through all of the results based directly on survival data
# and produce a readout containing:
# 
# - Regression summary tables
# - goodness of fit
# - extrapolations (short and long-term) for visual fit assessment
# 
# Each of these should have a separate section which at least states the identifiers
# (i.e., translating from numbers to text as in Section 3.4.2 above)
# 
# The best way to do this is with either Reduce or base for loops:
# 

# Produce all the KM, extrapolations and gof tables for decisions on the front-end

# Note whether or not the survival analysis report is run by the code is set in Excel as this takes a long time to produce
# This cannot be produced without access to PLD

if (i$dd_report_req_surv_reg=="Yes") {
  
  doc_surv <- f_surv_makeTSD14Report(
    fs_res = i$surv$reg,
    id     = i$id$ipd,
    lookup = i$lookup$ipd
  )
  print(doc_surv, target = "./4_Output/Survival_Analysis.docx")
  
  rm(doc_surv)
}


# 3.3.5 Comparative efficacy propagation (NMA) ---------------------------------------------------------------

# Pull in the data and calculate means by pop line mol endpoint reftrt and reftrial

# First read in RDS file containing the PH NMA coda samples


# 3.3.5.1.1 PH NMA data -----------------------------------------------------

# Option to read in PH NMA CODA from local machine, uncomment this and comment out the line below to use
RDS_path2 <- "./1_Data/PH_NMA_CODA.rds"
if (file.exists(RDS_path2)) {
  i$PHNMA <- readRDS(RDS_path2)
} else {
  i$PHNMA <- readRDS(rstudioapi::selectFile(
    caption = "Please select 'PH_NMA_CODA.rds'",
    label = "PH_NMA_CODA.rds",
    path = "./1_Data/",
    filter = "R Files (*.rds)",
    existing = TRUE
  ))
}


colnames(i$PHNMA$data) <- c("Run", "Population", "Line", "Molecule", "Endpoint", "Reference.treatment", "Reference.trial", "HR")
i$PHNMA$data$Reference.endpoint <- i$PHNMA$data$Endpoint


# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!

i$PHNMA$assume3L      <- i$PHNMA$data[Line==2,]
i$PHNMA$assume3L$Line <- 3
i$PHNMA$data         <- rbind(i$PHNMA$data,i$PHNMA$assume3L)

i$PHNMA$assumeTTD <- i$PHNMA$data[Endpoint==1,]
i$PHNMA$assumeTTD$Endpoint <- 2
i$PHNMA$data <- rbind(i$PHNMA$data,i$PHNMA$assumeTTD)

i$PHNMA$assumeTTP <- i$PHNMA$data[Endpoint==1,]
i$PHNMA$assumeTTP$Endpoint <- 3
i$PHNMA$data <- rbind(i$PHNMA$data,i$PHNMA$assumeTTP)

# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!

# Calculate the mean from the CODA samples for deterministic analysis

i$PHNMA$means <- i$PHNMA$data[,.(HR = mean(HR)),by=list(Population,Line,Molecule,Endpoint,Reference.treatment,Reference.trial)]


# 3.3.5.1.2 DETERMINISTIC CODA --------------------------------------------

# for the deterministic analysis we use the means. 
p$releff$CODA$PH <- i$PHNMA$means


# 3.3.5.2.1 FP NMA data -----------------------------------------------------

# Load in FP NMA data
i$FPNMA <- list()

#read in means for deterministic and PSA parameters for probabilistic

# option to read in from local machine, uncomment the below and comment out line 949 to use
RDS_path3 <- "./1_Data/FPNMA_means.rds"
if (file.exists(RDS_path3)) {
  i$FPNMA$means  <- readRDS(RDS_path3)
} else {
  i$FPNMA$means <- readRDS(rstudioapi::selectFile(
    caption = "Load in FP NMA CODA (FPNMA_means.rds)",
    label = "FPNMA_means.rds",
    path = "./1_Data/",
    filter = "R Files (*.rds)",
    existing = TRUE
  ))
}


#tidy means column names and timing
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "intervention_code"] <- "Molecule"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "reference_treatment_code"] <- "Reference.treatment"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "ref_trial_code"] <- "Reference.trial"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "population"] <- "Population"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "line"] <- "Line"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "endpoint"] <- "Endpoint"
colnames(i$FPNMA$means)[colnames(i$FPNMA$means) == "V1"] <- "HR"

i$FPNMA$means$time <- round(i$FPNMA$means$time * 52 / 12)

# means

# Rebasing to allow use of cabo as reference treatment in 2nd line
# repeats for means (stored in i which are later transferred to p) 
i$FPNMA$means <- f_rebase_for_cabo_as_ref_in_2L(FPNMAdata = i$FPNMA$means)

# Remove the now redundant objects we made in order to do this

# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
i$FPNMA$means <- f_3L_rel_effect_same_as_2L(FPNMAdata = i$FPNMA$means) 

# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!

# Create 1 row for each destination PLMTE, so that we know where to put the
# fp data without having to iterate much
i$FPNMA$destinations <- f_gen_destinations(fp_data = i$FPNMA$means)

# add in reference.trial 2
i$FPNMA$means <- f_add_reference_trial_2(fp_data = i$FPNMA$means)


# 3.3.5.2.2 DETERMINISTIC CODA --------------------------------------------

p$releff$CODA$FP <- i$FPNMA$means
p$releff$fp_dest <- i$FPNMA$destinations[!is.na(Molecule), ]

# limit to time horizon
p$releff$CODA$FP <- p$releff$CODA$FP[time <= p$basic$th, ]

#eliminate NAs in molecule
p$releff$fp_dest <- i$FPNMA$destinations[!is.na(Molecule), ]
# A note on i vs p ---------------------------------------------------------------

# P is for the parameters for one model scenario. the relative efficacy network is required
# in order to compute the s(t) for all the different PLMTEs we need to power the model with.
# Therefore, the samples which are used should go into p not i.
# 
# However, the full CODA samples are a different matter as particularly for the
# FPNMA these are large files and there's no need to copy paste this many times.
# 
# Instead when we get to the point of the releff network, THEN we can be putting 
# it into p. this is because for a particular probabilistic iteration,  
# scenario and so on we can pull through the right HRs to the right place!
# 

# 3.3.5.3 Empty relative efficacy network ---------------------------------

# Turn this into a list structure using the same naming convention as the rest of the model:
# 
# 
# For this population, population line molecule trial and endpoint, generate a list
# of spaces containing information on the relationship between other
# population line molecule trial and endpoint pairings and this one.
# 
# For example, we need to be able to do the following:
# 
# - HR applied to other subgroup for same line mol tr endpoint
# - HR applied to same subgroup for different line same mol tr endpoint
# - HR applied to same subgroup same line different mol different tr same endpoint
# - HR applied to same subgroup same line different mol same tr same endpoint
# - HR applied to same subgroup same line different mol tr endpoint
# 
# The best way to cope with all this is to basically list out where
# the extrapolation is coming from USING THE SAME NAMES AS IN i$surv$reg
# but with the addition of the selected distribution
# 
# So, we have a list with dest for destination (i.e. this extrapolation)
# origin (where it's coming from), and hr (what to apply to it)
# 
# The next step (a different function) populates orig
# 


p$releff$network <- f_NMA_generateNetwork(i$id$ipd,i$lookup$ipd)

# To visualize it a bit, the structure looks like tree roots. Like this:
# 1   Root                                    
# 2    ¦--pop_0                               
# 3    ¦   ¦--line_1                          
# 4    ¦   ¦   ¦--mol_0                       
# 5    ¦   ¦   ¦   ¦--trial_0                 
# 6    ¦   ¦   ¦   ¦   ¦--endpoint_0          
# 7    ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 8    ¦   ¦   ¦   ¦   ¦   °--orig            
# 9    ¦   ¦   ¦   ¦   ¦--endpoint_1          
# 10   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 11   ¦   ¦   ¦   ¦   ¦   °--orig            
# 12   ¦   ¦   ¦   ¦   ¦--endpoint_2          
# 13   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 14   ¦   ¦   ¦   ¦   ¦   °--orig            
# 15   ¦   ¦   ¦   ¦   ¦--endpoint_3          
# 16   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 17   ¦   ¦   ¦   ¦   ¦   °--orig            
# 18   ¦   ¦   ¦   ¦   ¦--endpoint_4          
# 19   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 20   ¦   ¦   ¦   ¦   ¦   °--orig            
# 21   ¦   ¦   ¦   ¦   ¦--endpoint_5          
# 22   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 23   ¦   ¦   ¦   ¦   ¦   °--orig            
# 24   ¦   ¦   ¦   ¦   ¦--endpoint_6          
# 25   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 26   ¦   ¦   ¦   ¦   ¦   °--orig            
# 27   ¦   ¦   ¦   ¦   °--endpoint_7          
# 28   ¦   ¦   ¦   ¦       ¦--dest            
# 29   ¦   ¦   ¦   ¦       °--orig            
# 80   ¦   ¦   ¦--mol_1                       
# 81   ¦   ¦   ¦   ¦--trial_0                 
# 82   ¦   ¦   ¦   ¦   ¦--endpoint_0          
# 83   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 84   ¦   ¦   ¦   ¦   ¦   °--orig            
# 85   ¦   ¦   ¦   ¦   ¦--endpoint_1          
# 86   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 87   ¦   ¦   ¦   ¦   ¦   °--orig            
# 88   ¦   ¦   ¦   ¦   ¦--endpoint_2          
# 89   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 90   ¦   ¦   ¦   ¦   ¦   °--orig            
# 91   ¦   ¦   ¦   ¦   ¦--endpoint_3          
# 92   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 93   ¦   ¦   ¦   ¦   ¦   °--orig            
# 94   ¦   ¦   ¦   ¦   ¦--endpoint_4          
# 95   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 96   ¦   ¦   ¦   ¦   ¦   °--orig            
# 97   ¦   ¦   ¦   ¦   ¦--endpoint_5          
# 98   ¦   ¦   ¦   ¦   ¦   ¦--dest            
# 99   ¦   ¦   ¦   ¦   ¦   °--orig            
# 100  ¦   ¦   ¦   ¦   °--... 2 nodes w/ 4 sub
# 101  ¦   ¦   ¦   °--... 2 nodes w/ 54 sub   
# 102  ¦   ¦   °--... 12 nodes w/ 956 sub     
# 103  ¦   °--... 5 nodes w/ 6288 sub         
# 104  °--... 3 nodes w/ 25463 sub             

# A function (or functions) is (are) required to do several things, IN THIS ORDER:
# 
#  1. Put the HRs from the PH NMA in the destinations, using the CODA sample identifiers to set the origin.
#  2. Put the time-varying HRs from the FP NMA in the destinations, using the identifiers to set the origins
#  3. Use the table R_table_eff_data_settings from Excel to apply any superseding & any assumption/ad-hoc stuff (HRs, assume equal to and so on)
#  4. Use the final network object to propagate relative efficacy throughout the network, producing a set of extrapolations with RE applied.
# 

# Generate DESTINATION trial number column if it doesn't already exist:
if(!"Trial" %in% colnames(i$R_table_eff_data_settings)) {
  i$R_table_eff_data_settings$Trial <- i$lookup$ipd$trial$Number[match(i$R_table_eff_data_settings$Trial.name.if.effectiveness.source.is.trial,i$lookup$ipd$trial$Description)]
}



# 3.3.5.4 Linking inputs from PH and FP NMAs --------------------------------------------


# Use the information we have from the PHNMA CODA sample to populate the corresponding
# places within p$releff$network
# 
p$releff$network <- f_NMA_linkPHNMA(
  network       = p$releff$network,
  hr_table      = p$releff$CODA$PH
)

# All treatments included in the network:
# unique(c(unique(p$releff$means$Molecule),unique(p$releff$means$Reference.treatment),unique(i$R_table_eff_data_settings$Origin.treatment),unique(i$R_table_eff_data_settings$Treatment)))

# Link in the fractional polynomial point estimate time-varying hazard ratios
p$releff$network <- f_NMA_linkFPNMA(
  network       = p$releff$network,
  destinations  = p$releff$fp_dest,
  hr_table      = p$releff$CODA$FP,
  time_horizon  = p$basic$th
)

# See e.g.
# p$releff$network$pop_0$line_1$mol_1$trial_0$endpoint_0
# 
# To see that the time-varying hazard has been passed along :)


# Remember that we only need to apply this one to the stuff that is NOT direct 
# survival analysis applied to the data. HOWEVER, we need those rows in the table
# as well as we need the "origin" distributional selections to fill in our output!



# 3.3.5.5 Assumptions from Excel ------------------------------------------

# Apply the setting dd_use_PHnma_for_FPnma from excel. This supersedes what's selected
# in efficacy settings, changing all FP NMA entries to PH NMA to force the model to use
# the PH NMA CODA sample over and above the FP NMA!

if (i$dd_use_PHnma_for_FPnma == "Yes") {
  i$which_fpnma <- which(i$R_table_eff_data_settings$Effectiveness.data.source == "FP_NMA")
  i$R_table_eff_data_settings$Effectiveness.data.source[i$which_fpnma] <- "PH_NMA"
}

# QC point: to test whether all entries to Effectiveness.data.source are in the
# lookup list, compare these two:
# 
# table(i$R_table_eff_data_settings$Effectiveness.data.source)
# 
# i$List_eff_datasources
# 
# i.e.,
# 
# setdiff(names(table(i$R_table_eff_data_settings$Effectiveness.data.source)),i$List_eff_datasources)
# 
# "0" is fine as that's empty cells, but there should be nothing else.
# 

# setdiff(names(table(i$R_table_eff_data_settings$Effectiveness.data.source)),i$List_eff_datasources)


p$releff$network <- f_NMA_AddAssumptionsToNetwork(
  network            = p$releff$network,
  phnma_table        = p$releff$CODA$PH,
  fpnma_table        = p$releff$CODA$FP,
  fpnma_destinations = p$releff$fp_dest,
  excel_table        = data.table(i$R_table_eff_data_settings),
  trial_flag         = i$List_eff_datasources[1],
  fpnma_flag         = i$List_eff_datasources[3],
  phnma_flag         = i$List_eff_datasources[2],
  et_flag            = i$List_eff_datasources[4],
  ahr_flag           = i$List_eff_datasources[5],
  verbose            = qc_mode
)


# Example to see all the HRs for population 0 in some nice tables 
# 
# f_releff_extract_all_HRs(p$releff$network)
# 
# assign to an object to browse through it


# 3.3.5.6 Propagating the network -----------------------------------------

# Now that we have collected together all of the survival extrapolations for all
# of the data we have, and all of the HRs, whether by assumption or derived
# via NMA, plus all of the underlying distributional selections, we can proceed
# to "proliferate" the evidence network, starting from the extrapolations we have already
# and working "outwards" from that point, arriving at ONLY the extrapolation for each
# possibility (rather than all distributions). From there, we can "cherry-pick"
# those population line molecule trial endpoint distribution HR application choices
# to "build-up" a given treatment pathway.


# These are still "raw" extraps, so should be in i as they're not used for 
# computation of pf. 

i$surv$extraps <- f_surv_getExtrapolations(regs = i$surv$reg)

# propagate the extrapolation and comparative efficacy data/assumptions:
# 
# Because this is for a particular iteration of the model, it should go into object p
# which is the input set for running one model.
# 
# 
p$surv$st <- f_releff_PropNetwork(
  network = p$releff$network,
  extraps = i$surv$extraps,
  dos     = 10,
  verbose = qc_mode,
  dist_lookups = p$basic$lookup$dist,
  excel_table = data.table(i$R_table_eff_data_settings)
)

# active trial depends on excel settings so it's hard to anticipate. 2 for all
# lines if rwe and 

# # pop 0 first line nivo mono trial 0
# f_qc_surv_gethead(p$surv$st,p = 0,l = 1,m = 1,t = 2,len = 3)
# 
# # pop 0 line 2 nivo mono trial 1
# f_qc_surv_gethead(p$surv$st,p = 0,l = 2,m = 0,t = 2,len = 3)
# 
# # follow pazo through the lines - so close! only OS in 4L but the others are all working
# # (just show me the first couple of rows from each for brevity). mol 5 selected
# # because it runs through all treatment lines
# f_qc_surv_gethead(p$surv$st,p = 0,l = 1,m = 5,t = 0,len = 3)
# f_qc_surv_gethead(p$surv$st,p = 0,l = 2,m = 5,t = 1,len = 3)
# f_qc_surv_gethead(p$surv$st,p = 0,l = 3,m = 5,t = 1,len = 3)
# f_qc_surv_gethead(p$surv$st,p = 0,l = 4,m = 5,t = 1,len = 3)
# 
# # Risk pop 1 first line nivo cabo from trial and then rwe:
# f_qc_surv_gethead(p$surv$st,1,1,1,2)


# 3.3.6 Extrapolation modifications ---------------------------------------

# This subsection focuses on applying the required modifications to the "raw" extrapolations
# now that we have propagated comparative efficacy. We need to first apply any adjustments for treatment effect waning 
# and then check the following:
# 
#  - mortality never falls below that of the general population
#  - The extrapolations which are NOT censored for death do not cross (in abs or hazard) OS

# During technical engagement functionality will be added here to adjust for the impact of prior adjuvant treatment as part of scenario analysis

if(i$dd_adjforprioradjuvant == "Yes") {
  p$surv$st <- f_surv_adjuvant_HR(
    st              = p$surv$st,
    adjuvant_impact = i$R_table_prior_IO_impact_eff,
    demo_table      = p$demo$table,
    lookup          = p$basic$lookup,
    verbose         = TRUE)
}



# Curve crossing should be fixed on the treatment sequence level, NOT the original PLMT level
# through E. This is because doing it that way round could lead to implausible
# extrapolations. For instance, if the PFS from somewhere else was being injected
# in (e.g. if trial 0 for OS and trial 1 from a different PLMT for PFS via equivalence
# assumption), it may be the case that the assumed equal PFS crosses OS.
# 
# Therefore, curve crossing fixes should be performed on the individual treatment
# pathways.
# 
# In conclusion, this section ONLY adjusts for gpop mortality, and curve crossing
# adjustment is performed AFTER f_seq_extrapCollector is used to pull through
# the correct extrapolations to the correct DESTINATION locations. 
# 
# The function is ONLY applied to OS lines.
# 
# Because we have pulled through the metatadata to our PLMTEs s(t) objects, 
# we can use the "dest" element to decide which row of Excel range R_table_ptchar
# to apply in the general population mortality
# 

# Firstly, treatment effect waning. This is done on unadjusted curves and is 
# not included in the model base-case.
#
# Make the information on patient characteristics from Excel into a data table

i$R_table_ptchar <- data.table(i$R_table_ptchar)





# BEFORE adjusting for general population, we must apply treatment effect waning.
#        Evidence has suggested that this predominantly affects 1L treatments,
#        and the table "R_table_TE_waning_settings" has been populated to capture
#        waning. This provides information on the "destination" and the nature
#        of the treatment effect waning. The "destination" information is then
#        used to link to the efficacy settings table "R_table_eff_data_settings"
#        in excel. This table contains the corresponding information on the "origin".
#        In the case of 1L treatments, this is the reference curve

##### Treatment effect waning application is defined in the main Excel workbook for input
##### For each treatment and outcome:
##### apply: yes / no
##### method: either absolute survival or hazards
##### when TE waning starts (years)
##### when TE waning is fully implemented (years)
##### linear graduation is used between the start and end time

# Note this function produces a warning which details when there is a danger of TE waning in the usual fashion (setting hazards the same)
# being implemented producing counterintuitive results. In this case we use the higher of the hazards to prevent this problem

if (sum(i$R_table_TE_waning_settings$apply.waning == "Yes") > 0) {
  p$surv$st <- f_surv_twaning_apply(
    st_list     = p$surv$st,
    tab_waning  = data.table(i$R_table_TE_waning_settings),
    tab_eff_set = data.table(i$R_table_eff_data_settings),
    verbose     = qc_mode
  )
}

# Generate all the gpop lines we can generate:
p$surv$gpop <- if (i$dd_age_sex_source == "Mean") f_surv_GenOSLines_det(
  R_table_ptchar         = i$R_table_ptchar,
  R_table_mort_lifeTable = i$R_table_mort_lifeTable,
  t_yr                   = p$basic$t_yr,
  lookups                = i$lookup
) else f_surv_GenOSLines_ipd(
  R_table_patientagesex  = i$R_table_patientagesex,
  R_table_mort_lifeTable = i$R_table_mort_lifeTable,
  t_yr                   = p$basic$t_yr,
  lookups                = i$lookup
)


# During QC model we produce a comparison between mean-based and IPD based general population OS lines:
if (qc_mode) {
  i$gpop <- list(
    means = f_surv_GenOSLines_det(
      R_table_ptchar         = i$R_table_ptchar,
      R_table_mort_lifeTable = i$R_table_mort_lifeTable,
      t_yr                   = p$basic$t_yr,
      lookups                = i$lookup
    ),
    ipd = f_surv_GenOSLines_ipd(
      R_table_patientagesex  = i$R_table_patientagesex,
      R_table_mort_lifeTable = i$R_table_mort_lifeTable,
      t_yr                   = p$basic$t_yr,
      lookups                = i$lookup
    )
  ) 
  
  # Compare all risk population first-line across the 2 methods
  
  i$gpop$plotdat <- data.table(
    t = rep(p$basic$t_yr,2),
    os = c(i$gpop$means$pop_0$line_1$os,i$gpop$ipd$pop_0$line_1$os),
    method = c(rep("Means",p$basic$th+1),rep("Patient data",p$basic$th+1))
  )
  
  i$gpop$comp_plot <- ggplot(i$gpop$plotdat, aes(x = t, y = os, colour = method)) + 
    geom_line() + 
    theme_classic() +
    theme(legend.position = "bottom", legend.title=element_blank()) + 
    labs(title = NULL, x = "Time (years)", y = "% Survival") + 
    scale_x_continuous(expand = expansion(mult = c(0,0.05))) + 
    scale_y_continuous(labels = scales::percent)
  
  if(qc_mode) {
    ggsave(
      filename = file.path("./4_Output/","gpop_1L_method_comparison.png"),
      plot = i$gpop$comp_plot,
      device = "png",
      units = "cm",
      width = 15
    )
  }
}

# adjust all OS lines in p$surv$st for gpop mortality
p$surv$st <- f_surv_gpopadjust(st      = p$surv$st,
                               gpop    = p$surv$gpop,
                               method  = "hazardmax",
                               verbose = qc_mode)

# Adjust PFS and TTD for OS - PFS and TTD cannot exceed OS
p$surv$st <- f_surv_PFSxOS(p$surv$st, if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})

p$surv$st <- f_surv_TTDxOS(p$surv$st, if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})

# Adjust TTP, PFS cannot go above TTP, this is done on absolute survival rather than allowing flexibility to look at hazards
p$surv$st <- f_surv_PFSxTTP(st = p$surv$st,method =  "abs")

# The below should produce a string of positive or 0s
# p$surv$st$pop_0$line_1$mol_7$trial_2$endpoint_3$st - p$surv$st$pop_0$line_1$mol_7$trial_2$endpoint_1$st

# Last assumption - 5L =4L. this moves over BSC when it comes after 4 active treatments.

p$surv$st <- lapply(p$surv$st, function(popu) {
  popu$line_5 <- popu$line_4
  return(popu)
})


# 3.3.7 Visual QC of final curves ---------------------------------------

# Here's an example of how to run the QC plots. These serve as a useful check
# later on too, because it reveals where some PLMTEs have not been populated.

# f_qc_surv_ExtrapPlot(
#   st   = p$surv$st,
#   popu = "pop_0",
#   li   = "line_1",
#   mo   = "mol_1",
#   tr   = "trial_2",
#   t_yr = p$basic$t_yr,
#   th = p$basic$th
# )
# 
# f_qc_surv_EstHazPlot(
#   st   = p$surv$st,
#   gpop = p$surv$gpop,
#   popu = "pop_0",
#   li   = "line_1",
#   mo   = "mol_1",
#   tr   = "trial_2",
#   t_yr = p$basic$t_yr,
#   th   = p$basic$th
# )


# Here is a QC method for visually having a look at the extrapolated survival
# after all relative efficacy has been applied. This creates a LOT of graphs. You can look through them by pressing the arrow on the plots window
if (qc_mode) {
  i$excel_destinations <- data.table(i$R_table_eff_data_settings)[Include.in.this.analysis.=="Yes",list(Population,Treatment.line,Molecule,Origin.trial,End.point)]
  i$surv$ExtrapSenseCheck <- Reduce(
    x = 1:nrow(i$excel_destinations),
    init = f_NMA_generateNetwork(i$id$ipd,i$lookup$ipd),
    accumulate = FALSE,
    f = function(prev, dest_row) {
      dat <- as.list(i$excel_destinations[dest_row,])
      
      d <- list(
        pop = paste0("pop_",dat$Population),
        line = paste0("line_",dat$Treatment.line),
        mol = paste0("mol_",dat$Molecule),
        trial = paste0("trial_",dat$Origin.trial),
        endpoint = paste0("endpoint_",dat$End.point)
      )
      
      cat(paste0(
        "Drawing plots: ",
        " row ", dest_row, " i$surv$ExtrapSenseCheck$",
        d$pop, "$",
        d$line, "$",
        d$mol, "$",
        d$trial, "$plots",
        "\n"
      ))
      
      p_extrap <- f_qc_surv_ExtrapPlot(
        st   = p$surv$st,
        popu = d$pop,
        li   = d$line,
        mo   = d$mol,
        tr   = d$trial,
        t_yr = p$basic$t_yr,
        th = p$basic$th
      )
      
      p_haz <- f_qc_surv_EstHazPlot(
        st   = p$surv$st,
        gpop = p$surv$gpop,
        popu = d$pop,
        li   = d$line,
        mo   = d$mol,
        tr   = d$trial,
        t_yr = p$basic$t_yr,
        th   = p$basic$th
      )
      
      plmt <- prev[[d$pop]][[d$line]][[d$mol]][[d$trial]]
      plmt$plots <- list(
        p_extrap = p_extrap,
        p_haz = p_haz
      )
      
      # Save 2 plots for each PLM available. Takes a while to run.
      if (qc_mode) {
        ggsave(
          filename = file.path("./4_Output",paste0("p_extrap_",paste(d[1:4],collapse="_"),".png")),
          plot     = plmt$plots$p_extrap,
          device = "png",
          units = "cm",
          width = 15
        )
        ggsave(
          filename = file.path("./4_Output",paste0("p_tp_",paste(d[1:4],collapse="_"),".png")),
          plot     = plmt$plots$p_haz,
          device = "png",
          units = "cm",
          width = 15
        )
      }
      
      prev[[d$pop]][[d$line]][[d$mol]][[d$trial]] <- plmt
      
      return(prev)
      
    }
  )
}




# 3.4 Preparation of p -------------------------------------

# !!!!!!!!!!!!!!!!!!!IMPORTANT, PLEASE READ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 
# Before the patient flow sheet can be computed, it is imperative to populate
# all of the elements of p which are required to compute it, such that the
# function f_pf_computePF() really only needs one argument, p, to run (with a 
# few extra arguments for verbosity and limiting which populations to run etc).
# 
# This is true of all analyses, including PSA. In the PSA p_psa will contain
# the components of p which require random number generation (RNG), whilst p 
# itself will be used for all the stuff that does not change.
# 
# Therefore, this section populates the "rest" of the model inputs (i.e. those
# that are not informing the sequences or disease models). This includes costs, 
# QALYs and AEs.
# 
# !!!!!!!!!!!!!!!!!!!IMPORTANT, PLEASE READ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 

# Passing along an input for the ps model in case it is needed. this assumes
# prop of pre-progression deaths.
p$surv$prop_deathinPFS <- i$dd_prop_deathinPFS


# 3.4.1 Demographics ------------------------------------------------------

# Demographics are simply processed from the tables in Excel.

p$demo$agg <- f_cleaning_ptchar(i$R_table_ptchar, i$lookup)

# Deterministic version is very easy.
p$demo$live <- p$demo$agg

# 3.4.2 QALYs -------------------------------------------------------------

# Utilities are applied to the disease model by treatment by line and whether the patient is on or off treatment
# Age adjustment is conducted multiplicatively in line with DSU guidance using earlier defined patient characteristics for age and sex

# Extracting from excel file 

p <- add_population_utility_params(p, psa = FALSE, .i = i)

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
    
    p$util$gpop <- lapply(ptc_L1, function(pop) adjust_utility(
      age            = pop$age,
      sex            = pop$sex,
      utilities      = base_utility,
      .patient_level = FALSE,
      .p             = p
    ))
    
  } else {
    # We will only include IPD from line 1, since the population
    # norm is applied according to absolute model time rather than
    # than time in state. We don't know which population they are
    # in, so we will replicate for pop_0, pop_1 and pop_2.
    ipd_L1 <- i$R_table_patientagesex$Line == 1
    p$util$gpop <- list()
    p$util$gpop$pop_0 <- adjust_utility(
      age            = i$R_table_patientagesex$Age[ipd_L1],
      sex            = if_else(i$R_table_patientagesex$Gender[ipd_L1] == "M", "male", "female"),
      utilities      = base_utility,
      .patient_level = TRUE,
      .p             = p
    )
    p$util$gpop$pop_1 <- p$util$gpop$pop_0
    p$util$gpop$pop_2 <- p$util$gpop$pop_0
  }
} else {
  p$util$gpop <- list(pop_0 = 1, pop_1 = 1, pop_2 = 1)
}                 

# Remove the object we made which is not used again
rm(base_utility)

i$QALYs <- list()

i$QALYs$utilities$means <- f_process_utilities(raw_utilities = i$R_table_util,
                                               PSA = FALSE,
                                               samples = FALSE)
# Sample code for PSA - temporariliy in file probabilistic_model_DEV.R

# For the deterministic analysis, pass the means into p for use in the model.
# We now have our population norms for calculating gpop utility multiplier
# compared to baseline as well as our HSUVs by PLM and treatment status. everything
# we need.
# 
p$util$mk <- data.table(i$QALYs$utilities$means)


# 3.4.3 AEs ---------------------------------------------------------------

# the method to apply AEs (one-off or per cycle) is defined in the Excel inputs file
# this applies to both costs and QALYs

# options related to the source of AE data are defined in the Excel workbook
# the functions below read in the information required to produce AE costs and utilities per cycle, this includes AE durations
# and the trial durations associated with the rates which are used to produce one off cost and QALY impacts

p$ae$aetype <- i$dd_apply_AE_options

# when one off we need trial duration:
p$ae$duration <- data.table(i$R_table_AE_rates)

# produce the per cycle impact for each of the AEs

p$ae$mk$per_cycle <- data.table(f_process_adverse_events(
  AE_costs = i$R_table_AE_costs,
  AE_disutil = i$R_table_AE_util,
  AE_duration = i$R_table_duration,
  AE_rate = i$R_table_AE_rates,
  comparators = i$lookup$trt,
  weeks_per_year = p$basic$cl_w / p$basic$cl_y,
  PSA = FALSE
))

# they all match up with column RCC_input_desc from Excel except BSC (mol_999). 
# Set BSC (assumed to have 0 AE impact):
p$ae$mk$per_cycle[trt == "BSC",]$trt <- i$lookup$ipd$mol[match(999,Number),]$RCC_input_desc

# Convert to numbers. Now it's ready for use in the patient flow sheet.
p$ae$mk$per_cycle$molecule <- i$lookup$ipd$mol[match(p$ae$mk$per_cycle$trt,RCC_input_desc),]$Number

# Add in the AE approach switch
p$ae$approach <- i$dd_apply_AE_options

# 3.4.4 Costs -------------------------------------------------------------

# Drug and admin costs and MRU costs are considered within this section
# Drug and admin costs are applied per treatment per line (as drug costs may differ depending on what line treatment is used at)
# The impact of stopping rules is consider as part of the calculation of drug and admin costs rather than in determining whether 
# patients in the on or off treatment health states
# MRU costs are applied per treatment and by on and off treatment status as the EAG was advised that MRU is the same across different lines
# When this model is expanded to a generic version flexibility to define per line will be added
# One off costs are included for treatment initiation at each line of treatment, terminal care at the end of life (applied on death) 
# and progression (radiotherapy and surgery costs)

i$cost <- list()

# Put the deterministic drug cost inputs into p:
# 
p$costs$mk <- f_process_cost_data(
  drug_and_admin  = i$R_table_drug_admin_costs,
  per_cycle_costs = i$R_table_MRU,
  time_horizon    = p$basic$th,
  max_trt_lines   = p$basic$R_maxlines,
  RDI_source      = i$dd_sc_RDI,
  verbose         = FALSE,
  PSA             = FALSE, 
  samples         = 1)

# For PSA, You can pull psa iteration like this (done with random as test):
# psa_it <- round(runif(1)*1000,0)
# ooc_psa <- one_off_costs[,c("Type.of.cost", "Type", "Apply.to",paste0("V",psa_it)),with = FALSE]
# setnames(ooc_psa, paste0("V",psa_it),"cost")
# print(ooc_psa)
# 
# Alternatively, passs the whole table into PSA version of model as it's small
# data. Can't do that for survival obviously, but we can for smaller data.
# 
p$costs$oneoff <-  f_process_other_cost_data(
  one_off_costs = i$R_table_MRU_oneoff,
  PSA = FALSE
)

# pull out the individual inputs required for the ST model. In a PSA run
# this would be the nth iteration.
p$costs$oneoff_mk <- p$costs$oneoff[,.(cost = sum(cost)),by=list(Apply.to)]
p$costs$oneoff_mk <- lapply(structure(
  p$costs$oneoff_mk$Apply.to,
  .Names = p$costs$oneoff_mk$Apply.to
), function(x) {
  p$costs$oneoff_mk[Apply.to == x, ]$cost
})

# These are then added to the first element of the cost vectors in the patient flow
# function.

#holding line for PSA - this will be replaced by function sampling from means and SEs

if (FALSE) {
  i$cost$drug_and_admin_cost_by_tunnel_state$PSA <- f_process_drug_and_admin_cost_data(
    raw_table = i$R_table_drug_admin_costs,
    PSA_samples = TRUE)
}  



# 3.4.4 Subsequent treatment -------------------------------------------------------------

# Read in cost and QALY consequences for subsequent treatments per first line option from Excel
# This information is only used whaen the PartSA model structure is selected

p$substrt$partsa  <- as.data.table(f_process_subs_txt_data(
  subs_txt = i$R_table_sub_txt_cum_costs,
  PSA             = FALSE
))


# 3.5 Population mapping --------------------------------------------------

# There are 6 populations. These are combinations of 
#   
#   - the risk populations (3 of them)
#   - whether or not patients have had prior adjuvant therapy with immune-oncology treatments (2 of them)
#   
#   This makes for a possible 6 combinations, which are mapped to each other
#   in the excel table i$r_overall_lookup_pop
#   
#   Currently the treatment sequences are sorted into 4 different populations. (combination)
#   Currently the final survival extrapolations are sorted into 3 different populations (risk)
#   Currently the costs, QALYs and AEs are sorted into 3 different populations (risk)
#   
#   Fortunately this table allows us to simply refer to the appropriate populations
#   whilst calculating the patient flow. As an example:
#   
#   for overall population 1:
#   
#   seq_pop  <- p$basic$lookup$pop_map[1,]$Sequencing.population.number
#   risk_pop <- p$basic$lookup$pop_map[1,]$Risk.population.number
#   
#   Then use seq_pop to pull sequences & risk_pop for everything else.
#   
#   

p$basic$lookup$pop_map <- data.table(i$r_overall_lookup_pop)

# SECOND NOTE: In later treatment lines the risk population is always 0.
# PLEASE READ  To prevent replicating more and more data, the model simply
#              pulls from (risk) population 0 for later lines within the function
#              f_seq_extrapCollector, which "collects" the correct extrapolations.
#              This simply pulls from (risk) population 0 if line is more than 1
#              
#              THIS IS A MODELLING ASSUMPTION AS WE WERE INFORMED THAT RISK IS NOT
#              MEASURED AT 2L AND PRIOR RISK STATUS DOES NOT IMPACT TREATMENT
#              OUTSIDE OF PRIOR THERAPY RECEIVED. A FULLY GENERIC MODEL WOULD NOT MAKE THIS ASSUMPTION
#              AND INSTEAD WOULD HAVE A LARGER TABLE IN THE EFFICACY SETTINGS
#              SHEET COVERING THE WHOLE MATRIX OF POP/LINE
#              



# 3.6 PATIENT FLOW ------------------------------------------------------
# Now that we have all of the disease evidence and modeled outcomes available,
# we can compute the disease model. In the deterministic case, this will simply 
# be a trace for both the PS and Markov models. 

i$R_table_eff_data_settings <- data.table(i$R_table_eff_data_settings)
p$releff$excel_table <- i$R_table_eff_data_settings


if(!str_trim(i$dd_model_struct) %in% c("State transition", "Partitioned survival")) stop(
  "The model structure is not one of the available options. These are 'State transition' and 'Partitioned survival'. Either code or the excel file are out of date or wrongly set up"
)


# The process of computing patient flow is as follows:
# 
#   - Create an empty object for housing the patient flow (a list)
#   - Populate the trace with either disease modelling approach
#   - Record metadata on the approach, as this determines what the trace looks like
#   - Compute costs based on time in state
#   - Compute QALYs based on time in state
#   - Compute AEs based on time in state
#     - Costs
#     - QALYs
#   - Return the populated pf list
# 


pf <- list()
res <- list()


# 3.6.1 Run the model -----------------------------------------------------

# Running the model can mean two different structures, a markov model which 
# we refer to as state transition or a partitioned survival approach. These
# both work using the same parameters object p. There is only one function
# to compute patient flow, f_pf_computePF. This takes several different arguments
# but the main one is p, or the parameters object to power the model with. 
# Everything within f_pf_computePF uses objects which are within p, such that
# a copy of p could be saved as a file and as long as all functions have been
# defined and packages loaded the model can be run directly from there. 

# The only structures allowed are state transition and partitioned survival. error
# if it is not one of those:
stopifnot(p$basic$structure %in% c("State transition","Partitioned survival"))

# Check the decision problem. If it's cabo+nivo only the first 3 overall
# populations are relevant as one cannot get cabo+nivo in pops 4-6
if(p$basic$decision_problem == "cabo+nivo") {
  p$basic$pops_to_run <- 1:3
} else {
  p$basic$pops_to_run <- NULL
}

# populate the pf object irrespective of model structure or overall populations
# to include.
# 
# Note that if you put n_cores as NULL or 1 then the model will run in sequence.

# For QC, you can easily just save i and p as a file so all you need to do is
# load libraries to repeat the results and QC things:

saveRDS(p,"./2_Scripts/standalone scripts/QC/p.rds")
saveRDS(i,"./2_Scripts/standalone scripts/QC/i.rds")
# 
# If you run the model using scenario 1 (rather than base case scenario 0)
# It will be a partitioned survival model. in that case it's a good idea to 
# save a backup of p and i under a different name so that you have an input
# set for the paritioned survival model to hand at any time.
# 
# saveRDS(p,"./2_Scripts/standalone scripts/QC/p_PS.rds")
# saveRDS(i,"./2_Scripts/standalone scripts/QC/i_PS.rds")
                          
# p <- readRDS("./2_Scripts/standalone scripts/QC/p.rds")
# i <- readRDS("./2_Scripts/standalone scripts/QC/i.rds")

# Make this NA to run single core:
tick <- Sys.time()
pf <- f_pf_computePF(
  p           = p,
  struct      = p$basic$structure,
  verbose     = FALSE,
  plots       = FALSE,
  just_pop    = p$basic$pops_to_run,
  just_nlines = NULL,
  just_seq    = NULL
)
print(Sys.time() - tick)


# if (is.na(keep_free_cores)) {
#   plan(sequential)
# } else {
#   plan(multisession(workers = max(availableCores()-keep_free_cores,1)))
# }




# 3.6.2 Compiling model results -------------------------------------------

# Depending on which model structure is used to compute patient flow, the results
# are processed in different ways. The below function simply routes the pf object
# to the appropriate function to process the results, and returns an object
# res containing the results:

Scenario_number <- i$R_Scenario_num

if(Scenario_number == 0) {detail_level <- 5} else {detail_level <- 4}
# Detail level guide:

# State transition
# 1 is top line results by sequence (costs, QALYs, LYs) and weighted average results
# 2 is just the top line table and non-dominated weighted average incremental results
# 3 includes weighted average pairwise comparisons
# 4 includes incremental analysis by sequence 
# 5 includes trace plots  

# PartSA analysis
# 1 is top line results for the PartSA analysis (costs, QALYs, LYs)
# 2 is incremental analysis as well
# 3 includes breakdown tables as well
# 4 includes state residency plots

res <- f_res_compute_results(
  pf              = pf,
  structure       = p$basic$structure,
  p               = p,
  detail_level    = detail_level,
  vs_mol          = 1,
  no_active_lines = i$R_max_trt_lines
)




# 3.6.3 Calculate severity modifier -----------------------------------

# The severity modifier for the weighted average first-line treatment comparison
# uses the best available treatment which is not nivo cabo (molecule 1) for the
# first 3 populations. 
# 
# This is because this is the best (i.e. most discounted QALYs) available 1L trt
# pathway set.
#
# Calculation is only provided for the state transition model


if (p$basic$structure == "State transition") {
  population_numbers <- if(sum(p$basic$pops_to_run == 1:3)>0){1:3} else{1:6}
  res$mk$qaly_shortfall_1_to_3 <- lapply(population_numbers, function(npa_pop) {
  
  lu_pop <- p$basic$lookup$pop_map
  lu_rpop <- p$basic$lookup$ipd$pop
  
  # npa_pop is overall population, we need to look up risk population from it:
  
  risk_pop_n <- lu_pop[match(npa_pop,lu_pop$Overall.population.number),]$Risk.population.number
  risk_pop <- lu_rpop[match(risk_pop_n,lu_rpop$Number),]$Description  
  
  i$R_table_ptchar <- as.data.table(i$R_table_ptchar)
  
  if (i$dd_age_sex_source == "Mean") {
    
    # So for this risk population, we need the baseline characteristics:
    bl_chars <- i$R_table_ptchar[Population == risk_pop & Treatment.line == 1,]
    bl_age  <- bl_chars$Starting.age..years..Mean
    bl_male <- 1-bl_chars$Starting...female.Mean
    
  } else {
    
    patient_sex_age_IPD <- as.data.table(i$R_table_patientagesex)
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
    .patient_level = if(i$dd_age_sex_source == "Mean") {FALSE} else {TRUE},  
    qalys = q_met,
    .i = i,
    .p = p
  )
  
  out <- cbind(out, SOC = comp_no_met)
  
  return(out)
  
})
}




# 3.6.4 Saving the results ------------------------------------------------------

# the results are in list format so should be saved as an R list file (.rds)
# The naming of the file should reflect the model structure used. The file
# produced has a time stamp to avoid overwriting previous results files.

Scenario_name <- i$R_Scenario_name    # Use ST for state transition, PS for Partitioned survival, LP for list price, cPAS for cPAS
Scenario_number <- i$R_Scenario_num
Run_date <- date()
if (p$basic$structure == "State transition") {
  saveRDS(res, paste0("./4_Output/ST_Scenario ",Scenario_number,"_",i$dd_drug_price_options,gsub(":","_",Run_date),".rds"))
} else {
  saveRDS(res, paste0("./4_Output/PartSA_Scenario ",Scenario_number,"_",i$dd_drug_price_options,gsub(":","_",Run_date),".rds"))
}












# 3.6.5 Outputting results to Word ------------------------------------------------------

# Outputting the results to word requires a series of functions to produce component
# parts to go in the report, and then consolidating them all into the requested output.

# Get the functions to produce the word document output:

# Produce an automatically generated word document with required results tables
# and automatically name it using the variables that are already inside of 
# i and p (including structure, which is located in p$basic$structure. Note
# that this comes from i$dd_model_struct, i.e. named range dd_model_struct from excel.
# please ensure that this is updated for PS model scenarios):

Word_width_inches      <-  29.7*0.3937

f_res_ProduceWordDoc(
  p                      = p,
  res                    = res,
  Scenario_name          = i$R_Scenario_name,
  Scenario_number        = i$R_Scenario_num,
  price_options          = i$dd_drug_price_options,
  Run_date               = Run_date,
  word_template_location = "./3_Functions/reporting/empty results doc.docx",
  Word_width_inches      = 29.7*0.3937,
  auto_save              = TRUE,
  verbose = TRUE
)




# END OF CODE -------------------------------------------------------


#### Additional changes were originally planned during Phase 2 of this pilot following use for the initial decision problem including
# - Addition of Shiny user interface
# - Genericisation of the code to allow wider use
# - Programming and analysis of model outputs related specifically to sequencing, this may include value of information analyses

# Unfortunately funding for this has not been confirmed currently. 
# If you are interested in discussing or funding further development please contact the PenTAG team at pentag@exeter.ac.uk


