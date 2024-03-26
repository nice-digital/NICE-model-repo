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

#### 2. Loading functions ###########


# This variable is used throughout the model to define whether to provide additional outputs useful for QC or not
# The model will take longer to run when this is set to TRUE
qc_mode <- FALSE


# This function allows parallel processing (similar to future apply)
plan(multisession)

options(crosstable_units="cm")

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




################# START REPORT GENERATION #######################

# note - change stem to wherever your sharepoint files are kept.  Subsequent lines should work ok
# make sure you have saved the input files into your data folder and changed the names of these in the code below
# search for path to find input files

stem <- "C:/Users/dl556/"

rds_file_path <- paste0(stem,"University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Results rds files")
output_file_path <- "./4_Output"
scenario_files_path <- paste0(stem,"University of Exeter/PATT Pathways - Renal cell carcinoma - Documents/04-2 ERG Model/Inputs front end file/Scenario set up files")

#get list of rds files
files <- list.files(rds_file_path)
files <- files[grepl("Scenario", files)]

for (f in files) {
  # for each rds - does the relevant word file exist? (i.e. already run).  If yes skip.
    # update list of output files
      output_files <- list.files(output_file_path)
      output_files <- output_files[grepl("Scenario", output_files)]
      
    # extract scenario number and PAS/list prices
      scenario <- substr(f, 9, regexpr("price", f)[1]-2)
      # check if word file already exists, if so skip
      if (length(output_files)>0) {
        if (sum(grepl(scenario, output_files)) > 0) {
          cat(f,": word output report already exists.  Skipping.\n")
          next
        }
      }
  #match to correct input file (scenario number and PAS/list)
  scenario_nr <- as.numeric(substr(scenario,1,regexpr("_", scenario)[1]-1))
  price_type <- substr(scenario,regexpr("_", scenario)[1] + 1, nchar(scenario))
  
  if (price_type == "List") {
    price_directory <- "/List price"
  } else if (price_type == "PAS") {
    price_directory <- "/PAS price"
  } else {
    stop("Failed to itentify price type for scenario")
  }
  
  excel_path <- paste0(scenario_files_path, price_directory,"/Scenario ", scenario_nr, ".xlsm")
  
  if (!file.exists(excel_path)) {
    stop("Unable to find scenario inputs file ", excel_path)
  }
    
  # The first part of this code pulls all of the named ranges from the excel workbook, expand the parameters table
  
  i <- f_excel_extract(excel_path, verbose = TRUE)
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
  
  # PSA is not yet included in this model, this will be added during technical engagement
  
  # For the future: this is what a PSA would look like. Note that each p in p_psa
  # needs to pull through one probabilistic iteration of all parameters.
  # 
  # future_lapply(1:n_psa, function(psa_iteration) {
  #   generate_pf(p_psa[[psa_iteration]], structure=i$...)
  # })
  # 
  
  
  p <- f_misc_param_generate_p(i)

  
  # Pass this into p so that p can be used to exclusively compute the model:
  p$basic$decision_problem <- i$decision_problem
  
  
  
  
  i$surv <- list()
  
  #### Read in survival data from Excel workbook 
  
  # Pull out the raw data from the IPD excel book - one named range per treatment at each line
  # Each reference curve is defined in Excel as time (weeks), event/censor (event coded as 1, censor as 0), patient group, line, molecule, trial and endpoint
  # Pull all of the named ranges from the excel workbook, expand the parameters table
  
  excel_path2 <- "./1_Data/IPD_Confidential _Rinput.xlsx"
  
  wb <- f_excel_extract(excel_path2, verbose = TRUE)
  
  i$surv$pld <- as.data.table(wb$`_xlnm._FilterDatabase`)
  
  rm(wb)
  
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
  
  
  p$basic$lookup$pop_map <- data.table(i$r_overall_lookup_pop)
  
  
  res <- readRDS(paste0(rds_file_path,"/",f))
  
  Scenario_number <- i$R_Scenario_num
  
  Scenario_name <- i$R_Scenario_name    # Use ST for state transition, PS for Partitioned survival, LP for list price, cPAS for cPAS
  
  Run_date <- date()
  
  # 3.4.2 QALYs -------------------------------------------------------------
  
  # Utilities are applied to the disease model by treatment by line and whether the patient is on or off treatment
  # Age adjustment is conducted multiplicatively in line with DSU guidance using earlier defined patient characteristics for age and sex
  
  # Extracting from excel file 
  
  p <- add_population_utility_params(p, psa = FALSE, .i = i)
  
  # Pre-calculate the population utility norms since they will be the same across
  # all sequences (though may vary across populations), and store in p
  
  i$R_table_ptchar <- as.data.table(i$R_table_ptchar)
  
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
  
  i$QALYs <- list()
  
  i$QALYs$utilities$means <- f_process_utilities(raw_utilities = i$R_table_util,
                                                 PSA = FALSE,
                                                 samples = FALSE)
  # Sample code for PSA
  
  if (FALSE) {
    i$QALYs$utilities$PSA <- f_process_utilities(i$R_table_util,
                                                 PSA = TRUE,
                                                 samples = 100)
  }
  
  # For the deterministic analysis, pass the means into p for use in the model.
  # We now have our population norms for calculating gpop utility multiplier
  # compared to baseline as well as our HSUVs by PLM and treatment status. everything
  # we need.
  # 
  # For probabilistic settings we'd need to pass along the nth iteration of i$QALYs$utilities$PSA
  # into the nth element in p_psa (i.e. copy of p). Due to size, this may need to be done
  # one iteration at a time, in which case p_psa <- p and then replacing elements, then
  # running the model saves a lot of memory.
  p$util$mk <- data.table(i$QALYs$utilities$means)
  
  f_res_ProduceWordDoc(
    p                      = p,
    res                    = res,
    Scenario_name          = i$R_Scenario_name,
    Scenario_number        = i$R_Scenario_num,
    price_options          = i$dd_drug_price_options,
    Run_date               = date(),
    word_template_location = "./3_Functions/reporting/empty results doc.docx",
    Word_width_inches      = 29.7*0.3937,
    auto_save              = TRUE
  )

}

# END OF CODE -------------------------------------------------------


#### A note on planned model changes
# This model is currently undergoing internal and external QC via the NICE DSU
# The following changes are planned to the code following technical engagement:
# - Incorporation of PSA
# - Incorporation of functionality to allow scenario analyses around the impact of prior adjuvant IO treatment


#### Additional changes will be made during Phase 2 of this pilot following use for the initial decision problem including
# - Addition of Shiny user interface
# - Genericisation of the code to allow wider use
# - Programming and analysis of model outputs related specifically to sequencing, this may include value of information analyses


