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

#note - change stem to wherever your sharepoint files are kept.  Subsequent lines should work ok
stem <- "E:/"

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
  
  # Set seed for PSA
  set.seed(1475)
  
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
  
  
  # 3.6.2 Calculate severity modifier -----------------------------------
  
  # The severity modifier for the weighted average first-line treatment comparison
  # uses the best available treatment which is not nivo cabo (molecule 1) for the
  # first 3 populations. 
  # 
  # This is because this is the best (i.e. most discounted QALYs) available 1L trt
  # pathway set.
  #
  # Calculation is only provided for the state transition model
  
  if(p$basic$decision_problem == "cabo+nivo") {
    p$basic$pops_to_run <- 1:3
  } else {
    p$basic$pops_to_run <- NULL
  }
  
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
  
  
  # 3.6.5 Outputting results to Word ------------------------------------------------------
  
  landscape <- prop_section(
    page_size = page_size(orient = "landscape")
  )
  
  portrait <- prop_section(
    page_size = page_size(orient = "landscape")
  )
  
  
  # Make a word document and add a table of contents to it with 4 levels.
  # This uses a template word document so that styles and so on are the same.
  
  i$dd_producewordreport <- "Yes"
  
  if (i$dd_producewordreport == "Yes") {word_results <- TRUE} else {word_results <- FALSE}


  if (word_results == TRUE) {
    
    doc_res <- read_docx("./3_Functions/reporting/empty results doc.docx")
    
    # Add a 1st level header for the overall population:
    doc_res <- doc_res %>% 
      body_add_par(paste0("Results of Model Run in R Scenario name: " , Scenario_name),style = "heading 1")  %>%
      body_add_par(paste0("Date and time run: ", Run_date))
    
    
    Word_width_inches <- 29.7*0.3937 # width of the side borders in the word_document output (in centimeters)
    
    # Producing report tables (state transition model) ------------------------------------------------------
    
    
    # Make a word document containing results tables using the object res
    
    # Produces a different format depending on model structure
    
    if(p$basic$structure=="State transition") {
      
      # make a table by overall population for the summary results
      ft_basic_bop <- do.call(rbind,lapply(structure(names(res$wa_summarised),.Names=names(res$wa_summarised)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- res$wa_summarised[[popu_txt]]
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Risk.population
        
        # popu$seq_pop <- seq_popu_lab
        popu$risk_pop <- rsk_popu_lab
        popu$L1 <- lu_mol[match(popu$L1,lu_mol$Number),]$Description
        
        return(popu)
        
      }))
      
      # Now do the same thing for the incremental analysis
      ft_wa_inc<- do.call(rbind,lapply(structure(names(res$weighted_incremental),.Names=names(res$weighted_incremental)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        if (is.null(res$weighted_incremental[[popu_txt]]$non_dominated)) {
          popu <- as.data.table(res$weighted_incremental[[popu_txt]])
          popu <-  data.table(popu, ic = 0, iq = 0, il = 0, ICER = "Dominant" )
          popu$str_dom <- NULL
          
        } else {
          popu <- as.data.table(res$weighted_incremental[[popu_txt]]$expanded_results)
          popu$ICER[popu$extdom==FALSE] <- as.character(paste0("£",round(popu$ICER[popu$extdom==FALSE] ,0))) 
          popu$ICER[popu$extdom==TRUE] <- "(ext dominated)" 
          popu$str_dom <- NULL
          popu$extdom <- NULL
          popu$r <- NULL
          
        }
        
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      
      # Create table combining pairwise and incremental ICERs 
      
      setDT(ft_basic_bop)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      ft_basic_bop <- ft_basic_bop[order( risk_pop, costs)] # order by increasing costs
      
      setDT(ft_wa_inc)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")] 
      ft_wa_inc <- ft_wa_inc[,c(1,3,2,4,5,6,7,8,9)]
      
      ft_pairwise <- do.call(rbind,lapply(structure(names(res$pairwise_vs_mol),.Names=names(res$pairwise_vs_mol)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- res$pairwise_vs_mol[[popu_txt]]
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Risk.population
        
        # popu$seq_pop <- seq_popu_lab
        popu$risk_pop <- rsk_popu_lab
        
        return(popu)
        
      }))
      
      
      
      setDT(ft_pairwise)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")] 
      
      ft_pairwise$Pairwise_ICER <- ft_pairwise$icer
      
      ft_pairwise$Pairwise_ICER[is.na(ft_pairwise$icer) != TRUE] <- as.character(paste0("£",round(ft_pairwise$icer[is.na(ft_pairwise$icer) != TRUE] ,0)))
      
      
      ft_pairwise[ft_pairwise$icer < 0 & ft_pairwise$iq < 0]$Pairwise_ICER <- "Cabo+nivo dominated"
      ft_pairwise[ft_pairwise$icer < 0 & ft_pairwise$iq > 0]$Pairwise_ICER <- "Cabo+nivo dominant"
      ft_pairwise[ft_pairwise$icer > 0 & ft_pairwise$iq < 0]$Pairwise_ICER <- paste0("SW quadrant ",ft_pairwise[ft_pairwise$icer > 0 & ft_pairwise$iq < 0]$Pairwise_ICER)
      
      ft_pairwise <- ft_pairwise[,.SD,.SDcols = c("L1", "costs", "qalys", "ly", "Pairwise_ICER","risk_pop")]
      
      ft_wa_inc <- merge(ft_pairwise, ft_wa_inc, all.x = TRUE)
      
      ft_wa_inc <- ft_wa_inc[,c(1,2,4,3,7,9,8,6,10,5)]
      
      ft_wa_inc <- ft_wa_inc[order( risk_pop, costs)] # order by increasing costs
      
      ft_wa_inc[is.na(ICER)]$ICER <- "(dominated)"
    
      
      # Pull out nearest comparators
      
      comparator_no_allrisk <-  ff_closest_comparator(res,"pop_1")
      comparator_no_favrisk <- ff_closest_comparator(res,"pop_2")
      comparator_no_IPrisk <- ff_closest_comparator(res,"pop_3")  
      
      
      # Create table for LY breakdown
      
      # All risk
      
      cabo_nivo_LY_allrisk <- res$weighted_model_undisc$pop_1[L1 == 1] %>%
        select(starts_with("ly"))
      
      comparator_LY_allrisk <- res$weighted_model_undisc$pop_1[L1 == comparator_no_allrisk] %>%
        select(starts_with("ly"))
      
      
      ft_LY_all_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_LY_allrisk,
        comparator_outcomes = comparator_LY_allrisk,
        comparator_no = comparator_no_allrisk,
        LYorQALY = "LY"
      )
      
      # Fav risk
      
      cabo_nivo_LY_favrisk <- res$weighted_model_undisc$pop_2[L1 == 1] %>% 
        select(starts_with("ly"))
      
      comparator_LY_favrisk <- res$weighted_model_undisc$pop_2[L1 == comparator_no_favrisk ] %>%     
        select(starts_with("ly"))
      
      
      ft_LY_fav_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_LY_favrisk, 
        comparator_outcomes = comparator_LY_favrisk, 
        comparator_no = comparator_no_favrisk, 
        LYorQALY = "LY") 
      
      
      # Int/poor risk
      cabo_nivo_LY_IPrisk <- res$weighted_model_undisc$pop_3[L1 == 1] %>% 
        select(starts_with("ly"))
      
      comparator_LY_IPrisk <- res$weighted_model_undisc$pop_3[L1 == comparator_no_IPrisk] %>%     
        select(starts_with("ly"))
      
      ft_LY_IP_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_LY_IPrisk, 
        comparator_outcomes = comparator_LY_IPrisk, 
        comparator_no = comparator_no_IPrisk, 
        LYorQALY = "LY") 
      
      
      
      
      # Create tables for QALY breakdown
      
      # All risk
      
      cabo_nivo_QALY_allrisk <- res$weighted_model_disc$pop_1[L1 == 1] %>% 
        select(starts_with("qaly"))
      
      cabo_nivo_AEQALY_allrisk <- cbind(res$weighted_model_disc$pop_1[L1 == 1] %>% 
                                          select(starts_with("ae_qaly")),BSC = 0)
      
      cabo_nivo_QALY_allrisk <-colSums(rbind(cabo_nivo_QALY_allrisk,  cabo_nivo_AEQALY_allrisk, use.names=FALSE))
      
      cabo_nivo_QALY_allrisk <- cabo_nivo_QALY_allrisk[c(1:2,9,3:8)]
      
      comparator_QALY_allrisk <- res$weighted_model_disc$pop_1[L1 == comparator_no_allrisk] %>%     
        select(starts_with("qaly"))
      
      comparator_AEQALY_allrisk <- cbind(res$weighted_model_disc$pop_1[L1 == comparator_no_allrisk] %>% 
                                           select(starts_with("ae_qaly")),BSC = 0)
      
      comparator_QALY_allrisk <-colSums(rbind(comparator_QALY_allrisk,  comparator_AEQALY_allrisk, use.names=FALSE))
      
      comparator_QALY_allrisk <- comparator_QALY_allrisk[c(1:2,9,3:8)]
      
      
      ft_QALY_all_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_QALY_allrisk, 
        comparator_outcomes = comparator_QALY_allrisk, 
        comparator_no = comparator_no_allrisk, 
        LYorQALY = "QALY") 
      
      
      # Fav risk
      
      cabo_nivo_QALY_favrisk <- res$weighted_model_disc$pop_2[L1 == 1] %>% 
        select(starts_with("qaly"))
      
      cabo_nivo_AEQALY_favrisk <- cbind(res$weighted_model_disc$pop_2[L1 == 1] %>% 
                                          select(starts_with("ae_qaly")),BSC = 0)
      
      cabo_nivo_QALY_favrisk <-colSums(rbind(cabo_nivo_QALY_favrisk,  cabo_nivo_AEQALY_favrisk, use.names=FALSE))
      
      cabo_nivo_QALY_favrisk <- cabo_nivo_QALY_favrisk[c(1:2,9,3:8)]
      
      comparator_QALY_favrisk <- res$weighted_model_disc$pop_2[L1 == comparator_no_favrisk] %>%     
        select(starts_with("qaly"))
      
      comparator_AEQALY_favrisk <- cbind(res$weighted_model_disc$pop_2[L1 == comparator_no_favrisk] %>% 
                                           select(starts_with("ae_qaly")),BSC = 0)
      
      comparator_QALY_favrisk <-colSums(rbind(comparator_QALY_favrisk,  comparator_AEQALY_favrisk, use.names=FALSE))
      
      comparator_QALY_favrisk <- comparator_QALY_favrisk[c(1:2,9,3:8)]
      
      
      ft_QALY_fav_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_QALY_favrisk, 
        comparator_outcomes = comparator_QALY_favrisk, 
        comparator_no = comparator_no_favrisk, 
        LYorQALY = "QALY") 
      
      
      # Int/poor risk
      
      cabo_nivo_QALY_IPrisk <- res$weighted_model_disc$pop_3[L1 == 1] %>% 
        select(starts_with("qaly"))
      
      cabo_nivo_AEQALY_IPrisk <- cbind(res$weighted_model_disc$pop_1[L1 == 1] %>% 
                                         select(starts_with("ae_qaly")),BSC = 0)
      
      cabo_nivo_QALY_IPrisk <-colSums(rbind(cabo_nivo_QALY_IPrisk,  cabo_nivo_AEQALY_IPrisk, use.names=FALSE))
      
      cabo_nivo_QALY_IPrisk <- cabo_nivo_QALY_IPrisk[c(1:2,9,3:8)]
      
      comparator_QALY_IPrisk <- res$weighted_model_disc$pop_3[L1 == comparator_no_IPrisk] %>%     
        select(starts_with("qaly"))
      
      comparator_AEQALY_IPrisk <- cbind(res$weighted_model_disc$pop_3[L1 == comparator_no_IPrisk] %>% 
                                          select(starts_with("ae_qaly")),BSC = 0)
      
      comparator_QALY_IPrisk <-colSums(rbind(comparator_QALY_IPrisk,  comparator_AEQALY_IPrisk, use.names=FALSE))
      
      comparator_QALY_IPrisk <- comparator_QALY_IPrisk[c(1:2,9,3:8)]
      
      
      ft_QALY_IP_tab <- ff_report_outcomes_breakdown(
        cabo_nivo_outcomes = cabo_nivo_QALY_IPrisk, 
        comparator_outcomes = comparator_QALY_IPrisk, 
        comparator_no = comparator_no_IPrisk, 
        LYorQALY = "QALY") 
      
      
      # Create tables for overall cost breakdown
      
      cost_type <- c("drug" , "admin" , "mru" , "eol" , "ae_cost")
      
      populations <- names(res$weighted_model_disc)
      
      summary_costs_table <- rbindlist(lapply(populations, function(popu) {
        treatments <- res$weighted_model_disc[[popu]]$L1
        rbindlist(lapply(treatments, function(mol) {
          ff_cost_table(
            disc_results = res$weighted_model_disc, 
            trt_no = mol, 
            pop = popu
          )
        }))
      }))
      
      summary_costs_table[, risk_pop := str_replace(Population, "Int/poor", "Intermediate / poor risk")]
      summary_costs_table[, Population := NULL]
      
      summary_costs_table <- summary_costs_table[order(risk_pop, Total)] # order by increasing costs
      
      ft_cost_tab <- summary_costs_table %>%
        rename(`Risk population` = risk_pop) %>%
        as_grouped_data(groups = "Risk population") %>% 
        as_flextable() %>% 
        width(., width = (Word_width_inches/(ncol(summary_costs_table)))) %>% 
        add_header_row(top = TRUE, values = c("","1L costs", "Subsequent treatment", "MRU","",""), colwidths = c(1,3,3,2,1,1)) %>%
        theme_box() |> 
        set_header_labels(
          values = list(
            Treatment = "Technologies",
            L1_drug = "Drug cost",
            L1_admin = "Admin cost",
            L1_ae = "AE cost",
            subs_drug = "Drug cost",
            subs_admin = "Admin cost",
            subs_ae = "AE cost",
            mru_1L = "1L",
            subs_mru = "Subsequent treatment",
            eol_cost = "EOL cost",
            Total = "Total cost"
          )
        ) %>%
        colformat_double(j=c(2:11), digits = 0, prefix = "£") %>%
        add_footer_lines("Abbreviations: admin, administration; AE, adverse event; EOL, end of life;  MRU, medical resource use") %>%
        # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
        bold( bold = TRUE, part="header") %>% 
        fontsize(i = NULL, size = 10, part = c("header")) %>%
        fontsize(i = NULL, size = 10, part = c("body")) %>%
        fontsize(i = NULL, size = 9, part = c("footer")) %>%
        align(i = ~ !is.na(`Risk population`), align = "left") %>% 
        align(i= NULL, align = "center", part = c("header")) %>%
        bold(i = ~ !is.na(`Risk population`))  %>%
        autofit() %>% 
        set_table_properties(layout = "autofit") 
      
      
      # produce break downs by population
      
      intervention_name <- p$basic$lookup$ipd$mol[Number == 1]$Description
      
      # all risk
      
      comparator_name <- p$basic$lookup$ipd$mol[Number == comparator_no_allrisk]$Description
      
      cost_breakdown_2 <- rbind(
        summary_costs_table[risk_pop == "All risk" & Treatment == intervention_name ], 
        summary_costs_table[risk_pop == "All risk" & Treatment == comparator_name]
      )
      
      # reshape the data:
      cb2 <- melt.data.table(cost_breakdown_2,id.vars = c("Treatment","risk_pop"))
      cb2$risk_pop <- NULL
      cb2 <- dcast.data.table(cb2, variable ~ Treatment)
      colnames(cb2) <- c("Type", "Int", "Comp")
      cb2$Inc <- cb2[,Int] - cb2[,Comp]
      cb2$abs <- abs(cb2$Inc)
      cb2$abs[10] <- sum(cb2$abs[1:9])
      cb2$abspercent <-  cb2$abs / cb2$abs[10] * 100
      
      cb2[,1] <- c("Drug acquisition cost (1L)", "Admin cost (1L)", "AE cost (1L)", "Drug acquisition cost (2L+)", "Admin cost (2L+)", "AE cost (2L+)", "MRU 1L", "MRU 2L+", "EOL","Total")
      
      cost_table_2_allrisk <- ff_cost_byrisk_table(cb2, comparator_no_allrisk)
      
      # favourable risk
      
      comparator_name <- p$basic$lookup$ipd$mol[Number == comparator_no_favrisk]$Description
      
      cost_breakdown_2 <- rbind(summary_costs_table[risk_pop == "Favourable risk" & Treatment == intervention_name ], summary_costs_table[risk_pop == "Favourable risk" & Treatment == comparator_name])
      cost_breakdown_2 <- as.data.table(x = t(cost_breakdown_2), stringsAsFactors = FALSE)
      cost_breakdown_2 <- cbind(colnames(summary_costs_table), cost_breakdown_2)
      cost_breakdown_2 <- cost_breakdown_2[2:11,]
      cost_breakdown_2[,2:3] <- lapply(cost_breakdown_2[,2:3], as.numeric)    
      colnames(cost_breakdown_2) <- c("Type", "Int", "Comp")
      cost_breakdown_2$Inc <- cost_breakdown_2[,Int] - cost_breakdown_2[,Comp]
      cost_breakdown_2$abs <- abs(cost_breakdown_2$Inc)
      cost_breakdown_2$abs[10] <- sum(cost_breakdown_2$abs[1:9])
      cost_breakdown_2$abspercent <-  cost_breakdown_2$abs / cost_breakdown_2$abs[10] * 100
      
      cost_breakdown_2[,1] <- c("Drug acquisition cost (1L)", "Admin cost (1L)", "AE cost (1L)", "Drug acquisition cost (2L+)", "Admin cost (2L+)", "AE cost (2L+)", "MRU 1L", "MRU 2L+", "EOL","Total")
      
      
      cost_table_2_favrisk <- ff_cost_byrisk_table(cost_breakdown_2, comparator_no_favrisk)
      
      
      # int / poor risk
      
      comparator_name <- p$basic$lookup$ipd$mol[Number == comparator_no_IPrisk]$Description
      
      cost_breakdown_2 <- rbind(summary_costs_table[risk_pop == "Intermediate / poor risk" & Treatment == intervention_name ], summary_costs_table[risk_pop == "Intermediate / poor risk" & Treatment == comparator_name])
      cost_breakdown_2 <- as.data.table(x = t(cost_breakdown_2), stringsAsFactors = FALSE)
      cost_breakdown_2 <- cbind(colnames(summary_costs_table), cost_breakdown_2)
      cost_breakdown_2 <- cost_breakdown_2[2:11,]
      cost_breakdown_2[,2:3] <- lapply(cost_breakdown_2[,2:3], as.numeric)    
      colnames(cost_breakdown_2) <- c("Type", "Int", "Comp")
      cost_breakdown_2$Inc <- cost_breakdown_2[,Int] - cost_breakdown_2[,Comp]
      cost_breakdown_2$abs <- abs(cost_breakdown_2$Inc)
      cost_breakdown_2$abs[10] <- sum(cost_breakdown_2$abs[1:9])
      cost_breakdown_2$abspercent <-  cost_breakdown_2$abs / cost_breakdown_2$abs[10] * 100
      
      cost_breakdown_2[,1] <- c("Drug acquisition cost (1L)", "Admin cost (1L)", "AE cost (1L)", "Drug acquisition cost (2L+)", "Admin cost (2L+)", "AE cost (2L+)", "MRU 1L", "MRU 2L+", "EOL","Total")
      
      
      cost_table_2_IPrisk <- ff_cost_byrisk_table(cost_breakdown_2, comparator_no_IPrisk)
      
      
      #### Scenario analysis tables
      
      
      # all risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_allrisk, "pop_1", p$basic$structure)
      
      Scenario_table_allrisk <- ff_scenario_table(Scenario_table)
      
      
      # favourable risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_favrisk, "pop_2", p$basic$structure)
      
      Scenario_table_favrisk <- ff_scenario_table(Scenario_table)
      
      # int/poor risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_IPrisk, "pop_3", p$basic$structure)
      
      Scenario_table_IPrisk <- ff_scenario_table(Scenario_table)
      
      # base case table
      
      ft_basecase <- ft_wa_inc %>%
        rename(`Risk population` = risk_pop) %>%
        as_grouped_data(groups = "Risk population") %>% 
        as_flextable() %>% 
        width(., width = (Word_width_inches/(ncol(ft_wa_inc)))) %>% 
        theme_box() |> 
        set_header_labels(
          values = list(
            L1 = "Technologies",
            costs = "Costs (£)",
            ly = "LYG",
            qalys = "QALYs",
            ic = "Inc. Costs",
            il = "Inc. LYG",
            iq = "Inc. QALYs",
            Pairwise_ICER = "ICER cabo + nivo vs comparator",
            ICER = "ICER incremental"
          )
        ) %>%
        flextable::colformat_double(j=c(2,5,8,9), digits = 0, prefix = "£") %>%
        flextable::colformat_double(j=c(3,4,6,7), digits = 2) %>%
        add_footer_lines("Abbreviations: ICER, incremental cost-effectiveness ratio; LYG, life-years gained; QALY, quality-adjusted life-year") %>%
        # add_header_row(colwidths = c(1,1, 2),values = c("","g1", "g2")) |> 
        bold( bold = TRUE, part="header") %>% 
        fontsize(i = NULL, size = 10, part = c("header")) %>%
        fontsize(i = NULL, size = 10, part = c("body")) %>%
        fontsize(i = NULL, size = 9, part = c("footer")) %>%
        align(i = ~ !is.na(`Risk population`), align = "left") %>% 
        bold(i = ~ !is.na(`Risk population`)) 
      
      
      # Severity modifier
      
      severity_table <- data.table(do.call(rbind, res$mk$qaly_shortfall_1_to_3)) 
      severity_table <- cbind(risk_pop = p$basic$lookup$pop_map$Risk.population[1:3], severity_table)
      
      severity_table <- rbind(severity_table, f_res_cabonivo_SevMod(
        res = res, 
        oo_pop_string = "Poor / intermediate risk", 
        pop_n = 3, 
        comp_numb = 5))
      
      severity_table <- rbind(severity_table, f_res_cabonivo_SevMod(
        res = res, 
        oo_pop_string = "Poor / intermediate risk", 
        pop_n = 3, 
        comp_numb = 8))
      
      setDT(severity_table)[, risk_pop := str_replace(risk_pop, "Favourable risk", "Fav")]
      setDT(severity_table)[, risk_pop := str_replace(risk_pop, "All risk", "All")]
      
      severity_table$SOC <- unlist(lapply(1:nrow(severity_table), function(mol) {
        p$basic$lookup$ipd$mol[Number == severity_table$SOC[mol]]$Description } ))
      
      
      
      
      
      ft_severity_mod <- ff_severity_table(severity_table)
      
      
      # Scenario analysis pairwise results
      
      ft_all_pairwise <- do.call(rbind,lapply(structure(names(res$pairwise_vs_mol),.Names=names(res$pairwise_vs_mol)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- res$pairwise_vs_mol[[popu_txt]]
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Risk.population
        
        # popu$seq_pop <- seq_popu_lab
        popu$risk_pop <- rsk_popu_lab
        
        return(popu)
        
      }))
      
      ft_all_pairwise$ICER[is.na(ft_all_pairwise$icer) != TRUE] <- as.character(paste0("£",round(ft_all_pairwise$icer[is.na(ft_all_pairwise$icer) != TRUE] ,0)))
      
      ft_all_pairwise[ft_all_pairwise$icer < 0 & ft_all_pairwise$iq < 0]$ICER <- "Cabo+nivo dominated"
      ft_all_pairwise[ft_all_pairwise$icer < 0 & ft_all_pairwise$iq > 0]$ICER <- "Cabo+nivo dominant"
      ft_all_pairwise[ft_all_pairwise$icer > 0 & ft_all_pairwise$iq < 0]$ICER <- paste0("SW quadrant ",ft_all_pairwise[ft_all_pairwise$icer > 0 & ft_all_pairwise$iq < 0]$ICER )
      
      ft_all_pairwise[,icer:=NULL]
      
      setDT(ft_all_pairwise)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      
      ft_all_pairwise_tab <- ff_scenario_pairwise_table(ft_all_pairwise )
      
      
      # Outputting report (state transition) ------------------------------------------------------
      
      # Add base case results.
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Base-case results (ordered in increasing cost)"),
                              bookmark = "tab1") %>%
        body_add_flextable(ft_basecase,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        
        body_add_break() 
      
      
      doc_res <- body_end_section_landscape(doc_res)
      
      doc_res <- doc_res %>% 
        body_add_par("Qualification for the severity modifier",style = "heading 2")  %>% 
        body_add_table_legend(paste0("Application of the severity modifier to the base case"),
                              bookmark = "tab2") %>%
        body_add_flextable(ft_severity_mod,
                           align = "left",
                           topcaption = TRUE,
                           split = TRUE) %>% 
        
        body_add_break() 
      
      
      doc_res <- doc_res %>% 
        body_add_par("Breakdowns by health state and cost category",style = "heading 2")  %>%
        body_add_table_legend(paste0("Summary of LY gain by health state (all risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_allrisk]$Description,")"),
                              bookmark = "tab3") %>%
        body_add_flextable(ft_LY_all_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of LY gain by health state (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_favrisk]$Description,")"),
                              bookmark = "tab4") %>%
        body_add_flextable(ft_LY_fav_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of LY gain by health state (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_IPrisk]$Description,")"),
                              bookmark = "tab5") %>%
        body_add_flextable(ft_LY_IP_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of QALY gain by health state (all risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_allrisk]$Description,")"),
                              bookmark = "tab1") %>%
        body_add_flextable(ft_QALY_all_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of QALY gain by health state (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_favrisk]$Description,")"),
                              bookmark = "tab6") %>%
        body_add_flextable(ft_QALY_fav_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of QALY gain by health state (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_IPrisk]$Description,")"),
                              bookmark = "tab7") %>%
        body_add_flextable(ft_QALY_IP_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      
      doc_res <- body_end_section_portrait(doc_res)
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of costs by health state"),
                              bookmark = "tab8") %>%
        body_add_flextable(ft_cost_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- body_end_section_landscape(doc_res)
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of predicted resource use by category of cost (all risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_allrisk]$Description,")"),
                              bookmark = "tab1") %>%
        body_add_flextable(cost_table_2_allrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of predicted resource use by category of cost (favourable risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_favrisk]$Description,")"),
                              bookmark = "tab9") %>%
        body_add_flextable(cost_table_2_favrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("Summary of predicted resource use by category of cost (intermediate / poor risk, cabo+nivo vs next best non-dominated comparator: " ,p$basic$lookup$ipd$mol[Number == comparator_no_IPrisk]$Description,")"),
                              bookmark = "tab10") %>%
        body_add_flextable(cost_table_2_IPrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break() 
      
      
      population_table <- as.data.table(cbind(risk_pop = p$basic$lookup$pop_map$Risk.population[1:3],pop_labels = c("pop_1","pop_2","pop_3")))
      setDT(population_table)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      
      
      
      for (popu in population_table$pop_labels) {
        for(mol in names(res$weighted_trace_plots[[popu]]$plots)) {
          print(paste0(popu,mol))
          doc_res <- doc_res %>% body_add_figure_legend(
            legend = paste0("Markov trace: ", population_table[pop_labels == popu]$risk_pop , ", ", p$basic$lookup$ipd$mol[Number == str_sub(mol, -1, -1)]$Description),
            bookmark = "fig1"
          ) %>%
            body_add_plot(print(res$weighted_trace_plots[[popu]]$plots[mol]), width = 6) %>%
            body_add_par(
              paste0(
                "Abbreviations: L1, 1st line; L2, 2nd line; L3, 3rd line; L4, 4th line; L5, 5th line"
              ), style = "Table footnote"
            ) %>% body_add_break() 
          
        }
      }
      
      
      doc_res <- doc_res  %>% body_add_break() %>% 
        body_add_par("Cost-effectiveness acceptability frontiers",style = "heading 2")  %>%  
        body_add_par(
          paste0("Cost-effectiveness acceptability frontiers are presented for all non-dominated treatments for each of the risk groups"))  %>%
        body_add_break() %>% 
        body_add_figure_legend(
          legend = paste0("Cost-effectiveness acceptability frontier – all risk"),
          bookmark = "fig2") %>% 
        body_add_plot(print(res$weighted_incremental$pop_1$p), height = 4) %>%
        body_add_par(
          paste0("Abbreviations: QALYs, quality-adjusted life-years"), style = "Table footnote") 
      
      doc_res <- doc_res %>% 
        body_add_figure_legend(
          legend = paste0("Cost-effectiveness acceptability frontier – favourable risk"),
          bookmark = "fig3") %>% 
        body_add_plot(print(res$weighted_incremental$pop_2$p), height = 4) %>%
        body_add_par(
          paste0("Abbreviations: QALYs, quality-adjusted life-years"), style = "Table footnote") 
      
      doc_res <- doc_res %>% 
        body_add_figure_legend(
          legend = paste0("Cost-effectiveness acceptability frontier – intermediate / poor risk"),
          bookmark = "fig4") %>% 
        body_add_plot(print(res$weighted_incremental$pop_3$p), height = 4) %>%
        body_add_par(
          paste0("Abbreviations: QALYs, quality-adjusted life-years"), style = "Table footnote")   %>%
        
        body_add_break() 
      
      doc_res <- body_end_section_portrait(doc_res)
      
      
      doc_res <- doc_res %>%
        body_add_par("Scenario analysis style tables",style = "heading 1")  %>%  
        body_add_table_legend(
          legend = paste0("Scenario analysis - all risk"),
          bookmark = "tab11") %>% 
        body_add_flextable(Scenario_table_allrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>%
        
        body_add_break() 
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis - favourable risk"),
          bookmark = "tab12") %>% 
        body_add_flextable(Scenario_table_favrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>%
        
        body_add_break() 
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis - intermediate / poor risk"),
          bookmark = "tab13") %>% 
        body_add_flextable(Scenario_table_IPrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>%
        
        body_add_break() 
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis pairwise comparison table"),
          bookmark = "tab14") %>% 
        body_add_flextable(ft_all_pairwise_tab ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>%
        
        body_add_break() 
      
    } else {
      
      # Producing report tables (PartSA) ------------------------------------------------------
      
      
      # Make LY table
      
      PartSA_Lys <- do.call(rbind,lapply(structure(names(res$ly),.Names=names(res$ly)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- as.data.table(res$ly[[popu_txt]])
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        popu <-  data.table(L1 = rownames(res$ly[[popu_txt]]), popu)
        popu$L1 <- lu_mol[match(popu$L1,lu_mol$RCC_input_desc),]$Description
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      setDT(PartSA_Lys)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      PartSA_Lys$Total <- rowSums(PartSA_Lys[,2:5])
      
      PartSA_LYs_table <- ff_PartSALY_table(PartSA_Lys)
      
      # Make QALYs table
      
      PartSA_QALYs <- do.call(rbind,lapply(structure(names(res$disc_qaly),.Names=names(res$disc_qaly)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- as.data.table(res$disc_qaly[[popu_txt]])
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        popu <-  data.table(L1 = rownames(res$disc_qaly[[popu_txt]]), popu)
        popu$L1 <- lu_mol[match(popu$L1,lu_mol$RCC_input_desc),]$Description
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      setDT(PartSA_QALYs)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      PartSA_QALYs$Total <- rowSums(PartSA_QALYs[,2:5])
      
      PartSA_QALYs_table <- ff_PartSAQALY_table(PartSA_QALYs)
      
      # Make costs table
      
      PartSA_costs <- do.call(rbind,lapply(structure(names(res$disc_cost),.Names=names(res$disc_cost)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- as.data.table(res$disc_cost[[popu_txt]])
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        popu <-  data.table(L1 = rownames(res$disc_cost[[popu_txt]]), popu)
        popu$L1 <- lu_mol[match(popu$L1,lu_mol$RCC_input_desc),]$Description
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      setDT(PartSA_costs)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      PartSA_costs$Total <- rowSums(PartSA_costs[,2:11])
      
      PartSA_costs_table <- ff_PartSAcost_table (PartSA_costs)
      
      # Make results table
      
      PartSA_wa <- do.call(rbind,lapply(structure(names(res$incremental),.Names=names(res$incremental)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        if (is.null(res$incremental[[popu_txt]]$non_dominated)) {
          popu <- as.data.table(res$incremental[[popu_txt]])
          popu <-  data.table(popu, ic = 0, iq = 0, il = 0, ICER = "Dominant" )
          popu$str_dom <- NULL
          
        } else {
          popu <- as.data.table(res$incremental[[popu_txt]]$expanded_results)
          popu$ICER[popu$extdom==FALSE] <- as.character(paste0("£",round(popu$ICER[popu$extdom==FALSE] ,0))) 
          popu$ICER[popu$extdom==TRUE] <- "(ext dominated)" 
          popu$str_dom <- NULL
          popu$extdom <- NULL
          popu$r <- NULL
          
        }
        
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      PartSA_wa <- PartSA_wa[,c(2,3,1,4,5,7,6,8,9)]
      
      
      
      PartSA_totals <- do.call(rbind,lapply(structure(names(res$tables$top_line),.Names=names(res$tables$top_line)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        lu_mol <- p$basic$lookup$ipd$mol
        
        popu <- as.data.table(res$tables$top_line[[popu_txt]])
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        popu$L1 <- lu_mol[match(popu$L1,lu_mol$Number),]$Description
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      PartSA_totals <- PartSA_totals[order( risk_pop, costs)] # order by increasing costs
      
      PartSA_totals$L1_risk  <- paste(PartSA_totals$L1, PartSA_totals$risk_pop)
      
      PartSA_wa$L1_risk  <- paste(PartSA_wa$L1, PartSA_wa$risk_pop)
      
      PartSA_results <- merge(PartSA_totals, PartSA_wa, all.x = TRUE)
      PartSA_results[is.na(ICER)]$ICER <- "(dominated)"
      
      PartSA_results <-  PartSA_results[,c(4,1,3,2,7,8,9,10,5)]
      
      PartSA_results <- PartSA_results[order( risk_pop, costs)] # order by increasing costs
      
      
      PartSA_Pairwise <- do.call(rbind,lapply(structure(names(res$tables$top_line),.Names=names(res$tables$top_line)), function(popu_txt) {
        
        lu_pop <- p$basic$lookup$pop_map
        
        popu <-  f_res_ICER_pairwiseVsoneTrt(res$tables$top_line[[popu_txt]], 1, p$basic$lookup$ipd$mol)
        popu_n <- as.numeric(gsub("pop_","",popu_txt))
        
        # seq_popu_lab <- lu_pop[Overall.population.number == popu_n,]$Sequencing.population
        rsk_popu_lab <- rep(lu_pop[Overall.population.number == popu_n,]$Risk.population,nrow(popu))
        
        # popu$seq_pop <- seq_popu_lab
        popu <- cbind(popu, risk_pop=rsk_popu_lab)
        
        return(popu)
        
      }))
      
      
      PartSA_Pairwise$Pairwise_ICER[is.na(PartSA_Pairwise$icer) != TRUE] <- as.character(paste0("£",round(PartSA_Pairwise$icer[is.na(PartSA_Pairwise$icer) != TRUE] ,0)))
      
      
      PartSA_Pairwise[PartSA_Pairwise$icer < 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER <- "Cabo+nivo dominated"
      PartSA_Pairwise[PartSA_Pairwise$icer < 0 & PartSA_Pairwise$iq > 0]$Pairwise_ICER <- "Cabo+nivo dominant"
      PartSA_Pairwise[PartSA_Pairwise$icer > 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER <- paste0("SW quadrant ",PartSA_Pairwise[PartSA_Pairwise$icer > 0 & PartSA_Pairwise$iq < 0]$Pairwise_ICER)
      
      PartSA_Pairwise_Scen <- PartSA_Pairwise
      
      PartSA_Pairwise <-  PartSA_Pairwise[,c(4,9,10)]
      
      
      PartSA_results <- merge(PartSA_results, PartSA_Pairwise, all.x = TRUE)
      
      PartSA_results <- PartSA_results[,c(1:8,10,9)]
      
      PartSA_results[ICER == 0]$ICER <- "(dominated)"
      PartSA_results[,6] <- as.numeric(unlist(PartSA_results[,6][[1]]))
      PartSA_results[,7] <- as.numeric(unlist(PartSA_results[,7][[1]]))
      PartSA_results[,8] <- as.numeric(unlist(PartSA_results[,8][[1]]))
      
      PartSA_results<- PartSA_results[order( risk_pop, costs)] # order by increasing costs
      setDT(PartSA_results)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      
      PartSA_results_tab <-ff_PartSAresults_table(PartSA_results)
      
      #### Scenario analysis tables
      
      comparator_no_allrisk <-  ff_closest_comparator_PartSA(res,"pop_1")
      comparator_no_favrisk <- ff_closest_comparator_PartSA(res,"pop_2")
      comparator_no_IPrisk <- ff_closest_comparator_PartSA(res,"pop_3")  
      
      
      # all risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_allrisk, "pop_1", p$basic$structure)
      
      Scenario_table_allrisk <- ff_scenario_table(Scenario_table)
      
      # favourable risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_favrisk, "pop_2", p$basic$structure)
      
      Scenario_table_favrisk <- ff_scenario_table(Scenario_table)
      
      # int/poor risk
      
      Scenario_table <- ff_scenario_output(res, Scenario_name, comparator_no_IPrisk, "pop_3", p$basic$structure)
      
      Scenario_table_IPrisk <- ff_scenario_table(Scenario_table)
      
      # Scenario analysis pairwise results
      
      setDT(PartSA_Pairwise_Scen)[, risk_pop := str_replace(risk_pop, "Int/poor", "Intermediate / poor risk")]
      
      PartSA_Pairwise_Scen <- PartSA_Pairwise_Scen[, c(4,1,2,3,5,6,7,10,9)]
      PartSA_Pairwise_Scen$ICER <- PartSA_Pairwise_Scen$Pairwise_ICER
      PartSA_Pairwise_Scen$Pairwise_ICER <- NULL
      
      ft_all_pairwise_tab <- ff_scenario_pairwise_table(PartSA_Pairwise_Scen )
      
      
      # Outputting report (PartSA) ------------------------------------------------------
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("PartSA life years"),
                              bookmark = "tab1") %>%
        body_add_flextable(PartSA_LYs_table,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break()
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("PartSA QALYs"),
                              bookmark = "tab2") %>%
        body_add_flextable(PartSA_QALYs_table,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break()
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("PartSA costs"),
                              bookmark = "tab3") %>%
        body_add_flextable(PartSA_costs_table,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break()
      
      doc_res <- doc_res %>% 
        body_add_table_legend(paste0("PartSA results (ordered in increasing cost)"),
                              bookmark = "tab4") %>%
        body_add_flextable(PartSA_results_tab,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) %>% 
        body_add_break()
      
      doc_res <- doc_res %>%
        body_add_par("Scenario analysis style table",style = "heading 1")  %>%  
        body_add_table_legend(
          legend = paste0("Scenario analysis style table - all risk"),
          bookmark = "tab5") %>% 
        body_add_flextable(Scenario_table_allrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>% 
        body_add_break()
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis style table - favourable risk"),
          bookmark = "tab6") %>% 
        body_add_flextable(Scenario_table_favrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>% 
        body_add_break()
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis style table - intermediate / poor risk"),
          bookmark = "tab7") %>% 
        body_add_flextable(Scenario_table_IPrisk ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE)  %>% 
        body_add_break()
      
      doc_res <- doc_res %>%
        body_add_table_legend(
          legend = paste0("Scenario analysis pairwise comparison table"),
          bookmark = "tab8") %>% 
        body_add_flextable(ft_all_pairwise_tab ,
                           align = "left", 
                           topcaption = TRUE,
                           split = TRUE) 
      
      
    }
    
    
    print(doc_res, target = paste0("./4_Output/Scenario ",Scenario_number,"_",i$dd_drug_price_options,gsub(":","_",Run_date),".docx"))
    rm(doc_res)
    
    
  }


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


