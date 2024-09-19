---
title: "Code walkthrough"
execute:
  keep-md: true
---



This is a single `.qmd` file that walks through the code in `Model_Structure.R`. It is rendered as a single file, then the output `.md` file is split into seperate `.md` files that can be incorporated into seperate `.qmd` pages.

This is done as running the code within seperate `.qmd` pages requires that the global environment is shared between them, which is difficult to achieve without long run times or very large files containing the environment.

SPLITMD_CODE1_START

This page contains some of the basic set-up steps like **loading the functions** and lots of the **model inputs**, as well as establishing the **possible treatment sequences**.

## Load required packages


::: {.cell}

```{.r .cell-code}
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
library(knitr, quiet = TRUE)
```
:::


## Set file paths


::: {.cell}

```{.r .cell-code}
# Set path to folders
d_path = "../../../1_Data"
f_path = "../../../3_Functions"
o_path = "../../../4_Output"
```
:::


## Define some of the model settings

The majority of this code block is dedicated to setting the model to **run sequentially or in parallel**. If running the base case using `Model_Structure.R` as provided (e.g. with pre-run survival analysis), you may find that the quickest option is to run the model sequentially. This is already set up by default, with `keep_free_cores <- NA`. For reference, the run time for this on an Intel Core i7-12700H with 32GB RAM running Ubuntu 22.04.4 Linux was 40 minutes.


::: {.cell}

```{.r .cell-code}
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
keep_free_cores <- NA
if (any(is.na(keep_free_cores), keep_free_cores<0)) {
  plan(sequential)
} else {
  plan(multisession(workers = max(availableCores()-keep_free_cores,1)))
}
```
:::


The three settings in this code block:

* `progressr::handlers("progress")` - one of the settings for how progress is reported whilst code is running
* `crosstable::options(crosstable_units="cm")` - the `crosstable` package generates descriptive statistics with function `crosstable()`, and this sets the unit for it, although it should be noted that this appears to be legacy code as it doesn't appear that `crosstable()` is used anywhere in the repository
* `qc_mode` - the input for `verbose` in `f_NMA_AddAssumptionsToNetwork()` which, if true, will mean that extra outputs are printed to the console


::: {.cell}

```{.r .cell-code}
# Other generic settings for the progress bar and units for table widths
handlers("progress")
options(crosstable_units="cm")

#### 2. Loading functions ###########


# This variable is used throughout the model to define whether to provide additional outputs useful for QC or not
# The model will take longer to run when this is set to TRUE
qc_mode <- FALSE
```
:::


## Import model functions

The functions for the analysis in `Model_Structure.R` are stored in the `3_Functions/` folder. Here, they are imported into our environment.


::: {.cell}

```{.r .cell-code}
# 2.1. Excel data extraction functions -----------------------------------------

#### These functions are used to extract parameters from the Excel input workbook for use in R
#### During Phase 2 a Shiny front-end will be added to the model which will allow an alternative mechanism to upload these types of inputs

source(file.path(f_path, "excel/extract.R"))

# 2.2. Treatment sequencing functions ----------------------------------------

#### Function: filter to active treatments and lines
##### Takes as an input the defined sequences, evaluation type and line to start the evaluation from 
##### Other input is % receiving each subs therapy at each line dependent on previous treatments received 
##### Reweights so that the % receiving each treatment sums to 100% within each arm / line being studied
##### Outputs a matrix that has the % receiving each possible combination

source(file.path(f_path, "sequencing/sequences.R"))

# 2.3. Survival analysis functions ---------------------------------------------

# Function: conduct survival analysis
##### by treatment, line, population and outcome fitted survival curves using Flexsurvreg (exp, Weibull, lognormal, loglog, Gompertz, gen gamma)
##### calculation of and adjustment for general population
##### adjustment for treatment effect waning

source(file.path(f_path, "survival/Survival_functions.R"))
source(file.path(f_path, "survival/other_cause_mortality.R"))
source(file.path(f_path, "survival/treatment_effect_waning.R"))

# 2.4 Misc functions ----------------------------------------------------------

### these functions enable smoother data cleaning and manipulation

source(file.path(f_path, "misc/other.R"))
source(file.path(f_path, "misc/shift_and_pad.R"))
source(file.path(f_path, "misc/cleaning.R"))

# 2.4.1 Functions imposing list structures -----------------------------------

source(file.path(f_path, "misc/nesting.R"))
source(file.path(f_path, "misc/discounting.R"))
source(file.path(f_path, "misc/qdirichlet.R"))
source(file.path(f_path, "misc/plotting.R"))
source(file.path(f_path, "misc/structure.R"))

# 2.4.2 Functions calculating HRs from FPNMA coefficients and other FPNMA manipulation ------

source(file.path(f_path, "misc/fpnma_fns.R"))


# 2.5 Utility functions -------------------------------------------------------

source(file.path(f_path, "utility/age_related.R"))
source(file.path(f_path, "costs_and_QALYs/utility_processing.R"))

# 2.6 AE functions --------------------------------------------------------

source(file.path(f_path, "adverse_events/AE_steps.R"))

# 2.7 Cost calculation functions --------------------------------------------

source(file.path(f_path, "costs_and_QALYs/cost_processing.R"))


# 2.8 State transition modelling functions --------------------------------

source(file.path(f_path, "markov/markov.R"))

# 2.9 Patient flow functions ----------------------------------------------

source(file.path(f_path, "patient_flow/overarching.R"))
source(file.path(f_path, "patient_flow/partitioned_survival.R"))
source(file.path(f_path, "patient_flow/markov.R"))
source(file.path(f_path, "patient_flow/drug_costs.R"))
source(file.path(f_path, "patient_flow/hcru_costs.R"))
source(file.path(f_path, "patient_flow/qalys.R"))
source(file.path(f_path, "patient_flow/ae.R"))



# 2.10 Results processing functions ---------------------------------------

source(file.path(f_path, "results/incremental_analysis.R"))
source(file.path(f_path, "results/model_averaging.R"))
source(file.path(f_path, "results/partitioned_survival.R"))
source(file.path(f_path, "misc/severity_modifier.R"))
source(file.path(f_path, "results/results_tables.R"))
source(file.path(f_path, "psa/psa functions.R"))



# 2.11 Office software outputs --------------------------------------------

source(file.path(f_path, "reporting/word_document_output.R"))
```
:::


## Get some of the model inputs

### Introductory comments

This section of the code is mostly comments that describe:

* The structure of `i` which contains model inputs
* That survival analysis will be conducted using state transition (markov) models and partitioned survival analysis (partSA)
* The five input files, which are detailed within the documentation on the page [Input data](../input_data.qmd)

There is a line of code to define `User_types`, but this appears to be legacy as it is not used anywhere.


::: {.cell}

```{.r .cell-code}
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
```
:::


### Get inputs from excel workbook and save as list `i`

Import the file at `excel_path` using `f_excel_extract()`, and then tidy `i$R_table_param` using `f_excel_cleanParams()`.

If the file doesn't exist, assuming the user is in RStudio, a dialog box will appear with the system files, and the user should then select a file from their directory. The dialog will open in `1_Data/`, show only `.xlsm` files, and the accept/ok button has the text `ID6184_RCC_model inputs....xlsm`.


::: {.cell}

```{.r .cell-code}
# The first part of this code pulls all of the named ranges from the excel workbook, expand the parameters table

#Option to define Excel path on local machine - comment in this and comment out the code below to select file
excel_path <- file.path(d_path, "ID6184_RCC_model inputs FAD version [UK RWE unredacted, ACIC redacted, cPAS redacted].xlsm")
#i <- f_excel_extract(excel_path, verbose = TRUE)

if (file.exists(excel_path)) {
  i <- f_excel_extract(excel_path, verbose = FALSE)
} else {
  i <- f_excel_extract(rstudioapi::selectFile(
    caption = "Select the Excel inputs file (ID6184_RCC_model inputs....xlsm)",
    label = "ID6184_RCC_model inputs....xlsm",
    path = "./1_Data/",
    filter = "Excel Files (*.xlsm)",
    existing = TRUE
  ), verbose = FALSE)
}

i <- c(i,f_excel_cleanParams(i$R_table_param))
```
:::


::: {.callout-note collapse="true"}

## What does `f_excel_extract()` do?

This function uses `openxlsx::getNamedRegions()` to find named regions in the workbook. These are all stored in a sheet called `named ranges`. As illustrated in [Input data](../input_data.qmd), these consist of two columns:

* A name for the region (e.g. `List_pop1_allowed`)
* The region, which consists of the: sheet, column/s and row/s (e.g. `=Lists!$BA$11:$BA$22`)

![](../../images/excel_named_range.png)

The function `f_excel_extract()` extracts each of these, saving them into a nested list called `i`. Each element in the list is a row from the `named ranges` sheet - for example:


::: {.cell}

```{.r .cell-code}
i$List_pop1_allowed
```

::: {.cell-output .cell-output-stdout}

```
 [1] "avelumab_plus_axitinib"      "axitinib"                   
 [3] "cabozantinib"                "everolimus"                 
 [5] "lenvatinib_plus_everolimus"  "cabozantinib_plus_nivolumab"
 [7] "nivolumab_monotherapy"       "pazopanib"                  
 [9] "sunitinib"                   "tivozanib"                  
```


:::
:::


The exception is the first element which is a copy of the named ranges sheet:


::: {.cell}

```{.r .cell-code}
kable(head(i[[1]]))
```

::: {.cell-output-display}


|Name                   |Cell.Range              |
|:----------------------|:-----------------------|
|apply_waning_to        |=Lists!$Y$10:$Y$11      |
|bc_settings_rng        |=Lists!$B$99:$B$174     |
|cabo_nivo_outcome_from |=Lists!$W$10            |
|cabo_nivo_outcomes     |=Lists!$X$10:$X$12      |
|count_bc_settings      |=Lists!$B$97            |
|dd_2ndline_NMA         |='Model settings'!$G$40 |


:::
:::


:::

::: {.callout-note collapse="true"}

## What does `f_excel_cleanParams()` do?

This function is applied to `i$R_table_param` which is the full table from the sheet `All parameters`.


::: {.cell}

```{.r .cell-code}
kable(head(i$R_table_param))
```

::: {.cell-output-display}


|Parameter.description                   |Parameter.name    |Mean.current.value |
|:---------------------------------------|:-----------------|:------------------|
|Include cabo nivo? (1=yes, 0=no)        |cabo_nivo_include |1                  |
|Include pem len? (1=yes, 0=no)          |pem_len_include   |1                  |
|Include panzopanib? (1=yes, 0=no)       |pazopanib_include |1                  |
|Include tivozanib? (1=yes, 0=no)        |tivozanib_include |1                  |
|Include sunitinib? (1=yes, 0=no)        |sunitinib_include |1                  |
|Include cabo monotherapy? (1=yes, 0=no) |cabo_include      |1                  |


:::
:::


The function `f_excel_cleanParams()`:

* Converts each row of the table into a list, with value from `Parameter.name` used as the list name
* Adds `Mean` which is simply `Mean.current.value` converted to numeric

These lists were then concatenated with `i`, so we can access them from `i` as follows:


::: {.cell}

```{.r .cell-code}
i$cabo_nivo_include
```

::: {.cell-output .cell-output-stdout}

```
$Parameter.description
[1] "Include cabo nivo? (1=yes, 0=no)"

$Parameter.name
[1] "cabo_nivo_include"

$Mean.current.value
[1] "1"

$Mean
[1] 1
```


:::
:::


:::

### Manually add some extra inputs to `i`


::: {.cell}

```{.r .cell-code}
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
```
:::


### Use `i` to make another list `p`

The list `p` is based on `i`, either copying over parameters or using them to calculate new parameters. The list is first created by `f_misc_param_generate_p()`. A few extra additions are then made to `p` in this code block, such as to set a maximum of 4 treatment lines before best supportive care.


::: {.cell}

```{.r .cell-code}
# The next step is to then "tidy up" i into another object, p. p doesn't necessarily
# have to house everything, only things that will change in PSA

p <- f_misc_param_generate_p(i)

# Set seed for PSA - note this is done in the script to run the PSA, not here!
# set.seed(1475)

# Max lines within the R model
p$basic$R_maxlines <- 4

# Pass this into p so that p can be used to exclusively compute the model:
p$basic$decision_problem <- i$decision_problem
```
:::


::: {.callout-note collapse="true"}

## What does `f_misc_param_generate_p()` do?

This function consists of relatively simple calculations using parameters from `i` to generate `p`. For the full overview of these calculations, check out the code for the function. To give an example though, the function code includes:

```
p <- list(
    basic = list(
      th   = ceiling(i$ui_time_horizon * 365.25 / 7), 
      th_y = i$ui_time_horizon, 
      ...
    )
    ...
)
```

In this example, we can see that is makes a copy of `i$ui_time_horizon` which has time horizon of model in years (40 years) (`p$basic$th_y`). It also converts the time horizon into weeks (`p$basic$th`).


::: {.cell}

```{.r .cell-code}
i$ui_time_horizon
```

::: {.cell-output .cell-output-stdout}

```
[1] 40
```


:::

```{.r .cell-code}
p$basic$th
```

::: {.cell-output .cell-output-stdout}

```
[1] 2088
```


:::

```{.r .cell-code}
p$basic$th_y
```

::: {.cell-output .cell-output-stdout}

```
[1] 40
```


:::
:::


:::

## Treatment sequences

### Find all possible sequences

The first step in determining the possible treatment sequences is to determine all possible combinations and orders of treatment, saving this as `i$sequences`.


::: {.cell}

```{.r .cell-code}
#### 3.2 Define sequences  ###########

#### This code produces a list of possible sequences per population based upon the rules defined for RCC
#### and the user input number of lines


# Add drug names to comparators vector extracted from inputs list.

i$sequences <- f_generate_sequences(
  comparators = i$List_comparators, 
  maxlines    = p$basic$R_maxlines
)
```
:::


::: {.callout-note collapse="true"}

## What does `f_generate_sequences()` do?

The input to this function is a list of all the possible treatments (14 options), and the maximum number of treatment lines (4).


::: {.cell}

```{.r .cell-code}
i$List_comparators
```

::: {.cell-output .cell-output-stdout}

```
 [1] "nivolumab_monotherapy"         "cabozantinib_plus_nivolumab"  
 [3] "nivolumab_plus_ipilimumab"     "lenvatinib_plus_pembrolizumab"
 [5] "avelumab_plus_axitinib"        "pazopanib"                    
 [7] "tivozanib"                     "sunitinib"                    
 [9] "cabozantinib"                  "lenvatinib_plus_everolimus"   
[11] "everolimus"                    "axitinib"                     
[13] "sorafenib"                     "placebo_BSC"                  
```


:::
:::


The function `f_generate_sequences()` outputs a table saved as `i$sequences`. This contains every possible order and combination of treatments.


::: {.cell}

```{.r .cell-code}
kable(head(i$sequences, 3))
```

::: {.cell-output-display}


|V1                     |V2       |V3           |V4                          |V5  |
|:----------------------|:--------|:------------|:---------------------------|:---|
|avelumab_plus_axitinib |axitinib |cabozantinib |cabozantinib_plus_nivolumab |BSC |
|avelumab_plus_axitinib |axitinib |cabozantinib |everolimus                  |BSC |
|avelumab_plus_axitinib |axitinib |cabozantinib |lenvatinib_plus_everolimus  |BSC |


:::

```{.r .cell-code}
dim(i$sequences)
```

::: {.cell-output .cell-output-stdout}

```
[1] 26404     5
```


:::
:::


These varied from a single treatment to up to four subsequent treatments, but always ended with best supportive case (BSC).


::: {.cell}

```{.r .cell-code}
kable(head(i$sequences[i$sequences$V3=="",], 3))
```

::: {.cell-output-display}


|V1                     |V2  |V3 |V4 |V5 |
|:----------------------|:---|:--|:--|:--|
|avelumab_plus_axitinib |BSC |   |   |   |
|axitinib               |BSC |   |   |   |
|cabozantinib           |BSC |   |   |   |


:::
:::


:::

### Filter to valid sequences

The table of all possible treatment sequences is then filtered down to valid sequences (for example, removing drugs if not allowed for a given population, or after another particular drug).

<!-- This section is set to eval FALSE as the progress statements are from cat() so cannot suppress output. Instead, actually run in dropdown below, where there, hide the code but can still see the output -->


::: {.cell}

```{.r .cell-code}
# restrict the pathways to those that are possible and permitted.
i$sequences <- as.data.frame(i$sequences)

populations <- i$i_nr_populations

seqs <- NULL
invisible(
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
)
rownames(seqs) <- NULL

i$sequences <- seqs

#### Uncomment this code to view the sequences and write the sequences defined to csv

# i$sequences
# write.csv(seqs, "4_Output/sequences.csv", row.names = F)
rm(s, seqs, populations)

# define number of cycles and a vector of the cycles 
```
:::


::: {.callout-note collapse="true"}

## About the populations being looped through

This code chunk restricted to valid sequences by population. In this analysis, there are four populations defined by time since an immuno-oncology (IO) treatment, and International Metastatic Renal Cell Carcinoma Database Consortium (IMDC) risk status. They are:

* pop1 >12m since IO, favourable risk
* pop2 >12m since IO, intermediate/poor risk
* pop3 <12m since IO, favourable risk
* pop4 <12m since IO, intermediate/poor risk

Each risk group has to be broken down by time since IO as there are five treatments that cannot be used within 12 months of an adjuvant IO treatment (4.3.5.6 in Assessment Report @lee_treatments_2023-1). An adjuvant treatment is one given alongside the primary treatment.

:::

::: {.callout-note collapse="true"}

## What criteria are there for valid treatments?

### By population and line

For **each population**, there are a list of valid treatments at **each line** of therapy (first-line through to fourth). For example, valid first-line treatments for population 1 are:


::: {.cell}

```{.r .cell-code}
# Pop1 1L treatments
f_get_L1_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
[1] "avelumab_plus_axitinib"      "cabozantinib_plus_nivolumab"
[3] "pazopanib"                   "sunitinib"                  
[5] "tivozanib"                  
```


:::
:::


### Only after

There are some treatments that can only come after other treatments. For example, for population 1:

* Axitinib can only be administered **after** a tyrosine kinase inhibitor (TKI) or cytokine treatment
* Everolimus can only be administered **after** a vascular endothelial growth factor (VEGF) treatment
* At 2L 3L or 4L, cabozantinib can only be administered **after** one of the listed treatments


::: {.cell}

```{.r .cell-code}
f_get_only_after_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$axitinib
[1] "avelumab_plus_axitinib"        "cabozantinib"                 
[3] "lenvatinib_plus_everolimus"    "cabozantinib_plus_nivolumab"  
[5] "pazopanib"                     "lenvatinib_plus_pembrolizumab"
[7] "sunitinib"                     "tivozanib"                    

$everolimus
[1] "avelumab_plus_axitinib"        "axitinib"                     
[3] "cabozantinib"                  "lenvatinib_plus_everolimus"   
[5] "cabozantinib_plus_nivolumab"   "pazopanib"                    
[7] "lenvatinib_plus_pembrolizumab" "sunitinib"                    
[9] "tivozanib"                    
```


:::

```{.r .cell-code}
f_get_2L_only_after_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$cabozantinib
[1] "avelumab_plus_axitinib"        "axitinib"                     
[3] "lenvatinib_plus_everolimus"    "cabozantinib_plus_nivolumab"  
[5] "pazopanib"                     "lenvatinib_plus_pembrolizumab"
[7] "sunitinib"                     "tivozanib"                    
```


:::
:::


### Not immediately after

Some treatments cannot come immediately after another treatment. For example, for population 1:

* Lenvatinib plus everolimus must **not come immediately after** nivolumab plus ipilimumab
* At 2L 3L or 4L, pazopanib and sunitinib and tivozanib must **not come immediately after** their respective listed treatments


::: {.cell}

```{.r .cell-code}
f_get_not_immediate_after_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$lenvatinib_plus_everolimus
[1] "nivolumab_plus_ipilimumab"
```


:::

```{.r .cell-code}
f_get_2L_only_immediate_after_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$pazopanib
[1] "avelumab_plus_axitinib"        "lenvatinib_plus_pembrolizumab"
[3] "cabozantinib_plus_nivolumab"   "nivolumab_plus_ipilimumab"    

$sunitinib
[1] "avelumab_plus_axitinib"        "lenvatinib_plus_pembrolizumab"
[3] "cabozantinib_plus_nivolumab"   "nivolumab_plus_ipilimumab"    

$tivozanib
[1] "avelumab_plus_axitinib"        "lenvatinib_plus_pembrolizumab"
[3] "cabozantinib_plus_nivolumab"   "nivolumab_plus_ipilimumab"    
```


:::
:::


### Only one from list

Some treatments are not allowed if another has already been given at any point prior.

For example, for population 1 there are five lists of treatments where **only one treatment from each list** is allowed.

There are no restrictions like this for population 1 specific to just 2L 3L or 4L treatments (hence, empty list).


::: {.cell}

```{.r .cell-code}
f_get_one_in_list_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$axitinib
[1] "avelumab_plus_axitinib" "axitinib"              

$cabozantinib
[1] "cabozantinib"                "cabozantinib_plus_nivolumab"

$everolimus
[1] "lenvatinib_plus_everolimus" "everolimus"                

$io
[1] "avelumab_plus_axitinib"        "nivolumab_plus_ipilimumab"    
[3] "cabozantinib_plus_nivolumab"   "nivolumab_monotherapy"        
[5] "lenvatinib_plus_pembrolizumab"

$nivolumab
[1] "nivolumab_plus_ipilimumab"   "cabozantinib_plus_nivolumab"
[3] "nivolumab_monotherapy"      

$TKIs
[1] "sunitinib" "pazopanib" "tivozanib"
```


:::

```{.r .cell-code}
f_get_2L_only_one_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
named list()
```


:::
:::


### Only one allowed before

In other cases, a treatment is not allowed if more than one of a particular category of treatment has been given. For example, for population 1:

* Of the listed treatments, **only one is allowed before** lenvatinib plus everolimus


::: {.cell}

```{.r .cell-code}
f_get_only_after_one_lists(i, 1)
```

::: {.cell-output .cell-output-stdout}

```
$lenvatinib_plus_everolimus
[1] "avelumab_plus_axitinib"        "axitinib"                     
[3] "cabozantinib"                  "cabozantinib_plus_nivolumab"  
[5] "pazopanib"                     "lenvatinib_plus_pembrolizumab"
[7] "sunitinib"                     "tivozanib"                    
```


:::
:::


:::

::: {.callout-note collapse="true"}

## What does `f_path_tx_restrict()` do?

The function `f_path_tx_restrict()` is defined in `sequences.R` (there is also a function of the same name in `rccFunctions.R` but this is not sourced).

Its purpose is to restrict the table of all possible sequences to just the valid sequences for each population, restricting it from **26404** rows with possible treatment sequences to just **744**.


::: {.cell}

```{.r .cell-code}
dim(i$sequences)
```

::: {.cell-output .cell-output-stdout}

```
[1] 26404     5
```


:::
:::


The inputs to this function are lists defining valid treatments by different criteria (e.g. line of therapy, subsequent treatments), as detailed in the note above. 

Within `f_path_tx_restrict()`, there are then several other functions which take these lists and use them to remove invalid sequences from the table.

For example, the allowed treatments identified using `f_get_allowed_lists()` are input to `f_path_tx_restrict()` as `allowed`. The function `f_path_allowed()` then uses that list to remove invalid drugs for a given population:

```
s <- f_path_allowed(s, allowed[[1]])
```

Looking at an excerpt of the code for `f_path_allowed()`, we can see that is adds "BSC" and no treatment ("") as valid options, and then only keeps rows if their treatments are in (`%in%`) the list of valid treatments.

```
rule <- c(rule, "BSC", "")

for (n in 1:ncol(perms)) {
  perms <- perms[perms[,n] %in% rule,]
}
```

:::

::: {.callout-note collapse="true"}

## View the sequence restrictions applied


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
Applying sequence restrictions to population 1 
Dropping drugs not allowed for this population.
applying rule: avelumab_plus_axitinib axitinib cabozantinib everolimus lenvatinib_plus_everolimus cabozantinib_plus_nivolumab nivolumab_monotherapy pazopanib sunitinib tivozanib are only allowed treatments.
Permutations before applying rule: 26404 
Permutations after applying rule : 5860 
applying rule: drug line restrictions.
Permutations before applying rule: 5860 
Permutations after applying rule : 628 
[1] "axitinib"
applying rule. axitinib is only allowed after avelumab_plus_axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 628 
Permutations after applying rule : 628 
[1] "everolimus"
applying rule. everolimus is only allowed after avelumab_plus_axitinib axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 628 
Permutations after applying rule : 628 
[1] "lenvatinib_plus_everolimus"
applying rule. lenvatinib_plus_everolimus is not allowed immediately after nivolumab_plus_ipilimumab 
Permutations before applying rule: 628 
Permutations after applying rule : 628 
[1] "axitinib"
applying rule axitinib : avelumab_plus_axitinib axitinib cannot be in one permutation
Permutations before applying rule: 628 
Permutations after applying rule : 559 
[1] "cabozantinib"
applying rule cabozantinib : cabozantinib cabozantinib_plus_nivolumab cannot be in one permutation
Permutations before applying rule: 559 
Permutations after applying rule : 520 
[1] "everolimus"
applying rule everolimus : lenvatinib_plus_everolimus everolimus cannot be in one permutation
Permutations before applying rule: 520 
Permutations after applying rule : 452 
[1] "io"
applying rule io : avelumab_plus_axitinib nivolumab_plus_ipilimumab cabozantinib_plus_nivolumab nivolumab_monotherapy lenvatinib_plus_pembrolizumab cannot be in one permutation
Permutations before applying rule: 452 
Permutations after applying rule : 400 
[1] "nivolumab"
applying rule nivolumab : nivolumab_plus_ipilimumab cabozantinib_plus_nivolumab nivolumab_monotherapy cannot be in one permutation
Permutations before applying rule: 400 
Permutations after applying rule : 400 
[1] "TKIs"
applying rule TKIs : sunitinib pazopanib tivozanib cannot be in one permutation
Permutations before applying rule: 400 
Permutations after applying rule : 202 
[1] "lenvatinib_plus_everolimus"
applying rule: lenvatinib_plus_everolimus can only be after ONE of avelumab_plus_axitinib axitinib cabozantinib cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 202 
Permutations after applying rule : 182 
[1] "cabozantinib"
applying rule: cabozantinib as 2L+ only allowed after avelumab_plus_axitinib axitinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 182 
Permutations after applying rule : 182 
[1] "pazopanib"
applying rule: pazopanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 182 
Permutations after applying rule : 172 
[1] "sunitinib"
applying rule: sunitinib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 172 
Permutations after applying rule : 162 
[1] "tivozanib"
applying rule: tivozanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 162 
Permutations after applying rule : 152 
Applying sequence restrictions to population 2 
Dropping drugs not allowed for this population.
applying rule: avelumab_plus_axitinib axitinib cabozantinib everolimus nivolumab_plus_ipilimumab lenvatinib_plus_everolimus cabozantinib_plus_nivolumab nivolumab_monotherapy pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib are only allowed treatments.
Permutations before applying rule: 26404 
Permutations after applying rule : 13344 
applying rule: drug line restrictions.
Permutations before applying rule: 13344 
Permutations after applying rule : 1036 
[1] "axitinib"
applying rule. axitinib is only allowed after avelumab_plus_axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 1036 
Permutations after applying rule : 1017 
[1] "everolimus"
applying rule. everolimus is only allowed after avelumab_plus_axitinib axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 1017 
Permutations after applying rule : 1004 
[1] "lenvatinib_plus_everolimus"
applying rule. lenvatinib_plus_everolimus is not allowed immediately after nivolumab_plus_ipilimumab 
Permutations before applying rule: 1004 
Permutations after applying rule : 984 
[1] "axitinib"
applying rule axitinib : avelumab_plus_axitinib axitinib cannot be in one permutation
Permutations before applying rule: 984 
Permutations after applying rule : 915 
[1] "cabozantinib"
applying rule cabozantinib : cabozantinib cabozantinib_plus_nivolumab cannot be in one permutation
Permutations before applying rule: 915 
Permutations after applying rule : 876 
[1] "everolimus"
applying rule everolimus : lenvatinib_plus_everolimus everolimus cannot be in one permutation
Permutations before applying rule: 876 
Permutations after applying rule : 773 
[1] "io"
applying rule io : avelumab_plus_axitinib nivolumab_plus_ipilimumab cabozantinib_plus_nivolumab nivolumab_monotherapy lenvatinib_plus_pembrolizumab cannot be in one permutation
Permutations before applying rule: 773 
Permutations after applying rule : 657 
[1] "lenvatinib_plus_everolimus"
applying rule lenvatinib_plus_everolimus : lenvatinib_plus_everolimus lenvatinib_plus_pembrolizumab cannot be in one permutation
Permutations before applying rule: 657 
Permutations after applying rule : 638 
[1] "nivolumab"
applying rule nivolumab : nivolumab_plus_ipilimumab cabozantinib_plus_nivolumab nivolumab_monotherapy cannot be in one permutation
Permutations before applying rule: 638 
Permutations after applying rule : 638 
[1] "TKIs"
applying rule TKIs : sunitinib pazopanib tivozanib cannot be in one permutation
Permutations before applying rule: 638 
Permutations after applying rule : 386 
[1] "lenvatinib_plus_everolimus"
applying rule: lenvatinib_plus_everolimus can only be after ONE of avelumab_plus_axitinib axitinib cabozantinib cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 386 
Permutations after applying rule : 359 
[1] "cabozantinib"
applying rule: cabozantinib as 2L+ only allowed after avelumab_plus_axitinib axitinib nivolumab_plus_ipilimumab lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 359 
Permutations after applying rule : 359 
[1] "pazopanib"
applying rule: pazopanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 359 
Permutations after applying rule : 339 
[1] "sunitinib"
applying rule: sunitinib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 339 
Permutations after applying rule : 319 
[1] "tivozanib"
applying rule: tivozanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 319 
Permutations after applying rule : 299 
Applying sequence restrictions to population 3 
Dropping drugs not allowed for this population.
applying rule: axitinib cabozantinib everolimus lenvatinib_plus_everolimus nivolumab_monotherapy pazopanib sunitinib tivozanib are only allowed treatments.
Permutations before applying rule: 26404 
Permutations after applying rule : 2080 
applying rule: drug line restrictions.
Permutations before applying rule: 2080 
Permutations after applying rule : 330 
[1] "axitinib"
applying rule. axitinib is only allowed after avelumab_plus_axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 330 
Permutations after applying rule : 330 
[1] "everolimus"
applying rule. everolimus is only allowed after avelumab_plus_axitinib axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 330 
Permutations after applying rule : 330 
[1] "lenvatinib_plus_everolimus"
applying rule. lenvatinib_plus_everolimus is not allowed immediately after nivolumab_plus_ipilimumab 
Permutations before applying rule: 330 
Permutations after applying rule : 330 
[1] "axitinib"
applying rule axitinib : avelumab_plus_axitinib axitinib cannot be in one permutation
Permutations before applying rule: 330 
Permutations after applying rule : 330 
[1] "everolimus"
applying rule everolimus : lenvatinib_plus_everolimus everolimus cannot be in one permutation
Permutations before applying rule: 330 
Permutations after applying rule : 288 
[1] "lenvatinib_plus_everolimus"
applying rule lenvatinib_plus_everolimus : lenvatinib_plus_everolimus lenvatinib_plus_pembrolizumab cannot be in one permutation
Permutations before applying rule: 288 
Permutations after applying rule : 288 
[1] "TKIs"
applying rule TKIs : sunitinib pazopanib tivozanib cannot be in one permutation
Permutations before applying rule: 288 
Permutations after applying rule : 120 
[1] "lenvatinib_plus_everolimus"
applying rule: lenvatinib_plus_everolimus can only be after ONE of axitinib cabozantinib pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 120 
Permutations after applying rule : 111 
[1] "cabozantinib"
applying rule: cabozantinib as 2L+ only allowed after avelumab_plus_axitinib axitinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 111 
Permutations after applying rule : 111 
[1] "pazopanib"
applying rule: pazopanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 111 
Permutations after applying rule : 111 
[1] "sunitinib"
applying rule: sunitinib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 111 
Permutations after applying rule : 111 
[1] "tivozanib"
applying rule: tivozanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 111 
Permutations after applying rule : 111 
Applying sequence restrictions to population 4 
Dropping drugs not allowed for this population.
applying rule: axitinib cabozantinib everolimus lenvatinib_plus_everolimus nivolumab_monotherapy pazopanib sunitinib tivozanib are only allowed treatments.
Permutations before applying rule: 26404 
Permutations after applying rule : 2080 
applying rule: drug line restrictions.
Permutations before applying rule: 2080 
Permutations after applying rule : 440 
[1] "axitinib"
applying rule. axitinib is only allowed after avelumab_plus_axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 440 
Permutations after applying rule : 440 
[1] "everolimus"
applying rule. everolimus is only allowed after avelumab_plus_axitinib axitinib cabozantinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 440 
Permutations after applying rule : 440 
[1] "lenvatinib_plus_everolimus"
applying rule. lenvatinib_plus_everolimus is not allowed immediately after nivolumab_plus_ipilimumab 
Permutations before applying rule: 440 
Permutations after applying rule : 440 
[1] "axitinib"
applying rule axitinib : avelumab_plus_axitinib axitinib cannot be in one permutation
Permutations before applying rule: 440 
Permutations after applying rule : 440 
[1] "everolimus"
applying rule everolimus : lenvatinib_plus_everolimus everolimus cannot be in one permutation
Permutations before applying rule: 440 
Permutations after applying rule : 384 
[1] "lenvatinib_plus_everolimus"
applying rule lenvatinib_plus_everolimus : lenvatinib_plus_everolimus lenvatinib_plus_pembrolizumab cannot be in one permutation
Permutations before applying rule: 384 
Permutations after applying rule : 384 
[1] "TKIs"
applying rule TKIs : sunitinib pazopanib tivozanib cannot be in one permutation
Permutations before applying rule: 384 
Permutations after applying rule : 198 
[1] "lenvatinib_plus_everolimus"
applying rule: lenvatinib_plus_everolimus can only be after ONE of axitinib cabozantinib pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 198 
Permutations after applying rule : 182 
[1] "cabozantinib"
applying rule: cabozantinib as 2L+ only allowed after avelumab_plus_axitinib axitinib lenvatinib_plus_everolimus cabozantinib_plus_nivolumab pazopanib lenvatinib_plus_pembrolizumab sunitinib tivozanib 
Permutations before applying rule: 182 
Permutations after applying rule : 182 
[1] "pazopanib"
applying rule: pazopanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 182 
Permutations after applying rule : 182 
[1] "sunitinib"
applying rule: sunitinib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 182 
Permutations after applying rule : 182 
[1] "tivozanib"
applying rule: tivozanib as 2L+ only allowed immediately after avelumab_plus_axitinib lenvatinib_plus_pembrolizumab cabozantinib_plus_nivolumab nivolumab_plus_ipilimumab 
Permutations before applying rule: 182 
Permutations after applying rule : 182 
```


:::
:::


:::

SPLITMD_CODE1_END

SPLITMD_CODE2_START

This page performs a **partitioned survival analysis** on the patient-level data, which is **real-world evidence (RWE)**. The analysis is performed in order to **extrapolate the survival curves** so they cover the full time horizon of the economic model (40 years).

## Import patient-level data

The patient-level data (or individual patient data (IPD)) is imported from `IPD_R_input_noACIC.xlsx`. This data has a row for each patient which states their population, line, treatment and trial, and then the time taken for them to experience an endpoint (e.g. overall survival) or be censored (i.e. stopped timing for some other reason). For more information, see the [Input data](../input_data.qmd) page.

The data is imported using `f_excel_extract()` which produces the object `wb`. This is a list with a single item: the table from the `IPD` sheet of the workbook. The code chunk converts this to a data table and filters it just the relevant columns.

As a survival time of 0 is not allowed, these are converted to 1 day (hence, `1/7` as the time unit of the analysis is weeks).


::: {.cell}

```{.r .cell-code}
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

excel_path2 <- file.path(d_path, "IPD_R_input_noACIC.xlsx")
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
```
:::


::: {.callout-note collapse="true"}

## Preview `i$surv$pld`


::: {.cell}

```{.r .cell-code}
kable(head(i$surv$pld))
```

::: {.cell-output-display}


| population| line| molecule| trial| endpoint|     timew| event_censor|
|----------:|----:|--------:|-----:|--------:|---------:|------------:|
|          0|    1|        1|     0|        0| 205.28572|            1|
|          0|    1|        1|     0|        0|  40.00357|            0|
|          0|    1|        1|     0|        0| 145.42857|            0|
|          0|    1|        1|     0|        0| 108.85714|            1|
|          0|    1|        1|     0|        0|  86.85714|            1|
|          0|    1|        1|     0|        0|  53.42857|            0|


:::
:::


:::

## Create look-ups to convert between numeric categories and labels

This section creates "look-ups" which enable us to convert between numeric and categorical variables.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the `i$id$ipd` look-up

The first step in creating these look-up tables was to find the unique values in each of the `i$r_pld_lookup_...` lists, which contain numeric categories.

This lookup simply contains the unique numeric values from the `i$r_pld_lookup_...` lists, either just as numeric values or appended with `pop`, `line`, `mol`, `trial` or `endpoint`.


::: {.cell}

```{.r .cell-code}
i$id$ipd
```

::: {.cell-output .cell-output-stdout}

```
$pop
pop_0 pop_1 pop_2 
    0     1     2 

$line
line_1 line_2 line_3 line_4 line_5 
     1      2      3      4      5 

$mol
  mol_0   mol_1   mol_2   mol_3   mol_4   mol_5   mol_6   mol_7   mol_8   mol_9 
      0       1       2       3       4       5       6       7       8       9 
 mol_10  mol_11  mol_12 mol_999 
     10      11      12     999 

$trial
trial_0 trial_1 trial_2 
      0       1       2 

$endpoint
endpoint_0 endpoint_1 endpoint_2 endpoint_3 endpoint_4 
         0          1          2          3          4 
```


:::
:::


:::


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the `i$lookup$ipd` look-up

This look-up consists of five data tables for population, treatment, molecule, trial and endpoint. Each converts between the numeric and two different categorical versions of each variable.


::: {.cell}

```{.r .cell-code}
kable(i$lookup$ipd)
```

::: {.cell-output-display}


<table class="kable_wrapper">
<tbody>
  <tr>
   <td> 

|Description              |RCC_input_desc                  | Number|
|:------------------------|:-------------------------------|------:|
|All                      |All risk groups                 |      0|
|Poor / intermediate risk |Poor or intermediate risk group |      1|
|Favourable risk          |Favourable risk group           |      2|

 </td>
   <td> 

|Description          |RCC_input_desc | Number|seq_col |R_id   |
|:--------------------|:--------------|------:|:-------|:------|
|Previously untreated |1L             |      1|V2      |line_1 |
|2nd line             |2L             |      2|V3      |line_2 |
|3rd line             |3L             |      3|V4      |line_3 |
|4th line             |4L             |      4|V5      |line_4 |
|BSC                  |5L             |      5|V6      |line_5 |

 </td>
   <td> 

|Description                   |RCC_input_desc                | Number|
|:-----------------------------|:-----------------------------|------:|
|Nivolumab monotherapy         |nivolumab_monotherapy         |      0|
|Cabozantinib plus nivolumab   |cabozantinib_plus_nivolumab   |      1|
|Nivolumab plus ipilimumab     |nivolumab_plus_ipilimumab     |      2|
|Lenvatinib plus pembrolizumab |lenvatinib_plus_pembrolizumab |      3|
|Avelumab plus axitinib        |avelumab_plus_axitinib        |      4|
|Pazopanib                     |pazopanib                     |      5|
|Tivozanib                     |tivozanib                     |      6|
|Sunitinib                     |sunitinib                     |      7|
|Cabozantinib                  |cabozantinib                  |      8|
|Lenvatinib plus everolimus    |lenvatinib_plus_everolimus    |      9|
|Everolimus                    |everolimus                    |     10|
|Axitinib                      |axitinib                      |     11|
|Sorafenib                     |sorafenib                     |     12|
|Placebo / BSC                 |placebo_BSC                   |    999|

 </td>
   <td> 

|Description         |RCC_input_desc | Number|
|:-------------------|:--------------|------:|
|CheckMate 9ER       |CheckMate_9ER  |      0|
|CheckMate 025       |CheckMate_025  |      1|
|Real world evidence |RWE            |      2|

 </td>
   <td> 

|Description |RCC_input_desc                    | Number|
|:-----------|:---------------------------------|------:|
|OS          |Overall survival                  |      0|
|PFS         |Progression free survival         |      1|
|TTD         |Time to treatment discontinuation |      2|
|TTP         |Time to progression               |      3|
|PPS         |Post progression survival         |      4|

 </td>
  </tr>
</tbody>
</table>


:::
:::


:::

This pre-existing lookup is simply copied but with a new name.


::: {.cell}

```{.r .cell-code}
i$lookup$dist <- i$r_pld_lookup_dist
```
:::


::: {.callout-note collapse="true"}

## View the `i$lookup$dist` look-up

This look-up is for the distributions used in the survival analyses.


::: {.cell}

```{.r .cell-code}
kable(i$lookup$dist)
```

::: {.cell-output-display}


|Description       |RCC_input_desc | Number|
|:-----------------|:--------------|------:|
|Generalised gamma |gengamma       |      0|
|Exponential       |exp            |      1|
|Weibull           |weibull        |      2|
|Log-normal        |lnorm          |      3|
|Gamma             |gamma          |      4|
|Gompertz          |gompertz       |      5|
|Log-logistic      |llogis         |      6|


:::
:::


:::


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the `i$lookup$trt` look-up

This look-up is for the treatments.


::: {.cell}

```{.r .cell-code}
i$lookup$trt
```

::: {.cell-output .cell-output-stdout}

```
        nivolumab_monotherapy   cabozantinib_plus_nivolumab 
                            0                             1 
    nivolumab_plus_ipilimumab lenvatinib_plus_pembrolizumab 
                            2                             3 
       avelumab_plus_axitinib                     pazopanib 
                            4                             5 
                    tivozanib                     sunitinib 
                            6                             7 
                 cabozantinib    lenvatinib_plus_everolimus 
                            8                             9 
                   everolimus                      axitinib 
                           10                            11 
                    sorafenib                           BSC 
                           12                           999 
```


:::
:::


:::

## Further pre-processing before survival analysis

### Copy items into `p`

Copy look-ups and IDs from `i` into `p`.


::: {.cell}

```{.r .cell-code}
# pass to p whenever i$lookup has been populated/updated.
p$basic$lookup <- i$lookup
p$basic$id <- i$id

# one can then simply i$lookup$trt["nivolumab"] or i$lookup$trt["sorafenib"] to 
# get the id numbers.
```
:::


### Tidy sequences for each population

Convert sequences into a data table, but making some amendments:

* Changing population to start from 0 - so `pop_0` to `pop_3` (rather than `pop_1` to `pop_4`) (as this aligns with populations elsewhere)
* Changing the columns to `line_1` to `line_5` (instead of `V1` to `V5`)
* Seperating the tables for each of the four populations


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$seq_clean`

As described above, we started with `i$sequences`:


::: {.cell}

```{.r .cell-code}
kable(head(i$sequences))
```

::: {.cell-output-display}


|V1   |V2                          |V3        |V4           |V5         |V6  |
|:----|:---------------------------|:---------|:------------|:----------|:---|
|pop1 |avelumab_plus_axitinib      |pazopanib |cabozantinib |everolimus |BSC |
|pop1 |avelumab_plus_axitinib      |sunitinib |cabozantinib |everolimus |BSC |
|pop1 |avelumab_plus_axitinib      |tivozanib |cabozantinib |everolimus |BSC |
|pop1 |cabozantinib_plus_nivolumab |pazopanib |axitinib     |everolimus |BSC |
|pop1 |cabozantinib_plus_nivolumab |pazopanib |everolimus   |axitinib   |BSC |
|pop1 |cabozantinib_plus_nivolumab |sunitinib |axitinib     |everolimus |BSC |


:::
:::


Which was tidied to create `i$seq_clean` with seperate dataframes for each population. For example, population 0:


::: {.cell}

```{.r .cell-code}
kable(head(i$seq_clean$pop_0))
```

::: {.cell-output-display}


|line_1                      |line_2    |line_3       |line_4     |line_5 |
|:---------------------------|:---------|:------------|:----------|:------|
|avelumab_plus_axitinib      |pazopanib |cabozantinib |everolimus |BSC    |
|avelumab_plus_axitinib      |sunitinib |cabozantinib |everolimus |BSC    |
|avelumab_plus_axitinib      |tivozanib |cabozantinib |everolimus |BSC    |
|cabozantinib_plus_nivolumab |pazopanib |axitinib     |everolimus |BSC    |
|cabozantinib_plus_nivolumab |pazopanib |everolimus   |axitinib   |BSC    |
|cabozantinib_plus_nivolumab |sunitinib |axitinib     |everolimus |BSC    |


:::
:::


:::

Numeric versions of the `i$seq_clean` tables are also created.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$seq_n` and `i$seq_ref`

These are numeric versions of the possible treatment sequences for each population.

They either represent each treatment just as a number (`i$seq_n`) or as mol_number (`i$seq_ref`).

For example, for population 0:


::: {.cell}

```{.r .cell-code}
kable(head(i$seq_n$pop_0))
```

::: {.cell-output-display}


| line_1| line_2| line_3| line_4| line_5|
|------:|------:|------:|------:|------:|
|      4|      5|      8|     10|    999|
|      4|      7|      8|     10|    999|
|      4|      6|      8|     10|    999|
|      1|      5|     11|     10|    999|
|      1|      5|     10|     11|    999|
|      1|      7|     11|     10|    999|


:::

```{.r .cell-code}
kable(head(i$seq_ref$pop_0))
```

::: {.cell-output-display}


|line_1 |line_2 |line_3 |line_4 |line_5  |
|:------|:------|:------|:------|:-------|
|mol_4  |mol_5  |mol_8  |mol_10 |mol_999 |
|mol_4  |mol_7  |mol_8  |mol_10 |mol_999 |
|mol_4  |mol_6  |mol_8  |mol_10 |mol_999 |
|mol_1  |mol_5  |mol_11 |mol_10 |mol_999 |
|mol_1  |mol_5  |mol_10 |mol_11 |mol_999 |
|mol_1  |mol_7  |mol_11 |mol_10 |mol_999 |


:::
:::


:::

These are then copied into `p`.


::: {.cell}

```{.r .cell-code}
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
```
:::


### Make a categorical version of the patient-level data

First, the look-ups are converted into lists.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$surv$lab_pld`

The look-up tables in `i$lookup$ipd` convert between the numeric and categorical variables, and we will use these to relabel the patient-level data. However, to do this, we need to convert them into lists, where the labels are the categories and the values are the numeric versions of each categories.

For example, the original population look-up table:


::: {.cell}

```{.r .cell-code}
kable(i$lookup$ipd$pop)
```

::: {.cell-output-display}


|Description              |RCC_input_desc                  | Number|
|:------------------------|:-------------------------------|------:|
|All                      |All risk groups                 |      0|
|Poor / intermediate risk |Poor or intermediate risk group |      1|
|Favourable risk          |Favourable risk group           |      2|


:::
:::


And the new list created from that:


::: {.cell}

```{.r .cell-code}
i$surv$lab_pld$population
```

::: {.cell-output .cell-output-stdout}

```
                     All Poor / intermediate risk          Favourable risk 
                       0                        1                        2 
```


:::
:::


:::

These look-up lists are then used to convert the patient-level data from numeric to categorical versions of each variable.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$surv$lab_pld$dat`

As mentioned, this code has converted the numeric patient-level data (`i$surv$pld`)...


::: {.cell}

```{.r .cell-code}
kable(head(i$surv$pld))
```

::: {.cell-output-display}


| population| line| molecule| trial| endpoint|     timew| event_censor|
|----------:|----:|--------:|-----:|--------:|---------:|------------:|
|          0|    1|        1|     0|        0| 205.28572|            1|
|          0|    1|        1|     0|        0|  40.00357|            0|
|          0|    1|        1|     0|        0| 145.42857|            0|
|          0|    1|        1|     0|        0| 108.85714|            1|
|          0|    1|        1|     0|        0|  86.85714|            1|
|          0|    1|        1|     0|        0|  53.42857|            0|


:::
:::


...Into categorical...


::: {.cell}

```{.r .cell-code}
kable(head(i$surv$lab_pld$dat))
```

::: {.cell-output-display}


|population |line                 |molecule                    |trial         |endpoint |     timew| event_censor|
|:----------|:--------------------|:---------------------------|:-------------|:--------|---------:|------------:|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       | 205.28572|            1|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       |  40.00357|            0|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       | 145.42857|            0|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       | 108.85714|            1|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       |  86.85714|            1|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER |OS       |  53.42857|            0|


:::
:::


:::

## Count the number of patients with each treatment, line, molecule, trial and endpoint in patient-level data

The first line of this code chunk counted the number of rows for each combination of population, line, molecule trial and endpoint (PLMTE - hence the name `n_by_plmte`).

The remaining lines convert the numeric categories into categorical versions


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$surv$n_by_plmte`


::: {.cell}

```{.r .cell-code}
kable(head(i$surv$n_by_plmte))
```

::: {.cell-output-display}


|population |line                 |molecule                    |trial               |endpoint |   N|
|:----------|:--------------------|:---------------------------|:-------------------|:--------|---:|
|All        |Previously untreated |Nivolumab monotherapy       |Real world evidence |OS       |  86|
|All        |Previously untreated |Nivolumab monotherapy       |Real world evidence |PFS      |  86|
|All        |Previously untreated |Nivolumab monotherapy       |Real world evidence |TTD      |  84|
|All        |Previously untreated |Nivolumab monotherapy       |Real world evidence |TTP      |  86|
|All        |Previously untreated |Nivolumab monotherapy       |Real world evidence |PPS      |  16|
|All        |Previously untreated |Cabozantinib plus nivolumab |CheckMate 9ER       |OS       | 323|


:::
:::


:::

## Run model

In `Model_Structure.R`, this section is not run as `i$dd_run_surv_reg` is set to "No", as set in the excel workbook.

However, if run, this function would:

**1. Run survival analysis**. Using `f_surv_runAllTSD14()` from `survival/` (see dropdown below). "TSD14" refers to technical support document 14 which is a methods guide from NICE for performing survival analysis, [available here](https://www.sheffield.ac.uk/nice-dsu/tsds/survival-analysis). The function;

* Runs through all possible populations, lines, molecules, trials and endpoints in `id` (`i$id$ipd`, the lookup of unique values for each of those items)
* If there is data available in `r_pld` (`i$surv$pld`, the patient level data), then it performs survival analysis using the function `flexsurv::flexsurvreg()`. This data will need to meet the threshold you set for the number of observations (default `28`).
* It repeats this for each with each of the distributions in `distnames`
* It saves the coefficients, variance covariance matrix, and goodness of fit statistics, and saves these as `fs_fits` (and also `gof`)
* The survival curves are then extrapolated using the fitted models over the specified time cycle (`t_cyc` (`p$basic$t_cyc`)) using the function `f_extrapolate()`. These are saved in the matrix `st` (survival at time t, or st for short)
* If creating plots, this is done using the function `f_extrap_plot()`
* Finally, the function returns the results for each combination

**2. Seperately, manually run survival analysis for best supportive care (BSC)**. This is done seperately as there is very little information available on BSC overall survival (OS). Hence, the best source of data is to use the pooled post-progression survival (PPS) data from 4L (fourth line) patients (since they are unlikely to receive something after that treatment). It has similar steps to `f_surv_runAllTSD14()`.


::: {.cell}

```{.r .cell-code}
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
  
  saveRDS(i$surv$reg, file = file.path(d_path, "Survival_analysis.rds"))
  
}
```
:::


## Load pre-run survival analysis

This section loads a provided pre-run survival analysis. It then limits the matrices with the extrapolated survival curves to the time horizon of the study (40 years).


::: {.cell}

```{.r .cell-code}
# to load in pre-run survival analysis select the RDS file here

# option to load from pre-specified file path on local machine, uncomment this and comment out the line below to use

RDS_path <- file.path(d_path, "survival_analysis_no_ipd_CompanyTTDTTPPPS_redacted.rds")
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
```
:::


::: {.callout-note collapse="true"}

## View `i$surv$reg`

The result of the survival analysis is a large nested list.

Below is an example of the result for population 2 (favourable risk) 1L treatment with avelumab plus axitinib, with an endpoint of progression-free survival.

For each distribution, it has `$coefs`, `$vcov` and `$fit` - for example, for weibull:


::: {.cell}

```{.r .cell-code}
i$surv$reg$pop_2$line_1$mol_4$trial_2$endpoint_1$fs_fits$weibull
```

::: {.cell-output .cell-output-stdout}

```
$coefs
    shape     scale 
0.1368876 4.9253396 

$vcov
            shape       scale
shape  0.02973812 -0.01789284
scale -0.01789284  0.04001588

$fit
      AIC       BIC    logLik 
 316.2418  320.5281 -156.1209 
```


:::
:::


The fit of each distribution is also summarised in a single table:


::: {.cell}

```{.r .cell-code}
kable(i$surv$reg$pop_2$line_1$mol_4$trial_2$endpoint_1$gof)
```

::: {.cell-output-display}


|         |      AIC|      BIC|    logLik|
|:--------|--------:|--------:|---------:|
|gengamma | 317.6759| 324.1053| -155.8380|
|exp      | 314.8361| 316.9792| -156.4180|
|weibull  | 316.2418| 320.5281| -156.1209|
|lnorm    | 316.5742| 320.8604| -156.2871|
|gamma    | 316.0581| 320.3443| -156.0290|
|gompertz | 316.8228| 321.1090| -156.4114|
|llogis   | 315.2343| 319.5206| -155.6172|


:::
:::


And finally, the survival times for each distribution are provided:


::: {.cell}

```{.r .cell-code}
kable(head(i$surv$reg$pop_2$line_1$mol_4$trial_2$endpoint_1$st))
```

::: {.cell-output-display}


|  gengamma|       exp|   weibull|     lnorm|     gamma|  gompertz|    llogis|
|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
| 1.0000000| 1.0000000| 1.0000000| 1.0000000| 1.0000000| 1.0000000| 1.0000000|
| 0.9985777| 0.9933915| 0.9964812| 0.9996291| 0.9972167| 0.9935903| 0.9980669|
| 0.9958140| 0.9868266| 0.9922259| 0.9979763| 0.9935048| 0.9872169| 0.9950622|
| 0.9922663| 0.9803051| 0.9876526| 0.9950948| 0.9893592| 0.9808795| 0.9914695|
| 0.9881616| 0.9738267| 0.9828687| 0.9912317| 0.9849185| 0.9745780| 0.9874457|
| 0.9836342| 0.9673911| 0.9779288| 0.9866047| 0.9802557| 0.9683123| 0.9830798|


:::
:::


If `draw_plots` was set to TRUE for `f_surv_runAllTSD14()` when this dataset was produced, then there is also a plot available - but in this case, it was not.


::: {.cell}

```{.r .cell-code}
i$surv$reg$pop_2$line_1$mol_4$trial_2$endpoint_1$plot
```

::: {.cell-output .cell-output-stdout}

```
NULL
```


:::
:::


:::

There are also some comments afterwards would provide some further details around the survival analysis.


::: {.cell}

```{.r .cell-code}
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
```
:::


## Make survival analysis report

If `i$dd_report_req_surv_reg=="Yes"` then a report is create - by default, this was set to "No".

This code will only work if you have run `f_surv_runAllTSD14()` with `draw_plots` set to TRUE. As the provided pre-run survival analysis does not include plots, this section will not run, and encounters an error if tried:

```
Error in optim(method = "BFGS", par = c(mu = 5.82584217367559, sigma = -0.208272130470619,  : 
  non-finite finite-difference value [1]
```


::: {.cell}

```{.r .cell-code}
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
  print(doc_surv, target = file.path(o_path, "Survival_Analysis.docx"))
  
  rm(doc_surv)
}
```
:::


SPLITMD_CODE2_END

SPLITMD_CODE3_START

On this page, the hazard ratios (HR) from some pre-run network meta-analyses (NMA) are applied to the extrapolated RWE survival curves generated on the previous page.

## Import and process results of proportional hazards NMA (PH NMA)

### Import results

Import results from the PH NMA, as described in [Input data](../input_data.qmd). The code chunk sets new column names for the result table, and creates a duplicate of the `Endpoint` column called `Reference.endpoint`.


::: {.cell}

```{.r .cell-code}
# 3.3.5 Comparative efficacy propagation (NMA) ---------------------------------------------------------------

# Pull in the data and calculate means by pop line mol endpoint reftrt and reftrial

# First read in RDS file containing the PH NMA coda samples


# 3.3.5.1.1 PH NMA data -----------------------------------------------------

# Option to read in PH NMA CODA from local machine, uncomment this and comment out the line below to use
RDS_path2 <- file.path(d_path, "PH_NMA_CODA.rds")
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
```
:::


::: {.callout-note collapse="true"}

## View `i$PHNMA$data`


::: {.cell}

```{.r .cell-code}
kable(head(i$PHNMA$data))
```

::: {.cell-output-display}


| Run| Population| Line| Molecule| Endpoint| Reference.treatment| Reference.trial|        HR| Reference.endpoint|
|---:|----------:|----:|--------:|--------:|-------------------:|---------------:|---------:|------------------:|
|   1|          0|    1|        4|        0|                   7|               0| 0.8958580|                  0|
|   1|          0|    1|        1|        0|                   7|               0| 0.7821468|                  0|
|   1|          0|    1|        8|        0|                   7|               0| 1.4324613|                  0|
|   1|          0|    1|        2|        0|                   7|               0| 0.7571763|                  0|
|   1|          0|    1|        5|        0|                   7|               0| 0.9492596|                  0|
|   1|          0|    1|        3|        0|                   7|               0| 0.8561836|                  0|


:::
:::


:::

### Set 3L to 2L, and TTD and TTP to PFS

This sections implements some assumptions where we the relative effectiveness of one group/endpoint is assumed to apply to another.

**2L and 3L**: The effectiveness of 3L treatments is assumed to be equal to their effectiveness at 2L, so it simply copies the data the effectiveness data from 2L (`i$PHNMA$data[Line==2,]`) but replaces the `Line` with `3` and appends it to the main data table (`i$PHNMA$data`).

**PFS, TTD and TTP**: The time to treatment discontinuation (TTD) and time to progression (TTP) HRs are assumed to be equal to the HR for progression-free survival (PFS). If we refer to the endpoint look-up table, we can see that PFS is endpoint 1, then TTD is 2 and TTP is 3. Hence, in the code below, the PFS results (`i$PHNMA$data[Endpoint==1,]`) are simply copied, but with the endpoint replaced with 2 or 3, before appending it back to the main data table.

::: {.callout-note collapse="true"}

## View endpoint lookup table


::: {.cell}

```{.r .cell-code}
i$lookup$ipd$endpoint
```

::: {.cell-output-display}

:::
:::


:::


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `i$PHNMA$data`

To view are implemented assumptions we can, for example, look at the results from 2L and 3L, and see that they are the same:


::: {.cell}

```{.r .cell-code}
kable(head(i$PHNMA$data[i$PHNMA$data$Line==2,]))
```

::: {.cell-output-display}


| Run| Population| Line| Molecule| Endpoint| Reference.treatment| Reference.trial|        HR| Reference.endpoint|
|---:|----------:|----:|--------:|--------:|-------------------:|---------------:|---------:|------------------:|
|   1|          0|    2|       11|        0|                  10|               1| 2.8940210|                  0|
|   1|          0|    2|        8|        0|                  10|               1| 0.6664276|                  0|
|   1|          0|    2|        9|        0|                  10|               1| 0.8515765|                  0|
|   1|          0|    2|        0|        0|                  10|               1| 0.6339945|                  0|
|   1|          0|    2|      999|        0|                  10|               1| 1.4503489|                  0|
|   1|          0|    2|       12|        0|                  10|               1| 2.9477019|                  0|


:::

```{.r .cell-code}
kable(head(i$PHNMA$data[i$PHNMA$data$Line==3,]))
```

::: {.cell-output-display}


| Run| Population| Line| Molecule| Endpoint| Reference.treatment| Reference.trial|        HR| Reference.endpoint|
|---:|----------:|----:|--------:|--------:|-------------------:|---------------:|---------:|------------------:|
|   1|          0|    3|       11|        0|                  10|               1| 2.8940210|                  0|
|   1|          0|    3|        8|        0|                  10|               1| 0.6664276|                  0|
|   1|          0|    3|        9|        0|                  10|               1| 0.8515765|                  0|
|   1|          0|    3|        0|        0|                  10|               1| 0.6339945|                  0|
|   1|          0|    3|      999|        0|                  10|               1| 1.4503489|                  0|
|   1|          0|    3|       12|        0|                  10|               1| 2.9477019|                  0|


:::
:::


:::

### Find mean HRs

The mean HR by population, line, molecule, endpoint, reference treatment and reference trial are calculated.

These are then add to `p`.


::: {.cell}

```{.r .cell-code}
# Calculate the mean from the CODA samples for deterministic analysis

i$PHNMA$means <- i$PHNMA$data[,.(HR = mean(HR)),by=list(Population,Line,Molecule,Endpoint,Reference.treatment,Reference.trial)]


# 3.3.5.1.2 DETERMINISTIC CODA --------------------------------------------

# for the deterministic analysis we use the means. 
p$releff$CODA$PH <- i$PHNMA$means
```
:::


::: {.callout-note collapse="true"}

## View `i$PHNMA$means`


::: {.cell}

```{.r .cell-code}
kable(head(i$PHNMA$means))
```

::: {.cell-output-display}


| Population| Line| Molecule| Endpoint| Reference.treatment| Reference.trial|        HR|
|----------:|----:|--------:|--------:|-------------------:|---------------:|---------:|
|          0|    1|        4|        0|                   7|               0| 0.7935352|
|          0|    1|        1|        0|                   7|               0| 0.7033510|
|          0|    1|        8|        0|                   7|               0| 0.8178626|
|          0|    1|        2|        0|                   7|               0| 0.7226217|
|          0|    1|        5|        0|                   7|               0| 0.9224925|
|          0|    1|        3|        0|                   7|               0| 0.7946571|


:::
:::


:::

## Import and process results of fractional polynomial NMA (FP NMA)

### Import results

Import results from the FP NMA, as described in [Input data](../input_data.qmd). The code chunk sets new column names for the result table, and converts the time from weeks to month (`*52` as 52 weeks per year, then `/12` to get months).


::: {.cell}

```{.r .cell-code}
# 3.3.5.2.1 FP NMA data -----------------------------------------------------

# Load in FP NMA data
i$FPNMA <- list()

#read in means for deterministic and PSA parameters for probabilistic

# option to read in from local machine, uncomment the below and comment out line 949 to use
RDS_path3 <- file.path(d_path, "FPNMA_means.rds")
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
```
:::


::: {.callout-note collapse="true"}

## View `i$FPNMA$means`


::: {.cell}

```{.r .cell-code}
kable(head(i$FPNMA$means))
```

::: {.cell-output-display}


| time| Molecule| Reference.treatment| Line| Endpoint| Population| Reference.trial|        HR|
|----:|--------:|-------------------:|----:|--------:|----------:|---------------:|---------:|
|    2|        4|                   7|    1|        1|          0|               0| 6.0866273|
|    3|        4|                   7|    1|        1|          0|               0| 0.7948374|
|    4|        4|                   7|    1|        1|          0|               0| 0.6307802|
|    5|        4|                   7|    1|        1|          0|               0| 0.6127038|
|    6|        4|                   7|    1|        1|          0|               0| 0.6141338|
|    7|        4|                   7|    1|        1|          0|               0| 0.6184885|


:::
:::


:::

### Set cabozantinib as reference


::: {.cell}

```{.r .cell-code}
# means

# Rebasing to allow use of cabo as reference treatment in 2nd line
# repeats for means (stored in i which are later transferred to p) 
i$FPNMA$means <- f_rebase_for_cabo_as_ref_in_2L(FPNMAdata = i$FPNMA$means)

# Remove the now redundant objects we made in order to do this
```
:::


### Set 3L to 2L

As implemented for PH NMA, the HRs for 3L are assumed to be the same as at 2L, and this is implemented here using `f_3L_rel_effect_same_as_2L()`.


::: {.cell}

```{.r .cell-code}
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
i$FPNMA$means <- f_3L_rel_effect_same_as_2L(FPNMAdata = i$FPNMA$means) 

# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
# IMPORTANT: 3L relative effectiveness is assumed the same as 2L!!!!
```
:::


### Find unique combinations

`f_gen_destinations()` finds unique combinations of population, line, molecule, endpoint, reference treatment and reference trial. A duplicate version of these results is made where the referebce.trial is set to 2 (real-world evidence (RWE)). These are combined and saved as `i$FPNMA$destinations`.


::: {.cell}

```{.r .cell-code}
# Create 1 row for each destination PLMTE, so that we know where to put the
# fp data without having to iterate much
i$FPNMA$destinations <- f_gen_destinations(fp_data = i$FPNMA$means)
```
:::


::: {.callout-note collapse="true"}

## View `i$FPNMA$destinations`


::: {.cell}

```{.r .cell-code}
kable(head(i$FPNMA$destinations))
```

::: {.cell-output-display}


| Population| Line| Molecule| Endpoint| Reference.treatment| Reference.trial|
|----------:|----:|--------:|--------:|-------------------:|---------------:|
|          0|    1|        4|        1|                   7|               0|
|          0|    1|        8|        1|                   7|               0|
|          0|    1|        3|        1|                   7|               0|
|          0|    1|        1|        1|                   7|               0|
|          0|    1|        2|        1|                   7|               0|
|          0|    1|        5|        1|                   7|               0|


:::
:::


:::

### Duplicate results with RWE as reference trial

`f_add_reference_trial_2()` creates a duplicate version of `i$FPNMA$means` where the trial is set to 2 (real-world evidence).


::: {.cell}

```{.r .cell-code}
# add in reference.trial 2
i$FPNMA$means <- f_add_reference_trial_2(fp_data = i$FPNMA$means)
```
:::


### Copy into p

Copy items into `p`, removing molecules with "NA" and limiting time horizon.


::: {.cell}

```{.r .cell-code}
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
```
:::


## Create nested list with PH NMA and FP NMA results

### Set-up empty list

This section uses `f_NMA_generateNetwork()` to create an empty list with the unique values for population > line > molecule > trial > endpoint (PLMTE). For each of these PLMTE, it then has:

* `$dest` - which holds the numbers for the PLMTE
* `$orig` - empty list of PLMTE plus dist and source
* `$hr` - set to 1
* `$fp` - empty list


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the empty `p$releff$network`

See example of PMTE 0 and L1:


::: {.cell}

```{.r .cell-code}
p$releff$network$pop_0$line_1$mol_0$trial_0$endpoint_0
```

::: {.cell-output .cell-output-stdout}

```
$dest
$dest$pop
[1] "pop_0"

$dest$line
[1] "line_1"

$dest$mol
[1] "mol_0"

$dest$trial
[1] "trial_0"

$dest$endpoint
[1] "endpoint_0"


$orig
$orig$pop
NULL

$orig$line
NULL

$orig$mol
NULL

$orig$trial
NULL

$orig$endpoint
NULL

$orig$dist
NULL

$orig$source
NULL


$hr
[1] 1

$fp
list()
```


:::
:::


:::

### Create numeric trial column in `i$R_table_eff_data_settings`

This generates a column `Trial` in the table `i$R_table_eff_data_settings` if it doesn't already exist. Whilst this  isn't used in making the empty nested list, it will be later used by`f_NMA_AddAssumptionsToNetwork()` and `f_releff_PropNetwork()`.

The new `Trial` column is a numeric version of the column `i$R_table_eff_data_settings$Trial.name.if.effectiveness.source.is.trial`. The conversion between the categorical and numeric versions is done using the look-up `i$lookup$ipd$trial`.


::: {.cell}

```{.r .cell-code}
# Generate DESTINATION trial number column if it doesn't already exist:
if(!"Trial" %in% colnames(i$R_table_eff_data_settings)) {
  i$R_table_eff_data_settings$Trial <- i$lookup$ipd$trial$Number[match(i$R_table_eff_data_settings$Trial.name.if.effectiveness.source.is.trial,i$lookup$ipd$trial$Description)]
}
```
:::


::: {.callout-note collapse="true"}

## View the new column

Filtering to the unique trial columns, we can see that this has either set Trial to 0 or 2.


::: {.cell}

```{.r .cell-code}
kable(unique(i$R_table_eff_data_settings %>% select(Trial.name.if.effectiveness.source.is.trial, Trial)))
```

::: {.cell-output-display}


|   |Trial.name.if.effectiveness.source.is.trial | Trial|
|:--|:-------------------------------------------|-----:|
|1  |0                                           |    NA|
|8  |Real world evidence                         |     2|


:::
:::


:::

### Add NMA HRs

These functions add the HRs from the PH NMA (`p$releff$CODA$PH`) and FP NMA (`p$releff$CODA$FP`) to the nested list.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the updated `p$releff$network`

See example:


::: {.cell}

```{.r .cell-code}
# PH NMA HR
p$releff$network$pop_0$line_1$mol_4$trial_0$endpoint_0$hr
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.7935352
```


:::

```{.r .cell-code}
# Time-varying FP NMA HRs
head(p$releff$network$pop_0$line_1$mol_4$trial_0$endpoint_0$fp$HR, 20)
```

::: {.cell-output .cell-output-stdout}

```
 [1] 0.9321210 0.7399991 0.7015206 0.6987356 0.7068145 0.7182071 0.7301949
 [8] 0.7417350 0.7524420 0.7622070 0.7710399 0.7790003 0.7861653 0.7926152
[15] 0.7984263 0.8036691 0.8084066 0.8126947 0.8165829 0.8201147
```


:::
:::


:::

## Add more "assumption-based" HRs

The column `i$R_table_eff_data_settings$Effectiveness.data.source` sets some of the sources of HRs.

::: {.callout-note collapse="true"}

## View `i$R_table_eff_data_settings$Effectiveness.data.source`


::: {.cell}

```{.r .cell-code}
table(i$R_table_eff_data_settings$Effectiveness.data.source)
```

::: {.cell-output .cell-output-stdout}

```

                      0             Apply HR to         Assume equal to 
                    247                      49                      44 
                 FP_NMA                  PH_NMA Trial survival analysis 
                     22                      39                      19 
```


:::
:::


:::

If the setting `dd_use_PHnma_for_FPnma` was set to "Yes" in the excel spreadsheet, then it will replace "FP_NMA" with "PH_NMA" in this list - thereby, forcing the model to use the PH NMA HRs instead of including some from FP NMA.

By default, this is currently set to "No".


::: {.cell}

```{.r .cell-code}
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
```
:::


The function `f_NMA_AddAssumptionsToNetwork()` adds the "assumption-based" HRs according to `i$R_table_eff_data_settings$Effectiveness.data.source`. There are five categories within this:

* **Trial** - `Trial survival analysis` - for these entries, they are a reference curve, so the function checks the origin and destination match, and then add some extra information from excel (e.g. dist, source, setting `fp` to empty list)
* **FP NMA or PH NMA** - `FP_NMA`, `PH_NMA` - as I’ve understood it, this populates the hazard ratio for a given comparison using the FP/PH NMA HR from the same comparison in a different trial - so, even though the HR data has come from a different context (different trial), we are saying it is still applicable to this treatment comparison.
* **ET flag** - `Assume equal to` - the treatment and reference treatment are considered to be equally effective, so the HR is set to 1
* **AHR flag** - `Apply HR to` - these entries have a HR provided in the Excel spreadsheet under a column `HR.to.apply`


::: {.cell}

```{.r .cell-code}
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
```
:::


## Apply HRs to the extrapolated RWE survival times

The function `f_surv_getExtrapolations()` creates a nested list with survival times `st` for each PLMTE, as copied from `i$surv$reg` (the RWE survival analysis).


::: {.cell}

```{.r .cell-code}
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
```
:::


The function `f_releff_PropNetwork()` then applies the HR in `p$releff$network` to the survival times, saving the resulting nested list as `p$surv$st`.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `p$surv$st`

This list is similar to `p$releff$network` but with the addition of survivals, as demonstrated below: 


::: {.cell}

```{.r .cell-code}
head(p$surv$st$pop_0$line_3$mol_0$trial_2$endpoint_2$st, 20)
```

::: {.cell-output .cell-output-stdout}

```
 [1] 1.0000000 0.9921156 0.9782498 0.9609062 0.9410968 0.9194740 0.8965183
 [8] 0.8726074 0.8480469 0.8230890 0.7979435 0.7727850 0.7477591 0.7229864
[15] 0.6985662 0.6745796 0.6510916 0.6281539 0.6058063 0.5840786
```


:::
:::


:::

## Apply modifications to the survival times

### Prior adjuvant treatments

By default, this is not run as `i$dd_adjforprioradjuvant` is set to "No".

However, if "Yes", it would adjust the survival times again for endpoint 0 (overall survival, OS) by applying another hazard ratio that is associated with prior adjuvant treatment. 

An adjuvant treatment is an additional treatment given alongside the primary treatment.


::: {.cell}

```{.r .cell-code}
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
```
:::


### Effect waning

Treatment effect waning is the loss of a treatment's benefits over time. It is only applied to some PLMTEs, as specified by `i$R_table_TE_waning_settings$apply.waning`.

This is implemented using the function `f_surv_twaning_apply()`, which reduces the effect to match the reference treatment once it reaches a certain cycle.


::: {.cell}

```{.r .cell-code}
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
```

::: {.cell-output .cell-output-stderr}

```
Warning in treatment_effect_waning_with_absolute(surv_active = plmte$st, :
Hazard rate in the active treatment is more than the hazard rate in the
reference treatment in 1890 cycles, the first of which is 198 . Using the
higher of the two hazards where they cross!
```


:::
:::


### General population mortality

In this section, the survival curves are calculated for the general population, and these are then used to adjust the treatment survival curves.

#### Generate the general population survival curves

The survival curves are for overall survival. They are calculated for each population and each line of therapy. They can be created using one of two functions:

* `f_surv_GenOSLines_det()` - calculates curves based on characteristics provided in a table
* `f_surv_GenOSLines_ipd()` - calculates curves based on individual patient data

The default option is to use patient-level data (`i$dd_age_sex_source == "PLD")`.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `p$surv$gpop`

As an example, the overall survival times for population 0 line 1 from the general population are:


::: {.cell}

```{.r .cell-code}
head(p$surv$gpop$pop_0$line_1$os, 10)
```

::: {.cell-output .cell-output-stdout}

```
 [1] 1.0000000 0.9996541 0.9993084 0.9989626 0.9986170 0.9982715 0.9979260
 [8] 0.9975807 0.9972355 0.9968905
```


:::
:::


:::

#### Compare curves generated by means v.s. patient-level data

If running the script in `qc_mode`, this code chunk will compare the general population survival curves produced based on average characteristics in a table to those produced from person-level data.


::: {.cell}

```{.r .cell-code}
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
```
:::


#### Adjust treatment curves based on general population

The function `f_surv_gpopadjust()` adjusts the survival curves based on the general population survival curves. This is to ensure that the advanced renal cell carcinoma (aRCC) patients are not projected to live longer than the general population.


::: {.cell}

```{.r .cell-code}
# adjust all OS lines in p$surv$st for gpop mortality
p$surv$st <- f_surv_gpopadjust(st      = p$surv$st,
                               gpop    = p$surv$gpop,
                               method  = "hazardmax",
                               verbose = qc_mode)
```
:::


### Curve overlap

A known limitation of partitioned survival analysis (which was used to extrapolate the RWE) is that it can produce curves where PFS lies above OS (which is impossible in real-life). Hence, in cases where this occurs, it was adjusted so that PFS <= OS (and also, TTD <= OS, and PFS <= TTP). This was implemented using the functions `f_surv_PFSxOS()`, `f_surv_TTDxOS()`, and `f_surv_PFSxTTP()`.


::: {.cell}

```{.r .cell-code}
# Adjust PFS and TTD for OS - PFS and TTD cannot exceed OS
p$surv$st <- f_surv_PFSxOS(p$surv$st, if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})

p$surv$st <- f_surv_TTDxOS(p$surv$st, if(i$dd_adj_cross_curves == "Use hazards"){"hazardmax"} else{"abs"})

# Adjust TTP, PFS cannot go above TTP, this is done on absolute survival rather than allowing flexibility to look at hazards
p$surv$st <- f_surv_PFSxTTP(st = p$surv$st,method =  "abs")

# The below should produce a string of positive or 0s
# p$surv$st$pop_0$line_1$mol_7$trial_2$endpoint_3$st - p$surv$st$pop_0$line_1$mol_7$trial_2$endpoint_1$st
```
:::


### Set 5L to 4L

Here, the survival times for 5L (ie. best supportive care, BSC) are set to those from 4L.


::: {.cell}

```{.r .cell-code}
# Last assumption - 5L =4L. this moves over BSC when it comes after 4 active treatments.

p$surv$st <- lapply(p$surv$st, function(popu) {
  popu$line_5 <- popu$line_4
  return(popu)
})
```
:::


## View the generated curves

If running with `qc_mode`, this would plot survival extrapolations and hazard curves for each PLMTE.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View some example plots

Although `qc_mode` is set to FALSE, we can still use the functions to view an example of a survival and hazard curves.


::: {.cell}

```{.r .cell-code}
f_qc_surv_ExtrapPlot(
  st   = p$surv$st,
  popu = "pop_0",
  li   = "line_1",
  mo   = "mol_1",
  tr   = "trial_2",
  t_yr = p$basic$t_yr,
  th = p$basic$th
)
```

::: {.cell-output-display}
![](code0a_walkthrough_files/figure-html/unnamed-chunk-75-1.png){width=672}
:::

```{.r .cell-code}
f_qc_surv_EstHazPlot(
  st   = p$surv$st,
  gpop = p$surv$gpop,
  popu = "pop_0",
  li   = "line_1",
  mo   = "mol_1",
  tr   = "trial_2",
  t_yr = p$basic$t_yr,
  th   = p$basic$th
)
```

::: {.cell-output-display}
![](code0a_walkthrough_files/figure-html/unnamed-chunk-75-2.png){width=672}
:::
:::


:::

SPLITMD_CODE3_END

SPLITMD_CODE4_START

This section processes data on demographics, QALYs, costs, AEs, costs, subsequent treatments, and population mappings.

SPLITMD_CODE4_END

SPLITMD_CODE5_START

SPLITMD_CODE5_END