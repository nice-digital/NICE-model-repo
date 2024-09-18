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

This page performs a **partitioned survival analysis** on the patient-level data, which is **real-world evidence (RWE)** sent by the data owners for that observational study. The analysis is performed in order to **extrapolate the survival curves** so they cover the full time horizon of the economic model (40 years).

## Import patient-level data

The patient-level data (or individual patient data (IPD)) is imported from `IPD_R_input_noACIC.xlsx`. This data has a row for each statement which states their population, line, treatment and trial, and then the time taken for them to experience an endpoint (e.g. overall survival) or be censored (i.e. stopped timing for some other reason). For more information, see the [Input data](../input_data.qmd) page.

The data is imported using `f_excel_extract()` which produces the object `wb` which is a list with a single item: the table from the `IPD` sheet of the workbook. The code chunk converts this to a data table and filters it just the relevant columns.

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
```
:::


## Further pre-processing before survival analysis


::: {.cell}

```{.r .cell-code}
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
```
:::


## Run model


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