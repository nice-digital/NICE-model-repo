
This section processes data on demographics, quality-adjusted life years (QALYs), costs, adverse events (AEs), costs, subsequent treatments, and population mappings.

## Pre-progression deaths

Add a parameter to `p` which provides the proportion of deaths in the progress-free survival group (so assuming there are some pre-progression deaths)

::: {.callout-note collapse="true"}

## View `i$dd_prop_deathinPFS`


::: {.cell}

```{.r .cell-code}
i$dd_prop_deathinPFS
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.1797068
```


:::
:::


:::


::: {.cell}

```{.r .cell-code}
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
```
:::


## Demographics

The population demographics are provided in a table which was from Excel, `i$R_table_ptchar`. The function `f_cleaning_ptchar()` converts this table into a nested list `p$demo$agg`.


::: {.cell}

```{.r .cell-code}
# 3.4.1 Demographics ------------------------------------------------------

# Demographics are simply processed from the tables in Excel.

p$demo$agg <- f_cleaning_ptchar(i$R_table_ptchar, i$lookup)

# Deterministic version is very easy.
p$demo$live <- p$demo$agg
```
:::


::: {.callout-note collapse="true"}

## View `p$demo$agg`

The original demographics table:


::: {.cell}

```{.r .cell-code}
kable(head(i$R_table_ptchar))
```

::: {.cell-output-display}


|Population               | Treatment.line| Starting.age..years..Mean| Starting...female.Mean| Starting...PorI.risk.Mean| Body.weight..kg..Mean| Prior.IO...in.12.months.Mean| Starting.age..years..n| Starting...female.n| Starting...PorI.risk.n| Body.weight..kg..n| Starting...PIR.n| Prior.IO...in.12.months.n| Starting.age..years..SE| Starting...female.SE| Starting...PorI.risk.SE| Body.weight..kg..SE| Prior.IO...in.12.months.SE|
|:------------------------|--------------:|-------------------------:|----------------------:|-------------------------:|---------------------:|----------------------------:|----------------------:|-------------------:|----------------------:|------------------:|----------------:|-------------------------:|-----------------------:|--------------------:|-----------------------:|-------------------:|--------------------------:|
|All                      |              1|                     64.40|              0.2903715|                 0.7755725|              83.38000|                            0|                   1311|                1319|                   1310|                114|             1310|                      1319|               0.2856284|            0.0124989|               0.0115269|            1.686593|                          0|
|Poor / intermediate risk |              1|                     64.20|              0.2952756|                 1.0000000|              81.26353|                            0|                   1011|                1016|                      0|                114|                0|                      1319|               0.3273973|            0.0143112|               0.0000000|            1.686593|                          0|
|Favourable risk          |              1|                     65.40|              0.2653061|                 0.0000000|              90.98037|                            0|                    293|                 294|                      0|                114|                0|                      1319|               0.5596696|            0.0257486|               0.0000000|            1.686593|                          0|
|All                      |              2|                     63.04|              0.2848101|                 0.0000000|              83.38000|                            0|                     NA|                  NA|                     NA|                 NA|               NA|                        NA|               0.4186000|            0.0179527|               0.0000000|            1.686593|                          0|
|All                      |              3|                     62.62|              0.2850467|                 0.0000000|              83.38000|                            0|                     NA|                  NA|                     NA|                 NA|               NA|                        NA|               0.7288700|            0.0308596|               0.0000000|            1.686593|                          0|
|All                      |              4|                     62.37|              0.2962963|                 0.0000000|              83.38000|                            0|                     NA|                  NA|                     NA|                 NA|               NA|                        NA|               1.2498500|            0.0621386|               0.0000000|            1.686593|                          0|


:::
:::


And a snippet from the nested list:


::: {.cell}

```{.r .cell-code}
p$demo$agg$pop_0$line_1
```

::: {.cell-output .cell-output-stdout}

```
$age
$age$mean
[1] 64.4

$age$se
[1] 0.2856284

$age$n
NULL


$pr_fem
$pr_fem$mean
[1] 0.2903715

$pr_fem$se
NULL

$pr_fem$n
NULL


$weight
$weight$mean
NULL

$weight$se
NULL

$weight$n
NULL


$prior_io
$prior_io$mean
NULL

$prior_io$se
NULL

$prior_io$se
NULL


$pr_i_rsk
$pr_i_rsk$mean
[1] 0.7755725

$pr_i_rsk$se
[1] 0.01152693

$pr_i_rsk$n
NULL
```


:::
:::


:::

## QALYs

The function `add_population_utility_params()` creates `p$util$pop_norms`, which are the EQ-5D (i.e. utility value) population norms for females and males. These are produced by `extract_utility_ageadjust_coefs()`, which generates age-adjusted utility values.


::: {.cell}

```{.r .cell-code}
# 3.4.2 QALYs -------------------------------------------------------------

# Utilities are applied to the disease model by treatment by line and whether the patient is on or off treatment
# Age adjustment is conducted multiplicatively in line with DSU guidance using earlier defined patient characteristics for age and sex

# Extracting from excel file 

p <- add_population_utility_params(p, psa = FALSE, .i = i)

# Pre-calculate the population utility norms since they will be the same across
# all sequences (though may vary across populations), and store in p
```
:::


::: {.callout-note collapse="true"}

## View `p$util$pop_norms`

As an example, the female values:


::: {.cell}

```{.r .cell-code}
p$util$pop_norms$female
```

::: {.cell-output .cell-output-stdout}

```
$l_coefficients
$l_coefficients[[1]]
    Age/10 (Age/10)^2  intercept 
   -0.0774     0.0064     0.2990 

$l_coefficients[[2]]
    Age/10 (Age/10)^2  intercept 
   -0.0147    -0.0003     0.8708 

$l_coefficients[[3]]
    Age/10 (Age/10)^2  intercept 
    0.2043    -0.0241     1.1659 


$l_mix_coef
$l_mix_coef[[1]]
   Age/10 intercept 
   0.4028   -4.4767 

$l_mix_coef[[2]]
   Age/10 intercept 
   0.1937   -1.3549 


$v_sigma
[1] 0.1282 0.0831 0.5230
```


:::
:::


:::

We start then with `base_utility`, which is a table with 2089 rows (which is the time horizon `TH`), where each row is a cycle. For each, the utility is set to 1. This basically means that every individual in the population starts with a full quality of life (utility value of 1).

However, if `i$dd_ageadjuutilities=="Yes"`, then it is adjusted based on the age and sex of either:

* The population (if based on means) (`i$dd_age_sex_source == "Mean"`)
* The individual patients (if based on individual patients)

The result is utilities for each population (pop_0, pop_1, pop_2) for each cycle.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `p$util$gpop`

In the base case, it produces age-adjusted utilities, with age from the patient-level data.


::: {.cell}

```{.r .cell-code}
i$dd_ageadjuutilities
```

::: {.cell-output .cell-output-stdout}

```
[1] "Yes"
```


:::

```{.r .cell-code}
i$dd_age_sex_source
```

::: {.cell-output .cell-output-stdout}

```
[1] "PLD"
```


:::
:::


The resulting utilities from, for example, population 0 are:


::: {.cell}

```{.r .cell-code}
head(p$util$gpop$pop_0, 20)
```

::: {.cell-output .cell-output-stdout}

```
 [1] 1.0000000 0.9999149 0.9998297 0.9997446 0.9996594 0.9995742 0.9994889
 [8] 0.9994037 0.9993184 0.9992331 0.9991477 0.9990623 0.9989769 0.9988915
[15] 0.9988061 0.9987206 0.9986351 0.9985496 0.9984640 0.9983784
```


:::
:::


:::

The function `f_process_utilities()` simply extracts and renames some of the columns from `i$R_table_util` to create `i$QALYs$utilities$means`. However, if `PSA = TRUE`, it would generate multiple samples of utility values for use in probabilistic sensitivity analysis.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

Original table of utilities:


::: {.cell}

```{.r .cell-code}
kable(head(i$R_table_util))
```

::: {.cell-output-display}


| Population|Population.name | Treatment.line| Molecule|Treatment.name                | OnTxtMean| OffTxtMean| PFSMean| PDMean| OnTxtSE| OffTxtSE|  PFSSE|   PDSE|
|----------:|:---------------|--------------:|--------:|:-----------------------------|---------:|----------:|-------:|------:|-------:|--------:|------:|------:|
|          0|All             |              1|        0|nivolumab_monotherapy         |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|
|          0|All             |              1|        1|cabozantinib_plus_nivolumab   |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|
|          0|All             |              1|        2|nivolumab_plus_ipilimumab     |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|
|          0|All             |              1|        3|lenvatinib_plus_pembrolizumab |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|
|          0|All             |              1|        4|avelumab_plus_axitinib        |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|
|          0|All             |              1|        5|pazopanib                     |     0.753|      0.753|   0.753|  0.683|  0.0753|   0.0753| 0.0753| 0.0683|


:::
:::


Table `f_process_utilities()`:


::: {.cell}

```{.r .cell-code}
kable(head(i$QALYs$utilities$means))
```

::: {.cell-output-display}


| Population|Population.name | Treatment.line| Molecule|Treatment.name                | OnTxt| OffTxt|   PFS|    PD|
|----------:|:---------------|--------------:|--------:|:-----------------------------|-----:|------:|-----:|-----:|
|          0|All             |              1|        0|nivolumab_monotherapy         | 0.753|  0.753| 0.753| 0.683|
|          0|All             |              1|        1|cabozantinib_plus_nivolumab   | 0.753|  0.753| 0.753| 0.683|
|          0|All             |              1|        2|nivolumab_plus_ipilimumab     | 0.753|  0.753| 0.753| 0.683|
|          0|All             |              1|        3|lenvatinib_plus_pembrolizumab | 0.753|  0.753| 0.753| 0.683|
|          0|All             |              1|        4|avelumab_plus_axitinib        | 0.753|  0.753| 0.753| 0.683|
|          0|All             |              1|        5|pazopanib                     | 0.753|  0.753| 0.753| 0.683|


:::
:::


:::

## AEs

This section processes data on adverse events. There are three outputs of this:

1. **AE type/approach** (`p$ae$aetype` and `p$ae$approach`) - which is defined by a single string (same value for each of those variables)
2. **AE duration** (`p$ae$duration`) - a table which, for particular molecules and treatment lines, the rate and duration of specific adverse events is listed
3. **AE impact** (`p$ae$mk$per_cycle`) - a table which, for particular molecules and treatment lines, the cost and QALY impact of adverse events is specified


::: {.cell}

```{.r .cell-code}
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
```

::: {.cell-output .cell-output-stdout}

```
Mean cost and QALY impact due to AEs per patient per week:
```


:::

```{.r .cell-code}
# they all match up with column RCC_input_desc from Excel except BSC (mol_999). 
# Set BSC (assumed to have 0 AE impact):
p$ae$mk$per_cycle[trt == "BSC",]$trt <- i$lookup$ipd$mol[match(999,Number),]$RCC_input_desc

# Convert to numbers. Now it's ready for use in the patient flow sheet.
p$ae$mk$per_cycle$molecule <- i$lookup$ipd$mol[match(p$ae$mk$per_cycle$trt,RCC_input_desc),]$Number

# Add in the AE approach switch
p$ae$approach <- i$dd_apply_AE_options
```
:::


::: {.callout-note collapse="true"}

## View the AE type/approach, and the duration and impact tables

AE type/approach:


::: {.cell}

```{.r .cell-code}
p$ae$aetype  # same as p$ae$approach
```

::: {.cell-output .cell-output-stdout}

```
[1] "one-off"
```


:::
:::


AE rate and duration:


::: {.cell}

```{.r .cell-code}
kable(head(p$ae$duration))
```

::: {.cell-output-display}


|Treatment.name        | Molecule| Treatment.line|AE                         | Rate.per.patient.per.week| duration_weeks|
|:---------------------|--------:|--------------:|:--------------------------|-------------------------:|--------------:|
|nivolumab_monotherapy |        0|              2|G3+ ALT increased          |                 0.0000000|       23.91518|
|nivolumab_monotherapy |        0|              2|G3+ Anaemia                |                 0.0028248|       23.91518|
|nivolumab_monotherapy |        0|              2|G3+ Decreased appetite     |                 0.0000000|       23.91518|
|nivolumab_monotherapy |        0|              2|G3+ Diarrhoea _TKI induced |                 0.0004556|       23.91518|
|nivolumab_monotherapy |        0|              2|G3+ Diarrhoea _IO induced  |                 0.0000000|       23.91518|
|nivolumab_monotherapy |        0|              2|G3+ Fatigue                |                 0.0010023|       23.91518|


:::
:::


AE cost and QALY impact:


::: {.cell}

```{.r .cell-code}
kable(head(p$ae$mk$per_cycle))
```

::: {.cell-output-display}


|trt                           | line|      cost|      QALYs| molecule|
|:-----------------------------|----:|---------:|----------:|--------:|
|nivolumab_monotherapy         |    1|  6.737719| -0.0000668|        0|
|cabozantinib_plus_nivolumab   |    1| 11.887276| -0.0001200|        1|
|nivolumab_plus_ipilimumab     |    1|  9.745281| -0.0000884|        2|
|lenvatinib_plus_pembrolizumab |    1| 14.363607| -0.0001422|        3|
|avelumab_plus_axitinib        |    1| 17.905744| -0.0001771|        4|
|pazopanib                     |    1| 14.707425| -0.0001450|        5|


:::
:::


:::

## Costs

This section processes cost data:

1. Costs of **drug administration and medical resource use** (MRU) (`p$costs$mk`)
2. **One-off** costs like end of life care and subsequent treatment (`p$costs$oneoff`)

There is some code within an `if (FALSE) {}` statement, meaning it is never run, but that would generate values for the probablistic sensitivity analysis (`i$cost$drug_and_admin_cost_by_tunnel_state$PSA`).


::: {.cell}

```{.r .cell-code}
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
```

::: {.cell-output .cell-output-stderr}

```
Warning in FUN(X[[i]], ...): NAs introduced by coercion
```


:::

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View the costs

`p$costs$mk` contains the costs per cycle for each molecule and line (hence length 2089, as 2089 weeks = 40 years). It has options of:

* Drug
* Admin
* MRU on
* MRU off


::: {.cell}

```{.r .cell-code}
# Length of costs
length(p$costs$mk$drug$line_1$mol_4)
```

::: {.cell-output .cell-output-stdout}

```
[1] 2089
```


:::

```{.r .cell-code}
# Example
head(p$costs$mk$drug$line_1$mol_4, 20)
```

::: {.cell-output .cell-output-stdout}

```
 [1] 5627.810    0.000 2666.496    0.000 5627.810    0.000 2666.496    0.000
 [9] 5627.810    0.000 2666.496    0.000 5627.810    0.000 2666.496    0.000
[17] 5627.810    0.000 2666.496    0.000
```


:::
:::


Whilst `p$costs$oneoff` and `p$costs$oneoff_mk` contain one-off costs.


::: {.cell}

```{.r .cell-code}
kable(p$costs$oneoff)
```

::: {.cell-output-display}


|Type.of.cost                                          |Type    |Apply.to |     cost|
|:-----------------------------------------------------|:-------|:--------|--------:|
|End of life (EOL) or terminal care 
(applied at death) |one-off |Death    | 8714.034|
|Subsequent radiotherapy (following progression)       |one-off |Prog     |  511.020|
|Subsequent surgery (following progression)            |one-off |Prog     | 5393.256|


:::

```{.r .cell-code}
p$costs$oneoff_mk
```

::: {.cell-output .cell-output-stdout}

```
$Death
[1] 8714.034

$Prog
[1] 5904.276
```


:::
:::


:::

## Subsequent treatments


::: {.cell}

```{.r .cell-code}
# 3.4.4 Subsequent treatment -------------------------------------------------------------

# Read in cost and QALY consequences for subsequent treatments per first line option from Excel
# This information is only used whaen the PartSA model structure is selected

p$substrt$partsa  <- as.data.table(f_process_subs_txt_data(
  subs_txt = i$R_table_sub_txt_cum_costs,
  PSA             = FALSE
))
```
:::


::: {.callout-note collapse="true"}

## View `p$substrt$parts`

This table contains information on the costs and AEs associated with subsequent treatment, as specified for each treatment and population.


::: {.cell}

```{.r .cell-code}
kable(head(p$substrt$partsa))
```

::: {.cell-output-display}


|Treatment                     |Population | drug_cost| admin_cost|  AE_cost| AE_QALY_impact|
|:-----------------------------|:----------|---------:|----------:|--------:|--------------:|
|nivolumab_monotherapy         |pop0       |      0.00|     0.0000|   0.0000|      0.0000000|
|cabozantinib_plus_nivolumab   |pop0       |  39268.59|   795.5377| 707.1494|     -0.0063590|
|nivolumab_plus_ipilimumab     |pop0       |      0.00|     0.0000|   0.0000|      0.0000000|
|lenvatinib_plus_pembrolizumab |pop0       |      0.00|     0.0000|   0.0000|      0.0000000|
|avelumab_plus_axitinib        |pop0       |  39608.96|   703.3352| 555.6036|     -0.0049325|
|pazopanib                     |pop0       |  54145.22|  4320.7039| 787.5639|     -0.0071626|


:::
:::


:::

## Population mappings

There are:

* Three risk populations (all risk, favourable risk, intermediate/poor risk)
* Two immuno-oncology (IO) populations (whether or not patients had prior adjuvant therapy with IO treatment within 12 months)

These leads to six possible combinations overall. We can use this table to map from those six, to those used in the analysis:

* Treatment sequences are sorted into 4 populations (`Sequencing.population.number`)
* Survival extrapolations, costs, QALYs and AEs are sorted into 3 populations (`Risk.population.number`)


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `p$basic$lookup$pop_map`


::: {.cell}

```{.r .cell-code}
kable(p$basic$lookup$pop_map)
```

::: {.cell-output-display}


|Sequencing.population                          |Risk.population | Sequencing.population.number| Risk.population.number| Overall.population.number|
|:----------------------------------------------|:---------------|----------------------------:|----------------------:|-------------------------:|
|All risk and favourable risk no prior adjuvant |All risk        |                            0|                      0|                         1|
|All risk and favourable risk no prior adjuvant |Favourable risk |                            0|                      2|                         2|
|Int / poor risk no prior adjuvant              |Int/poor        |                            1|                      1|                         3|
|All risk and favourable risk prior adjuvant    |All risk        |                            2|                      0|                         4|
|All risk and favourable risk prior adjuvant    |Favourable risk |                            2|                      2|                         5|
|Int / poor risk prior adjuvant                 |Int/poor        |                            3|                      1|                         6|


:::
:::


:::

