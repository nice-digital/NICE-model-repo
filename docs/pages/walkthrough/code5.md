
This page runs the **economic model**.

## Final pre-processing steps

A table with some settings is add to `p`.


::: {.cell}

```{.r .cell-code}
# 3.6 PATIENT FLOW ------------------------------------------------------
# Now that we have all of the disease evidence and modeled outcomes available,
# we can compute the disease model. In the deterministic case, this will simply 
# be a trace for both the PS and Markov models. 

i$R_table_eff_data_settings <- data.table(i$R_table_eff_data_settings)
p$releff$excel_table <- i$R_table_eff_data_settings
```
:::


::: {.callout-note collapse="true"}

## View `p$releff$excel_table`


::: {.cell}

```{.r .cell-code}
kable(head(p$releff$excel_table))
```

::: {.cell-output-display}


| Population|Population.name        | Treatment.line| Molecule|Treatment.name                |Include.in.this.analysis. | End.point|End.point.name   |Effectiveness.data.source |Trial.name.if.effectiveness.source.is.trial | Origin.population| Origin.line| Origin.treatment| Origin.trial| Origin.endpoint|Curve.fit..for.survival.analysis. | HR.to.apply| HR.95..CI..LCL.| HR.95..CI..UCL.| Trial|
|----------:|:----------------------|--------------:|--------:|:-----------------------------|:-------------------------|---------:|:----------------|:-------------------------|:-------------------------------------------|-----------------:|-----------:|----------------:|------------:|---------------:|:---------------------------------|-----------:|---------------:|---------------:|-----:|
|          1|Poor/intermediate risk |              1|        0|nivolumab_monotherapy         |No                        |         0|Overall survival |0                         |0                                           |                 1|           1|                7|            2|               0|Exponential                       |       0.000|             0.0|           0.000|    NA|
|          1|Poor/intermediate risk |              1|        1|cabozantinib_plus_nivolumab   |Yes                       |         0|Overall survival |FP_NMA                    |0                                           |                 1|           1|                7|            2|               0|Exponential                       |       0.000|             0.0|           0.000|    NA|
|          1|Poor/intermediate risk |              1|        2|nivolumab_plus_ipilimumab     |Yes                       |         0|Overall survival |FP_NMA                    |0                                           |                 1|           1|                7|            2|               0|Exponential                       |       0.000|             0.0|           0.000|    NA|
|          1|Poor/intermediate risk |              1|        3|lenvatinib_plus_pembrolizumab |Yes                       |         0|Overall survival |Apply HR to               |0                                           |                 1|           1|                1|            2|               0|NA                                |       1.134|             0.8|           1.619|    NA|
|          1|Poor/intermediate risk |              1|        4|avelumab_plus_axitinib        |No                        |         0|Overall survival |PH_NMA                    |0                                           |                 1|           1|                7|            2|               0|Exponential                       |       0.000|             0.0|           0.000|    NA|
|          1|Poor/intermediate risk |              1|        5|pazopanib                     |Yes                       |         0|Overall survival |Assume equal to           |0                                           |                 1|           1|                7|            2|               0|Exponential                       |       0.000|             0.0|           0.000|    NA|


:::
:::


:::

The script also confirms that the specified model structure is a valid option (from `i$dd_model_struct`, and again a bit later from `p$basic$structure`).

::: {.callout-note collapse="true"}

## View default model structure


::: {.cell}

```{.r .cell-code}
i$dd_model_struct
```

::: {.cell-output .cell-output-stdout}

```
[1] "State transition"
```


:::

```{.r .cell-code}
p$basic$structure
```

::: {.cell-output .cell-output-stdout}

```
[1] "State transition"
```


:::
:::


:::

The comments in the code provide some description of the modelling, which is a survival analysis using either a state transition (markov) model or partitioned state survival model.


::: {.cell}

```{.r .cell-code}
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
```
:::


This code then checks whether the decision problem is "cabozantinib plus nivolumab" - in which case, the first three populations are relevant, as the final three populations can't receive this treatment.

::: {.callout-note collapse="true"}

## View default decision problem


::: {.cell}

```{.r .cell-code}
p$basic$decision_problem
```

::: {.cell-output .cell-output-stdout}

```
[1] "cabo+nivo"
```


:::
:::


:::


::: {.cell}

```{.r .cell-code}
# Check the decision problem. If it's cabo+nivo only the first 3 overall
# populations are relevant as one cannot get cabo+nivo in pops 4-6
if(p$basic$decision_problem == "cabo+nivo") {
  p$basic$pops_to_run <- 1:3
} else {
  p$basic$pops_to_run <- NULL
}
```
:::


The objects `p` and `i`, which contain the parameters used to run the model, can be save as `.rds` files to make it easier to re-run the model.


::: {.cell}

```{.r .cell-code}
# populate the pf object irrespective of model structure or overall populations
# to include.
# 
# Note that if you put n_cores as NULL or 1 then the model will run in sequence.

# For QC, you can easily just save i and p as a file so all you need to do is
# load libraries to repeat the results and QC things:

saveRDS(p, file.path(s_path,"standalone scripts/QC/p.rds"))
saveRDS(i, file.path(s_path, "standalone scripts/QC/i.rds"))
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
```
:::


## Run the model

The survival analysis run using `f_pf_computePF()`. This function can take a while to run. For the purposes of this walkthrough, the code has been modified using `if(FALSE){}` to prevent the function from running and, instead, a set of pre-run model results are loaded.

::: {.callout-note collapse="true"}

## Further info about pre-run model results

The model result file is very large and, as such, cannot be easily synced with GitHub. Therefore, in order to run this script on a new machine (as described in the walkthrough README), you will need to run `Model_Structure.R` (or the walkthrough `.qmd` file) to generate the file.

:::

The function `f_pf_compute_PF()` will use either:

* `f_pf_computePF_mk()` (if state transition model, i.e. Markov model)
* `f_pf_computePF_ps()` (if partitioned survival model)

The purpose of this is to model how patients move through different health states/treatments, from which we can find outputs like costs and utilities.

**State transition model**. Steps include:

* Getting relevant population data and creating `sequence_list` with the allowed treatment sequences for a given population
* Creating `util_gpop_mat` - a matrix of health state utility values that account for ageing
* Transition probabilities (`TPs`) are calculated for each sequence, then used to construct Markov traces (`TRACE`) (ie. matrix of state vectors across all model cycles)
* These are used to calculate costs through the treatment sequence (`pf_costs`) (e.g. cost per cycle `cpc`, costs on initiation `coi`, and costs on death `cod`), QALYs `pf_qalys`) and adverse events (`pf_ae`)

**Partitioned survival model**. Steps include:

* Getting relevant population data and creating the utility matrix `util_gpop_mat`
* Using `f_pf_partSA_state_res()`, simulates patient states (e.g. PFS) across time horizon for each treatment sequence, returns the costs (`costs`) and QALYs (`qalys`)


::: {.cell}

```{.r .cell-code}
# Make this NA to run single core:
if(FALSE){
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
}

# if (is.na(keep_free_cores)) {
#   plan(sequential)
# } else {
#   plan(multisession(workers = max(availableCores()-keep_free_cores,1)))
# }

# Save result
# computepf_file = file.path(d_path, "computepf_example.rds")
# saveRDS(pf, file=computepf_file)
# Load result
computepf_file = file.path(d_path, "computepf_example.rds")
pf <- readRDS(computepf_file)
```
:::


## Analyse the model results

The model results are compiled using `f_res_compute_results()` which will then use either:

* `f_res_compute_results_mk()` for the state transition model
* `f_res_compute_results_ps()` for the partitioned survival model

These each in turn utilities various other functions to summarise the model results. There are various functions, depending on the level of detail desired in the plots. This is set to either 4 or 5, depending on whether it is scenario 0.

The summary functions include:

* `f_pf_mk_summary()` to get `undisc` and `disc` results (undiscounted and discounted)
* `f_res_mk_incremental()` to get `incremental`
* `f_res_wa_model()` to get `weighted_model` discounted or undiscounted, summarised using `f_res_sum_weighted_model()`
* `f_pf_wa_sr_plots()` to get `weighted_trace_plots`

And more! As summarised in code comments below, the outcomes of this from each detail level are:

1. Top line results by sequence and weighted average results (costs, QALYs, LYs)
2. Include incremental analysis (from `f_res_mk_incremental()`)
3. Include weighted average pairwise comparisons
4. Include incremental analysis by sequence
5. Include trace plots

Except for the partitioned survival, 4 is breakdown tables and 4 is state residency plots.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## View `res`

A few examples of results from `res` are shown below...


::: {.cell}

```{.r .cell-code}
kable(head(res$undisc$pop_1$ly$breakdown, 3))
```

::: {.cell-output-display}


|trt_n |trt                                       |    L1_on|    L1_off|       BSC| L2_on| L2_off| L3_on| L3_off| L4_on| L4_off|
|:-----|:-----------------------------------------|--------:|---------:|---------:|-----:|------:|-----:|------:|-----:|------:|
|1→999 |Cabozantinib plus nivolumab→Placebo / BSC | 1.945363| 0.1088617| 0.3625033|    NA|     NA|    NA|     NA|    NA|     NA|
|5→999 |Pazopanib→Placebo / BSC                   | 1.144445| 0.1148029| 0.3673995|    NA|     NA|    NA|     NA|    NA|     NA|
|7→999 |Sunitinib→Placebo / BSC                   | 1.144445| 0.1148029| 0.3673995|    NA|     NA|    NA|     NA|    NA|     NA|


:::

```{.r .cell-code}
kable(head(res$incremental$pop_1$expanded_results, 3))
```

::: {.cell-output-display}


|trt_n    |trt                                |    costs|    qalys|       ly|L1                 |str_dom |extdom |       ic|        iq|       il|     ICER|  r|
|:--------|:----------------------------------|--------:|--------:|--------:|:------------------|:-------|:------|--------:|---------:|--------:|--------:|--:|
|7→999    |Sunitinib→Placebo / BSC            | 22826.06| 1.085297| 1.626648|7&#124;999         |FALSE   |FALSE  |    0.000| 0.0000000| 0.000000|     0.00|  1|
|5→999    |Pazopanib→Placebo / BSC            | 24690.10| 1.085678| 1.626648|5&#124;999         |FALSE   |TRUE   |       NA|        NA|       NA|       NA|  2|
|7→10→999 |Sunitinib→Everolimus→Placebo / BSC | 31562.18| 1.399302| 2.169939|7&#124;10&#124;999 |FALSE   |FALSE  | 8736.118| 0.3140051| 0.543291| 27821.58|  3|


:::

```{.r .cell-code}
res$weighted_trace_plots$pop_1$plots$L1_1
```

::: {.cell-output-display}
![](code0a_walkthrough_files/figure-html/unnamed-chunk-109-1.png){width=672}
:::

```{.r .cell-code}
res$wa_summarised$pop_1$costs
```

::: {.cell-output .cell-output-stdout}

```
[1] 268538.83  80398.95 100004.76  77049.92
```


:::

```{.r .cell-code}
res$pairwise_vs_mol$pop_1$qalys
```

::: {.cell-output .cell-output-stdout}

```
[1] 2.253814 1.712196 1.674692 1.684174
```


:::
:::


:::

## Severity modifier

A severity modifier is implemented for the state transition model (and not for the partitioned survival model). These are used to weight QALYs for conditions with a high degree of severity.

The result of this section is `res$mk$qaly_shortfall_1_to_3` which is a list of three lists, each containing the QALY adjustments for each population.


::: {.cell}

```{.r .cell-code}
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
```
:::


::: {.callout-note collapse="true"}

## How is the severity modifier calculated?

There are first a few processing steps before the function to calculate severity, `calc_severity_modifier()`. Looking at an example of `npa_pop=1`, we start with getting the look-up for the population and the population mappings.


::: {.cell}

```{.r .cell-code}
npa_pop = 1

lu_pop <- p$basic$lookup$pop_map
lu_rpop <- p$basic$lookup$ipd$pop

lu_pop
```

::: {.cell-output-display}

:::

```{.r .cell-code}
lu_rpop
```

::: {.cell-output-display}

:::
:::


The risk population can then be identified - in this case "0" or "All risk".


::: {.cell}

```{.r .cell-code}
# npa_pop is overall population, we need to look up risk population from it:

risk_pop_n <- lu_pop[match(npa_pop,lu_pop$Overall.population.number),]$Risk.population.number
risk_pop <- lu_rpop[match(risk_pop_n,lu_rpop$Number),]$Description  

risk_pop_n
```

::: {.cell-output .cell-output-stdout}

```
[1] 0
```


:::

```{.r .cell-code}
risk_pop
```

::: {.cell-output .cell-output-stdout}

```
[1] "All"
```


:::
:::


The patient characteristics table from excel is converted into a data table (if not already).


::: {.cell}

```{.r .cell-code}
i$R_table_ptchar <- as.data.table(i$R_table_ptchar)
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


With our source of IPD (as from `i$dd_age_sex_source`), this code then created lists with the ages and genders of each person in the IPD.


::: {.cell}

```{.r .cell-code}
patient_sex_age_IPD <- as.data.table(i$R_table_patientagesex)
patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="M","male")
patient_sex_age_IPD$Gender <- replace(patient_sex_age_IPD$Gender, patient_sex_age_IPD$Gender=="F","female")

bl_age <- patient_sex_age_IPD[Line ==1]$Age 
bl_male <- patient_sex_age_IPD[Line ==1]$Gender

head(bl_age)
```

::: {.cell-output .cell-output-stdout}

```
[1] 46.36414 74.09583 55.18070 69.14784 78.15880 79.32854
```


:::

```{.r .cell-code}
head(bl_male)
```

::: {.cell-output .cell-output-stdout}

```
[1] "male"   "male"   "female" "male"   "male"   "male"  
```


:::
:::


The name of the population is created based on `npa_pop`.


::: {.cell}

```{.r .cell-code}
pna_txt <- names(res$wa_summarised)[npa_pop]
pna_txt
```

::: {.cell-output .cell-output-stdout}

```
[1] "pop_1"
```


:::
:::


The economic model results for that population are identified for treatments except L1 treatment of molecule 1.


::: {.cell}

```{.r .cell-code}
tab <- res$wa_summarised[[pna_txt]][L1 != 1,]
kable(tab)
```

::: {.cell-output-display}


| L1|     costs|    qalys|       ly|
|--:|---------:|--------:|--------:|
|  5|  80398.95| 1.712196| 2.836778|
|  6| 100004.76| 1.674692| 2.764879|
|  7|  77049.92| 1.684174| 2.780869|


:::
:::

From that, the treatment with the maximum QALYs is identified, and the QALYs and treatment were extracted


::: {.cell}

```{.r .cell-code}
met <- tab[which.max(qalys),]
q_met <- met$qalys
comp_no_met <- met$L1

kable(met)
```

::: {.cell-output-display}


| L1|    costs|    qalys|       ly|
|--:|--------:|--------:|--------:|
|  5| 80398.95| 1.712196| 2.836778|


:::
:::


These results were then used in the function `calc_severity_modifier()`. That function uses `get_severity_modifier()` to find the appropriate severity modifier according to the discounted QALYs for those with the condition and those without.

This was repeated for the three populations, generating three results lists, which have been combined into the table below:


::: {.cell}

```{.r .cell-code}
# The three lists of results, combined into a single table
kable(bind_rows(res$mk$qaly_shortfall_1_to_3))
```

::: {.cell-output-display}


| qaly_soc| qaly_gpop|   abs_sf|   prop_sf| modifier| SOC|
|--------:|---------:|--------:|---------:|--------:|---:|
| 1.712196|  10.38197| 8.669771| 0.8350798|        1|   5|
| 2.243521|  10.38197| 8.138445| 0.7839021|        1|   5|
| 2.013113|  10.38197| 8.368854| 0.8060952|        1|   3|


:::
:::


:::

## Save the results

This code chunks saves the `res` list as an `.rds` file, although that has been disabled for the purpose of this walkthrough.


::: {.cell}

```{.r .cell-code}
# 3.6.4 Saving the results ------------------------------------------------------

# the results are in list format so should be saved as an R list file (.rds)
# The naming of the file should reflect the model structure used. The file
# produced has a time stamp to avoid overwriting previous results files.

if(FALSE){
  Scenario_name <- i$R_Scenario_name    # Use ST for state transition, PS for Partitioned survival, LP for list price, cPAS for cPAS
  Scenario_number <- i$R_Scenario_num
  Run_date <- date()
  if (p$basic$structure == "State transition") {
    saveRDS(res, file.path(o_path, paste0(
      "ST_Scenario ", Scenario_number,"_",i$dd_drug_price_options,gsub(":","_",Run_date),".rds")))
  } else {
    saveRDS(res, file.path(o_path, paste0(
      "PartSA_Scenario ",Scenario_number,"_",i$dd_drug_price_options,gsub(":","_",Run_date),".rds")))
  }
}
```
:::


## Creating the report

The function `f_res_ProduceWordDoc()` produces a word document report. Depending on `p$basic$structure`, will either use `f_res_ProduceWordDoc_ST()` (st - state transition) or `f_res_ProduceWordDoc_PS()` (ps - partitioned survival) to produce the report content. This applies results tables specific to the cabo+nivo population to the word document, in the format required by NICE.

This report is explored in detail on the [Report](../report.qmd) page.


::: {.cell}

```{.r .cell-code}
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

Run_date <- date()
f_res_ProduceWordDoc(
  p                      = p,
  res                    = res,
  Scenario_name          = i$R_Scenario_name,
  Scenario_number        = i$R_Scenario_num,
  price_options          = i$dd_drug_price_options,
  Run_date               = Run_date,
  word_template_location = file.path(f_path, "reporting/empty results doc.docx"),
  Word_width_inches      = 29.7*0.3937,
  auto_save              = FALSE,
  verbose = TRUE
)
```

::: {.cell-output .cell-output-stdout}

```
[0;32mWord output for Base case. (scen #0): [0m
[0;31mWord output for Base case. (scen #0): LY breakdown[0m
[0;33mWord output for Base case. (scen #0): QALY breakdown[0m
[0;34mWord output for Base case. (scen #0): Cost breakdown[0m
[0;35mWord output for Base case. (scen #0): Scenario tables[0m
[0;36mWord output for Base case. (scen #0): Severity modifier[0m
[0;37mWord output for Base case. (scen #0): Pairwise results[0m
[0;40mWord output for Base case. (scen #0): Generating word document from results...[0m
$L1_1
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_5
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_7
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_6
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_1
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_5
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_7
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_6
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_8
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_1
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_3
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_2
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_5
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_7
```


:::

::: {.cell-output .cell-output-stdout}

```

$L1_6
```


:::


::: {.cell-output .cell-output-stdout}

```
rdocx document with 147 element(s)

* styles:
                         Normal                       heading 1 
                    "paragraph"                     "paragraph" 
                      heading 2                       heading 3 
                    "paragraph"                     "paragraph" 
                      heading 4                       heading 5 
                    "paragraph"                     "paragraph" 
                      heading 6                       heading 7 
                    "paragraph"                     "paragraph" 
                      heading 8                       heading 9 
                    "paragraph"                     "paragraph" 
         Default Paragraph Font                    Normal Table 
                    "character"                         "table" 
                        No List            annotation reference 
                    "numbering"                     "character" 
                annotation text               Comment Text Char 
                    "paragraph"                     "character" 
             annotation subject            Comment Subject Char 
                    "paragraph"                     "character" 
                         Strong                    Normal (Web) 
                    "character"                     "paragraph" 
                 List Paragraph                      Table Text 
                    "paragraph"                     "paragraph" 
               Table Text Char1                NICEReport-table 
                    "character"                         "table" 
                 Heading 3 Char              Numbered heading 1 
                    "character"                     "paragraph" 
                Bullet indent 1            Bullet indent 1 Char 
                    "paragraph"                     "character" 
           Bullet indent 1 last              Numbered heading 2 
                    "paragraph"                     "paragraph" 
        Numbered heading 2 Char           Numbered level 3 text 
                    "character"                     "paragraph" 
                 Heading 1 Char                  Heading 2 Char 
                    "character"                     "character" 
                        Mention                   footnote text 
                    "character"                     "paragraph" 
             Footnote Text Char             ACIC-Clear-MainText 
                    "character"                     "character" 
           ACIC-Clear-TableText                             AIC 
                    "character"                     "character" 
                      app:head1                       app:head2 
                    "paragraph"                     "paragraph" 
                      app:head3                    Balloon Text 
                    "paragraph"                     "paragraph" 
              Balloon Text Char                       list:bull 
                    "character"                     "paragraph" 
                 list:bull Char                            blue 
                    "character"                     "paragraph" 
           Borderless (Relaxed)                   Bullet left 1 
                        "table"                     "paragraph" 
             Bullet left 1 Char              Bullet left 1 last 
                    "character"                     "paragraph" 
        Bullet left 1 last Char                   Bullet left 2 
                    "character"                     "paragraph" 
             Bullet left 2 last         Bullet left 2 last Char 
                    "paragraph"                     "character" 
                  Bullet left 3                Bullet list left 
                    "paragraph"                     "numbering" 
                        caption                    Caption Char 
                    "paragraph"                     "character" 
                     centhead12                            cf01 
                    "paragraph"                     "character" 
                            CIC             c-mrkdwn__highlight 
                    "character"                     "character" 
                        Default                    Document Map 
                    "paragraph"                     "paragraph" 
              Document Map Char            EndNote Bibliography 
                    "character"                     "paragraph" 
      EndNote Bibliography Char      EndNote Bibliography Title 
                    "character"                     "paragraph" 
EndNote Bibliography Title Char                             eop 
                    "character"                     "character" 
                        findhit               FollowedHyperlink 
                    "character"                     "character" 
                         footer                     Footer Char 
                    "paragraph"                     "character" 
             footnote reference              Grid Table 1 Light 
                    "character"                         "table" 
    Grid Table 1 Light Accent 3                    Grid Table 2 
                        "table"                         "table" 
                   Grid Table 3                    Grid Table 4 
                        "table"                         "table" 
                         header                     Header Char 
                    "paragraph"                     "character" 
                 Heading 4 Char                  Heading 5 Char 
                    "character"                     "character" 
                 Heading 6 Char                  Heading 7 Char 
                    "character"                     "character" 
                 Heading 8 Char                  Heading 9 Char 
                    "character"                     "character" 
                    Hidden Text                Hidden Text Char 
                    "paragraph"                     "character" 
                       HTAtable                       Hyperlink 
                        "table"                     "character" 
                      left head                      lefthead12 
                    "paragraph"                     "paragraph" 
                      lh:NonTOC                     lh:NonTOC12 
                    "paragraph"                     "paragraph" 
               lh:NonTOC12 Char                     list:indent 
                    "character"                     "paragraph" 
               list:indent Char                list:indent bull 
                    "character"                     "paragraph" 
          list:indent bull Char                        list:num 
                    "character"                     "paragraph" 
                  list:num Char                        Mention1 
                    "character"                     "character" 
            NICE figure caption        NICE figure caption Char 
                    "paragraph"                     "character" 
                    NICE normal                NICE normal Char 
                    "paragraph"                     "character" 
                    NoNum:Head1                     NoNum:Head2 
                    "paragraph"                     "paragraph" 
                    NoNum:Head3                NoNum:Head3 Char 
                    "paragraph"                     "character" 
                    NoNum:Head4                NoNum:Head4 Char 
                    "paragraph"                     "character" 
                    NoNum:Head5                   normaltextrun 
                    "paragraph"                     "character" 
                    page number                       paragraph 
                    "character"                     "paragraph" 
               Placeholder Text                  Plain Table 21 
                    "character"                         "table" 
                    quote:block                     quote:short 
                    "paragraph"                     "character" 
                       Subtitle                   Subtitle Char 
                    "paragraph"                     "character" 
                    superscript                      Table text 
                    "character"                     "paragraph" 
                Table text Char                  Table bullet 1 
                    "character"                     "paragraph" 
            Table bullet 1 Char                  Table bullet 2 
                    "character"                     "paragraph" 
            Table bullet 2 Char                  Table footnote 
                    "character"                     "paragraph" 
            Table footnote Char                      Table Grid 
                    "character"                         "table" 
           table of authorities                table of figures 
                    "paragraph"                     "paragraph" 
                     table:bull                           Title 
                    "paragraph"                     "paragraph" 
                     Title Char                           toc 1 
                    "character"                     "paragraph" 
                          toc 2                           toc 3 
                    "paragraph"                     "paragraph" 
                          toc 4                           toc 5 
                    "paragraph"                     "paragraph" 
                          toc 6                           toc 7 
                    "paragraph"                     "paragraph" 
                          toc 8                           toc 9 
                    "paragraph"                     "paragraph" 
                    TOC Heading                      TOC_Header 
                    "paragraph"                     "paragraph" 
                       TOC_Page              Unresolved Mention 
                    "paragraph"                     "character" 
            Unresolved Mention1             Unresolved Mention2 
                    "character"                     "character" 
            Unresolved Mention3             Unresolved Mention4 
                    "character"                     "character" 
            Unresolved Mention5  List Table 7 Colorful Accent 6 
                    "character"                         "table" 
                       Revision                    Style1lh:TOC 
                    "paragraph"                     "paragraph" 
                         lh:TOC                       msonormal 
                    "paragraph"                     "paragraph" 
                           xl65                            xl66 
                    "paragraph"                     "paragraph" 
                           xl67                            xl68 
                    "paragraph"                     "paragraph" 
                           xl69                            xl70 
                    "paragraph"                     "paragraph" 
                           xl71                            xl72 
                    "paragraph"                     "paragraph" 
                           xl73                            xl74 
                    "paragraph"                     "paragraph" 
                           xl75                            xl76 
                    "paragraph"                     "paragraph" 
                           xl77                     gnd-iwgdh3b 
                    "paragraph"                     "character" 
              HTML Preformatted          HTML Preformatted Char 
                    "paragraph"                     "character" 

* Content at cursor location:
  level num_id text style_name content_type
1    NA     NA              NA    paragraph
```


:::
:::


## Future plans

The final section of code outlines some of the future plans for this model.


::: {.cell}

```{.r .cell-code}
# END OF CODE -------------------------------------------------------


#### Additional changes were originally planned during Phase 2 of this pilot following use for the initial decision problem including
# - Addition of Shiny user interface
# - Genericisation of the code to allow wider use
# - Programming and analysis of model outputs related specifically to sequencing, this may include value of information analyses

# Unfortunately funding for this has not been confirmed currently. 
# If you are interested in discussing or funding further development please contact the PenTAG team at pentag@exeter.ac.uk
```
:::


