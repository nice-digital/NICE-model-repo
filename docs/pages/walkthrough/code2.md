
This page performs a **survival analysis** on the patient-level data, which uses **real-world evidence (RWE)** in the **base case**. The analysis is performed in order to **extrapolate the survival curves** so they cover the full time horizon of the economic model (40 years).

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


