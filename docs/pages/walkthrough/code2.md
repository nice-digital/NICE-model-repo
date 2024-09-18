
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


