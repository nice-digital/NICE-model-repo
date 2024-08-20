
# This script is an example of how the patient flow sheet calculations will take place
# There will be one for each treatment pathway (as these are the states of the world
# to be compared to each other). The result will be one 2-dimensional numeric matrix
# per treatment pathway. Selected values from the colSums() of this matrix will be
# the results that feed into results calculations, and also all the necessary breakdowns
# of costs, LYs, QALYs, doses received and so on, so this matrix should be quite 
# exhaustive and should think about the future (i.e. what tables will be needed,
# yes, but also what tables will be needed in the future too)

# A THROWAWAY EXAMPLE of ONE p - one iteration of one scenario of the parameters
# list. obviously this will all come from Excel extraction
require(data.table)
p <- list(
  basic = list(
    TH_yr      = 40,
    cl_yr      = 28 / 365.25,
    cl_mo      = 1  / 12,
    cl_we      = 28 / 7,
    cl_da      = 28,
    disc_c     = 0.035,
    disc_q     = 0.035,
    TH         = ceiling((365.25/28) * 40),
    t_yr       = (seq(0,ceiling((365.25/28) * 40)) * (28 / 365.25))[1:ceiling((365.25/28) * 40)],
    drug_names = c("drug 1", "drug 2")
  ),
  path = list(
    pathways = 1:2,
    labels   = c("Drug 1 only", "Drug 1 followed by Drug 2")
  ),
  drug = list(
    rel_by_pathway = list("drug 1", c("drug 1", "drug 2")),
    dosing_by_drug       = list(
      list(
        name  = "drug 1",
        form  = "IV",
        basis = "mg/kg",
        dose  = 5,
        dur   = NULL,
        titr  = NULL,
        wean  = NULL
      ),
      list(
        name  = "drug 2",
        form  = "oral",
        basis = NULL,
        dose  = 3,
        dur   = 26,
        titr  = list(
          c(dose = 1, dur = 4),
          c(dose = 2, dur = 8)
        ),
        wean  = NULL
      )
    ),
    final_cost = list(
      50,
      100
    )
  )
)

names(p$basic$drug_names)    <- p$basic$drug_names
names(p$drug$dosing_by_drug) <- p$basic$drug_names
names(p$drug$final_cost)     <- p$basic$drug_names

# Discount factors. Discounting is done after 
p$basic$dfac_c <- 1/((1+p$basic$disc_c) ^ p$basic$t_yr)
p$basic$dfac_q <- 1/((1+p$basic$disc_q) ^ p$basic$t_yr)



# make a list of results r to house...the results

r <- list(
  det  = list(),
  psa  = list(),
  scen = list(),
  owsa = list(),
  evppi = list()
)

# FOR AN EXAMPLE - calculate the patient flow sheet of the deterministic base-case
#                  analysis using the input list p. This list should contain EVERYTHING
#                  needed to run the results. DO NOT refer to things outside of 
#                  p whilst inside this Reduce() call. This is the most central
#                  aspect of the CE model and MUST be as transparent as possible.

r$det$pf <- lapply(p$path$pathways, function(pathway) {
  
  
  # Make the top row of the patient flow sheet for this treatment pathway. This
  # varies only in some of the top row values. Drug costs and such will have
  # different top row values per relevant drugs per line etc. In this simplified
  # example, we've got drug 1 and drug 2. only drug 1 in pathway 1, drug 1 and
  # drug 2 in pathway 2.
  
  top_row <- list(
    cycle  = 0,
    t_yr   = 0,
    t_mo   = 0,
    t_we   = 0,
    t_da   = 0,
    dfac_c = 1,
    dfac_q = 1
  )
  
  # Drugs that are relevant in this treatment pathway
  dr_nam         <- p$basic$drug_names
  
  # retrieve a vector of the relevant drugs for this treatment pathway
  relevant_drugs <- p$drug$rel_by_pathway[[pathway]]
  
  # To check which pathway we're currently in:
  # which_pathway <- p$path$labels[pathway]
  
  # retrieve the dosing information for the relevant drugs for this treatment pathway
  # Note one square bracket allows subsetting of a list for multiple entries
  dosing_info <- p$drug$dosing_by_drug[relevant_drugs]
  
  # Add in 0s for ALL POSSIBLE DRUGS. This then ensures that all patient flow sheets
  # for all treatment pathways have the same columns, allowing one function
  # to summarise them all.
  top_row[paste0("d_cost_undisc_",p$basic$drug_names)] <- 0
  top_row[paste0("d_cost_disc_"  ,p$basic$drug_names)] <- 0
  
  # now, only for the relevant drugs, add in their drug cost in the first cycle:
  top_row[paste0("d_cost_undisc_",relevant_drugs)] <- unlist(p$drug$final_cost[relevant_drugs])
  
  pf_cycle_list <- Reduce(
    accumulate = TRUE,
    x = 1:p$basic$TH,
    init = as.list(top_row),
    f = function(previous_cycle, cycle_number) {
      
      # calculate the values for the next cycle using
      #  - cycle_number: the current cycle number - note x goes from 1:TH not 0:TH as top_row is cycle 0
      #  - previous_cycle: the result of reduce for the previous cycle
      #  - p: a complete parameters list for this scenario of the model.
      #   - This should contain ALL inputs required to calculate the patient flow
      #   - p is specific to this iteration for this structural scenario
      #   - p should be substituted for e.g. PSA[[1]] or OWSA$ub[[this_parameter]] to run new PF sheets
      #   - NOTHING AT ALL from outside should be used here, except for functions.
      #   - All elements in top_row are named so that it is clear what data is used
      #   - All elements in top_row AND previous_cycle must have length 1
      #  - functions to calculate this_cycle given p, previous_cycle, and cycle_number
      
      # start from the values from the previous cycle
      prev <- as.list(previous_cycle)
      this <- prev
      
      # Calculate values for this cycle - start from the basics:
      this$cycle <- prev$cycle + 1
      this$t_yr  <- prev$t_yr  + p$basic$cl_yr
      this$t_mo  <- prev$t_mo  + p$basic$cl_mo
      this$t_we  <- prev$t_we  + p$basic$cl_we
      this$t_da  <- prev$t_da  + p$basic$cl_da
      
      # We can use dosing_info as an argument to a function whcih calculates
      # consumption of each drug in our cycle length period given the point in
      # time and the potential rules, RDIs and so on.
      
      # assign the relevant drug costs to the relevant columns
      this[paste0("d_cost_undisc_",relevant_drugs)] <- p$drug$final_cost[relevant_drugs]
      
      # Note we do not perform discounting here as it's inefficient. instead
      # take the discount factors e.g., p$basic$dfac_c and apply them at the end!
      
      # note, when using reduce you can refer to objects in previous_cycle, which
      # gets around the biggest pitfall of lapply vs for, which is referring to
      # previous iterations.
      
      # ... you get the idea, repeat this stuff for every column in the PFS sheet
      
      # once the cylce is calculated, unlist the results. this puts it back into
      # a named vector, so that at the end, they can all be rbound together into
      # one numeric matrix!
      as.data.table(this)
    }
  )
  
  # This sticks them together row wise and then separates them column wise!
  res <- as.list(rbindlist(pf_cycle_list)[1:p$basic$TH])
  
  res$dfac_c <- p$basic$dfac_c
  res$dfac_q <- p$basic$dfac_q
  
  res[paste0("d_cost_disc_",p$basic$drug_names)]   <- lapply(res[paste0("d_cost_undisc_",p$basic$drug_names)], function(x) x * res$dfac_c)
  
  # Before returning this treatment pathway, bind the results together row wise, or
  # get a huge computational gain by defining it as a matrix. This results
  # in a 2d matrix for each treatment pathway which is not restricted by size, labelling,
  # number of treatment pathways, approaches to calculating columns, underlying settings
  # in p, whether the analysis is probabilistic or deterministic, whether it's a
  # scenario analysis or a base-case analysis.
  
  return(as.data.table(res))
  
})

r$det$pf
