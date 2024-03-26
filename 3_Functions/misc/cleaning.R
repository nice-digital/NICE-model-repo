# Script dedicated to functions which clean up inputs from excel for entry
# into p (or to be used to generate e.g. probabilsitic draws, bounds and so on).



# Patient characteristics table -------------------------------------------


#' Function to clean the patient characteristic table up into a neat list for
#' entry into p.
#' 
#' @param R_table_ptchar named range "R_table_ptchar" from excel inputs
#' @param lookups i$lookup as generated in the model_structure script.
#' 
f_cleaning_ptchar <- function(R_table_ptchar, lookups) {
  tab <- data.table(i$R_table_ptchar)
  
  # Translate population into numbers:
  tab$Population <- lookups$ipd$pop$Number[match(tab$Population,lookups$ipd$pop$Description)]
  
  u_pop  <- structure(unique(tab$Population),.Names = paste0("pop_",unique(tab$Population)))
  u_line <- structure(unique(tab$Treatment.line),.Names = paste0("line_",unique(tab$Treatment.line)))
  
  # produce a list by population and line which presents cleaned inputs for
  # each input. This can then go into p for a deterministic case or be used to
  # generate probabilistic draws.
  lapply(u_pop, function(popu) {
    lapply(u_line, function(line) {
      ta <- tab[Population == popu & Treatment.line == line,]
      if (nrow(ta) == 0) {
        ta <- tab[Population == 0 & Treatment.line == line,]
      }
      
      # Make a tidy list:
      list(
        age    = list(
          mean = ta$Starting.age..years..Mean,
          se = ta$Starting.age..years..SE,
          n = ta$Starting.age..years....10.n
        ),
        pr_fem = list(
          mean = ta$Starting...female.Mean,
          se = ta$Starting...female...10.SE,
          n = ta$Starting...female...10.n
        ),
        weight = list(
          mean = ta$Body.weight..kg...10.Mean,
          se = ta$Body.weight..kg...10.SE,
          n = ta$Body.weight..kg...10.n
        ),
        prior_io = list(
          mean = ta$Prior.IO...in.12.months..10.Mean,
          se = ta$Prior.IO...in.12.months..10.SE,
          se = ta$Prior.IO...in.12.months..10.n
        ),
        pr_i_rsk = list(
          mean = ta$Starting...PorI.risk.Mean,
          se = ta$Starting...PorI.risk.SE,
          n = ta$Starting...PorI.risk...10.n
        )
      )
      
    })
  })
}

