

#' Function f_NestedApply: a function to impose the nesting structure required
#'                         for the model. This function makes sure that
#'                         within the structure, an arbitrary set of code can
#'                         be automatically applied to it without the need
#'                         to repeatedly type out all of the nesting. This
#'                         keeps all of the bits of data (which may differ in
#'                         terms of shape, size, class, themselves being nested,
#'                         a list or even a regression object) separate in their
#'                         own named spaces, whilst also allowing the convenience
#'                         of universally applying a function (one could feasibly 
#'                         add conditional logic into the arbitrary code here
#'                         to allow additional flexibility)
#'
#'@param mylist a list which strictly follows the set list structure. if it differs this will error!
#'@param f any R function, OR a series of lines of code wrapped in "{}" just like a lapply 
#'
f_NestedApply <- function(mylist,f) {
  lapply(mylist, function(this_pop) {
    lapply(this_pop, function(this_line) {
      lapply(this_line, function(this_mol) {
        lapply(this_mol, function(this_trial) {
          lapply(this_trial, function(this_endpoint) {
            f(this_endpoint)
          })
        })
      })
    })
  })
}


# list_surv_survivaldata in the model is defined early and uses this structure.
# to do something like see what the unique values are in certain columns
# of this dataset:


# # Testing that each of he identifier columns for the dataset matchess the place
# # it is in the list:
# f_NestedApply(list_surv_survivaldata, f = function(obj) {
#   unlist(sapply(obj,unique)[c("population", "line", "molecule", "trial", "endpoint")])
# })
# 
# # maximum follow up:
# f_NestedApply(list_surv_survivaldata, f = function(obj) max(obj$timew))
# 
# # Earliest recorded non-censor event
# f_NestedApply(list_surv_survivaldata, f = function(obj) {
#   c(
#     event_censor_0 = min(obj[event_censor == 0,]$timew),
#     event_censor_1 = min(obj[event_censor == 1,]$timew)
#   )
# })
# 
# # Basic survival analysis
# all_TSD14 <- f_NestedApply(list_surv_survivaldata, f = function(obj) {
#   if (nrow(obj) > 0) {
#     # Do survival analysis if the data has any rows!
#     
#     lapply (c(
#       gengamma  = "gengamma",
#       exp       = "exp",
#       weibull   = "weibull",
#       lnorm     = "lnorm",
#       gamma     = "gamma",
#       gompertz  = "gompertz",
#       llogis    = "llogis",
#       lognormal = "lognormal"
#     ), function(distr) {
#       reg <-
#         flexsurvreg(
#           formula = Surv(timew, event_censor) ~ 1,
#           data = obj,
#           dist = "gengamma"
#         )
#       
#       return(list(
#         coef = coefficients(reg),
#         vcov = vcov(reg),
#         fit = c(AIC = AIC(reg), BIC = BIC(reg), ll = logLik(reg))
#       ))
#     })
#   }
# })
# 


#' Version of `f_NestedApply` with definable levels. Allows applying on PLMTE level
#' 
#' @param mylist list following PLMTE nesting format
#' @param f a function to apply like sum, length, head etc. can use `function(x) ...` format if useful
#' @param levels an indicator of levels
#' 
function_NestedApply <- function(mylist,f,levels = 5) {
  lapply(mylist, function(this_pop) {
    if(levels>1) {lapply(this_pop, function(this_line) {
      if(levels>2){lapply(this_line, function(this_mol) {
        if(levels>3){lapply(this_mol, function(this_trial) {
          if(levels>4){lapply(this_trial, function(this_endpoint) {
            f(this_endpoint)
          })}
        })}
      })}
    })}
  })
}
