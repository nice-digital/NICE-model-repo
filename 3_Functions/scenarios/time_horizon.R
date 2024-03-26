
#' function to update i$surv$reg so that the extrapolations in there for
#' deterministic analysis are reduced down to the time horizon
#' 
#' @param regs `i$surv$reg` from model, after loading in or running them
#' @param TH time horizon `+1` cycles
#' 
f_scen_TH_impose_surv <- function(regs, TH) {
  lapply(regs, function(popu) {
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
}