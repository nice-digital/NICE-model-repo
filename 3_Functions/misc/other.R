
#' Test the lengths of multiple objects are the same
f_misc_all_same_length = function(...){
  length(unique(lengths(list(...)))) == 1
}


#' Very useful function to extract a PLMTE from a nested list. Used within
#' the relative efficacy and throughout the model.
#' 
#' @param dat nested list following the PLMTE format (population line molecule trial endpoint)
#' @param l list containing location information. MUST have pop line mol trial endpoint entries
#' @tr_or allows alternative entry for `l$trial` in case it's needed to look at another trial's PLM&E
#' 
#' 
f_misc_get_plmte <- function(dat,l, tr_or = NULL) {
  list2env(l,envir = environment())
  if (!is.null(tr_or)) {
    dat[[pop]][[line]][[mol]][[tr_or]][[endpoint]]
  } else {
    dat[[pop]][[line]][[mol]][[trial]][[endpoint]]
  }
}


#' fun function to print to console with some colour
f_misc_colcat <- function(txt,col_num = 32) {
  cat(
    paste0(
      "\033[0;",
      col_num,
      "m",
      txt,
      "\033[0m",
      "\n"
    )
  )
}