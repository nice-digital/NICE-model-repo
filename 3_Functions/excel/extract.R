# This script contains the functions required to extract data from excel.
# The primary function f_excel_extract extracts all named ranges from an excel file,
# and puts all of those objects into an object of the same name in a flat list.

#' Function to automatically extract all named ranges from an excel workbook
#' 
#' 
#' @param path_to_excel_file the path to the excel file. Ensure that this is normalized for the operating system you are using
#' @param verbose if FALSE (the default), no printing will occur. If TRUE the name and file path of each named range will be printed to trace errors.
#' 
f_excel_extract <- function(path_to_excel_file, verbose = FALSE) {
  
  # require(openxlsx)
  
  # Note that pre-loading the workbook to an object vastly improves computational
  # efficiency, because it's in RAM and not on disk. Reading from disk repeatedly
  # is extremely slow and even a small task like this can take time unnecessarily.
  # Load workbook, get the named regions from it ready for lapply.
  
  wb          <- loadWorkbook(path_to_excel_file)
  sheets      <- wb$sheet_names
  named       <- getNamedRegions(wb)
  nams        <- unlist(as.list(named))
  names(nams) <- nams
  
  # Cycle through all the named ranges, making a decision on cleaning as we go
  output <- lapply(nams, function(named_range) {
    
    if (verbose) cat(paste0("Extracting named range ", named_range, " from ", path_to_excel_file, "\n"))
    
    # It is impossible to presume the inclusion or exclusion of headers without
    # further information. This is because if the headers are included in the columns
    # the columns become character data, and tables may contain text in rows
    # so we can't try to convert to numbers to identify whether the header
    # row is erroneously included in the named range. Therefore we assume
    # that all tables for R DO include the column names, and these are converted
    # to R-safe column names (i.e. replacing special characters with dots)
    
    dat <- read.xlsx(wb,namedRegion = named_range, check.names = TRUE,colNames = TRUE)
    
    if (nrow(dat) == 0) {
      # if it has 0 rows, then the input must be a 1 cell one, so we can safely assume
      # no row names
      dat <- read.xlsx(wb,namedRegion = named_range,colNames = FALSE,rowNames = FALSE)
      
      if (all(length(dat) == 1, names(dat)[1] == "X1")) {
        dat <- unlist(dat, use.names = FALSE)
      } else if(nrow(dat) == 1) {
        dat <- unlist(read.xlsx(wb,namedRegion = named_range,colNames = FALSE,rowNames = FALSE), use.names = F)
      } 
      
      return(dat)
    } else if(ncol(dat) == 1) {
      # if there's 1 column, then it's actually just a vector:
      dat <- unlist(read.xlsx(wb,namedRegion = named_range,colNames = FALSE,rowNames = FALSE), use.names = F)
    } else {
      # logically we can't identify whether a table has appropriate names or not
      # at this point. this will have to be left to the cleaning functions.
      return(dat)
    }
  })
  
  # Drop the workbook at the end (I think in {} this would happen anyway, but want to ensure!)
  rm(wb)
  return(output)
}



#' Function to take the table "R_table_param" from the "PATT RCC_model inputs.xlsx workbook 
#' and expand it out into a list of lists which can then be expanded during cleaning (for
#' things such as adding in a vcov and such)
#' 
#' @param tab tbe named range "R_table_param" from "PATT RCC_model inputs.xlsx"
#' 
#' @examples 
#' 
#' i <- f_excel_extract(path_to_excel_file, verbose = TRUE)
#' f_excel_cleanParams(i$R_table_param)
#' 
f_excel_cleanParams <- function(tab, verbose = FALSE) {
  require(data.table)
  tab <- as.data.table(tab)
  params <- tab$Parameter.name
  names(params) <- params
  
  list_for_i <- lapply(params, function(param) {
    
    if (verbose) cat(paste0("Excel extraction - tidying parameters table: ",param,"\n"))
    
    # for each param, grab the row as a list so we can add to it with things
    # that don't have to be 1x1
    out_list <- as.list(tab[Parameter.name == param,])
    
    # there are several columns that should be numeric but aren't due to other
    # row values being text
    suppressWarnings(
      if (!is.na(as.numeric(out_list$Mean))) out_list$Mean <- as.numeric(out_list$Mean)
    )
    
    return(out_list)
  })
  
  return(list_for_i)
}
