# Functions for RCC model

# created 08/02/2023 by Ed Wilson

#' Function to generate treatment sequences through generating all permutations of comparators
#' 
#' 
#' @param comparators A vector of comparator names. Note that the 
#' @param maxlines The maximum number of treatment lines (such that each row of permutations output has that many columns)  
#' 
#' 
generate_sequences <- function(comparators, maxlines = 5) {
  
  strategies <- list()
  
  # The maximum number of living states is the number of lines (BSC is always the last line) 
  # The number of model health states is number of lines +2 to include death too
  max_living_states <- maxlines 
  
  # List out all the possible treatment pathways by line, entering blank values
  # when looking at lines past the first. So for 4th line patients there's only their
  # treatment plus BSC to follow.
  # 
  # lapply makes the list per line, rbindlist binds data.tables efficiently
  # 
  rbindlist(lapply(1:(maxlines -  1), function(line) {
    
    # For this line, the number of remaining living states is the difference between
    # line and the max of living states:
    living_states_left <- max_living_states-line
    
    # use gtools' function to efficiently generate all possible permutations
    # of treatments given the possible comparators, given starting from this
    # line of therapy (i.e. so that at later lines there are less columns.)
    out <- data.table(permutations(n = length(comparators), r = living_states_left, v = comparators))
    
    # Adding in BSC as the rightmost column (data.table naming convention is VN
    # where N is number. take the number of rows, and if it's less than the max
    # then add empty columns to ensure there are maxlines+1 columnns at the end
    out[[paste0("V",dim(out)[2]+1)]] <- "BSC"
    nco <- dim(out)[2]
    
    # Make sure the output always has columns equal to the number of treatment lines
    if (nco < max_living_states) {
      ncol_to_add <- max_living_states - nco
      for (col_to_add in 1:ncol_to_add) {
        out[[paste0("V",nco+col_to_add)]] <- ""
      }
    }
    
    return(out)
  }))
}







#' Utility function to apply a filter sequentially to selected columns of a data.table
#' based on a text string. There is a data.table way to do this using .SDCols but it's
#' hard to read, so this is more transparent. Assumes the columns are named V per
#' data.table defaults.
#' 
#' @param dat_tab a `data.table` containing the sequenceis
#' @param lines numeric vector of the treatment lines to remove the treatment from
#' @param trt character entry for the treatment to remove from the designated lines
#' 
f_path_RemoveFromLines <- function(dat_tab, lines, trt) {
  if (length(lines) == 1) {
    return(dat_tab[eval(as.name(paste0("V",lines))) != trt,])
  } else {
    Reduce(
      x          = paste0("V",lines),
      accumulate = FALSE,
      init       = dat_tab,
      f = function(tab_last_round, this_col) {
        return(tab_last_round[eval(as.name(this_col)) != trt,])
      }
    )
  }
  
}




#' Utiltiy function to apply the "only after" ruleset to a treatment using subset
#' filtering.
#' 
#' Note that this is a reverse filter (note the ! in front of the conditions). Handy
#' for only returning those results which don't match the conditions. in this case
#' it's essentially like performing a filter for those conditions.
#' 
#' It's very efficient and easy to read. 
#' 
#' @description At each line, the treatments in rule aren't allowed after trt. filter out.
#' 
#' @param perms permutatins of possible treatments
#' @param trt the treatment to apply this filtering alg to
#' @param rule a character vector of those treatments that aren't allowed after trt
#' 
f_path_notAllowedAfter <- function(perms, trt, rule) {
  perms <- subset(perms,!(V1 == trt & (V2 %in% rule | V3 %in% rule | V4 %in% rule | V5 %in% rule)))
  perms <- subset(perms,!(V2 == trt & (V3 %in% rule | V4 %in% rule | V5 %in% rule)))
  perms <- subset(perms,!(V3 == trt & (V4 %in% rule | V5 %in% rule)))
  perms <- subset(perms,!(V4 == trt & (V5 %in% rule)))
  return(perms)
}


#' Function as above but the opposite - only allowing the rule treatments before
#' the given treatment (trt). Filters out starting from 2L as it's pointless to 
#' filter first line based on previous treatments!
f_path_onlyAllowedAfter <- function(perms, trt, rule) {
  perms <- subset(perms,!(V2 == trt & !V1 %in% rule))
  perms <- subset(perms,!(V3 == trt & (!V1 %in% rule & !V2 %in% rule)))
  perms <- subset(perms,!(V4 == trt & (!V1 %in% rule & !V2 %in% rule & !V3 %in% rule)))
  perms <- subset(perms,!(V5 == trt & (!V1 %in% rule & !V2 %in% rule & !V3 %in% rule & !V4 %in% rule)))
  return(perms)
}

#' another inversion of the above functions. this one doesn't allow the rule treatments
#' to come before trt, such that if you say no vegf treatment before trt it will
#' filter out prior vegf
f_path_notAllowedBefore <- function(perms, trt, rule) {
  perms <- subset(perms,!(V2 == trt & V1 %in% rule))
  perms <- subset(perms,!(V3 == trt & (V1 %in% rule | V2 %in% rule)))
  perms <- subset(perms,!(V4 == trt & (V1 %in% rule | V2 %in% rule | V3 %in% rule)))
  perms <- subset(perms,!(V5 == trt & (V1 %in% rule | V2 %in% rule | V3 %in% rule | V4 %in% rule)))
  return(perms)
}


# List of the treatments that other treatments are only allowed subsequently to
# take these from excel so they're always up to date!

apply_tx_restrictions <- function(sequences) {
  #Rules as per Dawn's spreadsheet "Treatment Sequences"
  # move to 3_functions once written
  
  # 1. Ave_axi 1L tx with advanced RCC
  for (i in 2:ncol(sequences)) {
    sequences <- sequences[sequences[,i]!="ave_axi",]
  }
  
  # 2. Axitinib only 2L+ and after failure of sunitinib or a cytokine
  # remove first line axitinib
  sequences <- sequences[sequences[,1]!="axitinib",]
  onlyafters <- c("sunitinib") #note - what are the cytokine drugs?  Need to add
  
  sequences[,ncol(sequences)+1] <- TRUE
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no axitinib
    if (!("axitinib" %in% sequences[i,])) next
    #remove sequences with axitinib but no of the 'onlyafters' drugs
    if (!(onlyafters %in% sequences[i,])) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    #remove sequences where onlyafters are before axitinib
    for (j in 1:length(onlyafters)) {
      if (match(onlyafters[j], sequences[i,]) > match("axitinib", sequences[i,])){
        sequences[i,ncol(sequences)] <- FALSE
      } 
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] #drop violating sequences
  sequences[,ncol(sequences)] <- NULL
  rm(onlyafters)
  
  # 3. Belzutifan only after PDI (ipi_nivo, nivo_cabo, nivolumab, pem_len) 1L
  #    and VEGF TKI (ave_axi, axitinib, cabozantinib, len_evro, nivo_cabo, pazopinib, pem_len, sunitinib, tivozanib)
  onlyafters1 <- c("ipi_nivo", "nivo_cabo", "nivolumab", "pem_len")
  onlyafters2 <- c("ave_axi", "axitinib", "cabozantinib", "len_evro",
                   "nivo_cabo", "pazopinib", "pem_len", "sunitinib", "tivozanib")
  
  #remove any 1L belzutifan
  sequences <- sequences[sequences[,1]!="belzutifan",]
  
  sequences[,ncol(sequences)+1] <- TRUE
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no belzutifan
    if (!("belzutifan" %in% sequences[i,])) next
    
    #remove sequences without PD1/PD-L1 inhibitor 
    if (!(sequences[i,1] %in% onlyafters1)) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences with belzutifan but no of the 'onlyafters' drugs
    if (sum(onlyafters2 %in% sequences[i,])==0) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences without VEGF TKI before belzutifan
    if (sum(sequences[i,] %in% onlyafters2) > 0) {
      if(match("belzutifan", sequences[i,]) < min(match(onlyafters2,sequences[i,], nomatch = 99))) {
        sequences[i,ncol(sequences)] <- FALSE
        next
      }
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] #drop violating sequences
  sequences[,ncol(sequences)] <- NULL
  rm(onlyafters1, onlyafters2)
  
  # 4. cabozantinib monotherapy for advanced renal cell carcinoma
  #- as first-line treatment of adult patients with intermediate or poor risk 
  #- in adults following prior vascular endothelial growth factor (VEGF)-targeted therapy "
  onlyafters <- c("ave_axi", "axitinib", "cabozantinib",
                  "len_evro", "nivo_cabo", "pazopinib", 
                  "pem_len", "sunitinib", "tivozanib")
  sequences[,ncol(sequences)+1] <- TRUE
  
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no cabozantinib
    if (!("cabozantinib" %in% sequences[i,])) next
    
    #keep sequences with cabozantinib first line
    if (sequences[i,1] == "cabozantinib") next
    
    #remove sequences without VEGF
    if (sum(sequences[i,] %in% onlyafters) == 0) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences with cabozantinib before VEGF
    if(match("cabozantinib", sequences[i,]) < min(match(onlyafters,sequences[i,], nomatch = 99))) {
      sequences[i,ncol(sequences)] <- FALSE
      print(i)
      next
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] 
  sequences[,ncol(sequences)] <- NULL
  
  #Note - rule 4 doesn't drop any sequences - CHECK
  
  
  # 5. Evrolimus only after VEGF-targeted therapy
  onlyafters <- c("ave_axi", "axitinib", "cabozantinib", "len_evro",
                  "nivo_cabo", "pazopinib", "pem_len", "sunitinib", "tivozanib")
  
  #remove any 1L evrolimus
  sequences <- sequences[sequences[,1]!="evrolimus",]
  
  sequences[,ncol(sequences)+1] <- TRUE
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no evrolimus
    if (!("evrolimus" %in% sequences[i,])) next
    
    #remove sequences with evrolimus but no of the 'onlyafters' drugs
    if (sum(onlyafters %in% sequences[i,])==0) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences with evrolimus before VEGF 
    if (sum(sequences[i,] %in% onlyafters) > 0) {
      if(match("evrolimus", sequences[i,]) < min(match(onlyafters,sequences[i,], nomatch = 99))) {
        sequences[i,ncol(sequences)] <- FALSE
        next
      }
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] #drop violating sequences
  sequences[,ncol(sequences)] <- NULL
  rm(onlyafters)
  
  # 6. ipi_nivo only first line
  for (i in 2:ncol(sequences)) {
    sequences <- sequences[!(sequences[,i] == "ipi_nivo"),]
  }
  
  # 7. len_evro following 1 prior VEGF-targeted therapy
  onlyafters <- c("ave_axi", "axitinib", "cabozantinib", "len_evro",
                  "nivo_cabo", "pazopinib", "pem_len", "sunitinib", "tivozanib")
  
  #remove any 1L len_evro
  sequences <- sequences[sequences[,1]!="len_evro",]
  
  sequences[,ncol(sequences)+1] <- TRUE
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no len_evro
    if (!("len_evro" %in% sequences[i,])) next
    
    #remove sequences with len_evro but none of the 'onlyafters' drugs
    if (sum(onlyafters %in% sequences[i,])==0) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences with len_evro before VEGF and more than 1 VEGF 
    if (sum(sequences[i,] %in% onlyafters) > 0) {
      #reject if len_evro is before a VEGF
      if (match("len_evro", sequences[i,]) < min(match(onlyafters,sequences[i,], nomatch = 99))) {
        sequences[i,ncol(sequences)] <- FALSE
      }
      #reject if nr VEGF therapies prior to len_evro > 1
      if (sum(match(onlyafters,sequences[i,], nomatch = 99) < match("len_evro", sequences[i,])) > 1) {
        sequences[i,ncol(sequences)] <- FALSE
      }
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] #drop violating sequences
  sequences[,ncol(sequences)] <- NULL
  rm(onlyafters)
  
  # 8. nivo_cabo first line only
  for (i in 2:ncol(sequences)) {
    sequences <- sequences[!(sequences[,i] == "nivo_cabo"),]
  }
  
  # 9. nivolumab monotherapy only 2L+
  sequences <- sequences[!(sequences[,1] == "nivolumab"),]
  
  # 10. pazopinib 1L abd for patients who have received prior cytokine therapy
  # (CHECK - what is cytokine therapy?)
  sequences <- sequences[!(sequences[,1] == "nivolumab"),]
  
  # 11. pem_len 1L only
  for (i in 2:ncol(sequences)) {
    sequences <- sequences[!(sequences[,i] == "pem_len"),]
  }
  
  # 12. sunitinib no restrictions
  
  # 13. tivozanib 1L and no prior VEG-F and mTOR
  notafters <- c("ave_axi", "axitinib", "cabozantinib", "evrolimus",
                 "len_evro", "nivo_cabo", "pazopinib", 
                 "pem_len", "sunitinib", "tivozanib")
  sequences[,ncol(sequences)+1] <- TRUE
  
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no tivozanib
    if (!("tivozanib" %in% sequences[i,])) next
    
    #keep sequences with tivozanib first line
    if (sequences[i,1] == "tivozanib") next
    
    #remove sequences with tivozanib before VEGF or mTOR
    if(match("tivozanib", sequences[i,]) > min(match(notafters,sequences[i,], nomatch = 99))) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] 
  sequences[,ncol(sequences)] <- NULL
  
  # 14. cannot use >1 IO
  IO <- c("ave_axi", "ipi_nivo", "nivo_cabo", "nivolumab", "pem_len") #Note - CHECK (which are IOs?)
  sequences[,ncol(sequences)+1] <- TRUE
  
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no one or less IO
    if (sum(IO %in% sequences[i,]) < 2) next
    
    sequences[i,ncol(sequences)] <- FALSE
  }
  sequences <- sequences[sequences[,ncol(sequences)],] 
  sequences[,ncol(sequences)] <- NULL  
  
  # 15. tivozanib, pazopinib, sunitinib can only be used after ipi_nivo
  drugs <- c("tivozanib", "pazopinib", "sunitinib")
  sequences[,ncol(sequences)+1] <- TRUE
  for (i in 1:nrow(sequences)) {
    #ignore sequences with no tivozanib, pazopinib, sunitinib
    if (sum(drugs %in% sequences[i,]) == 0) next
    
    #remove sequences with tivozanib, pazopinib, sunitinib and no ipi_nivo  
    if (!("ipi_nivo" %in% sequences[i,])) {
      sequences[i,ncol(sequences)] <- FALSE
      next
    }
    
    #remove sequences with tivozanib, pazopinib, sunitinib before ipi_nivo
    if(match("ipi_nivo", sequences[i,]) > min(match(drugs,sequences[i,], nomatch = 99))) {
      sequences[i,ncol(sequences)] <- FALSE
      print(i)
      next
    }
  }
  sequences <- sequences[sequences[,ncol(sequences)],] 
  sequences[,ncol(sequences)] <- NULL
}


f_path_tx_restrict <- function(sequences, subs_to, prev_to, one_in_class) {
  
  s <- sequences
  
  #  ave_axi is only permitted in first line therapy. therefore, filter out
  #  all subsequent uses of this treatment. Determine the columns and then apply the filter.
  #  Use the function f_path_RemoveFromLines to remove "ave_axi" from 2L+
  #  
  #  Showing arguments once here, then will be a one-liner from here:
  
  s <- f_path_RemoveFromLines(
    dat_tab = s,
    lines   = 2:ncol(s),
    trt     = "ave_axi"
  )
  
  # Next, axitinib 2L+
  s <- f_path_RemoveFromLines(s,1,"axitinib")
  
  # filter out unacceptable treatments previous to axitinib at all lines 2L+
  # sub_axi <- subs_to$axitinib
  s <- f_path_notAllowedAfter(perms = s,trt = "axitinib",rule = subs_to$axitinib)
  s <- f_path_onlyAllowedAfter(perms = s,trt = "axitinib",rule = subs_to$axitinib)
  
  # To QC this function e.g. sub_axi %in% unique(unlist(unique(s[V2 == "axitinib",list(V3,V4,V5,V6)]),use.names = F))
  
  # cabo is allowed 1L as mono, OR 2L+ for those that have had a VEGF before
  
  s <- f_path_notAllowedAfter(s,"cabozantinib",subs_to$cabozantinib)
  s <- f_path_onlyAllowedAfter(s,"cabozantinib",subs_to$cabozantinib)
  
  # evero post vegf, 2L+
  s <- f_path_RemoveFromLines(s,1,"evero")
  s <- f_path_notAllowedAfter( s,"everolimus",subs_to$everolimus)
  s <- f_path_onlyAllowedAfter(s,"everolimus",subs_to$everolimus)
  
  # ipi_nivo only 1L
  s <- f_path_RemoveFromLines(s,2:ncol(s),"ipi_nivo")
  
  # len_evero is only allowed after 1 prior VEGF
  # I do not think the current functions handle this - to discuss
  s <- f_path_RemoveFromLines(s,1,"len_evero")
  s <- f_path_notAllowedAfter( s,"len_evero",subs_to$len_evero)
  s <- f_path_onlyAllowedAfter(s,"len_evero",subs_to$len_evero)
  
  # nivo_cabo 1L
  s <- f_path_RemoveFromLines(s,2:ncol(s),"nivo_cabo")
  
  # nivo mono only 2L+
  s <- f_path_RemoveFromLines(s,1,"nivolumab")
  s <- f_path_notAllowedAfter( s,"len_evero",subs_to$nivolumab)
  s <- f_path_onlyAllowedAfter(s,"len_evero",subs_to$nivolumab)
  
# pazo is 1L or no prior VEGF
  # check if this handles or correctly
  # also applies to sunitinib and tivozanib
  
  s <- f_path_notAllowedBefore( s,"tivozanib",subs_to$tivozanib)
  s <- f_path_onlyAllowedAfter(s,"tivozanib",subs_to$tivozanib)
 
  # pem_len 1L
  s <- f_path_RemoveFromLines(s,2:ncol(s),"pem_len")
  
  # sunitinib is 1L or no prior VEGF
  
  s <- f_path_notAllowedBefore( s,"sunitinib",subs_to$sunitinib)
  s <- f_path_onlyAllowedAfter(s,"sunitinib",subs_to$sunitinib)
  
  # tivozanib 1L or no prior VEGF
  s <- f_path_notAllowedBefore( s,"pazopinib",subs_to$pazopinib)
  s <- f_path_onlyAllowedAfter(s,"pazopinib",subs_to$pazopinib)
  
  # Cannot repeat treat with IOs. for each io, impose not allowed after and not
  # allowed before, for the other io's in the vector. this should make it impossible
  # for a double treat to slip through
  s <- Reduce(
    x = one_in_class$io,
    accumulate = FALSE,
    init = s,
    f = function(prev, io) {
      not_allowed_after_these <- one_in_class$io[which(one_in_class$io != io)]
      prev <- f_path_notAllowedAfter(prev,io,not_allowed_after_these)
      out  <- f_path_notAllowedBefore(prev,io,not_allowed_after_these)
    }
  )
    
  # return the possible sequences now that treatment rules have been applied
  return(s)
}

