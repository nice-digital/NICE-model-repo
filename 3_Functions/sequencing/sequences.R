# Sequencing functions

#' Function to generate treatment sequences through generating all permutations of comparators
#' 
#' 
#' @param comparators A vector of comparator names. Note that the 
#' @param maxlines The maximum number of treatment lines (such that each row of permutations output has that many columns)  
#' 
#' 
f_generate_sequences <- function(comparators, maxlines = 4) {
  
  strategies <- list()
  
  # The maximum number of living states is the number of lines +1 (for BSC). 
  # The number of model health states is number of lines +2 to include death too
  max_living_states <- maxlines + 1
  
  # List out all the possible treatment pathways by line, entering blank values
  # when looking at lines past the first. So for 5th line patients there's only their
  # treatment plus BSC to follow.
  # 
  # lapply makes the list per line, rbindlist binds data.tables efficiently
  # 
  rbindlist(lapply(1:maxlines, function(line) {
    
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

# Functions to implement sequencing rules based on lists of what treatments are allowed after each other
f_get_only_after_lists <- function(i, population) {
  output <- i[grep("only_after", names(i))]
  output <- output[!(grepl("only_after_one", names(output)))]
  output <- output[!(grepl("_2L_only_after", names(output)))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_only_after", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_not_immediate_after_lists <- function(i, population) {
  output <- i[grep("not_immediate_after", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_not_immediate_after", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_one_in_list_lists <- function(i, population) {
  output <- i[grep("one_allowed", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_only_one_allowed", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_only_after_one_lists <- function(i, population) {
  output <- i[grep("only_after_one", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_only_after_one", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_2L_only_after_lists <- function(i, population) {
  output <- i[grep("2L_only_after", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_2L_only_after", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_2L_only_immediate_after_lists <- function(i, population) {
  output <- i[grep("_2L_only_immediate_after", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_2L_only_immediate_after", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_2L_only_one_lists <- function(i, population) {
  output <- i[grep("_2L_only_one", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  names(output) <- sub("_2L_only_one", "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_allowed_lists <- function(i, population) {
  output <- i[grep("allowed", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  output <- output[!(grepl("only_one_allowed", names(output)))]
  names(output) <- sub(paste0("List_","pop",population,"_"), "", names(output))
  f_check_drugnames(i, output)
  return(output)
}

f_get_L1_lists <- function (i, population){
  output <- i[grep("_1L", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  output <- unlist(output[[1]])
  f_check_drugnames(i, output)
  return(output)
}
  
f_get_L2_lists <- function (i, population){
  output <- i[grep("_2L", names(i))]
  output <- output[!(grepl("_2L_only_after", names(output)))]
  output <- output[grep(paste0("pop", population), names(output))]
  output <- unlist(output[[1]])
  f_check_drugnames(i, output)
  return(output) 
}

f_get_L3_lists <- function (i, population){
  output <- i[grep("_3L", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  output <- unlist(output[[1]])
  f_check_drugnames(i, output)
  return(output) 
}

f_get_L4_lists <- function (i, population){
  output <- i[grep("_4L", names(i))]
  output <- output[grep(paste0("pop", population), names(output))]
  output <- unlist(output[[1]])
  f_check_drugnames(i, output)
  return(output)
}

#' Check for typos in drug names
#' Verifies that all names - name of list and names of drugs entered are correct
f_check_drugnames <- function(i, output) {
  if(sum(unlist(output) %in% i$List_comparators) != length(unlist(output))) {
    message("One or more drugs not found in following list(s).")
    print(output)
    stop("Check source file for typos")
  }
}


#' Utility function to apply the "only after" ruleset to a treatment using subset
#' filtering.
#' 
#' Note that this is a reverse filter (note the ! in front of the conditions). Handy
#' for only returning those results which don't match the conditions. in this case
#' it's essentially like performing a filter for those conditions.
#' 
#' @description At each line, the treatments in rule aren't allowed after trt. filter out.
#' 
#' @param perms permutatins of possible treatments
#' @param trt the treatment to apply this filtering alg to
#' @param rule a character vector of those treatments that aren't allowed after trt
#' 

f_path_onlyAllowedAfter <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  #so flag perms that have trt before drugs in rule
  cat("applying rule.", trt, "is only allowed after", rule,"\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  #first off delete perms with trt first line
  perms <- perms[perms[,1] != trt,]

  for (n in 2:ncol(perms)) {
    # as.data.frame needed to ensure code works for when n=ncol(perms) 
    # (as it reduces to a character vector and crashes the code)
    
    # flag perms with any of the drugs in 'rule' before line n 
    rule_drugs_before_line_n <- as.data.frame(apply(as.data.frame(perms[,1:(n-1)]), 1, function(x) x %in% rule))
    if (ncol(rule_drugs_before_line_n)>nrow(rule_drugs_before_line_n)) rule_drugs_before_line_n <- t(rule_drugs_before_line_n) #needed to ensure format of temp for same reason as.data.frame is needed
    rule_drugs_before_line_n <- as.logical(apply(rule_drugs_before_line_n, 1, sum))
    
    #flag perms with trt in line n
    trt_in_line_n <- perms[,n] == trt   
    
    violators <- trt_in_line_n & !rule_drugs_before_line_n
    
    
    # #temp lines for testing only
    # perms[,6] <- rule_drugs_before_line_n
    # perms[,7] <- trt_in_line_n
    # perms[,8] <- violators

    #remove violating perms
    perms <- perms[!violators,]
    
  }

  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
}

f_path_notAllowedImmediateAfter <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  #so flag perms that have trt before drugs in rule
  cat("applying rule.", trt, "is not allowed immediately after", rule,"\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  for (n in 2:ncol(perms)) {
    # as.data.frame needed to ensure code works for when n=ncol(perms) 
    # (as it reduces to a character vector and crashes the code)
    
    # flag perms with any of the drugs in 'rule' in line n-1 
    rule_drugs_before_line_n <- as.data.frame(apply(as.data.frame(perms[,(n-1)]), 1, function(x) x %in% rule))
    if (ncol(rule_drugs_before_line_n)>nrow(rule_drugs_before_line_n)) rule_drugs_before_line_n <- t(rule_drugs_before_line_n) #needed to ensure format of temp for same reason as.data.frame is needed
    rule_drugs_before_line_n <- as.logical(apply(rule_drugs_before_line_n, 1, sum))
    
    #flag perms with trt in line n
    trt_in_line_n <- perms[,n] == trt   
    
    violators <- trt_in_line_n & rule_drugs_before_line_n
    
    
    # #temp lines for testing only
    # perms[,6] <- rule_drugs_before_line_n
    # perms[,7] <- trt_in_line_n
    # perms[,8] <- violators
    # 
    #remove violating perms
    perms <- perms[!violators,]
    
  }
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
}

f_path_one_in_list <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  #so flag perms that have trt before drugs in rule
  cat("applying rule", trt, ":", rule, "cannot be in one permutation\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  temp <- apply(perms, 1, function(x) sum(x %in% rule))
  
  ##temp lines for testing only
  #perms[,6] <- temp
  
  violators <- temp > 1
  
  #remove violating perms
  perms <- perms[!violators,]
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
}

f_path_only_after_one <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  cat("applying rule:", trt, "can only be after ONE of", rule, "\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  #first off delete perms with trt first line
  perms <- perms[perms[,1] != trt,]
  
  #if drug is 2L then can only have 1 of the tx in rule so count from col 3
  #exclude last column as this is always BSC
  for (n in 3:(ncol(perms)-1)) {
    #if statement below ensures code only runs when there are at least 2 different
    #entries in the perms column.  If is only one then is only BSC or ""
    #code crashes if empty column of "" so this stops it
    if (length(unique(perms[,n]))>1) {
      # flag perms with any of the drugs in 'rule' before line n 
      rule_drugs_before_line_n <- as.data.frame(
          apply(as.data.frame(perms[,1:(n-1)]), 1, function(x) x %in% rule)
        )
      if (ncol(rule_drugs_before_line_n) > nrow(rule_drugs_before_line_n)) {
        rule_drugs_before_line_n <- t(rule_drugs_before_line_n)
      } #this line needed to ensure format of temp for same reason as.data.frame is needed
      
      rule_drugs_before_line_n <- apply(rule_drugs_before_line_n, 1, sum)
      rule_drugs_before_line_n <- rule_drugs_before_line_n >= 2
      
      #flag perms with trt in line n
      trt_in_line_n <- perms[,n] == trt   
      
      violators <- trt_in_line_n & rule_drugs_before_line_n
      
      
      # #temp lines for testing only
      # perms[,6] <- rule_drugs_before_line_n
      # perms[,7] <- trt_in_line_n
      # perms[,8] <- violators
   
      #remove violating perms
      perms <- perms[!violators,]
    }
  }
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)

}

f_path_allowed <- function(perms, rule) {
  #objective of function is to identify perms that violate rule and exclude
  cat("applying rule:", rule, "are only allowed treatments.\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  rule <- c(rule, "BSC", "")
  
  for (n in 1:ncol(perms)) {
    perms <- perms[perms[,n] %in% rule,]
  }
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
  
}

f_path_drug_lines <- function(perms, L1, L2, L3, L4){
  cat("applying rule: drug line restrictions.\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  n <- 0
  for (line in list(L1,L2,L3,L4)){
    n <- n + 1
    line <- c(line, "BSC", "")
    #print(line)
    perms <- perms[(perms[,n] %in% line),]
  }
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
}

f_path_2L_only_after <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  cat("applying rule:", trt, "as 2L+ only allowed after", rule, "\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  #exclude last column as this is always BSC
  for (n in 2:(ncol(perms)-1)) {

    # flag perms with any of the drugs in 'rule' before line n 
    rule_drugs_before_line_n <- as.data.frame(
      apply(as.data.frame(perms[,1:(n-1)]), 1, function(x) x %in% rule)
    )
    if (ncol(rule_drugs_before_line_n) > nrow(rule_drugs_before_line_n)) {
      rule_drugs_before_line_n <- t(rule_drugs_before_line_n)
    } #this line needed to ensure format of temp for same reason as.data.frame is needed
    
    rule_drugs_before_line_n <- apply(rule_drugs_before_line_n, 1, sum)
    rule_drugs_before_line_n <- rule_drugs_before_line_n >= 1
    
    #flag perms with trt in line n
    trt_in_line_n <- perms[,n] == trt   
    
    violators <- trt_in_line_n & !rule_drugs_before_line_n
    
    
    # #temp lines for testing only
    # perms[,6] <- rule_drugs_before_line_n
    # perms[,7] <- trt_in_line_n
    # perms[,8] <- violators
    
    #remove violating perms
    perms <- perms[!violators,]
  }
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
  
}

f_path_2L_only_immediate_after <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  cat("applying rule:", trt, "as 2L+ only allowed immediately after", rule, "\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")
  
  #exclude last column as this is always BSC
  for (n in 2:(ncol(perms)-1)) {
    #flag perms with trt AND (rule[1] or rule[..] or rule[n])
    perms_with_rule_drugs <- apply(perms,2, function(i){i %in% rule})
    perms_with_rule_drugs <- apply(perms_with_rule_drugs,1,function(i){sum(i)>0})
    
    perms_with_trt <- apply(perms,2, function(i){i %in% trt})   
    perms_with_trt <- apply(perms_with_trt, 1, function(i){sum(i)>0})   
    
    perms_with_trt_AND_rule <- perms_with_rule_drugs & perms_with_trt
    
    trt_at_line_n <- perms[,n] == trt
    rule_not_at_line_n_minus_1 <- !(perms[,(n-1)] %in% rule)
  
    violators <- perms_with_trt_AND_rule & trt_at_line_n & rule_not_at_line_n_minus_1 
    
    # #temp lines for testing only
    # perms[,6] <- perms_with_trt_AND_rule
    # perms[,7] <- trt_at_line_n
    # perms[,8] <- rule_not_at_line_n_minus_1
    # perms[,9] <- violators
    
    #remove violating perms
    perms <- perms[!violators,]
  }
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
  
}


f_path_2L_only_one <- function(perms, trt, rule) {
  #objective of function is to identify perms that violate rule and exclude
  cat("applying rule:", trt, "as 2L+ only allowed one of", rule, "\n")
  cat("Permutations before applying rule:", nrow(perms), "\n")

    
  rule_drugs_2Lplus <- as.data.frame(
    apply(as.data.frame(perms[,2:ncol(perms)]), 1, function(x) x %in% rule)
  )
  
  if (ncol(rule_drugs_2Lplus) > nrow(rule_drugs_2Lplus)) {
    rule_drugs_2Lplus <- t(rule_drugs_2Lplus)
  } 
    
  rule_drugs_2Lplus <- apply(rule_drugs_2Lplus, 1, sum)
  violators <- rule_drugs_2Lplus > 1
    
  #remove violating perms
  perms <- perms[!violators,]
  
  cat("Permutations after applying rule :", nrow(perms),"\n")
  return(perms)
  
}

f_path_tx_restrict <- function(sequences, 
                               allowed, L1, L2, L3, L4,
                               only_after, not_immediate_after, 
                               one_in_list, only_after_one,
                               L2_only_after, L2_only_immediate_after, 
                               L2_only_one) {
  
  s <- sequences
  
  cat("Dropping drugs not allowed for this population.\n")
  s <- f_path_allowed(s, allowed[[1]])
  
  s <- f_path_drug_lines(s, L1, L2, L3, L4)
  
  if (length(only_after) > 0) {
    for (n in 1:length(only_after)) {
      print(names(only_after)[n])
      s <- f_path_onlyAllowedAfter(s, names(only_after)[n], only_after[[n]])
    }
  }
  
  if (length(not_immediate_after) > 0) {
    for (n in 1:length(not_immediate_after)) {
      print(names(not_immediate_after)[n])
      s <- f_path_notAllowedImmediateAfter(s, names(not_immediate_after)[n], not_immediate_after[[n]])
    }
  }   
    
  if (length(one_in_list) > 0) {
    for (n in 1:length(one_in_list)) {
      print(names(one_in_list)[n])
      s <- f_path_one_in_list(s, names(one_in_list)[n], one_in_list[[n]])
    }
  }
  
  if (length(only_after_one) > 0) {
    for (n in 1:length(only_after_one)) {
      print(names(only_after_one)[n])
      s <- f_path_only_after_one(s, names(only_after_one)[n], only_after_one[[n]])
    }
  }
  
  if (length(L2_only_after) > 0) {
    for (n in 1:length(L2_only_after)) {
        print(names(L2_only_after)[n])
      s <- f_path_2L_only_after(s, names(L2_only_after)[n], L2_only_after[[n]])
    }
  }

  if (length(L2_only_immediate_after) > 0) {
    for (n in 1:length(L2_only_immediate_after)) {
      print(names(L2_only_immediate_after)[n])
      s <- f_path_2L_only_immediate_after(s, names(L2_only_immediate_after)[n], L2_only_immediate_after[[n]])
    }
  }
  
  if (length(L2_only_one) > 0) {
    for (n in 1:length(L2_only_one)) {
      print(names(L2_only_one)[n])
      s <- f_path_2L_only_one(s, names(L2_only_one)[n], L2_only_one[[n]])
    }
  }
  
    # return the possible sequences now that treatment rules have been applied
  return(s)
}
