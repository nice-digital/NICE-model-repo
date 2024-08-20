# This script defines the functions for extrapolating the treatment line residency
# by line, on/off treatment, BSC and dead. This includes several functions:
# 
# - A function specific to this project to take the 2d numeric matrix of transition probabilities and convert it to the list structure required for the M compiler
#   - Separates columns out into separate matrices by treatment line
#   - Calculates any 1-sum(x) 
#   - splits the matrices into numeric vectors for efficiency later
# - A function to compile the sparse matrix, M, containing the transitions within and across treatment lines per model cycle, which are contained in matrix TP
#   - Input is a NAMED LIST of lists of numeric vectors, and a debug flag which checks probabilities and final trace all sum to 1 at the cost of speed
#     - Lines of therapy MUST be called Lx where x is treatment line
#       - lines of therapy contain 2 elements: on and off relating to on or off treatment
#         - on treatment contains 4 numeric vectors named stay, disc, next, die
#         - off treatment contains 3 numeric vectors named stay, next, die
#         - Summing the nth element from these vectors should always sum to 1 and this is checked if debug==TRUE
#     - BSC MUST be called "bsc"
#       - bsc only contains stay and dead probabilities
#   - Expects each line object to contain 3 or 4 numeric vectors called "stay", "disc", "next", "die"
# - A function to compute a Markov trace from the initial population vector p (summing to 1) as efficiently as possible
#   - The top 2 rows of M change each model cycle (transitions for first-line) according to TP
#   - the rest of M is static
#   - Extrapolation will therefore use block matrix multiplications for efficiency. Before extrapolating:
#     - break off M into M1 and M2 where M1 is the first 2 rows of M and M2 is the rest
#     - Place first set of values for M1 from TP$L1 into M1
#   - Extrapolation then simplifies to this each model cycle:
#     - p_tplus1 = (p_t[1:2] %*% M1) + (p_t[3:N] %*% M2)
#     - both of these produce a 1xN vector where N is the number of columns in M
#     - vectorised sum is near-instant in R, so this is efficient and reduces the size of the matrix multiplication



# Data manipulation function ----------------------------------------------

# This function manipulates the 2d numeric matrix TP into the list format needed
# for the M compiler

requireNamespace("data.table")
requireNamespace("tidyverse")
requireNamespace("collapse")
requireNamespace("Matrix")
requireNamespace("openxlsx")



#' Function to prepare a 2d numeric matrix containing transition probabilities
#' describing a treatment sequencing model for compilation into the large sparse matrix
#' required to extrapolate such a model's treatment line residency. 
#' 
#' NO LONGER USED IN THE MODEL. Legacy code during formation of this strategy.
#' Includes unit testing so remains in the repository as could be used in future.
#'  
#' Note that the naming of columns is strict, and the function makes the following assumptions: 
#'  - There is a column for treatment cycle called "c" in the first column
#'  - There MUST be a prefix Ln where n is treatment line
#'  - The TPs MUST be in the following order for each treatment line:
#'    - Discontinuation probability per cycle
#'    - Next treatment probability per cycle
#'    - Death probability per cycle
#'  - The last line of treatment must have a next line probability (i.e. to go onto BSC or no treatment)
#'  - There MUST be a column for death probability whilst on no treatment
#'  
#' @param tp 2d numeric matrix containing cycle number, transition probabilities. must follow rules and ordering!
#' @param n_lines the number of treatment lines NOT INCLUDING NO TREATMENT AT THE END
#' @param debug double check TPs sum to one (they always do, but more features could be added in future)
#'
#'
#' @examples
#' # using the toy example in the NICE PATT pathways project, produce TP vectors per line and state
#' TP <- f_markov_M_prep(
#'   tp = read.xlsx("./1_Data/TP matrix expansion toy example.xlsm",sheet = "Sheet3"),
#'   n_lines = 4,
#'   debug = FALSE
#' )
#'
f_markov_M_prep <- function(tp,n_lines, debug = FALSE) {
  
  # Convert to data.table for efficiency, as.list and good subsetting later:
  
  if (!all(c("data.table", "data.frame") %in% class(tp))) {
    tp <- as.data.table(tp)
  }
  
  # Calculate target number of columns in data and check it's correct, error if not.
  target_cols <- 1 + (n_lines * 3) + 1
  if(dim(tp)[2] != target_cols) stop(paste0(
    "The transition probability matrix should have ", target_cols, 
    " columns, but instead it has ", dim(tp)[2], ". Each treatment line must have ",
    " a column for discontinuation, next line, death, column 1 must be cycle number ",
    " and the last column must be death probability per cycle whilst on no treatment"
  ))
  
  # Rename the columns to standardize. NOTE that this assumes the correct ordering.
  # We can't possibly know ex ante what these columns will be called, so we have to do this.
  colnames(tp) <- c("c", unlist(lapply(1:n_lines, function(x){paste0("L",x,c("_disc", "_next", "_death"))})),"NT_death")
  
  # next step - separate these out into our list:
  
  line_names        <- c(paste0("L",1:n_lines),"NT")
  names(line_names) <- line_names
  
  # Break the TPs up into separate matrices per treatment line, adding in the required
  # columns:
  line_list <- lapply(line_names, function(treatment_line) {
    which_nams   <- which(substr(colnames(tp),1,2)== treatment_line)
    tp_this_line <- as.list(tp[,which_nams,with=FALSE])
    
    # add up the other elements in a vectorised way to produce the probability of leaving
    # and then p(stay) = 1-p(leave)
    tp_this_line[[paste0(treatment_line,"_stay")]] <- 1-Reduce(`+`,tp_this_line)
    
    # Finally, since we've already separated them via list depth, we can drop the 
    # prefixes from their names:
    names(tp_this_line) <- str_replace(names(tp_this_line),paste0(treatment_line,"_"),"")
    
    return(tp_this_line)
  })
  
  # Check all probabilities sum to 1 at all times
  if (debug) {
    cat("Checking that TPs always sum to 1 for every line:\n")
    unlist(lapply(lapply(line_list, function(tline) Reduce(`+`,tline)) , function(trt_line) all(trt_line == 1)))
  }
  
  # so now we have a nicely organised list of transition probabilities that we can
  # efficiently pull out ([] subsetting is slower than .subset2() and $ within lists)
  
  return(line_list)
  
}


# N calculator ------------------------------------------------------------

#' Function to calculate the value of N for the sequencing modelling function.
#' Includes columns for first-line therapy, then tunnel columns for all subsequent lines
#' plus one for dead patients.
#' 
#' @param n_lines number of treatment lines NOT INCLUDING no treatment at the end
#' @param TH time horizon of the model in cycles
#' 
f_markov_calcN <- function(n_lines, TH) {
  
  # Breadkown of N:
  # 
  #   - 2 for first line on and off treatment
  #   - TH*2 for each subsequent treatment which is not BSC/NT
  #   - TH for BSC/NT
  #   - 1 for dead
  
  return(2 + (2 * TH * (n_lines-1)) + TH + 1)
  
}



# Topleft finder ----------------------------------------------------------

#' Function to locate the top left cell in the matrix M to enter values diagonally
#' from that point for each treatment line. For example 2L on treatment tunnel
#' starts from 3 3 and goes south easterly for time horizon cells. 2L off treatment
#' tunnel starts from 3 TH+3 and goes south easterly etc
#' 
#' @param tun_n The number of tunnels
#' @param TH the time horizon in cycles, including cycle 0
#' @param nt the "no tunnel" states before the tunnels start. Usually 2 for 1L on treatment, 1L off treatment.
#' 
f_markov_topleftFinder <- function(tun_n, TH, nt) {
  
  # Calculate the indices
  tun_r_on       <- sum(nt    , ((2*TH) * (tun_n-1)), 1)
  tun_c_on_stay  <- sum(nt + 1, ((2*TH) * (tun_n-1)), 1)
  tun_c_on_stop  <- tun_c_on_stay + TH
  tun_c_on_next  <- tun_c_on_stop + TH - 1
  tun_r_off      <- tun_c_on_stop - 1
  tun_c_off_stay <- tun_c_on_stop
  tun_c_off_next <- tun_c_on_next
  
  # Put them in a matrix that can be used for co-ordinates in f_vec_in_mat_diag
  outMat <- matrix(
    c(
      tun_r_on , tun_c_on_stay,
      tun_r_on , tun_c_on_stop,
      tun_r_on , tun_c_on_next,
      tun_r_off, tun_c_off_stay,
      tun_r_off, tun_c_off_next
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      c(
        paste0(paste0("L",tun_n + 1,"_"),c("stay","disc","next"),"_on"),
        paste0(paste0("L",tun_n + 1,"_"),c("stay","next"),"_off")
      ),
      c("row", "col")
    )
  )
  return(outMat)
}



# M compiler --------------------------------------------------------------

#' Function to compile a transition probability matrix defining a treatment sequence,
#' where patients move from one treatment to the next, or straight to death, including
#' on and off treatment tunnels. To be used as "block matrices" with a variable and static component.
#' 
#' NOTE, this has been superseded by `f_markov_M_compiler_TPupdate` which functions exactly the same
#' except for taking tp to be pre-processed rather than the result of `f_markov_M_prep`. To exapand,
#' In the patient flow `f_pf_computePF_mk` function, the function `f_pf_mk_ComputeTPs` is now
#' used to compute the TPs in the final form required by `M`. This changed some code inside 
#' this function (simply what to do with the values of `TPs` and how to manipulate them).
#' 
#' WARNING: This function assumes the cycle number is 1 (i.e. M_1 is the result). This is because you should not
#'          be replacing the TPM M including the whole tunnel each cycle, Instead
#'          use this function just once at baseline to create the structure, and then
#'          replace the required values for transitions relating to first-line therapy
#' 
#' @param tp table of TPs. MUST have 1st column cycle, then disc, next, death sequentially for each line, plus no treatment death (e.g. t L1_disc L1_next L1_death L2_disc...NT_death)
#' @param n_lines Number of treatment lines total, including first line and BSC (for 4 tunnels this is 6 with BSC)
#' @param TH time horizon in model cycles (the number of rows in the original tp before using `f_markov_M_prep` on it)
#' @param N The required dimensions of M - use `f_markov_calcN` to calculate it
#' @param nt non-tunnel rows - almost always 2, but if for example second line is not a tunnel (e.g. relapse and remission model) then it would be 4 and tunnels would start at 3L ONLY 2 has been tested
#' 
#' 
f_markov_M_compiler <- function(tp, n_lines, TH, N, nt = 2) {
  
  require(data.table)
  require(collapse)
  require(Matrix)
  
  # Assume cycle = 1, which it should be every time this function is used!
  c <- 1
  
  # Now we have probability to stay in state as well as probability of 
  # discontinuation, next line, death, the rest is putting these columns
  # or derivitives of them into M
  
  cat("Generating matrix M\n")
  # Part 2: define M (WARNING this uses a lot of RAM)
  
  M      <- Matrix(data = 0,nrow = N, ncol = N,sparse = TRUE)
  M[N,N] <- 1 # dead is dead, you can't be come undead.
  
  cat("Finding diagonal start points\n")
  # Part 3: compute top left elements
  tl_elements <- do.call(
    rbind,
    lapply(1:(n_lines-1), function(i_tun) {
      f_markov_topleftFinder(
        tun_n = i_tun,
        TH    = TH,
        nt    = 2 
      )
    })
  )
  
  cat("Adding in BSC/NT tunnel M\n")
  # Add in BSC. patients coming out of the last line of therapy go in here if they're not dead,
  # so it's the same row, next column as the end of the lad
  tl_elements <- rbind(
    tl_elements, 
    matrix(
      data     = c(max(tl_elements[, "col"]), max(tl_elements[, "col"]) + 1),
      nrow     = 1,
      dimnames = list("BSC_stay", c("row", "col"))
    )
  )
  
  # Part 4: populate the matrix:
  # 
  # Right, now we have all the data, all the starting points, and the matrix
  # to put it all into. We can proceed to populate. Let's start with first line
  # as it's different from the others.
  cyc_prob <- lapply(.subset2(tp,"L1"),function(x) x[1])
  vals_1L  <- matrix(
    c(
      1, 1, .subset2(cyc_prob,"stay"),
      1, 2, .subset2(cyc_prob,"disc"),
      1, 3, .subset2(cyc_prob,"next"),
      1, N, .subset2(cyc_prob,"death"),
      2, 2, .subset2(cyc_prob,"stay") + .subset2(cyc_prob,"disc"),
      2, 3, .subset2(cyc_prob,"next"),
      2, N, .subset2(cyc_prob,"death")
    ),
    nrow  = 7,
    byrow = TRUE
  )
  M[vals_1L[,1:2]] <- vals_1L[,3]
  
  
  # Note that we made all of our replacements in one go by making a matrix
  # of matrix co-ordinates and values. If we do the same for all of our values
  # we can enter all of the data needed in one fell swoop. We're therefore going
  # to take the same approach as above, but on a bigger scale!
  # 
  # To do that, we're going to do what we did in the definition of our function
  # f_vec_in_mat_diag - make an index, and then compile our co-ordinate matrix
  # using one row of tl_elements at a time along with the paired data from tp.
  # 
  # The index will always be 0 to TH-1, and will simply be added to the start
  # row and start column to identify where to put the data.
  
  indx <- 0:(TH - 1)
  
  co_ordinate_matrix <- lapply(1:(n_lines-1), function(line_n) {
    # Identify the starting points for this line:
    tl_tun <- f_markov_topleftFinder(
      tun_n = line_n,
      TH    = TH,
      nt    = nt
    )
    
    # Make an id for this line which is used to put the right numbers in the right place
    tline   <- paste0("L",line_n)
    
    # Pull out the columns to put as the data to enter in M
    line_tp <- .subset2(tp,tline)
    p_stay  <- .subset2(line_tp, "stay")
    p_disc  <- .subset2(line_tp, "disc")
    p_next  <- .subset2(line_tp, "next")
    p_death <- .subset2(line_tp, "death")
    
    # make co-ordinate values and value values for the on and off treatment tunnels including
    # death for this treatment line. Note that this assumes equivalency for on/off treatment
    # states for death probability and next probability. Off treatment stay needs recalculating
    # as there's no discontinue, so it's 1-p_death-p_next
    co_list <- list(
      # on-treatment tunnel:
      co_on_stay   = matrix(c(tl_tun[1,"row"] + indx, tl_tun[1,"col"] + indx,p_stay) ,ncol=3),
      co_on_disc   = matrix(c(tl_tun[2,"row"] + indx, tl_tun[2,"col"] + indx,p_disc) ,ncol=3),
      co_on_next   = matrix(c(tl_tun[3,"row"] + indx, tl_tun[3,"col"] + indx,p_next) ,ncol=3),
      co_on_death  = matrix(c(tl_tun[1,"row"] + indx, rep(N,TH)             ,p_death),ncol=3),
      
      # Off treatment tunnel for this line:
      co_off_stay  = matrix(c(tl_tun[4,"row"] + indx, tl_tun[4,"col"] + indx,1 - p_death - p_next),ncol=3),
      co_off_next  = matrix(c(tl_tun[5,"row"] + indx, tl_tun[5,"col"] + indx,p_next              ),ncol=3),
      co_off_death = matrix(c(tl_tun[4,"row"] + indx, rep(N,TH)             ,p_death             ),ncol=3)
    )
    
    # Bind it together into one larger matrix and return that
    return(do.call(rbind,co_list))
  })
  
  # co_ordinate_matrix now has all of the data organised except for the BSC/NT tunnel.
  # We add that here:
  
  bsc_start_point <- f_markov_topleftFinder(
    tun_n = n_lines - 1,
    TH    = TH,
    nt    = nt
  )[5,"col"]
  co_bsc <- rbind(
    matrix(c(bsc_start_point + indx, bsc_start_point + 1 + indx,tp$NT$stay) ,ncol=3),
    matrix(c(bsc_start_point + indx, rep(N,TH)                 ,tp$NT$death),ncol=3)
  )
  
  # bind it all together into one massive matrix of co-ordinates and values:  
  co_ordinate_matrix <- rbind(do.call(rbind,co_ordinate_matrix),co_bsc)
  
  cat("Inserting values\n")
  # Put all values into the matrix in one vectorised command:
  M[co_ordinate_matrix[,1:2]] <- co_ordinate_matrix[,3]
  
  # The last cycle of BSC is an issue, stay probability isn't entered.
  M[N-1,N-1] <- 1-M[N-1,N]
  
  # Well all rows sum to 1 so that's promising:
  # all(round(rowSums(M),10) == 1) # reads true:
  return(M)
  
}

# Markov trace calculator -------------------------------------------------

# Simple function to take the results of f_markov_M_compiler and use them to 
# extrapolate the population using a block-matrix approach. Requires the first-line 
# transitions from TP (i.e. TP$L1 after applying f_markov_M_prep)
# 
# We try to do this as efficiently as possible as it's a big number crunch, hence
# the sparse matrix approach, and splitting to blocks:


#' Function to extrapolate a treatment sequence using a block sparse matrix multiplication
#' method. THIS FUNCTION IS NOT USED AND HAS BEEN SUPERSEDED BY `f_markov_M_extrapolator`
#' 
#' @param L1_tp transition probabilities for first-line treatment, resulting from `f_markov_M_prep`. will change within M each model cycle (hence block matrix approach)
#' @param TH time horizon in model cycles. Not inferred so that a shorter time horizon can be entered if desired.
#' @param M TP matrix compiled using `f_markov_M_compiler`
#' @param N dimension of M, calculate either as `dim(M)[1]`, `dim(M)[2]`, or use `f_markov_calcN`
#' @param p (optional) initial population vector. Usually sparse matrix 1xN with value of 1 in 1,1. Assumed so if not entered, but does allow mixed lines at baseline if desired.
#'
#' @details superseded by `f_markov_M_extrapolator`. `f_markov_M_extrapolator`
#'          does the same job as f_markov_sequencExtrapolator, but expects the 
#'          `L1_tp` argument to be the result of `f_pf_mk_ComputeTPs`, rather
#'          than computing those TPs "in-house". This function is tested via
#'          unit test so remains.
#' 
f_markov_sequencExtrapolator <- function(L1_tp,TH,M,N,p=NULL) {
  
  require(Matrix)
  require(utils)
  
  # If the user doesn't provide p, make it and put everyone in first-line therapy
  # at baseline (90+% of cases this is true, though sometimes reasonably people choose BSC/NT
  # from baseline, or skip into subsequent treatment lines)
  if (is.null(p)) {
    p      <- Matrix(0,nrow = 1,ncol = N,sparse = TRUE)
    p[1,1] <- 1
  }
  
  # Reduce is quite efficient here, I imagine it's hard to beat efficiency wise.
  # Reduce is a kind of advanced for loop which lets you iterate over an entire list
  # each round in a safe environment. This is useful here as without committing things
  # to .GlobalEnv repeatedly (slow) we can iterate over p extrapolating the population
  # in every health state in a very efficient way:
  
  # Remember that every row of M after the 2nd one is STATIC, it never changes!
  
  M1 <- M[1:2,]
  M2 <- M[3:N,]
  
  cat(paste0("Extrapolating a model with ", N, " health states for ",TH," model cycles, give me a minute...\n"))
  
  pb <- txtProgressBar(
    min = 1,
    max = TH-1,
    initial = 1,
    width = 50,
    style = 3,
    char = "="
  )
  
  trace_separate <- Reduce(
    x = 1:(TH-1),
    init = p,
    accumulate = TRUE,
    f = function(prev, t) {
      
      setTxtProgressBar(pb, t)
      
      # get the 1L probabilities, put them in a co-ordinate matrix and replace
      # those elements of M. Note that t starts from 1 but model cycles start from 0.
      # Therefore the first element of this output is p_1 (not p_0). p_0 can be added
      # at the end as we already know it, so why calculate it again :)
      
      cyc_prob <- lapply(L1_tp,\(x) x[t])
      vals_1L  <- matrix(
        c(
          1, 1, .subset2(cyc_prob,"stay"),
          1, 2, .subset2(cyc_prob,"disc"),
          1, 3, .subset2(cyc_prob,"next"),
          1, N, .subset2(cyc_prob,"death"),
          2, 2, .subset2(cyc_prob,"stay") + .subset2(cyc_prob,"disc"),
          2, 3, .subset2(cyc_prob,"next"),
          2, N, .subset2(cyc_prob,"death")
        ),
        nrow  = 7,
        byrow = TRUE
      )
      M1[vals_1L[,1:2]] <- vals_1L[,3]
      
      # Now that M has been updated, we can apply the matrix multiplication as
      # a block matrix multiplication by doing the variable and static parts separately
      # and adding them up:
      
      (prev[1:2] %*% M1) + (prev[3:N] %*% M2)
      
    }
  )
  
  cat(paste0("\nDone! Sticking the rows together to give you your trace:\n"))
  
  # Stick them together efficiently by defining a numeric matrix (much faster than do.call(rbind,x)).
  # This is the result. Note that all rows sum to 1 to at least 12 decimal places :)
  # all(round(rowSums(trace),12)==1)
  # 
  # Make this sparse as well as there's loads of zeros in there
  # 
  return(Matrix(
    matrix(
      data = unlist(lapply(trace_separate,as.numeric)),
      nrow = TH,
      ncol = N,
      byrow = TRUE
    ),
    sparse = TRUE
  ))
  
}


# Adding up selected columns ----------------------------------------------

#' Function to "consolidate" the expanded trace (whcih may have 18k+ columns!) into
#' the descrete health states that HEOR typically use to comprehend a cost-effectiveness
#' model setting. 
#' 
#' @param full_trace This is `M`. the fully expanded block-diagonal sparse matrix
#' @param L1_onoff_adjustment First-line on and off treatment adjustment. `TPs$L1$death_on + TPs$L1$next_on` where `f_pf_mk_ComputeTPs` was used to get `TPs`
#' @param split_list a list which can be used to add up columns before and after a time point (e.g. 10 cycles in that state)
#' @param TH Time horizon in cycles accounting for cycle 0
#' @param n_lines nubmer of treatment lines
#' @param discFacQ discount FACTOR for QALYs (i.e. discount per cycle going from time t to TH)
#' @param discFacC discount FACTOR for costs (i.e. discount per cycle going from time t to TH)
#' 
#' 
#' @details We can't show all the health states as there are many thousands. It
#'          uses too much data and is not informative. therefore consolidated
#'          is more useful because we can observe patient flow through the
#'          treatment lines. 
#'          
#'          For example for a model with `TH=2080` and 4 lines there are `f_markov_calcN(4,2080)=14563`
#'          health states. This is too many to plot or use informatively to 
#'          characterize a decision problem.
#'          
#'          Instead we add them up by discrete health state so that we can
#'          generate plots and compute all non-cost elements of the model 
#'          (which do not currently require per time in state level calculations)
#' 
f_markov_traceConsolidator <-
  function(full_trace,
           L1_onoff_adjustment,
           split_list = NULL,
           TH,
           n_lines,
           discFacQ,
           discFacC) {
    
  
  # Validation first:
  if(n_lines > 1) {
    target_names <- c(
      "L1_on",
      "L1_off",
      unlist(lapply(paste0("L",2:n_lines),function(lam) {
        paste0(lam,c("_on","_off"))
      })),
      "BSC"
    )
  } else {
    target_names <- c(
      "L1_on",
      "L1_off",
      "BSC"
    )
  }
  
  if (!is.null(split_list)) {
    if(!all(names(split_list) %in% target_names)) stop("The names in your split_list don't match what they should be!")
  }
  
  N <- ncol(full_trace)
  
  # Good, so we have everything we need in split list.
  # Now we work out some co-ordinates. it's simply a series of multiplications of TH
  
  # Figure out which columns of the trace correspond to what treatment lines:
  col_starts <- c(1:3,(TH * 1:((N-3)/TH))+3)
  names(col_starts) <- c(target_names,"dead")
  
  # Calculate entrants entrants for all lines 2+.
  L1_on <- full_trace[,1]
  L1_exit <- L1_on - shift(L1_on,type = "lead",fill = 0)
  L1_off_entrants <- L1_exit - (L1_on*L1_onoff_adjustment)
  
  entrants <- full_trace[,col_starts[3:length(col_starts)]]
  entrants[,ncol(entrants)] <- c(0,diff(entrants[,ncol(entrants)]))
  colnames(entrants) <- c(target_names[3:length(target_names)],"dead")
  
  entrants <- cbind(L1_off = L1_off_entrants,entrants)
  
  
  line_labs <- names(col_starts)
  names(line_labs) <- line_labs
  
  # We can generate some useful stuff right from the off. for example, OS & first-line-PFS (s for survival):
  
  s <- list()
  
  # Overall and 1L population
  s$OS <- rowSums(full_trace[,1:(N-1)])
  s$L1 <- rowSums(full_trace[,1:2])
  
  # total line populations over time by treatment status (i.e. with no splits):
  s$full_lines <- lapply(1:(length(col_starts) - 1), function(living_state) {
    
    col_index <- col_starts[living_state]:(col_starts[living_state+1]-1)
    
    if (length(col_index) == 1) {
      full_trace[,col_index]
    } else {
      rowSums(full_trace[,col_index])
    }
  })
  s$full_lines[[length(s$full_lines) + 1]] <- full_trace[,N]
  s$full_lines <- matrix(unlist(s$full_lines),ncol = length(line_labs))
  colnames(s$full_lines) <- line_labs
  
  s$entrants <- entrants
  s$col_starts <- col_starts
  
  # discFacC <- diag(discFacC)
  # discFacQ <- diag(discFacQ)
  
  # Discounted versions:
  s$disc <- list(
    C = list(
      full_lines = discFacC * s$full_lines,
      entrants   = discFacC * s$entrants
    ),
    Q = list(
      full_lines = discFacQ * s$full_lines,
      entrants   = discFacQ * s$entrants
    )
  )
  
  # split list is optional
  if (!is.null(split_list)) {
    # cycle through the line labs to help us select the correct columns. Note, no point
    # in splitting the dead population so we're not doing it.
    s$split_pop <- lapply(3:(length(line_labs) - 1), function(iline) {
      tline          <- line_labs[iline]
      split_settings <- split_list[[tline]]
      
      if (split_settings$any_split == FALSE) {
        # not splitting this line, don't return anything as we've already calculated 
        # it in full_lines
        return(NULL)
      } else {
        # We're going to split this population by before and after the multiple 
        # split times we've been given:
        
        col_index <- col_starts[iline]:(col_starts[iline+1]-1)
        
        split_pops <- do.call(
          cbind,
          lapply(split_settings$t, function(t_split) {
            matrix(
              data = c(
                rowSums(full_trace[,col_index[1:(t_split-1)]]),
                rowSums(full_trace[,col_index[t_split:length(col_index)]])
              ),
              ncol = 2,
              byrow = FALSE,
              dimnames = list(
                NULL,
                c(paste0(tline,"_split",t_split,"_before"),paste0(tline,"_split",t_split,"_after"))
              )
            )
          })
        )
      }
    })
    names(s$split_pop) <- line_labs[3:(length(line_labs) - 1)]
  }
  
  return(s)
}



# Updated M related functions ---------------------------------------------

#' Second version of f_markov_M_compiler. Functionality is the same as f_markov_M_compiler,
#' but the expected values of `tp` are different (the result of function `f_pf_mk_ComputeTPs` instead of `f_markov_M_prep`)
#' 
#' @param tp the result of `f_pf_mk_ComputeTPs` to compute the TPs from the extrapolated survival
#' @param n_lines the number of treatment lines
#' @param nt Number of non-tunnel states before the tunnels start (default is 2 for 1L_on and 1L_off)
#' @param verbose additional console output on progress of compiling `M`
#' 
#' 
#' @details See function `f_markov_M_compiler`. Generates a block-diagonal sparse
#'          matrix, with diagonal elements starting from entry to the first tunnel
#'          state and ending at death. Matrix-multiplication of a vector of population
#'          against this `M` computes tunnel states because the only ways out of
#'          the tunnels are off treatment, death, next line for on-treatment states
#'          and death and next line for off-treatment states.
#'          
#'          this approach also allows different complex vectors of costs and so on
#'          to be applied given time in state, such that those going to a line
#'          will first get the initiation cost of that line then follow the dosing
#'          schedule and any stopping rules, titration, stopping rules within.
#'          
#'          This `M` approach allows a great deal of flexibility in a markov sequencing
#'          setting, which would not be possible in Excel.
#' 
f_markov_M_compiler_TPupdate <- function(tp, n_lines, TH, nt = 2, verbose = FALSE) {
  
  # Compute the size of the required matrix M
  N <- f_markov_calcN(n_lines = n_lines, TH = TH)
  
  if (verbose) cat(paste0("Generating matrix M (",N,"x",N,")\n"))
  if (verbose) cat(paste0("state count: 2 for 1L on/off, ", TH*2, " for each active line (",(TH*2)*(n_lines-1),"), ", TH, " for BSC, and 1 for dead\n"))
  # Part 2: define M
  
  M      <- Matrix(data = 0,nrow = N, ncol = N,sparse = TRUE)
  M[N,N] <- 1 # dead is dead, you can't be come undead.
  
  if (verbose) cat("Finding diagonal start points\n")
  # Part 3: compute top left elements
  if (n_lines > 1) {
    line_vec <- 1:(n_lines-1)
    tl_elements <- do.call(
      rbind,
      lapply(line_vec, function(i_tun) {
        f_markov_topleftFinder(
          tun_n = i_tun,
          TH    = TH,
          nt    = 2 
        )
      })
    )
  } 
  
  if (verbose) cat("Adding in BSC/NT tunnel M\n")
  # Add in BSC. patients coming out of the last line of therapy go in here if they're not dead,
  # so it's the same row, next column as the end of the lad
  
  if (n_lines == 1) {
    tl_elements <- f_markov_topleftFinder(
      tun_n = 1,
      TH    = TH,
      nt    = 2 
    )[1,]
  } else {
    tl_elements <- rbind(
      tl_elements, 
      matrix(
        data     = c(max(tl_elements[, "col"]), max(tl_elements[, "col"]) + 1),
        nrow     = 1,
        dimnames = list("BSC_stay", c("row", "col"))
      )
    )
  }
  
  
  # Part 4: populate the matrix:
  # 
  # Right, now we have all the data, all the starting points, and the matrix
  # to put it all into. We can proceed to populate. Let's start with first line
  # as it's different from the others.
  # Assume cycle = 1, which it should be every time this function is used!
  cyc_prob <- lapply(.subset2(tp,"L1"),function(x) x[1])
  vals_1L  <- matrix(
    c(
      1, 1, .subset2(cyc_prob,"stay_on"),
      1, 2, .subset2(cyc_prob,"disc_on"),
      1, 3, .subset2(cyc_prob,"next_on"),
      1, N, .subset2(cyc_prob,"death_on"),
      2, 2, .subset2(cyc_prob,"stay_off"),
      2, 3, .subset2(cyc_prob,"next_off"),
      2, N, .subset2(cyc_prob,"death_off")
    ),
    nrow  = 7,
    byrow = TRUE
  )
  M[vals_1L[,1:2]] <- vals_1L[,3]
  
  
  # Note that we made all of our replacements in one go by making a matrix
  # of matrix co-ordinates and values. If we do the same for all of our values
  # we can enter all of the data needed in one fell swoop. We're therefore going
  # to take the same approach as above, but on a bigger scale!
  # 
  # To do that, we're going to do what we did in the definition of our function
  # f_vec_in_mat_diag - make an index, and then compile our co-ordinate matrix
  # using one row of tl_elements at a time along with the paired data from tp.
  # 
  # The index will always be 0 to TH-1, and will simply be added to the start
  # row and start column to identify where to put the data.
  
  indx <- 0:(TH - 1)
  
  if (n_lines > 1) {
    # this is for 2L+, 1L does not have a tunnel, it's just the first 2 columns (on and off trt)
    co_ordinate_matrix <- lapply(1:(n_lines-1), function(line_n) {
      
      # Identify the starting points for this line:
      tl_tun <- f_markov_topleftFinder(
        tun_n = line_n,
        TH    = TH,
        nt    = nt
      )
      
      # Make an id for this line which is used to put the right numbers in the right place
      tline   <- paste0("L",line_n+1)
      
      # Pull out the columns to put as the data to enter in M
      line_tp <- .subset2(tp,line_n+1)
      
      # make co-ordinate values and value values for the on and off treatment tunnels including
      # death for this treatment line. Note that this assumes equivalency for on/off treatment
      # states for death probability and next probability. Off treatment stay needs recalculating
      # as there's no discontinue, so it's 1-p_death-p_next
      co_list <- list(
        # on-treatment tunnel:
        co_on_stay   = matrix(c(tl_tun[1,"row"] + indx[-TH], tl_tun[1,"col"] + indx[-TH],.subset2(line_tp, "stay_on")[-TH]) ,ncol=3),
        co_on_disc   = matrix(c(tl_tun[2,"row"] + indx[-TH], tl_tun[2,"col"] + indx[-TH],.subset2(line_tp, "disc_on")[-TH]) ,ncol=3),
        co_on_next   = matrix(c(tl_tun[3,"row"] + indx     , rep(tl_tun[3,"col"], TH)   ,.subset2(line_tp, "next_on")     ) ,ncol=3),
        co_on_death  = matrix(c(tl_tun[1,"row"] + indx     , rep(N,TH)                  ,.subset2(line_tp, "death_on")    ),ncol=3),
        
        # Off treatment tunnel for this line:
        co_off_stay  = matrix(c(tl_tun[4,"row"] + indx[-TH], tl_tun[4,"col"] + indx[-TH],.subset2(line_tp, "stay_off")[-TH]),ncol=3),
        co_off_next  = matrix(c(tl_tun[5,"row"] + indx     , rep(tl_tun[5,"col"], TH)   ,.subset2(line_tp, "next_off") ),ncol=3),
        co_off_death = matrix(c(tl_tun[4,"row"] + indx     , rep(N,TH)                  ,.subset2(line_tp, "death_off")),ncol=3)
      )
      
      # Bind it together into one larger matrix and return that
      return(do.call(rbind,co_list))
    })
  }
  
  
  # co_ordinate_matrix now has all of the data organised except for the BSC/NT tunnel.
  # We add that here. Note that the start point is still correct if n_lines is 1
  # so we're looking at "tunnel 0"
  
  bsc_start_point <- f_markov_topleftFinder(
    tun_n = n_lines - 1,
    TH    = TH,
    nt    = nt
  )[5,"col"]
  
  
  co_bsc <- rbind(
    matrix(c(bsc_start_point + indx, bsc_start_point + 1 + indx,.subset2(tp$NT,"stay_on")) ,ncol=3),
    matrix(c(bsc_start_point + indx, rep(N,TH)                 ,.subset2(tp$NT,"death_on")),ncol=3)
  )
  
  # bind it all together into one massive matrix of co-ordinates and values:  
  if (n_lines > 1) {
    co_ordinate_matrix <- rbind(do.call(rbind,co_ordinate_matrix),co_bsc)
  } else {
    co_ordinate_matrix <- co_bsc
  }
  
  
  if (verbose) cat("Inserting values\n")
  # Put all values into the matrix in one vectorised command:
  M[co_ordinate_matrix[,1:2]] <- co_ordinate_matrix[,3]
  
  # The last cycle of BSC is an issue, stay probability isn't entered.
  M[N-1,N-1] <- 1-M[N-1,N]
  
  # Well all rows sum to 1 so that's promising:
  # all(round(rowSums(M),10) == 1) # reads true:
  return(M)
  
}

#' Updated version of the sequence extrapolator `f_markov_sequencExtrapolator` which
#' takes results from `f_pf_mk_ComputeTPs` to populate `L1_tp`. Otherwise the same
#' as `f_markov_sequencExtrapolator`.
#' 
#' @param L1_tp transition probabilities for first-line treatment, resulting from `f_pf_mk_ComputeTPs`. will change within M each model cycle (hence block matrix approach)
#' @param TH time horizon in model cycles. Not inferred so that a shorter time horizon can be entered if desired.
#' @param M TP matrix compiled using `f_markov_M_compiler`
#' @param N dimension of M, calculate either as `dim(M)[1]`, `dim(M)[2]`, or use `f_markov_calcN`
#' @param p (optional) initial population vector. Usually sparse matrix 1xN with value of 1 in 1,1. Assumed so if not entered, but does allow mixed lines at baseline if desired.
#'
#' @details See `f_markov_sequencExtrapolator`
#' 
f_markov_M_extrapolator <- function(L1_tp,TH,M,N,p=NULL,prog=FALSE,verbose=FALSE) {
  
  # If the user doesn't provide p, make it and put everyone in first-line therapy
  # at baseline (90+% of cases this is true, though sometimes reasonably people choose BSC/NT
  # from baseline, or skip into subsequent treatment lines)
  if (is.null(p)) {
    p      <- Matrix(0,nrow = 1,ncol = N,sparse = FALSE)
    p[1,1] <- 1
  }
  
  # Reduce is quite efficient here, I imagine it's hard to beat efficiency wise.
  # Reduce is a kind of advanced for loop which lets you iterate over an entire list
  # each round in a safe environment. This is useful here as without committing things
  # to .GlobalEnv repeatedly (slow) we can iterate over p extrapolating the population
  # in every health state in a very efficient way:
  
  # Remember that every row of M after the 2nd one is STATIC, it never changes!
  
  M1 <- M[1:2,]
  M2 <- M[3:N,]
  
  if (verbose) cat(paste0("Extrapolating a model with ", N, " health states for ",TH," model cycles, give me a minute...\n"))
  
  if (prog) {
    pb <- txtProgressBar(
      min = 1,
      max = TH-1,
      initial = 1,
      width = 50,
      style = 3,
      char = "="
    )
  }
  
  tp1 <- as.list(t(as.data.table(L1_tp)))
  
  
  # Update, make all the top slices of M for all cycles before extrapolating
  m1_list <- lapply(1:length(L1_tp$disc_on), function(model_cycle) {
    coo <- matrix(
      c(1,1,.subset2(.subset2(L1_tp, "stay_on"), model_cycle),
        1,2,.subset2(.subset2(L1_tp, "disc_on"), model_cycle),
        1,3,.subset2(.subset2(L1_tp, "next_on"), model_cycle),
        1,N,.subset2(.subset2(L1_tp, "death_on"), model_cycle),
        2,2,.subset2(.subset2(L1_tp, "stay_off"), model_cycle),
        2,3,.subset2(.subset2(L1_tp, "next_off"), model_cycle),
        2,N,.subset2(.subset2(L1_tp, "death_off"), model_cycle)
      ),
      nrow  = 7,
      byrow = TRUE
    )
    M1[coo[,1:2]] <- coo[,3]
    return(M1)
  })
  
  trace_separate <- Reduce(
    x = 1:(TH-1),
    init = p,
    accumulate = TRUE,
    f = function(prev, t) {
      if(prog) setTxtProgressBar(pb, t)
      (prev[1:2] %*% .subset2(m1_list,t)) + (prev[3:N] %*% M2)
    }
  )
  
  # free up RAM
  rm(m1_list)
  rm(M2)
  rm(M)
  
  if (verbose) cat(paste0("\nDone! Sticking the rows together to give you your trace:\n"))
  
  # Stick them together efficiently by defining a numeric matrix (much faster than do.call(rbind,x)).
  # This is the result. Note that all rows sum to 1 to at least 12 decimal places :)
  # all(round(rowSums(trace),12)==1)
  # 
  # Make this sparse as well as there's loads of zeros in there
  # 
  return(Matrix(
    unlist(lapply(trace_separate, function(x) x@x)),
    nrow = TH,
    byrow = TRUE,
    sparse = TRUE
  ))
  
}



