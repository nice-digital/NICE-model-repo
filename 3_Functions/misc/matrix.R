# Extremely useful functions for matrix operations



# diagonal vector assignment ----------------------------------------------

# Function to allow diagonal insertion of a vector into a matrix in an arbitrary
# direction and position. Extremely useful for the tunnel state stuff central
# to this disease model.

# Courtesy of https://stackoverflow.com/questions/28746023/how-to-insert-values-from-a-vector-diagonally-into-a-matrix-in-r
# there is in fact a way to stick elements diagonally into an existing matrix:


#'Function to insert a vector diagonally into a matrix at an arbitrary location.
#'Assumes south easterly direction of insertion unless option is FALSE
#'
#'@param A A matrix - doesn't have to be square, but the bounds cannot be breached by b
#'@param b A vector - obviously mustn't go outside of the matrix when populating
#'@param start_row the row to "start" from: for nw and ne direction, it will start from `start_row + length(b) - 1`
#'@param start_col the column to "start" from: for sw and nw direction, it will start from `start_col + length(b) - 1`
#'
#'@examples
#'
#'M <- matrix(nrow=10,ncol=10)
#'
#'f_vec_in_mat_diag(
#' A         = M,
#' b         = runif(5),
#' start_row = 3,
#' start_col = 4,
#' direction = "se"
#')
#'
#'
#'f_vec_in_mat_diag(
#' A         = M,
#' b         = runif(5),
#' start_row = 3,
#' start_col = 4,
#' direction = "ne"
#')
#'
#'
f_vec_in_mat_diag <- function(A, b, start_row = 1, start_col = 1, direction = "se") {
  if(direction == "se") {
    # Starting from top left moving south east (default)
    indx <- 0:(length(b) - 1)
    ind_mat <- matrix(c(start_row + indx, start_col + indx),ncol=2)
  } else if (direction == "nw") {
    # starting from bottom right moving north west
    indx <- (length(b) - 1):0
    ind_mat <- matrix(c(start_row + indx, start_col + indx),ncol=2)
  } else if (direction == "sw") {
    # starting from top right going south west
    indxr <- 0:(length(b) - 1)
    indxc <- (length(b) - 1):0
    ind_mat <- matrix(c(start_row + indxr, start_col + indxc),ncol=2)
    stop("south west (sw) and north west (nw) are not supported...yet!")
  } else if (direction == "ne") {
    # starting from bottom left going north east
    indxr <- (length(b) - 1):0
    indxc <- 0:(length(b) - 1)
    ind_mat <- matrix(c(start_row + indxr, start_col + indxc),ncol=2)
  } else {
    stop('direction must be one of "se" (default), "nw", "sw", or "ne"')
  }
  # put the values in the matrix using the co-ordinate matrix & return the result
  A[ind_mat] <- b
  return(A)
}



# So, we now have a Markov trace for our treatment sequence.
# We want to apply the cost vector that we have to diagonals
# starting from 1,1 to TH,TH, then 2,1 to TH,TH-1 then 3,1 to TH,TH-2
# and so on with that same vector each time, but 1 shorter!

# The MOST efficient way to populate a matrix in R is to define a co-ordinate
# matrix, which populates a bunch of elements all at once. This is a matrix
# of 3 columns, row, column and value. 

# As we know what our values are each time, and we also know our rows and columns,
# we can define one of these. To make it easier, I define some functions which
# make the value column


# Only used within this script
f_mat_diagDownCM <- function(TH,vals) {
  
  # V is for values - this is the matrix which is going to get
  # populated in the end
  V      <- Matrix(data = 0,nrow = TH, ncol = TH,sparse = TRUE)
  V[TH,TH] <- 1 
  
  # Make the co-ordinate matrix by first defining the 3 columns, then
  # collapsing them into a matrix:
  one_to_TH <- 1:TH
  row_numbers <- Reduce(
    x = 1:(TH-1),
    init = one_to_TH,
    accumulate = TRUE,
    f = function(prev, cyc) {
      (one_to_TH + cyc)[1:(TH-cyc)]
    }
  )
  
  # To see what this looks like, look at this:
  # row_numbers[(length(row_numbers)-3):length(row_numbers)]
  # 
  # This shows the last 4 runs, with the first element being the cycle that someone
  # enters the state in. Collapsing this with unlist() gives us our rows column:
  row_numbers <- unlist(row_numbers,use.names = FALSE)
  
  
  # column number
  col_numbers <- Reduce(
    x = 1:(TH-1),
    init = one_to_TH,
    accumulate = TRUE,
    f = function(prev, cyc) {
      1:(TH-cyc)
    }
  )
  
  # To see this:
  # col_numbers[(length(col_numbers)-3):length(col_numbers)]
  # 
  # It moves down into the bottom left such that the people coming into this line 
  # in cycle TH-1 are in the TH-1 th row in the 1st column: 
  col_numbers <- unlist(col_numbers,use.names = FALSE)
  
  # So the co-ordinates for the last 10 elements to put in our matrix V
  # are:
  # cbind(tail(row_numbers,10),tail(col_numbers,10))
  
  # Now all we need to do is the same thing again for the values:
  vals_to_enter <- Reduce(
    x = 1:(TH-1),
    init = vals,
    accumulate = TRUE,
    f = function(prev, cyc) {
      vals[1:(TH-cyc)]
    }
  )
  vals_to_enter <- unlist(vals_to_enter,use.names = FALSE)
  
  # vals_to_enter[(length(vals_to_enter)-3):length(vals_to_enter)]
  # cbind(tail(row_numbers,10),tail(col_numbers,10),tail(vals_to_enter,10))
  
  # So now we have all 3 of our co-ordinates!:
  
  co_ord <- matrix(
    c(row_numbers,col_numbers, vals_to_enter),
    ncol = 3
  )
  
  V[co_ord[,1:2]] <- co_ord[,3]
  
  return(V)
  
}

if (FALSE) {
  # Function to populate a matrix which should be element-wise multiplied by a vertical slice (columns between two points) of TRACE
  library(Matrix)
  
  TH <- 2089
  
  # example dose schedule, dose increases over time to 100, patients are dosed every other week throughout.
  dose_over_time <- c(c(10,0),c(20,0),c(35,0),rep(c(50,0),2),rep(c(100,0),ceiling((TH-10)/2)))[1:TH]
  
  # Apply stopping rule to dose schedule
  stopping_rule <- 104
  dose_over_time[stopping_rule:TH] <- 0
  
  # Apply RDI
  rdi <- 0.8
  dose_over_time <- dose_over_time * rdi
  
  # Not going to do wastage in this example, just cost per mg for simplicity:
  cost_per_mg <- 5
  
  # Cost schedule including dose changes over time and stopping rule.
  cost_over_time <- dose_over_time * cost_per_mg
  
  # This is what drug cost in this hypothetical looks like over time, thus correctly allowing people to discontinue treatment between doses.
  plot(cost_over_time[1:(52*3)],type="l")
  
  # Load in some fitted survival data from the model (you have this rds file in your branch so it should just run)
  st_list <- readRDS("./1_Data/st_list_for_dawn.rds")
  
  # Get our compiling and extrapolating functions
  source(file.path("3_Functions/markov/markov.R"))
  source(file.path("3_Functions/patient_flow/markov.R"))
  source(file.path("./3_Functions/markov/TEMP M compiler.R"))
  
  
  # Work out transition probabilities
  TPs <- f_pf_mk_ComputeTPs(st_list)
  
  # Calculate the trace - note that n_lines does NOT include BSC (it's active treatment lines), take 1 of the length of st_list. This fixes empty columns issue I think!:
  M <- f_markov_M_compiler_TPupdate(
    tp      = TPs, 
    n_lines = length(st_list)-1, 
    TH      = TH, 
    nt      = 2, 
    verbose = TRUE
  )
  
  TRACE <- f_markov_M_extrapolator(
    L1_tp = TPs$L1,
    TH = TH,
    M = M,
    N = f_markov_calcN(n_lines = length(st_list)-1, TH = TH),
    verbose = TRUE,
    prog = TRUE
  )
  
  
  
  cost_matrix <- f_mat_diagDownCM(TH,cost_over_time)
  
  # Let's say that we're going to apply these costs to 2L on treatment: 
  line_2_starts_where <- 3
  line_2_cols <- line_2_starts_where:(TH+line_2_starts_where-1)
  
  line_2_on_pop <-  TRACE[,line_2_cols]
  
  # Note one can't get into 2L on until cycle 3 start, so 
  cost_per_cycle <- line_2_on_pop[3:TH,] * cost_matrix[1:(TH-2),]
  
  # Results:
  plot(rowSums(cost_per_cycle), type="l")
  sum(cost_per_cycle)
  
  # So for pre-104 week population:
  rowSums(cost_per_cycle[,1:104])
  
  # Rest should be 0
  rowSums(cost_per_cycle[,105:TH])
  
  
  # Now, using the consolidator:
  split_list <- lapply(1:(((length(st_list))*2)-1), function(line_n) {
    list(any_split = FALSE, t = NULL)
  })
  names(split_list) <- c(
    unlist(lapply(1:(length(st_list)-1), function(line_n) {
      c(paste0("L",line_n,"_on"),paste0("L",line_n,"_off"))
    })),
    "BSC"
  )
  
  # Add in the details for our hypothetical treatment:
  # 
  #   - Dose loading lasts for 10 cycles
  #   - Stopping rule at 104 cycles
  # 
  split_list$L2_on$any_split <- TRUE
  split_list$L2_on$t         <- c(10,104)
  
  # Get the consolidated trace:
  consolidated_trace <- f_markov_traceConsolidator(
    full_trace = TRACE,
    split_list = split_list,
    TH         = TH,
    n_lines    = length(st_list)-1
  )
  
  # Just to check our populations are right:
  sum(rowSums(line_2_on_pop) - consolidated_trace$full_lines[,"L2_on"])
  
  # Combinate 1 has dose loading of 10 cycles, full dose up to cycle 104 then 0 cost after:
  cost_during_loading <- consolidated_trace$split_pop$L2_on[,"L2_on_split10_before"] * mean(cost_over_time[1:10])
  full_dose_cost      <- (consolidated_trace$split_pop$L2_on[,c("L2_on_split104_before")] - consolidated_trace$split_pop$L2_on[,c("L2_on_split10_before")]) * mean(cost_over_time[11:104])
  
  # Add them up:
  combinate_1_cost <- cost_during_loading + full_dose_cost
  
  # Compare:
  sum(combinate_1_cost)
  sum(cost_per_cycle)
  
  # Very similar, but not the same. Wonder if it's a reasonable approximation?
  # 
  # 
  # % error:
  ((sum(combinate_1_cost)-sum(cost_per_cycle))/sum(cost_per_cycle))*100
  
  # 0.823% difference between doing it the simple way and fully multiplying out.
  # 
  # Is it worth the massive increase in computations?
  
  
}
