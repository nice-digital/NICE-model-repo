# Shift a vector and pad with zeroes

# takes a vector of n values and turn into a matrix
# eg x = {1,2,3,4,5}
# output =
# 1 0 0 0 0 
# 2 1 0 0 0
# 3 2 1 0 0
# 4 3 2 1 0
# 5 4 3 2 1 
f_shift_and_pad <- function(x) {
  #little trick to allow use of apply function:
  #add column number as first row of the matrix
  shifted        <- matrix(0,nrow = length(x) + 1, ncol = length(x))
  shifted[1,]  <- seq(1:ncol(shifted))
  shifted[-1,] <- x
  
  shifted <- apply(shifted, 2, function(x) shift(x[-1], x[1]-1, fill=0))
  
  shifted
}