
# this script is designed to run through the progression of more and more complicated
# Markov models, to ultimately arrive at what we need for this project (the most complex case).
# This final case is:
# 
#  Extended Markov model - time-varying transitions and multiple time-varying tunnel states of time-varying duration
#  
#  - transition probability matrix for each model cycle to represent MSM results (???)
#  - Including multiple tunnel states between the main health states (treatment lines)
#  - Tunnel states change length per model cycle so sometimes its 2 cycles, sometimes its 3 kind of thing
#   - depends on the difference between time to next treatment (TTNTx) and time to discontinuation (TTDisc)
#   - both of these curves change over time per extrapolations, so tunnel changes length every cycle
#   - possibility that mortality inside tunnel is dependent on multiple things:
#     - The time that the cohort slice entered the tunnel (i.e., t_enter)
#     - The time inside the tunnel given t_enter (t_tun | t_enter)
#   - Exiting tunnel is only possible for 2 reasons
#     - Death per (OS | t_tun, t_enter)
#     - Made it to the end of the tunnel and go to the next treatment line
#   
# 


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(gtools)

# Simplest Markov model ---------------------------------------------------

# Simply repeatedly multiplying vector by matrix:

p <- c(1,0,0)

tpm <- matrix(
  c(0.99, 0.005 , 0.005,
    0.001, 0.9   , 0.099,
    0    , 0     , 1),
  nrow = 3,
  byrow = TRUE
)

res_as_list <- Reduce(
  x = 1:1000,
  init = p,
  accumulate = TRUE,
  f = function(previous_cycle_result, cycle_number) {
    previous_cycle_result %*% tpm
  }
)

res <- do.call(rbind,res_as_list)

# This is the trace, showing the proportion of people in each state at each time:
res

colnames(res) <- c("sick", "sicker", "dead")


res <- as.data.table(res)

res$cycle <- 0:(dim(res)[1]-1)

plot_res <- melt(res,id.vars = "cycle",variable.name = "state",value.name = "pop")



p1 <- ggplot(plot_res,aes(x = cycle, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")


# Time-varying tpms -------------------------------------------------------

# Just for an example, make some probabilistic draws

row_1 <- gtools::rdirichlet(1000,round(c(0.99 , 0.005, 0.005)*1000))
row_2 <- gtools::rdirichlet(1000,round(c(0.001, 0.9   , 0.099)*1000))
row_3 <- c(0,0,1)

tpm <- lapply(1:1000, function(cyc) {
  matrix(
    c(
      row_1[cyc,],
      row_2[cyc,],
      row_3
    ),
    nrow = 3,
    byrow = TRUE
  )
})


res_as_list <- Reduce(
  x = 1:1000,
  init = p,
  accumulate = TRUE,
  f = function(previous_cycle_result, cycle_number) {
    previous_cycle_result %*% tpm[[cycle_number]]
  }
)

res <- do.call(rbind,res_as_list)

# This is the trace, showing the proportion of people in each state at each time:
res

colnames(res) <- c("sick", "sicker", "dead")


res <- as.data.table(res)

res$cycle <- 0:(dim(res)[1]-1)

plot_res <- melt(res,id.vars = "cycle",variable.name = "state",value.name = "pop")



p2 <- ggplot(plot_res,aes(x = cycle, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")



# Fixed-time TPMs ---------------------------------------------------------

tpm <- list(
  matrix(
    c(0.99, 0.005 , 0.005,
      0.001, 0.9   , 0.099,
      0    , 0     , 1),
    nrow = 3,
    byrow = TRUE
  ),
  matrix(
    c(0.98, 0.015 , 0.005,
      0.001, 0.95 , 0.049,
      0    , 0     , 1),
    nrow = 3,
    byrow = TRUE
  )
)

switch_time <- 25


res_as_list <- Reduce(
  x = 1:1000,
  init = p,
  accumulate = TRUE,
  f = function(previous_cycle_result, cycle_number) {
    
    if(cycle_number < switch_time) {
      previous_cycle_result %*% tpm[[1]]
    } else {
      previous_cycle_result %*% tpm[[2]]
    }
  }
)

res <- do.call(rbind,res_as_list)

# This is the trace, showing the proportion of people in each state at each time:
res

colnames(res) <- c("sick", "sicker", "dead")


res <- as.data.table(res)

res$cycle <- 0:(dim(res)[1]-1)

plot_res <- melt(res,id.vars = "cycle",variable.name = "state",value.name = "pop")



p3 <- ggplot(plot_res,aes(x = cycle, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")


p1
p2
p3
