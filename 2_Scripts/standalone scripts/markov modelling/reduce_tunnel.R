
# Stand-alone script which was used during the model development stages. Read 
# at your leisure. This shows some of the progression towards using sparse matrix
# block diagonal multiplications to compute the markov model.


# Libraries ---------------------------------------------------------------

library(gtools) # for dirichlet
library(data.table)
library(tidyverse)
library(collapse)
library(bigalgebra)

# Example data ------------------------------------------------------------

tp_top    <- rdirichlet(99,c(95,2,1))
tp_mid    <- rdirichlet(99,c(0,93,7))
tp_bottom <- c(0,0,1)

# List of TPMs with some random variation cycle to cycle
tp <- lapply(1:dim(tp_top)[1], function(cycle) {
  matrix(
    c(
      tp_top[cycle,],
      tp_mid[cycle,],
      tp_bottom
    ),
    nrow = 3,
    byrow = TRUE
  )
})


# Reduce through TPMs

prev_cycle <- c(1,0,0)
tun_length <- 2
tun_death  <- c(0,0.15,0.1)

mk_top_row <- matrix(c(1,0,0),nrow=1)
tun_top_row <- matrix(c(0,0,0),nrow=1)



# Attempts at cracking it -------------------------------------------------



# ~ What a tunnel state does ----------------------------------------------

# A tunnel state pulls patients out of the state residency into a special state
# with a fixed duration. This state can only be exited via death or making it to the end.
# 
# In a Markov model, this means that simply repeatedly applying matrix multiplication
# is incorrect, since transition probabilities are proportional, so you can't work out
# the trace with a tunnel without building the tunnel into the multiplication steps. 

# Let's go ultra simple and try to work through a tunnel state that interrupts 1->2 for 
# THREE CYCLES, and kills off 10% each cycle in there, then feeding into state 2 afterwards:

tpm <- tp[[1]]

p <- c(1,0,0)

tun_dur <- 3
tun_qx <- c(0.1,0.1,0.15)

# Ignoring the tunnel
tr_noTun <- Reduce(
  x = 1:100,
  init = p,
  accumulate = TRUE,
  f = function(previous_result, cycle_number) {
    previous_result %*% tpm
  }
)

tr_noTun <- do.call(rbind,tr_noTun)

head(tr_noTun)

# Including the tunnel but fixing the duration at 3 cycles:

TH <- 100

p <- matrix(nrow = TH,
            ncol = 3 + 3)    # Note the extra state for "inside the tunnel" allowing each row in p to sum to 1 is not here!!!
p[1,] <- c(1,rep(0,3 + 3 - 1))       # Note the extra state for "inside the tunnel" allowing each row in p to sum to 1 is not here!!!

p[is.na(p)] <- 0

colnames(p) <- c("1L", "tun1", "tun2", "tun3", "2L", "dead")


# Note that tun_trace lags behind p by 1 cycle. To adjust for that, the code starts
# in the 2nd row

# OS within the tunnel - if tunnel qx is time-varying simply replace with S(t) vector.
# done that here as an example.
tun_os <- cumprod(c(1,(1-tun_qx)))

# cumulative mort in the tunnel:
tun_mort <- 1-tun_os

# so for each batch of patients that go into the tunnel, we have an overall survival
# and cumulative mortality line (just multiply the above by entry pop).
# 
# With that we can track populations inside the tunnel, and calculate those cohorts
# all at once


# Next, we can't just do the above Reduce approach of iteratively multiplying the 
# population in cycle-1 by the tpm to get the new population, since there are patients
# getting taken out of it. Instead, what we could do is expand the TPM:
tpm <- as.data.table(tpm)

tpm$V4 <- 0
tpm$V5 <- tpm[,V2]
tpm$V6 <- tpm[,V3]
tpm$V2 <- 0
tpm$V3 <- 0

# Set up the tunnel
tpm[1,"V2"] <- tpm[1,"V5"]
tpm[1,"V5"] <- 0

tpm2 <- do.call(rbind,lapply(1:nrow(tpm), function(rw) {
  if (rw == 1 | rw == 3) {
    return(as.numeric(tpm[rw,]))
  } else if (rw == 2) {
    do.call(rbind,list(
      rep(0,ncol(tpm)),
      rep(0,ncol(tpm)),
      rep(0,ncol(tpm)),
      as.numeric(tpm[rw,])
    ))
  }
}))


dimnames(tpm2) <- list(c("1L", "tun1", "tun2", "tun3", "2L", "dead"),c("1L", "tun1", "tun2", "tun3", "2L", "dead"))

tpm2["tun1","dead"] <- tun_qx[1]
tpm2["tun1","tun2"] <- 1-tun_qx[1]
tpm2["tun2","dead"] <- tun_qx[2]
tpm2["tun2","tun3"] <- 1-tun_qx[2]
tpm2["tun3","dead"] <- tun_qx[3]
tpm2["tun3","2L"]   <- 1-tun_qx[3]

# Now that we have expanded the transition probability matrix to incorporate the
# tunnel state directly, it should be the case that we can now simply multiply it out!

trace_tun <- matrix(
  unlist(Reduce(
    x          = 2:TH,
    init       = p[1,],
    accumulate = TRUE,
    f = function(prev_pop, cyc) {
      prev_pop %*% tpm2
    }
  )),
  nrow = TH,
  byrow = TRUE,
  dimnames = list(NULL,c("1L", "tun1", "tun2", "tun3", "2L", "dead"))
)

# Success! this is a 3 cycle tunnel state. Let's write it into a function where we can defined a FIXED tunnel duration.


# Making fixed tunnel into a function -------------------------------------

TH <- 100

state_names <- c(paste0("line_",1:5),"bsc","dead")

bl_pop        <- c(1,rep(0,length(state_names)-1))
names(bl_pop) <- state_names

tp_1L   <- rdirichlet(TH,c(95,2 ,0 ,0 ,0 ,3 ,4))
tp_2L   <- rdirichlet(TH,c(0 ,93,6 ,0 ,0 ,3 ,5))
tp_3L   <- rdirichlet(TH,c(0 ,0 ,92,5 ,0 ,3 ,6))
tp_4L   <- rdirichlet(TH,c(0 ,0 ,0 ,91,3 ,3 ,7))
tp_5L   <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,90,3 ,8))
tp_BSC  <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,0 ,85,9))
tp_dead <- c(0, 0, 0, 0, 0, 0, 1)

# List of TPMs with some random variation cycle to cycle to simulate time-varying
# transition probabilities between states. Patients can't go backwards in line, and
# can't skip lines except for going straight to bsc (i.e. choosing no treatment/palliative care)
tp <- lapply(1:TH, function(cycle) {
  matrix(
    c(tp_1L[cycle,], 
      tp_2L[cycle,], 
      tp_3L[cycle,], 
      tp_4L[cycle,],
      tp_5L[cycle,],
      tp_BSC[cycle,],
      tp_dead
    ),
    nrow = 7,
    byrow = TRUE,
    dimnames = list(state_names,state_names)
  )
})


# Other parameters for our function:

tun_pos  <- 1
tun_len  <- 3
tun_mort <- c(0.1,0.08,0.15)
tpm      <- tp[[1]]

# Now for the function! we're going to make it possible to expoand the transition probability
# matrix at any position and insert a tunnel state there. This function will
# expand ONE TPM. you can then use lapply to apply it to all of them in a list
# or apply to time invariant as you wish!

f_mk_FixedTunnelMarkovExpand <- function(tpm, tun_pos, tun_len, tun_mort) {
  
  # always start with validation:
  # Determine whether tun_mort is a single value or a vector:
  if (!length(tun_mort) %in% c(1,tun_len)) stop("tunnel mort must be either 1 number, or a vector of the same length as tun_len!")
  if (tun_len == 1) stop("tun_len should be at least 2")
  
  # Establish the entrance to the tunnel. For example if tun_pos is 1 then
  # it's the transition tpm[1,2], if it's 3 then it's tpm[3,4]. Move that TP from its 
  # place in aft to a separate object, and then set it to 0 in aft.
  
  pr_enter_tun           <- tpm[tun_pos,tun_pos+1]
  tpm[tun_pos,tun_pos+1] <- 0
  
  # Let's make a way to appropriately name this tunnel (i.e. tun_ prefix and numeric positoins)
  
  if (!is.null(dimnames(tpm)[[1]])) {
    tunnel_start_state <- dimnames(tpm)[[1]][tun_pos]
    tun_lab <- paste0("tun_",tunnel_start_state,"_",1:tun_len)
  } else {
    tun_lab <- paste0("tun_",tun_pos,"_",tun_pos+1,"_",1:tun_len)
  }
  
  # So now we need to slice up the matrix and add some zeros. We're going to 
  # insert a tunnel state to tpm that STARTS WHEN EXITING tun_pos
  # state (which allows you to enter the number 1). 
  # This tunnel is tun_len cycles long and has tun_mort associated
  # mortality whilst inside the tunnel
  # 
  # This means that we should first add the zero rows in the right place, then
  # the zero columns (it's less faff that way)
  
  bef <- tpm[1:tun_pos,,drop=FALSE]
  aft <- tpm[(tun_pos+1):dim(tpm)[1],,drop=FALSE]
  
  empty_row <- rep(0,dim(tpm)[2])
  names(empty_row) <- dimnames(bef)[[2]]
  empty_rows <- do.call(rbind,lapply(1:tun_len, function(tun_cyc) empty_row))
  rownames(empty_rows) <- tun_lab
  
  tpm <- do.call(rbind,list(
    bef,
    empty_rows,
    aft
  ))
  
  # Now do the same thing for columns too:
  bef <- tpm[,1:tun_pos,drop=FALSE]
  aft <- tpm[,(tun_pos+1):dim(tpm)[2],drop=FALSE]
  
  empty_col <- rep(0,dim(tpm)[1])
  names(empty_col) <- dimnames(bef)[[1]]
  empty_cols <- do.call(cbind,lapply(1:tun_len, function(tun_cyc) empty_col))
  colnames(empty_cols) <- tun_lab
  
  tpm <- do.call(cbind,list(
    bef,
    empty_cols,
    aft
  ))
  
  # Great! now we have a larger transition probability matrix which has a space
  # for all of the possible transitions, including those for the tunnel state. 
  # From here, we can systematically populate the transition probabilities within
  # the tunnel, using tun_mort, depending on whether it's a single number or
  # vector of numbers
  
  # First off, put the transition probability we pulled out at the start of the
  # tunnel:
  tpm[tun_pos,tun_pos+1] <- pr_enter_tun
  
  if (length(tun_mort) == 1) {
    tun_mort <- rep(tun_mort,tun_len)
  }
  
  # Now cycle through the cycles in the tunnel, informing the correct cell in
  # the tpm with the probability of staying in the tunnel or death, such that the
  # final tpm returned from this function has rows that sum to 1 and a tunnel in it!
  return(
    Reduce(
      x = 1:tun_len,
      accumulate = FALSE,
      init = tpm,
      f = function(prev, index) {
        # the index'th cycle within the tunnel - patient can die or transition into
        # the next cycle of the tunnel (or out of the tunnel)
        prev[tun_pos + index, ncol(prev)]          <- tun_mort[index]
        prev[tun_pos + index, tun_pos + index + 1] <- 1-tun_mort[index]
        return(prev)
      }
    )
  )
  
}


# now let's test this function:
tpm2 <- f_mk_FixedTunnelMarkovExpand(
  tpm      = tpm,
  tun_pos  = 2,
  tun_len  = 2,
  tun_mort = c(0.15,0.35)
)


bl_pop <-
  c(
    line_1 = 1,
    line_2 = 0,
    tun_line_2_1 = 0, 
    tun_line_2_2 = 0,
    line_3 = 0,
    line_4 = 0,
    line_5 = 0,
    bsc = 0,
    dead = 0
  )

test_trace <- do.call(
  rbind,
  Reduce(
    x = 1:(100-1),
    init = bl_pop,
    accumulate = TRUE,
    f = function(prev, cycle) {
      prev %*% tpm2
    }
  )
)


rowSums(test_trace)

test_trace <- as.data.table(test_trace)
test_trace$t <- 0:(nrow(test_trace)-1)


test_trace <- melt(test_trace,id.vars = "t",variable.name = "state", value.name = "pop")

ggplot(test_trace, aes(x = t, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")

# WOOP! all rows sum to 1 and it has a tunnel state in it which we can place anywhere
# and change the length of.

# Problem is, the same cohort should consistently be exposed to those TPs, so this
# approach has the weakness that tun_mort can't be dependent on absolute time since
# baseline, just on time in tunnel given line before entering the tunnel. 


# Now let's try adding another tunnel to this tunnel-state matrix!

tpm3 <- f_mk_FixedTunnelMarkovExpand(
  tpm = tpm2,
  tun_pos = 5,
  tun_len = 3,
  tun_mort = c(0.1,0.2,0.15)
)
tpm4 <- f_mk_FixedTunnelMarkovExpand(
  tpm = tpm3,
  tun_pos = 9,
  tun_len = 2,
  tun_mort = c(0.05,0.1)
)
tpm5 <- f_mk_FixedTunnelMarkovExpand(
  tpm = tpm4,
  tun_pos = 12,
  tun_len = 4,
  tun_mort = c(0.05,0.09,0.1,0.15)
)

# Now that we've got a fully fleshed out transition probability matrix:

# this will create the vector of column names and print it to the console for your
# convenience. Instead we do it programatically:
# dput(colnames(tpm5))

st_nam        <- colnames(tpm5)
bl_pop        <- rep(0,length(st_nam))
names(bl_pop) <- st_nam
bl_pop[1]     <- 1

test_trace <- do.call(
  rbind,
  Reduce(
    x = 1:(100-1),
    init = bl_pop,
    accumulate = TRUE,
    f = function(prev, cycle) {
      prev %*% tpm5
    }
  )
)

# All rows still sum to 1 :)
rowSums(test_trace)

test_trace   <- as.data.table(test_trace)
test_trace$t <- 0:(nrow(test_trace)-1)
test_trace   <- melt(test_trace,id.vars = "t",variable.name = "state", value.name = "pop")
ggplot(test_trace, aes(x = t, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")


# As there are so many individual time points within tunnels, we might just want
# to see the full tunnel population for each tunnel. This is a bit fiddly to do
# programatically but worth it


test_trace <- do.call(
  rbind,
  Reduce(
    x = 1:(100-1),
    init = bl_pop,
    accumulate = TRUE,
    f = function(prev, cycle) {
      prev %*% tpm5
    }
  )
)

orig_nams <- colnames(tpm)
names(orig_nams) <- orig_nams

full_nam <- colnames(tpm5)
tun_nams <- full_nam[grep("tun_",full_nam)]

tun_pops <- lapply(orig_nams, function(orig_state) {
  
  red_nam <- full_nam[grep(orig_state,full_nam)]
  
  # If there's no tunnel state following this state just return a NULL value
  if(length(red_nam) == 1) return(NULL)
  
  # If there is a tunnel, we now know the columns to add up to get our total
  # tunnel state residency by time t
  col_sel <- red_nam[which(red_nam != orig_state)]
  
  rowSums(test_trace[,col_sel])
  
})

# Reduce down to just the tunnel populations
tun_pops <- tun_pops[which(!sapply(tun_pops,is.null))]

# rename them appropriately
names(tun_pops) <- paste0("tun_",names(tun_pops))

# turn them into a data.table
tun_pops <- as.data.table(tun_pops)

# remove all tunnels breakdowns from the trace, and add the summed ones back in:
trace2 <- cbind(test_trace[,orig_nams],tun_pops)

# All still sum to 1 in rows :)
rowSums(trace2)


# Now we can plot:
trace2   <- as.data.table(trace2)
trace2$t <- 0:(nrow(trace2)-1)
trace2   <- melt(trace2,id.vars = "t",variable.name = "state", value.name = "pop")
ggplot(trace2, aes(x = t, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")

# So, we can now include any number of tunnel states of any length to any TPM.
# We're getting somewhere now!


# the big one -------------------------------------------------------------

# ok now we want to do time varying transition probabilities with different exit
# time for tunnel per cycle, so that the mortality in the tunnel stays the same, but 
# its length can move around from cycle to cycle

vec_tun_len <- sample(1:4,TH,TRUE)
max_tun_len <- max(vec_tun_len)

# To make it even more complicated, the tunnel mortality can vary over time IF 
# mortality in the tunnel is dependent on both t and t_tun (i.e., time inside
# tunnel):

# Just for an example I'm going to generate a bunch of beta draws, separated by cycle.
# The reason it's a list is because each element is a different length, so you can't
# just bind them together. 
# 
# Programming wise, it's just a case of l_tun_os[[cohort]][cohort_cycle,] instead of l_tun_os[cycle,]
# so it's trivial to adapt to it. in other words, for each cohort there's a totally 
# different set of TPs to the death state whilst within the tunnel, given
# absolute time AND time since entering the tunnel (i.e., p(death | t, t_tun)).
# I cannot see why any further level of branching would be required beyond this,
# and in most CE modelling cases I can see probability of death being only a function 
# of time in tunnel OR absolute time. In both these cases, one would simply adapt
# the below lapply to appropriately generate these death probability matrices.
# 
# To clarify - here we're saying that mortality whilst within tunnel depends on
# both the cycle you went into the tunnel AND how long you've been in the tunnel given
# the time that you entered the tunnel. This accounts for "counting back"

l_tun_os <- lapply(1:length(vec_tun_len), function(cyc) {
  # Make a matrix which contains our (example) random numbers for probability of death
  # in tunnel. obviously these will be different when we're using real data. Transpose
  # it so we can split it by row afterwards. Speeds things up later so we may as well
  # do it now.
  
  tun_len_this_cohort <- vec_tun_len[cyc]
  
  tp_i <- t(matrix(
    data = rbeta((TH - cyc + 1) * tun_len_this_cohort,10,100),
    nrow = TH - cyc + 1,
    ncol = tun_len_this_cohort,
    byrow = TRUE,
    dimnames = list(NULL,NULL)
  ))
  tp_i <- split(tp_i, rep(1:ncol(tp_i), each = nrow(tp_i)))
  names(tp_i) <- NULL
  
  # extend the values with 100% death probability to capture the impossible people
  # who shouldn't be there anyway. This way will force rowsums to be 
  
  if (max_tun_len - tun_len_this_cohort > 0) {
    tun_ext_len <- rep(1,max_tun_len - tun_len_this_cohort)
    return(lapply(tp_i, function(this_cycle_this_cohort) c(this_cycle_this_cohort, tun_ext_len)))
  } else {
    return(tp_i)
  }
})

# We now refer to each death vector in a very simple way: l_tun_os[[cohort]][[model_cycle]]
# This gets us the p(death | t, t_tun) where this time t_tun is the time at which
# that cohort ENTERED the tunnel.
# 
# This allows rule-based extrapolations (e.g. if there's an OS line for tunnel
# before and after 6 months from baseline) or relative rule-based extrapolations
# (e.g. if 6 months after entering the previous line get this OS, otherwise that one)
# allowing for really complicated interactions between lines and tunnels to be
# accomplished whilst assuring that no patients get lost.

# To recap, we have a tunnel state that changes in duration every cycle, and in
# mortality rate as well (so that each cohort has its own unique OS line!)

# So let's start with tp, the set of TPs we made earlier that goes to the TH

# For each model cycle we need to make a vector of EXPANDED TPMs that has TH elements in it.
# This should have a tunnel which is max_tun_len in length, but has TP of 0 to the
# tunnel cycles which lie beyond the extent for that cohort's tunnel. To illustrate,
# let's imagine the following:
# 
# - Patients in cycle 1 have a tunnel between 1 and 2 that is 2 cycles in length
# - patients in cycle 2 have a tunnel between 1 and 2 that is 3 cycles in length
# - the value of max_tun_len is 4 cycles (i.e. the longest tunnel in the time horizon)
# 
# So, for our cycle 1 people, we have to make a tunnel with 4 cycles, but make it
# impossible to go from tun_line_1_2 to tun_line_1_3.
# 
# We actually do not NEED to set the probabilities in subsequent cycles to 0, we
# just need to make it impossible to continue down the tunnel. this means simply
# replacing the 1-death probability of going to the next tunnel cycle with 0 and
# put that number in the probability to go to the next treatment line!
# 
# 

# As we have different lengths of tunnel and different mortality within tunnel,
# we can now generate our huge blob of expanded TPMs!
# 
# Here, let's put in our tunnel between lines 3 and 4. this lets enough patients
# slip through other states so that they're not ALL going in the tunnel.

ex_tun_pos <- 3

tp_tun_vary <- lapply(1:TH, function(this_cycle_cohort) {
  
  # remember that here we're going to generate TH matrices for EACH cycle cohort!
  
  # Extend the length of the tunnel to match the max and kill everyone that's doing
  # the impossible - this will mean that if it's wrong, the rows won't sum to 1 at
  # the end!
   
  # The length of the tunnel for this cohort in cycles. This lets us know where
  # to move the exit probability in the matrix from. Where to move it to is a function
  # of this value AND the positioning of the tunnel (i.e. tun_pos when expanding)
  this_tun_mort <- l_tun_os[[this_cycle_cohort]]
  
  # So, we have a list of vectors. Each one of these is the conditional death probability
  # patients inside the tunnel face at each model cycle. This can count back 
  # because when this list was made it can be made to count back, so it's all
  # led by the data by this point.
    
  lapply(this_cycle_cohort:TH, function(this_cohort_this_cycle) {
    
    # For this cycle for this cohort, Expand the tp matrix (which ONLY depends on t)
    # and insert the correct tunnel mortality estimates (which depend on t and t|t_tun)
    
    which_tun_mort <- this_cohort_this_cycle - this_cycle_cohort + 1
    
    tpm <- f_mk_FixedTunnelMarkovExpand(
      tpm      = tp[[this_cohort_this_cycle]], # for this model cycle (for this absolute t)
      tun_pos  = ex_tun_pos,
      tun_len  = max_tun_len,
      tun_mort = this_tun_mort[[which_tun_mort]]
    )
    
    # Right, now we have a TPM which kills everyone when they get to the nth cycle
    # inside the tunnel, if the tunnel length is less than max_tun_len. Therefore
    # we simply need to make the probability of going into that next tunnel cycle
    # 0, and the probability of going out of the tunnel to the next treatment line
    # equal to the probability that we did have for going to the next tunnel cycle.
    # 
    # To simplify to mechanics 
    #   - the row number is ex_tun_pos 
    #   - the origin column is ex_tun_pos + length(this_tun_mort[[this_cohort_this_cycle]]) 
    #   - the destination column is ex_tun_pos + max_tun_len + 1
    
    tun_len_this_cohort <- vec_tun_len[this_cycle_cohort]
    
    if (max_tun_len - tun_len_this_cohort > 0) {
      row      <- ex_tun_pos + tun_len_this_cohort
      col_from <- ex_tun_pos + tun_len_this_cohort + 1 
      col_to   <- ex_tun_pos + max_tun_len         + 1
      
      p_ft              <- tpm[row,col_from]
      tpm[row,col_from] <- 0
      tpm[row,col_to]   <- p_ft
    }
    
    
    # Great, so it's now impossible to go too far into the tunnel as the transition
    # into that tunnel cycle (and therefore any that follow it) has 0% probability.
    # Patients that would have gone into that next tunnel (i.e. haven't died) are 
    # now going out of the tunnel directly to the next treatment line.
    
    return(tpm)
    
  })
  
})


# Simplified version for programming --------------------------------------

# Version with cycle after line 1 for varying cycles. easier to run through
# stuff if it's from the first state:

ex_tun_pos <- 1
tp_tun_vary <- lapply(1:TH, function(this_cycle_cohort) {
  this_tun_mort <- l_tun_os[[this_cycle_cohort]]
  lapply(this_cycle_cohort:TH, function(this_cohort_this_cycle) {
    which_tun_mort <- this_cohort_this_cycle - this_cycle_cohort + 1
    tpm <- f_mk_FixedTunnelMarkovExpand(
      tpm      = tp[[this_cohort_this_cycle]], # for this model cycle (for this absolute t)
      tun_pos  = ex_tun_pos,
      tun_len  = max_tun_len,
      tun_mort = this_tun_mort[[which_tun_mort]]
    )
    tun_len_this_cohort <- vec_tun_len[this_cycle_cohort]
    if (max_tun_len - tun_len_this_cohort > 0) {
      row      <- ex_tun_pos + tun_len_this_cohort
      col_from <- ex_tun_pos + tun_len_this_cohort + 1 
      col_to   <- ex_tun_pos + max_tun_len         + 1
      
      p_ft              <- tpm[row,col_from]
      tpm[row,col_from] <- 0
      tpm[row,col_to]   <- p_ft
    }
    return(tpm)
  })
})

tun_from_state  <- "line_1"
tun_entry_state <- "tun_line_1_1"

empty_row <- tp_tun_vary[[1]][[1]][1,]
empty_row <- sapply(empty_row,function(x) 0)

bl <- tp_tun_vary[[1]][[1]][1,]
bl <- sapply(bl,function(x) 0)
bl[1] <- 1

# calculate the population that never goes in the tunnel.
trace_no_tun <- Reduce(
  x = 1:(TH-1),
  accumulate = TRUE,
  init = bl,
  f = function(prev,cyc) {
    this_tp <- tp_tun_vary[[1]][[cyc]]
    this_tp[tun_from_state,tun_entry_state] <- 0
    prev %*% this_tp
  }
)
trace_no_tun <- matrix(
  unlist(trace_no_tun),
  nrow = TH,
  byrow = TRUE,
  dimnames = list(NULL,dimnames(tp_tun_vary[[1]][[1]])[[2]])
)

non_tun_pop <- 1-rowSums(trace_no_tun)
tun_entry   <- non_tun_pop - shift(non_tun_pop,1,0)

empty_trace <- matrix(
  rep(0,TH*length(bl)),
  nrow = TH,
  dimnames = list(NULL,names(bl))
)

# This is the tricky-dicky bit - we create TH traces, one for each cohort
# entering the tunnel, and then we extrapolate the state residency of each one
# of them one by one
tun_traces <- lapply(1:TH, function(cyc) {
  pop <- empty_trace
  pop[cyc,tun_entry_state] <- tun_entry[cyc]
  pop
})

# Now, we have the trace without the tunnel state in it (losing people), and 
# individual traces for every entry point to the tunnel. We can now apply the
# correct TPMs to each of these traces to calculate them. the sum of all the matrices
# is then the full trace with varying-length-time-varying tunnel state!

tun_list <- lapply(1:TH, function(tunnel_entry_time) {
  
  # pull out the cohort that entered the tunnel at cycle tunnel_entry_time
  cohort <- tun_traces[[tunnel_entry_time]]
  
  if (tunnel_entry_time == TH) return(cohort)
  
  # Calculate a trace for this population, using tp_tun_vary[[tunnel_entry_time]][1:(TH-tunnel_entry_time)]
  Reduce(
    x = 1:(TH-tunnel_entry_time),
    init = cohort,
    accumulate = FALSE,
    f = function(p, c) {
      p[tunnel_entry_time + c,] <- p[tunnel_entry_time + c - 1,] %*% tp_tun_vary[[tunnel_entry_time]][[c]]
      p
    }
  )
})

# Unit test - all summing to 1????
rowSums(Reduce(`+`,tun_list)) + rowSums(trace_no_tun)

# I....I think I did it....I think I just solved time varying duration varying tunnel states...

trace_list <- list(
  trace_no_tun,
  tun_list
)

full_trace <- Reduce(`+`,tun_list) + trace_no_tun

rowSums(full_trace)

full_trace   <- as.data.table(full_trace)
full_trace$t <- 0:(nrow(full_trace)-1)
full_trace   <- melt(full_trace,id.vars = "t",variable.name = "state", value.name = "pop")
ggplot(full_trace, aes(x = t, y = pop, colour = state)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "bottom")




# ~~ Two tunnels? ---------------------------------------

# ok this is all well and good for 1 tunnel, but what about when there are further tunnels, which also vary by 
# duration depending on cycle of entry, and also have OS depending on cycle of entry? What do we do then?
# 
# We might have to go row by row:


# Let's use tp again so we have time varying TPMs

str(tp,list.len = 2)

# Let's imagine we have 2 tunnels, one varying between 1 and 4 cycles, and the other
# varying between 2 and 8 cycles in length, depending on the cycle number of entry, 
# which I'm going to call e forthwith.
tun_len <- list(
  sample(1:4,TH,TRUE),
  sample(2:8,TH,TRUE)
)

# We want the ORIGINAL position of each tunnel. To test robustness, we're going to say
# that the tunnels are upon exiting 1L and 3L:
tun_pos <- list(1,3)


# 2 layers of nesting here. first is the tunnel, t, and the other is cycle of entry e
# h_tun[[t]][[e]] is the hazard to be applied to the cohort that entered tunnel t
# at cycle e, for each cycle inside the tunnel, d. To recap:
# 
#  - c is model cycle
#  - t is tunnel
#  - e is cycle of entry into tunnel (e.g. cohort entered tunnel t in c = 2, therefore e = 2)
#  - d is time inside tunnel (e.g. cohort entered tunnel in cycle 2, and it is now cycle 5, d = 3)
# 
# The 1's are added to keep them all the same length - this helps when creating the TPMs, otherwise
# the population vectors produced at each model cycle will have different lengths, which will make 
# it more annoying to bind it all together. Better to keep it tidy in this way instead!
# 
# Obviously, these are just beta draws and instead will have to come from some data
# or analysis for the CE model. This could potentially be slices from an OS extrapolation
# for being between treatment lines (censoring for going onto next treatment to isolate
# the hazard of just death whilst waiting for next treatment line)
# 
# if, for example the expected time between discontinuing 1L and starting 2L is
# 5 weeks for those discontinuing in cycle 1, then use the hazard for those 5 weeks starting
# from c=5? Hazard for those with e=6 should be h(death | t = 1, c = 6, e = 6) then
# h(death | t = 1, c = 7, e = 6) and so on to make the OS line for the cohort with
# that value of e, correct? We would then get (TH-1) OS lines per t, like the below:
# 
h_tun <- lapply(
  X = tun_len,
  FUN = function(D) {
    lapply(1:length(D), function(c) {
      d_given_te <- D[c]
      tp_i       <- rbeta(d_given_te,floor(rnorm(1,10,2)),floor(rnorm(1,100,2)))
      return(c(tp_i, rep(1,max(D) - d_given_te)))
    })
  }
)


# Now we perform the transition probability matrix expansion. To walk through it:
# 
#  - the cohort with e=1 has TH TPMs, each applying the mortality in h_tun for e=1 (i.e. lapply(h_tun, function(x) x[[1]]))
#  - the cohort with e=2 has TH-1 TPMs, each applying the mortality in h_tun for e=2 (i.e. lapply(h_tun, function(x) x[[2]]))
#  - ...
#  - the cohort with e=TH-1 has 1 TPM, applying the mortality in h_tun for e=TH-1 (i.e. lapply(h_tun, function(x) x[[TH]]))
# 
# In that sense, it's kind of like a pyramid or a right-angled triangle facing point down, with the
# number of TPMs required reducing as we go further into the model
# 

n_tun <- length(tun_pos)

# Work out what the final positions in the TPMs will be for our 2 tunnels
tun_positions <- Reduce(
  x = 1:length(tun_pos),
  init = tun_pos,
  accumulate = FALSE,
  f = function(prev, tun) {
    if (tun == 1) {
      return(prev)
    } else {
      max_len <- max(tun_len[[tun-1]])
      prev[[tun]] <- prev[[tun - 1]] + max_len + tun_pos[[tun]] - 1
      return(prev)
    }
  }
)


tp_triangle <- lapply(1:TH, function(e) {
  
  # Get mortality and tunnel length for cohorts entering tunnels at cycle c (i.e. where e = c)
  this_tun_mort <- lapply(h_tun, function(t) t[[e]])
  this_tun_len  <- lapply(tun_len, function(t) t[[e]])
  
  # Now, d is between e and the time horizon TH. This is our iterator, as it tells us
  # which TPM to use (i.e. transitions for non tunnel states for cycle c)
  d <- e + 0:(TH - e - 1)
  
  # Now, iterating along d, take the correct TPM, and expand it to include all the
  # tunnels, applying this_tun_mort and making the tunnel exit adjustment for this_tun_len.
  # This will result in length(d) TPMs, all for the correct c and e for all t. Note, The
  # index is actually c here as it's that absolute time step in the model!
  lapply(d, function(c) {
    Reduce(
      x = 1:length(this_tun_mort),
      init = tp[[c]],
      accumulate = FALSE,
      f = function(prev_tpm, tun) {
        
        # expand the tpm for model cycle c using mortality for entry cycle e along d
        tpm_expanded <- f_mk_FixedTunnelMarkovExpand(
          tpm      = prev_tpm,
          tun_pos  = tun_positions[[tun]],
          tun_len  = max(tun_len[[tun]]),
          tun_mort = this_tun_mort[[tun]]
        )
        
        # Correct the matrix to make patients exit the tunnel corresponding to
        # d given e. Identify the row we're changing & which column is going where.
        row_from_to  <- tun_positions[[tun]] + this_tun_len[[tun]]
        col_from     <- tun_positions[[tun]] + this_tun_len[[tun]] + 1
        col_to       <- tun_positions[[tun]] + max(tun_len[[tun]]) + 1
        
        # Make the changes to the TPM
        tpm_expanded[row_from_to,col_to]   <- tpm_expanded[row_from_to,col_from]
        tpm_expanded[row_from_to,col_from] <- 0
        
        # it is now impossible to go beyond the end of the tunnel, and instead
        # patients that survive will go to the next treatment line. the Reduce
        # will then add the next tunnel state if there are any
        return(tpm_expanded)
      }
    )
  })
})

# To illustrate - the length of the list for each e goes down by 1
unlist(lapply(tp_triangle,length))


# tp_triangle has TH-1 TPMs for the first e, TH-2 for the 2nd and so on.

empty_pop <- unlist(lapply(as.list(tp_triangle[[1]][[1]][1,]),function(x) 0))
bl_pop    <- empty_pop
bl_pop[1] <- 1

empty_trace <- matrix(
  data = rep(empty_pop,TH),
  nrow = TH,
  byrow = TRUE,
  dimnames = list(NULL,names(empty_pop))
)

bl_trace      <- empty_trace
bl_trace[1,1] <- 1

# Right, let's try to apply these tpms now. The way you do this is you separate out
# tunnel entrants from everyone else, and put the tunnel entrants into their 
# e-based slot. To calculate a few of them through:

# apply the TPM for e=1 c=1
pop_c2 <- bl_pop %*% tp_triangle[[1]][[1]]

# These are people going into the tunnels at cycle 2 so e is 2 (i.e. c + 1)
e <- 2

# Feed the tunnel entrants for this e into their space for trace:
pop_e2  <- empty_trace
pop_e2[e,unlist(tun_positions)+1] <- pop_c2[unlist(tun_positions)+1]

# Patients entering either tunnel with e of 2 can now have a trace calculated using Reduce:

tr_e2 <- rbind(
  do.call(rbind, lapply(1:(e - 1), function(c) empty_pop)),
  do.call(
    rbind,
    Reduce(
      x = (e + 1):TH,
      init = pop_e2[e, ],
      accumulate = TRUE,
      f = function(tr, cyc) {
        tr %*% tp_triangle[[1]][[cyc - e]]
      }
    )
  )
)

# Right, that's the trace for people entering tunnels in cycle 2. let's check we're not
# losing anyone:

rowSums(tr_e2)

# Great, so we're not losing anyone and we're tracing through cohort e=2, through
# all tunnels.

# Lovely. so essentially all we have to do now is repeat the above for c = 1:TH to get
# TH Markov traces. Each time, we get our new cohort, "split them off" from the rest,
# calculate the trace for that cohort, and move onto the next cycle c. We should end up with
# TH traces, one for those that never go in a tunnel and TH-1 for people going in tunnels
# at different times
# 









# So, as before we've got a nested list with the following levels:
# 
#   1. Cycle that the cohort ENTERS the tunnel, cyc
#   2. Cycles since that cohort that ENTERED the tunnel at cyc have BEEN in the tunnel
# 
# l_tun_os2[[1]][[1]] Entered the tunnel in cycle 1, first cycle in the tunnel
# l_tun_os2[[2]][[1]] Entered the tunnel in cycle 2, first cycle in the tunnel
# l_tun_os2[[5]][[3]] Entered the tunnel in cycle 5, 3rd cycle in the tunnel
# 
# 




# Turning the example into a function -------------------------------------

# So we've cracked it I think withouth making any assumptions about where the tunnel 
# state is or how long it is. now we have two challenges remaining:
# 
#   1. turn it into a function
#   2. make it cope with multiple tunnel states
#   
#   



#' Function to compute transition probability matrices including a varying-duration
#' tunnel state which can then be used with another function to calculate the required traces 
#' 
#' @param list_tun_os a list of lists of mortality vectors. The top level is tunnel number, the 2nd level is entry cycle, the 3rd level is cycles since tunnel entry
#' @param tpm_list transition probability matrix list BEFORE tunnels inserted - one for each model cycle, or list of repeated if time invariant
#' @param init_pop (optional) initial population. will be assumed 100% in first cycle if left blank
#' 
f_mk_varyTunneltpms <- function(list_tun_os, tpm_list, tun_position, vec_tun_len, init_pop = NULL) {
  
  # Step 1 - generate the set of TPMs taking tunnel mort given time and time in tun
  
  TH <- length(vec_tun_len)
  
  if(is.null(init_pop)) {
    init_pop <- c(1,rep(0,dim(tpm_list[[1]])[2]-1))
    names(init_pop) <- dimnames(tpm_list[[1]])[[2]]
  }
  
  # The first layer of nesting is the time of entry to the tunnel, and the second layer
  # is tinme since entering the tunnel.
  
  tp_tun_vary <- lapply(1:TH, function(this_cycle_cohort) {
    
    # For this entry time to the tunnel, repeatedly expand the TPM to take into
    # account all tunnel states:
    this_tun_mort <- l_tun_os[[this_cycle_cohort]]
    lapply(this_cycle_cohort:TH, function(this_cohort_this_cycle) {
      which_tun_mort <- this_cohort_this_cycle - this_cycle_cohort + 1
      tpm <- f_mk_FixedTunnelMarkovExpand(
        tpm      = tp[[this_cohort_this_cycle]], # for this model cycle (for this absolute t)
        tun_pos  = ex_tun_pos,
        tun_len  = max_tun_len,
        tun_mort = this_tun_mort[[which_tun_mort]]
      )
      tun_len_this_cohort <- vec_tun_len[this_cycle_cohort]
      if (max_tun_len - tun_len_this_cohort > 0) {
        row      <- ex_tun_pos + tun_len_this_cohort
        col_from <- ex_tun_pos + tun_len_this_cohort + 1 
        col_to   <- ex_tun_pos + max_tun_len         + 1
        
        p_ft              <- tpm[row,col_from]
        tpm[row,col_from] <- 0
        tpm[row,col_to]   <- p_ft
      }
      return(tpm)
    })
  })
  
  return(tp_tun_vary)
  
}



# ~ Testing area ----------------------------------------------------------

TH <- 100

state_names <- c(paste0("line_",1:5),"bsc","dead")

bl_pop        <- c(1,rep(0,length(state_names)-1))
names(bl_pop) <- state_names

tp_1L   <- rdirichlet(TH,c(95,2 ,0 ,0 ,0 ,3 ,4))
tp_2L   <- rdirichlet(TH,c(0 ,93,6 ,0 ,0 ,3 ,5))
tp_3L   <- rdirichlet(TH,c(0 ,0 ,92,5 ,0 ,3 ,6))
tp_4L   <- rdirichlet(TH,c(0 ,0 ,0 ,91,3 ,3 ,7))
tp_5L   <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,90,3 ,8))
tp_BSC  <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,0 ,85,9))
tp_dead <- c(0, 0, 0, 0, 0, 0, 1)

# List of TPMs with some random variation cycle to cycle to simulate time-varying
# transition probabilities between states. Patients can't go backwards in line, and
# can't skip lines except for going straight to bsc (i.e. choosing no treatment/palliative care)
tp <- lapply(1:TH, function(cycle) {
  matrix(
    c(tp_1L[cycle,], 
      tp_2L[cycle,], 
      tp_3L[cycle,], 
      tp_4L[cycle,],
      tp_5L[cycle,],
      tp_BSC[cycle,],
      tp_dead
    ),
    nrow = 7,
    byrow = TRUE,
    dimnames = list(state_names,state_names)
  )
})

# some example of the tunnel length vector. This could be derived as floor(TTNT - TTDisc)?
vec_tun_len <- vec_tun_len

tun_positions <- 1


# Some random numbers for the tunnel death probs to show it varying over time. Note that
# the length goes down as the cycle goes up because there's no point in having numbers
# that go past the time horizon.
# 
# The reason it's in an extra layer of list is to illustrate that you will have data
# like this PER TUNNEL STATE. This example is just 1 tunnel.

l_tun_os <- lapply(1:TH, function(cyc) {
    # Make a matrix which contains our (example) random numbers for probability of death
    # in tunnel. obviously these will be different when we're using real data. Transpose
    # it so we can split it by row afterwards. Speeds things up later so we may as well
    # do it now.
    
    tun_len_this_cohort <- vec_tun_len[cyc]
    
    tp_i <- t(matrix(
      data = rbeta((TH - cyc + 1) * tun_len_this_cohort,10,100),
      nrow = TH - cyc + 1,
      ncol = tun_len_this_cohort,
      byrow = TRUE,
      dimnames = list(NULL,NULL)
    ))
    tp_i <- split(tp_i, rep(1:ncol(tp_i), each = nrow(tp_i)))
    names(tp_i) <- NULL
    
    # extend the values with 100% death probability to capture the impossible people
    # who shouldn't be there anyway. This way will force rowsums to be 
    
    if (max_tun_len - tun_len_this_cohort > 0) {
      tun_ext_len <- rep(1,max_tun_len - tun_len_this_cohort)
      return(lapply(tp_i, function(this_cycle_this_cohort) c(this_cycle_this_cohort, tun_ext_len)))
    } else {
      return(tp_i)
    }
  })



tp_tun1 <- f_mk_varyTunneltpms(
  list_tun_os  = l_tun_os,
  tpm_list     = tp,
  tun_position = 1,
  vec_tun_len  = vec_tun_len
)

# Right - so to add multiple tunnels like this, you need:
# 
# A list of mortality for each tunnel state
# 
# ...I think that's it?
# 
# 

# Just applying the same mortality for now to avoid extra lines for no reason
l_tun_os2 <- l_tun_os




# Change of tac -----------------------------------------------------------

# ok, so upon discussion with Dawn Lee, it's going to be a massive matrix
# which takes the TPs for each cycle given cycles in tunnel and places them in 
# a square matrix with dimensions (Ntun * TH_in_cycles) + Nstates. So, if there's
# 4 treatment lines, 3 tunnels between them, BSC+death, weekly CL and 40 year TH, 
# the dimensions of the matrix are (3 * (40*52)) + 1 + 1 = (6242X6242) FOR EACH CYCLE,
# so for each treatment pathway we end up with TH 6242x6242 matrices, which then
# go into a Reduce to calculate the trace, and then get discarded.


# size_of_mat_in_MB <- format(object.size(matrix(
#   rep(runif(6242*6242)),
#   nrow=6242,
#   ncol=6242
# )), units="MB", standard="SI")

# Which results in 311.7 MB of RAM per matrix. 311.7 * 52 * 40 = 648,336 MB (648 GB) of RAM
# requirement to compute the matrices. Obviously this is assuming that there are no 0s,
# which might reduce burden. Let's try it again with a diagonal matrix

# size_of_diag_in_MB <- format(object.size(diag(runif(6242))), units="MB", standard="SI")

# Sadly it's the same amount of RAM...Let's check the size of it actually IN RAM:

# M <- diag(runif(6242))

# Yep, 311MB of RAM per matrix.
# 
# library(biganalytics)
# 
# M <- as.big.matrix(M)
# 
# rep(1,6242) %*% M
# 
# Pt <- t(as.matrix(rep(1,6242)))
# 
# # Calculate the next row using M and :
# PM  <- bigalgebra::dgemm(A = Pt,B = M)
# PM2 <- bigalgebra::dgemm(A = PM,B = M)
# PM3 <- bigalgebra::dgemm(A = PM2,B = M)

# This would work but would obviously slow things down.


state_names <- state_names[c(1:4,6:7)]


tp_1L   <- rdirichlet(TH,c(95,2 ,0 ,0 ,1 ,4))
tp_2L   <- rdirichlet(TH,c(0 ,93,6 ,0 ,1 ,5))
tp_3L   <- rdirichlet(TH,c(0 ,0 ,92,5 ,1 ,6))
tp_4L   <- rdirichlet(TH,c(0 ,0 ,0 ,91,1 ,7))
tp_5L   <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,90,8))
tp_BSC  <- rdirichlet(TH,c(0 ,0 ,0 ,0 ,88,9))
tp_dead <- c(0, 0, 0, 0, 0, 1)

# List of TPMs with some random variation cycle to cycle to simulate time-varying
# transition probabilities between states. Patients can't go backwards in line, and
# can't skip lines except for going straight to bsc (i.e. choosing no treatment/palliative care)
tp <- lapply(1:TH, function(cycle) {
  matrix(
    c(tp_1L[cycle,], 
      tp_2L[cycle,], 
      tp_3L[cycle,], 
      tp_4L[cycle,],
      tp_BSC[cycle,],
      tp_dead
    ),
    nrow = length(state_names),
    byrow = TRUE,
    dimnames = list(state_names,state_names)
  )
})


# Load in the TPs and tidy columns up
TPs <- openxlsx::read.xlsx("./1_Data/TP matrix expansion toy example.xlsm",sheet = "Sheet3")
TPs <- as.data.table(TPs)
colnames(TPs) <- c("c", unlist(lapply(1:4, function(x){paste0("L",x,c("_disc", "_next", "_death"))})),"NT_death")


TP1 <- matrix(nrow=10,ncol=10)
TP1[is.na(TP1)] <- 0

cyc <- as.numeric(TPs[1,])
nc <- dim(TP1)[2]

TP1[1,1]  <- 1-sum(cyc[2:4])
TP1[1,c(2:3,nc)] <- cyc[2:4]

TP1[2,2]  <- TP1[1,1] + TP1[1,2]
TP1[2,c(3,nc)]  <- TP1[1,c(3,nc)]

c3 <- c(4:5,nc)
TP1[3,3]   <- 1-sum(cyc[5:7])
TP1[3,c3]  <- cyc[5:7]

TP1[4,4]  <- TP1[3,3] + TP1[3,4]
TP1[4,c(5,nc)]  <- TP1[3,c(5,nc)]
  
c5 <- c(6:7,nc)
TP1[5,5]   <- 1-sum(cyc[8:10])
TP1[5,c5]  <- cyc[8:10]

TP1[6,6]  <- TP1[5,5] + TP1[5,6]
TP1[6,c(7,nc)]  <- TP1[5,c(7,nc)]

c7 <- c(8:9,nc)
TP1[7,7]   <- 1-sum(cyc[11:13])
TP1[7,c7]  <- cyc[11:13]

TP1[8,8]  <- TP1[7,7] + TP1[7,8]
TP1[8,c(9,nc)]  <- TP1[7,c(9,nc)]

TP1[9,9]   <- 1-cyc[14]
TP1[9,10]  <- cyc[14]

TP1[nc,nc] <- 1

# Test:
all(round(rowSums(TP1),12)==1)



# ~ repeating the VBA provided --------------------------------------------

# Three macros are provided, one to create the empty matrix, one to populate it
# and one to test it. let's replicate that action.
# 
# The action is essentially to insert a set of diagonal matrices into the full matrix.
# These should have NO overlap with each other as it should be impossible to transition
# between them
# 




# Create matrix outline

TH <- nrow(TPs)

# Define all the labels for the rows and columns
labs <- c(
  "L1_on",
  "L1_off",
  paste0("L2_on_c" ,1:TH),
  paste0("L2_off_c",1:TH),
  paste0("L3_on_c" ,1:TH),
  paste0("L3_off_c",1:TH),
  paste0("L4_on_c" ,1:TH),
  paste0("L4_off_c",1:TH),
  paste0("BSC_c"   ,1:TH),
  "dead"
)

len_m <- length(labs)

# Let's start with a diagonal matrix - 100% chance of staying in the same state for 
# all states.
M                   <- diag(nrow=len_m)
M[row(M) == col(M)] <- 0
dimnames(M)         <- list(labs,labs)
# all(rowSums(M) == 0) # test that it's full of zeros - it is a lot faster to do this
# than it is to make a matrix full of NAs and replace all NAs (the default) with 0s
M[len_m,len_m]      <- 1

# Size of this matrix in MB - not actually that bad, it's only 1.6GB of RAM!
# Make sure to only use this method when there's lots of RAM to spare, or it could error out.
format(object.size(M),units = "GB")

# Rows 3 onward in this matrix are FIXED. once this matrix is populated, the values
# in the following cells can simply be replaced and the vector of state residency
# can then be multiplied by the updated matrix to produce the next state residency.
# 
# This is a huge efficiency gain and ensures that all transitions are taken into
# account. If the matrix gets really really big, we can use the "bigmemory" package
# which can perform matrix multiplication of files using C++ directly. This performs
# the operation on disk so you're only limited by your HDD. 
# 
# Rows 1 and 2 control line 1, and are simply replaced with movement probabilities
# from the associated row of TPs above each cycle. the rest of the matrix
# is populated systematically once.

gc()

# Let's have a bit of a sneak peak before we get started. a big square of zero values.
f_mat_firstFewCells <- function(mat, n=10) {mat[1:n,1:n]}

f_mat_firstFewCells(M)

# I programmed a function based on an excellent stack exchange answer, credit to David Arenburg at
# https://stackoverflow.com/questions/28746023/how-to-insert-values-from-a-vector-diagonally-into-a-matrix-in-r
# which taught me that one can insert values into an existing matrix in a vectorised manner
# by supplying a matrix of co-ordinate values to the square brackets, like so:
# 
# A[co_ord_matrix] <- b
# 
# where A is a matrix (not necessarily square), and b is a vector whcih must of course
# fit inside the matrix. the co_ord_matrix must be 2 columns and must contain
# integer values.
# 
# Function is called f_vec_in_mat_diag and it is a huge time and effort saver here
# because we only need to know the top left element location for each vector
# to insert to the big matrix!
# 
# 

source("./3_Functions/misc/matrix.R")

# For example, in our case, the "stay in tunnel" probabilities for line 2 on treatment
# tunnel start from element 3,4 in M. to demonstrate with some dummy data:

f_vec_in_mat_diag(
  A = matrix(nrow = 10,ncol = 10),
  b = 1:5,
  start_row = 3,
  start_col = 4,
  direction = "se"
)

# As the values inside of the matrix all come from or are derived from the
# table TPs, all I need to know is positions within the matrix (top left cells only).


# So these are the top left cells for each tunnel

n_tun_excluding_bsc <- 3
lab_co_ord <- unlist(lapply(1:n_tun_excluding_bsc, function(tun) {
  line_on <- paste0("tun ",tun+1,"L on")
  line_on <- paste0(line_on,c(" stay", " stop", " next"))
  line_off <- paste0("tun ",tun+1,"L off")
  line_off <- paste0(line_off,c(" stay", " next"))
  c(line_on, line_off)
}))

lab_co_ord <- c(lab_co_ord, "tun BSC survive")

# the rows and columns in the workbook and found these are numbered manually. There
# are gaps of 1 row or column in there (row sums to 0 in excel). 
tl_cells <- matrix(
  c(
    3    , 4    , # tun 2L on  stay
    3    , 2084 , # tun 2L on  stop
    3    , 4163 , # tun 2L on  next
    2083 , 2084 , # tun 2L off stay
    2083 , 4163 , # tun 2L off next
    
    4163 , 4164 , # tun 3L on  stay
    4163 , 6244 , # tun 3L on  stop
    4163 , 8323 , # tun 3L on  next
    6243 , 6244 , # tun 3L off stay
    6243 , 8323 , # tun 3L off next
    
    8323 , 8324 , # tun 4L on  stay
    8323 , 10404, # tun 4L on  stop
    8323 , 12483, # tun 4L on  next
    10403, 10404, # tun 4L off stay
    10403, 12483, # tun 4L off next
    
    12483, 12484  # tun BSC survive
  ),
  nrow = (n_tun_excluding_bsc*5 + 1),
  byrow = TRUE,
  dimnames = list(lab_co_ord,c("row", "col"))
)

# now, we could systematically derive these using the time horizon TH and the 
# non-time-horizon bounded rows (i.e. the first 2). something like:

nt <- 2

# For each tunnel except the last one which has no adjoining tunnel, there are
# 4 transitions, stay, discontinue, next line and death whilst in the on-treatment tunnel.
# there are 3 transitions (stay, next, death) when in the off-treatment tunnel.
# Death transitions are always in the final column of the matrix, so can be
# handled separately. we're only concerned with the diagonal bits for now!
# 
# As a toy example, let's calculate the rows and columns above using parameters:
# 
# tun_n, the tunnel number (starting from 1)
# nt non tunnel rows (always the number of rows of non-tunnel states that come 
# before the first tunnel in the series)
# TH time horizon in cycles
# TPs our table of transition probabilities

# First tunnel:
tun_n          <- 1

tun_r_on       <- sum(nt    , ((2*TH) * (tun_n-1)), 1)
tun_c_on_stay  <- sum(nt + 1, ((2*TH) * (tun_n-1)), 1)
tun_c_on_stop  <- tun_c_on_stay + TH
tun_c_on_next  <- tun_c_on_stop + TH - 1
tun_r_off      <- tun_c_on_stop - 1
tun_c_off_stay <- tun_c_on_stop
tun_c_off_next <- tun_c_on_next

# Checked and these match


# second tunnel:
tun_n          <- 2

tun_r_on       <- sum(nt    , ((2*TH) * (tun_n-1)), 1)
tun_c_on_stay  <- sum(nt + 1, ((2*TH) * (tun_n-1)), 1)
tun_c_on_stop  <- tun_c_on_stay + TH
tun_c_on_next  <- tun_c_on_stop + TH - 1
tun_r_off      <- tun_c_on_stop - 1
tun_c_off_stay <- tun_c_on_stop
tun_c_off_next <- tun_c_on_next


# Looks good. Let's formalise into a second function which gets our upper left
# element indices for us, given a tunnel and a time horizon:
# 
# ASSUMES THE TUN N IS FOR THE NEXT LINE, tun_n 1 is 2L etc
# 

f_diag_mat_topleftFinder <- function(tun_n, TH, nt) {
  
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
        paste0(paste0("tun ",tun_n + 1,"L on "),c("stay","stop","next")),
        paste0(paste0("tun ",tun_n + 1,"L off "),c("stay","next"))
      ),
      c("row", "col")
    )
  )
  return(outMat)
}

# updated version with naming of the rows that lines up with the columns in tp 
# defined during f_seq_tpm_compiler
f_diag_mat_topleftFinder2 <- function(tun_n, TH, nt) {
  
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


# Now we can much more easily get our top lefts for all our diagonals:

tl_tun1 <- f_diag_mat_topleftFinder(
  tun_n = 1,
  TH = TH,
  nt = 2
)



# We can then cycle down this list putting values into M using our function.
# For instance, let's do the "stay" probability for the first tunnel, which is
# derived by whats left in TPs:

TPs$L1_stay   <- 1-rowSums(TPs[,list(L1_disc, L1_next, L1_death)])
TPs$L2_stay   <- 1-rowSums(TPs[,list(L2_disc, L2_next, L2_death)])
TPs$L3_stay   <- 1-rowSums(TPs[,list(L3_disc, L3_next, L3_death)])
TPs$L4_stay   <- 1-rowSums(TPs[,list(L4_disc, L4_next, L4_death)])
TPs$BSC_death <- 1-TPs$BSC_stay

# CHECK: sum to 1
all(round(TPs$L1_stay + TPs$L1_disc + TPs$L1_next + TPs$L1_death,12) == 1)
all(round(TPs$L2_stay + TPs$L2_disc + TPs$L2_next + TPs$L2_death,12) == 1)
all(round(TPs$L3_stay + TPs$L3_disc + TPs$L3_next + TPs$L3_death,12) == 1)
all(round(TPs$L4_stay + TPs$L4_disc + TPs$L4_next + TPs$L4_death,12) == 1)
all(round(TPs$BSC_stay + TPs$BSC_death,12) == 1)


# Let's try and do the whole first tunnel state with the first cycle 
# probabilities for 1L treatment as well, and with mortality

# 1L TPM for cycle 1. Let's just calculate L1 stay to avoid having to do it again and again


# Now pull out the right TPs for the first cycle, and taking the same strategy as
# for the diagonal values, apply the values all at the same time:
# 
# Remember, this is for 1L for each cycle
cyc_prob <- as.list(TPs[1,])
vals_1L  <- matrix(
  c(
    1, 1    , .subset2(cyc_prob,"L1_stay"),
    1, 2    , .subset2(cyc_prob,"L1_disc"),
    1, 3    , .subset2(cyc_prob,"L1_next"),
    1, len_m, .subset2(cyc_prob,"L1_death"),
    2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
    2, 3    , .subset2(cyc_prob,"L1_next"),
    2, len_m, .subset2(cyc_prob,"L1_death")
  ),
  nrow  = 7,
  byrow = TRUE
)
M[vals_1L[,1:2]] <- vals_1L[,3]


# Check on first 2 rows:
all(rowSums(M[1:2,])==1)

# Great, so the above is our process once per cycle to update the matrix before
# running it again. All we now need to do is correctly build the rest of the matrix
# just once


# 2L tunnel, staying in the tunnel:
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L2_stay,
  start_row = tl_tun1["tun 2L on stay","row"],
  start_col = tl_tun1["tun 2L on stay","col"],
  direction = "se"
)

# 2L tunnel: Discontinuation (period of no treatment between lines i.e., 2L-->off trt | 2L)
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L2_disc,
  start_row = tl_tun1["tun 2L on stop","row"],
  start_col = tl_tun1["tun 2L on stop","col"],
  direction = "se"
)
# Move directly to next line (i.e. 2L-->3L | 2L)
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L2_next,
  start_row = tl_tun1["tun 2L on next","row"],
  start_col = tl_tun1["tun 2L on next","col"],
  direction = "se"
)

# Die within the tunnel
M[tl_tun1["tun 2L on stop","row"]:(tl_tun1["tun 2L on stop","row"] + TH - 1),len_m] <- TPs$L2_death

f_mat_firstFewCells(M)

# ok so does it add to 1?

all(rowSums(M[1:50,])==1)


# Brilliant, we've successfully added the first on-treatment tunnel. Now let's add the off-treatment
# tunnel state
M <- f_vec_in_mat_diag(
  A         = M,
  b         = 1 - TPs$L2_next - TPs$L2_death,
  start_row = tl_tun1["tun 2L off stay","row"],
  start_col = tl_tun1["tun 2L off stay","col"],
  direction = "se"
)
M <- f_vec_in_mat_diag(
  A         = M,
  b         = TPs$L2_next,
  start_row = tl_tun1["tun 2L off next","row"],
  start_col = tl_tun1["tun 2L off next","col"],
  direction = "se"
)
M[tl_tun1["tun 2L off next","row"]:(tl_tun1["tun 2L off next","row"] + TH - 1),len_m] <- TPs$L2_death


# Still testing well to 12 dp
all(round(rowSums(M[tl_tun1["tun 2L off next","row"]:(tl_tun1["tun 2L off next","row"]+TH),]),12)==1)


# So, we've done it! we've added the first tunnel state representing second-line treatment


# Now, let's add the second tunnel state (3L treatment)


tl_tun2 <- f_diag_mat_topleftFinder(
  tun_n = 2,
  TH    = TH,
  nt    = 2
)

# 3L on-treatment tunnel:
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L3_stay,
  start_row = tl_tun2["tun 3L on stay","row"],
  start_col = tl_tun2["tun 3L on stay","col"],
  direction = "se"
)
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L3_disc,
  start_row = tl_tun2["tun 3L on stop","row"],
  start_col = tl_tun2["tun 3L on stop","col"],
  direction = "se"
)
M <- f_vec_in_mat_diag(
  A = M,
  b = TPs$L3_next,
  start_row = tl_tun2["tun 3L on next","row"],
  start_col = tl_tun2["tun 3L on next","col"],
  direction = "se"
)
M[tl_tun2["tun 3L on stop","row"]:(tl_tun2["tun 3L on stop","row"] + TH - 1),len_m] <- TPs$L3_death

# Off treatment tunnel:
M <- f_vec_in_mat_diag(
  A         = M,
  b         = 1 - TPs$L3_next - TPs$L3_death,
  start_row = tl_tun2["tun 3L off stay","row"],
  start_col = tl_tun2["tun 3L off stay","col"],
  direction = "se"
)
M <- f_vec_in_mat_diag(
  A         = M,
  b         = TPs$L3_next,
  start_row = tl_tun2["tun 3L off next","row"],
  start_col = tl_tun2["tun 3L off next","col"],
  direction = "se"
)
M[tl_tun2["tun 3L off next","row"]:(tl_tun2["tun 3L off next","row"] + TH - 1),len_m] <- TPs$L3_death


all(round(rowSums(M[1:8000,]),12) == 1)

# Amazing, we've done it!




# Formalising into a function ---------------------------------------------


# The function will have 4 steps:
# 
#   1. Calculate the extra columns needed from TPs
#   2. Pre-define the matrix M (with the labs)
#   3. Compute the top left element in M associated with each tunnel
#   4. populate the tunnels for on treatment, off treatment, and BSC
# 


#' Function to compile a transition probability matrix defining a treatment sequence,
#' where patients move from one treatment to the next, or straight to death, including
#' on and off treatment tunnels
#' 
#' WARNING: This function assumes the cycle number is 1. This is because you should not
#'          be replacing the TPM M including the whole tunnel each cycle, Instead
#'          use this function just once at baseline to create the structure, and then
#'          replace the required values for transitions relating to first-line therapy
#' 
#' @param tp table of TPs. MUST have 1st column cycle, then disc, next, death sequentially for each line, plus no treatment death (e.g. t L1_disc L1_next L1_death L2_disc...NT_death)
#' @param n_lines Number of treatment lines total, including first line and BSC (for 4 tunnels this is 6 with BSC)
#' @param include_bsc whether to include a tunnel for BSC or skip BSC. If FALSE tp has one less column
#' @param nt non-tunnel rows - almost always 2, but if for example second line is not a tunnel (e.g. relapse and remission model) then it would be 4 and tunnels would start at 3L
#' 
#' 
f_seq_tpm_compiler <- function(tp, n_lines, include_bsc = TRUE, nt = 2) {
  
  require(data.table)
  require(collapse)
  
  # Assume cycle = 1, which it should be every time this function is used!
  c <- 1
  
  # calculate time horizon based on TP table
  TH <- nrow(tp)
  
  target_col <- 1 + (n_lines * 3)
  if(include_bsc) target_col <- target_col + 1
  
  if(include_bsc==FALSE) stop("Excluding BSC is not supported yet, please contact Darren Burns or Dawn Lee to build it!")
  
  # Validation - the number of columns of tp should adhere to the documentation above the function
  if(ncol(tp) != target_col) stop(
    paste0(
      "tp should have ",
      target_col,
      " columns, but it has ",
      ncol(tp),
      ". Please ensure that the columns are cycle number, 1L discontinuation, 1L next treatment, 1L death, 2L discontinuation ... no-treatment death"
    )
  )
  
  cat("re-labelling columns in tp\n")
  # Right, now we can get on with it. Start by auto-labelling the columns:
  colnames(tp) <- c(
    "c",
    unlist(lapply(1:n_lines, function(tline) {
      paste0(paste0("L",tline,"_"),c("disc","next","death"))
    })),
    "BSC_death"
  )
  
  # data.tables are easier to work with when lots of rows:
  if(class(tp)[1] != "data.table") tp <- as.data.table(tp)
  
  cat("calculating extra columns in tp\n")
  # now we have standardised naming and standardised format, we can calculate
  # our new columns:
  tp$L1_stay  <- 1-rowSums(tp[,list(L1_disc, L1_next, L1_death)])
  tp$L2_stay  <- 1-rowSums(tp[,list(L2_disc, L2_next, L2_death)])
  tp$L3_stay  <- 1-rowSums(tp[,list(L3_disc, L3_next, L3_death)])
  tp$L4_stay  <- 1-rowSums(tp[,list(L4_disc, L4_next, L4_death)])
  tp$BSC_stay <- 1-tp$BSC_death
  
  # Now we have probability to stay in state as well as probability of 
  # discontinuation, next line, death, the rest is putting these columns
  # or derivitives of them into M
  
  cat("Generating matrix M\n")
  # Part 2: define M (WARNING this uses a lot of RAM)
  labs <- c(
    "L1_on",
    "L1_off",
    unlist(lapply(2:n_lines, function(tline) {
      c(paste0("L",tline,"_on_c" ,1:TH),paste0("L",tline,"_off_c" ,1:TH))
    })),
    paste0("BSC_c"   ,1:TH),
    "dead"
  )
  len_m               <- length(labs)
  M                   <- diag(nrow=len_m)
  M                   <- recode_num(M,`1`=0)
  M[len_m,len_m]      <- 1
  
  cat("Finding diagonal start points\n")
  # Part 3: compute top left elements
  tl_elements <- do.call(
    rbind,
    lapply(1:(n_lines-1), function(i_tun) {
      f_diag_mat_topleftFinder2(
        tun_n = i_tun,
        TH    = TH,
        nt    = 2 
      )
    })
  )
  
  cat("Generating co-ordinates for entry into M\n")
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
  cyc_prob <- as.list(tp[1,])
  vals_1L  <- matrix(
    c(
      1, 1    , .subset2(cyc_prob,"L1_stay"),
      1, 2    , .subset2(cyc_prob,"L1_disc"),
      1, 3    , .subset2(cyc_prob,"L1_next"),
      1, len_m, .subset2(cyc_prob,"L1_death"),
      2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
      2, 3    , .subset2(cyc_prob,"L1_next"),
      2, len_m, .subset2(cyc_prob,"L1_death")
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
    tl_tun <- f_diag_mat_topleftFinder(
      tun_n = line_n,
      TH    = TH,
      nt    = nt
    )
    
    # Make an id for this line which is used to put the right numbers in the right place
    tline   <- paste0("L",line_n)
    
    # Pull out the columns to put as the data to enter in M
    line_tp <- as.matrix(tp)[,which(grepl(tline,colnames(tp)))]
    p_stay  <- line_tp[,which(grepl("stay",colnames(line_tp)))]
    p_disc  <- line_tp[,which(grepl("disc",colnames(line_tp)))]
    p_next  <- line_tp[,which(grepl("next",colnames(line_tp)))]
    p_death <- line_tp[,which(grepl("death",colnames(line_tp)))]
    
    # make co-ordinate values and value values for the on and off treatment tunnels including
    # death for this treatment line:
    co_list <- list(
      # on-treatment tunnel:
      co_on_stay   = matrix(c(tl_tun[1,"row"] + indx, tl_tun[1,"col"] + indx,p_stay) ,ncol=3),
      co_on_disc   = matrix(c(tl_tun[2,"row"] + indx, tl_tun[2,"col"] + indx,p_disc) ,ncol=3),
      co_on_next   = matrix(c(tl_tun[3,"row"] + indx, tl_tun[3,"col"] + indx,p_next) ,ncol=3),
      co_on_death  = matrix(c(tl_tun[1,"row"] + indx, rep(len_m,TH)         ,p_death),ncol=3),
      
      # Off treatment tunnel for this line:
      co_off_stay  = matrix(c(tl_tun[4,"row"] + indx, tl_tun[4,"col"] + indx,1 - p_death - p_next),ncol=3),
      co_off_next  = matrix(c(tl_tun[5,"row"] + indx, tl_tun[5,"col"] + indx,p_next              ),ncol=3),
      co_off_death = matrix(c(tl_tun[4,"row"] + indx, rep(len_m,TH)         ,p_death             ),ncol=3)
    )
    
    # Bind it together into one larger matrix and return that
    return(do.call(rbind,co_list))
  })
  
  # co_ordinate_matrix now has all of the data organised except for the BSC tunnel.
  # We add that here:
  
  bsc_start_point <- f_diag_mat_topleftFinder(
    tun_n = n_lines - 1,
    TH    = TH,
    nt    = nt
  )[5,"col"]
  co_bsc <- rbind(
    matrix(c(bsc_start_point + indx, bsc_start_point + 1 + indx,tp$BSC_stay) ,ncol=3),
    matrix(c(bsc_start_point + indx, rep(len_m,TH)             ,tp$BSC_death) ,ncol=3)
  )
  
  # bind it all together into one massive matrix of co-ordinates and values:  
  co_ordinate_matrix <- rbind(do.call(rbind,co_ordinate_matrix),co_bsc)
  
  cat("Inserting values\n")
  # Put all values into the matrix in one vectorised command:
  M[co_ordinate_matrix[,1:2]] <- co_ordinate_matrix[,3]
  
  # The last cycle of BSC is an issue, stay probability isn't entered.
  M[len_m-1,len_m-1] <- 1-M[len_m-1,len_m]
  
  # Well all rows sum to 1 so that's promising:
  # all(round(rowSums(M),10) == 1) # reads true:
  return(M)
  
}

TPs <- openxlsx::read.xlsx("./1_Data/TP matrix expansion toy example.xlsm",sheet = "Sheet3")
TPs <- as.data.table(TPs)
colnames(TPs) <- c("c", unlist(lapply(1:4, function(x){paste0("L",x,c("_disc", "_next", "_death"))})),"NT_death")



M <- f_seq_tpm_compiler(
  tp          = TPs,
  n_lines     = 4,
  include_bsc = TRUE,
  nt          = 2
)

# Test that all rows sum to 1 each:
all(round(rowSums(M),10) == 1)

# Generating a trace ------------------------------------------------------

# Now that we have M for cycle 1 we can follow this process:
# 
#   1. update 1L transitions in M 
#   2. p_(t-1) %*% M
#   
# That's it!

tr <- matrix(
  ncol = dim(M)[2],
  nrow = TH
)
tr[is.na(tr)] <- 0
tr[1,1]       <- 1


# names(bl_pop) <- dimnames(M)[[2]]

TPs$L1_stay  <- 1-rowSums(TPs[,list(L1_disc, L1_next, L1_death)])


TRACE <- Reduce(
  x = 1:10,
  init = tr,
  accumulate = FALSE,
  f = function(prev, c) {
    if (c == 1) {
      # If it's the first cycle, we've calculated M already so just multiply it
      cat("cycle 1\n")
      prev[2,] <- as.numeric(prev[1,] %*% M)
      return(prev)
    } else {
      # It's a subsequent cycle, replace the necessary cells for first line, everything else stays the same:
      cat(paste0("cycle ",c,"\n"))
      cyc_prob <- as.list(TPs[c,])
      vals_1L  <- matrix(
        c(
          1, 1    , .subset2(cyc_prob,"L1_stay"),
          1, 2    , .subset2(cyc_prob,"L1_disc"),
          1, 3    , .subset2(cyc_prob,"L1_next"),
          1, len_m, .subset2(cyc_prob,"L1_death"),
          2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
          2, 3    , .subset2(cyc_prob,"L1_next"),
          2, len_m, .subset2(cyc_prob,"L1_death")
        ),
        nrow  = 7,
        byrow = TRUE
      )
      M[vals_1L[,1:2]] <- vals_1L[,3]
      
      prev[c+1,] <- as.numeric(prev[c,] %*% M)
      return(prev)
    }
  }
)

# HAHAHAHA we did it Dawn! It's slow but it's correct I think, at least for the first few cycles!
rowSums(TRACE[1:10,])
f_mat_firstFewCells(TRACE)


# now, this is unacceptably slow, so we need to figure out how to optimise it a bit:

library(profvis)

# The assignment of values to a large matrix is the problem:

microbenchmark::microbenchmark(
  get_value = {M[vals_1L[,1:2]]},
  set_value = {M[vals_1L[,1:2]] <- vals_1L[,3]},
  times = 1000
)

# To prove it:


profvis({
  TRACE <- Reduce(
    x = 1:10,
    init = tr,
    accumulate = FALSE,
    f = function(prev, c) {
      if (c == 1) {
        # If it's the first cycle, we've calculated M already so just multiply it
        cat("cycle 1\n")
        prev[2,] <- as.numeric(prev[1,] %*% M)
        return(prev)
      } else {
        # It's a subsequent cycle, replace the necessary cells for first line, everything else stays the same:
        cat(paste0("cycle ",c,"\n"))
        cyc_prob <- as.list(TPs[c,])
        vals_1L  <- matrix(
          c(
            1, 1    , .subset2(cyc_prob,"L1_stay"),
            1, 2    , .subset2(cyc_prob,"L1_disc"),
            1, 3    , .subset2(cyc_prob,"L1_next"),
            1, len_m, .subset2(cyc_prob,"L1_death"),
            2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
            2, 3    , .subset2(cyc_prob,"L1_next"),
            2, len_m, .subset2(cyc_prob,"L1_death")
          ),
          nrow  = 7,
          byrow = TRUE
        )
        M[vals_1L[,1:2]] <- vals_1L[,3]
        
        prev[c+1,] <- as.numeric(prev[c,] %*% M)
        return(prev)
      }
    }
  )
})

# So, 200/280 ms spent assigning values to M, even though it's only a few. This is because
# R copy and pastes the ENTIRE 1.6GB matrix to do this (See R documentation for proof of this)
# 
# link: https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Subset-assignment
# also: https://github.com/r-lib/R6/issues/201
# 
# Now, there is a package called bigalgebra, which stores the matrix in a reference
# style object, meaning making changes to that object will not copy paste it, but will
# change the values directly on the disc or in RAM. that's what we need here. 
# 
# This package also has the ability to perform matrix multiplication, but ONLY on
# to 
# 

M <- as.big.matrix(M)

# tr[[1]] %*% M now doesn't work, but bigalgebra has a function for it

tr2 <- list(as.big.matrix(matrix(tr[1,],nrow = 1)))

# We would do bigalgebra::dgemm(A=tr[[1]],B=M) to matrix multiply the vector
# of state occupancy in t-1 by the matrix M which has been updated to incorporate
# this cycle's first-line transitions

i  <- list(tr=tr2,M=M)

profvis({
  TRACE <- Reduce(
    x = 1:10,
    init = i,
    accumulate = FALSE,
    f = function(prev, c) {
      if (c == 1) {
        # If it's the first cycle, we've calculated M already so just multiply it
        cat("cycle 1\n")
        prev$tr[[2]] <- bigalgebra::dgemm(A = prev$tr[[1]],B = prev$M)
        return(prev)
      } else {
        # It's a subsequent cycle, replace the necessary cells for first line, everything else stays the same:
        cat(paste0("cycle ",c,"\n"))
        cyc_prob <- as.list(TPs[c,])
        vals_1L  <- matrix(
          c(
            1, 1    , .subset2(cyc_prob,"L1_stay"),
            1, 2    , .subset2(cyc_prob,"L1_disc"),
            1, 3    , .subset2(cyc_prob,"L1_next"),
            1, len_m, .subset2(cyc_prob,"L1_death"),
            2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
            2, 3    , .subset2(cyc_prob,"L1_next"),
            2, len_m, .subset2(cyc_prob,"L1_death")
          ),
          nrow  = 7,
          byrow = TRUE
        )
        
        # Assign values to the big.matrix which ONLY changes those values and 
        # doesn't copy paste the matrix to do it.
        
        prev$M[vals_1L[,1:2]] <- vals_1L[,3]
        
        prev$tr[[c+1]] <- bigalgebra::dgemm(A = prev$tr[[c]],B = prev$M)
        return(prev)
      }
    }
  )
})
 

TRACE$tr <- do.call(rbind,lapply(TRACE$tr,bigmemory::as.matrix))

rowSums(TRACE$tr)
f_mat_firstFewCells(TRACE$tr)

# Great, so there's a big computational gain, I think about 5x quicker on my computer.
# 
# Now let's have a look at just how long this is going to take to compute using 50 cycles
# as a starting point:
# 

t1 <- Sys.time()

tr3 <- tr2[[1]]

# reset M 
cyc_prob <- as.list(TPs[1,])
vals_1L  <- matrix(
  c(
    1, 1    , .subset2(cyc_prob,"L1_stay"),
    1, 2    , .subset2(cyc_prob,"L1_disc"),
    1, 3    , .subset2(cyc_prob,"L1_next"),
    1, len_m, .subset2(cyc_prob,"L1_death"),
    2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
    2, 3    , .subset2(cyc_prob,"L1_next"),
    2, len_m, .subset2(cyc_prob,"L1_death")
  ),
  nrow  = 7,
  byrow = TRUE
)
M[vals_1L[,1:2]] <- vals_1L[,3]

f_mat_firstFewCells(M)

# See what 
TRACE <- Reduce(
  x = 1:50,
  init = tr3,
  accumulate = TRUE,
  f = function(prev, c) {
    cat(paste0("cycle ",c,"\n"))
    cyc_prob <- as.list(TPs[c,])
    vals_1L  <- matrix(
      c(
        1, 1    , .subset2(cyc_prob,"L1_stay"),
        1, 2    , .subset2(cyc_prob,"L1_disc"),
        1, 3    , .subset2(cyc_prob,"L1_next"),
        1, len_m, .subset2(cyc_prob,"L1_death"),
        2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
        2, 3    , .subset2(cyc_prob,"L1_next"),
        2, len_m, .subset2(cyc_prob,"L1_death")
      ),
      nrow  = 7,
      byrow = TRUE
    )
    M[vals_1L[,1:2]] <- vals_1L[,3]
    return(bigalgebra::dgemm(A = prev,B = M))
  }
)


t2 <- Sys.time()

print(t2 - t1)

t50 <- as.numeric(t2 - t1)

# Running one treatment line with this method will take approximately this much time:
t_th_mins <- ((t50 / 50) * TH) / 60
print(paste0("Matrix method - one trace for one treatment line will take approximately ", round(t_th_mins,2), " minutes to run..."))

TRACE <- do.call(rbind,lapply(TRACE,bigmemory::as.matrix))
f_mat_firstFewCells(TRACE)
round(rowSums(TRACE),12) == 1


# Hmm...that's a long time! This might be one of those cases where it's actually
# better to use a for loop:

tr4 <- tr2

t3 <- Sys.time()
for (c in 1:50) {
  cat(paste0("cycle ",c,"\n"))
  cyc_prob <- as.list(TPs[c,])
  vals_1L  <- matrix(
    c(
      1, 1    , .subset2(cyc_prob,"L1_stay"),
      1, 2    , .subset2(cyc_prob,"L1_disc"),
      1, 3    , .subset2(cyc_prob,"L1_next"),
      1, len_m, .subset2(cyc_prob,"L1_death"),
      2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
      2, 3    , .subset2(cyc_prob,"L1_next"),
      2, len_m, .subset2(cyc_prob,"L1_death")
    ),
    nrow  = 7,
    byrow = TRUE
  )
  M[vals_1L[,1:2]] <- vals_1L[,3]
  tr4[[c+1]] <- bigalgebra::dgemm(A = tr4[[c]],B = M)
}
t4 <- Sys.time()

t_th_mins <- ((as.numeric(t4 - t3) / 50) * TH) / 60
print(paste0("Matrix method for loop - one trace for one treatment line will take approximately ", round(t_th_mins,2), " minutes to run..."))

tr4 <- do.call(rbind,lapply(tr4,bigmemory::as.matrix))
f_mat_firstFewCells(tr4)
round(rowSums(tr4),12) == 1



# Sparse matrices ---------------------------------------------------------

# I have since learned that in R you can have sparse matrices if you use the Matrix package

library(Matrix)

M <- Matrix(
  as.matrix(M),
  sparse = TRUE
)


pb <- txtProgressBar(
  min = 1,
  max = TH-1,
  initial = 1,
  width = 50,
  style = 3,
  char = "="
)

out_list <- Reduce(
  x = 1:(TH-1),
  init = list(
    tpm = M,
    p   = list(Matrix(tr[1,],nrow = 1,sparse = T))
  ),
  accumulate = FALSE,
  f = function(prev,c) {
    setTxtProgressBar(pb, c)
    cyc_prob <- as.list(TPs[c,])
    vals_1L  <- matrix(
      c(
        1, 1    , .subset2(cyc_prob,"L1_stay"),
        1, 2    , .subset2(cyc_prob,"L1_disc"),
        1, 3    , .subset2(cyc_prob,"L1_next"),
        1, len_m, .subset2(cyc_prob,"L1_death"),
        2, 2    , .subset2(cyc_prob,"L1_stay") + .subset2(cyc_prob,"L1_disc"),
        2, 3    , .subset2(cyc_prob,"L1_next"),
        2, len_m, .subset2(cyc_prob,"L1_death")
      ),
      nrow  = 7,
      byrow = TRUE
    )
    prev$tpm[vals_1L[,1:2]] <- vals_1L[,3]
    
    prev$p[[c+1]] <- prev$p[[c]] %*% prev$tpm
    
    return(prev)
    
  }
)

close(pb)

# To get the trace:
sparse_trace <- matrix(unlist(lapply(out_list$p,as.numeric),use.names = F),nrow = TH,byrow = TRUE)
f_mat_firstFewCells(sparse_trace)
rowSums(sparse_trace)

