#' Apply treatment effect waning
#' 
#' @param surv_active  A survivor function as a vector giving the probability of
#'                     having not experienced the event of interest at the start
#'                     of each cycle
#' @param surv_ref     A reference survivor function - as the treatment effect
#'                     wanes, the hazard function will switch to the one given
#'                     by this curve
#' @param start_cycle  The cycle in which the treatment effect begins to wane
#' @param finish_cycle The cycle in which the treatment effect finishes waning
#' @param apply_waning Change this to FALSE to bypass waning altogether
#' @param wcc          With the default value (0.5), the proportion of treatment
#'                     effect which remains is based on the mid-point of the
#'                     cycle. Set `wcc = 0` to base it on the start of the cycle
#'                     or `wcc = 1` for the end of the cycle.
#' 
#' @details The hazard rate in the survivor function adjusted for treatment
#'   effect waning is assumed to be $r h_1(t_i) + (1-r) h_0(t_i)$, where $r$ is
#'   proportion of treatment effect remaining, $h_1$ is the active treatment
#'   hazard rate and $h_0$ is the reference treatment hazard rate. $r = 1$ for
#'   any cycle earlier than `start_cycle` and $r = 0$ for `finish_cycle` and
#'   subsequent cycles. For cycles from `start_cycle` and earlier than
#'   `finish_cycle`, linear interpolation is used based on the mid-point of the
#'   cycle assuming that $r = 1$ at the start of `start_cycle` and $r = 0$ at
#'   the start of `finish_cycle`. If `start_cycle` and `finish_cycle` are the
#'   same, then there is no interpolation as the effect wanes immediately.
#'   NOTE: The first cycle is numbered 0.
treatment_effect_waning <- function(
    surv_active, surv_ref, start_cycle, finish_cycle, apply_waning = TRUE,
    wcc = 0.5
  ) {
  
  if (apply_waning == FALSE) return(surv_active)
  
  n <- length(surv_active)
  
  if (n < start_cycle) return(surv_active)
  
  # Check inputs are valid
  stopifnot(start_cycle <= finish_cycle)
  stopifnot(length(surv_active) == length(surv_ref))
  
  res <- numeric(n)
  
  # Copy survival curve exactly before `start_cycle`
  res[1:(start_cycle-1)] <- surv_active[1:(start_cycle-1)]
  
  # Calculate per cycle hazard
  h_active <- log(surv_active[1:(n-1)]) - log(surv_active[2:n])
  h_ref    <- log(surv_ref[1:(n-1)]) - log(surv_ref[2:n])
  
  if (any(h_ref < h_active)) {
    total_warnings <- sum(h_ref < h_active)
    cycle_first_warning <- which.max(h_ref < h_active)
    warning_msg <- paste(
      "Hazard rate in the active treatment is more than the hazard rate in the",
      "reference treatment in", total_warnings, "cycles, the first of which is",
      (cycle_first_warning - 1)
    )
    warning(warning_msg)
  }
  
  # Calculate waning
  t <- seq_along(surv_active) - 1 + wcc
  r <- if (start_cycle < finish_cycle)
    approx(
      x    = c(0, start_cycle, finish_cycle, Inf),
      y    = c(1, 1, 0, 0),
      xout = t[1:(n-1)]
    )$y
  else
    rep(c(1, 0), times = c(start_cycle, n - start_cycle - 1))
  
  # Calculate adjusted hazards
  h_adj <- r * h_active + (1 - r) * h_ref
  p_adj <- exp(-h_adj)
  
  # Produce adjusted survivor function
  s_adj <- Reduce(
    f          = `*`,
    x          = p_adj,
    init       = 1,
    accumulate = TRUE
  )
  
  res[start_cycle:n] <- s_adj[start_cycle:n]
  
  res
  
}


#' Apply treatment effect waning MODIFIED TO INCLUDE ABSOLUTE WANING
#' 
#' @param surv_active  A survivor function as a vector giving the probability of
#'                     having not experienced the event of interest at the start
#'                     of each cycle
#' @param surv_ref     A reference survivor function - as the treatment effect
#'                     wanes, the hazard function will switch to the one given
#'                     by this curve
#' @param start_cycle  The cycle in which the treatment effect begins to wane
#' @param finish_cycle The cycle in which the treatment effect finishes waning
#' @param apply_waning Change this to FALSE to bypass waning altogether
#' @param wcc          With the default value (0.5), the proportion of treatment
#'                     effect which remains is based on the mid-point of the
#'                     cycle. Set `wcc = 0` to base it on the start of the cycle
#'                     or `wcc = 1` for the end of the cycle.
#' @param method       Apply treatment effect waning to the hazards (default)
#'                     or the absolute survival at time t (sometimes used in 
#'                     cost-effectiveness models)
#' 
#' @details The hazard rate in the survivor function adjusted for treatment
#'   effect waning is assumed to be $r h_1(t_i) + (1-r) h_0(t_i)$, where $r$ is
#'   proportion of treatment effect remaining, $h_1$ is the active treatment
#'   hazard rate and $h_0$ is the reference treatment hazard rate. $r = 1$ for
#'   any cycle earlier than `start_cycle` and $r = 0$ for `finish_cycle` and
#'   subsequent cycles. For cycles from `start_cycle` and earlier than
#'   `finish_cycle`, linear interpolation is used based on the mid-point of the
#'   cycle assuming that $r = 1$ at the start of `start_cycle` and $r = 0$ at
#'   the start of `finish_cycle`. If `start_cycle` and `finish_cycle` are the
#'   same, then there is no interpolation as the effect wanes immediately.
#'   NOTE: The first cycle is numbered 0.
treatment_effect_waning_with_absolute <- function(
    surv_active, surv_ref, start_cycle, finish_cycle, apply_waning = TRUE,
    wcc = 0.5, method = "h", if_cross_use_worst = TRUE
  ) {
  
  # method must be h or a (hazard or absolute)
  stopifnot(method %in% c("h", "a"))
  
  if (apply_waning == FALSE) return(surv_active)
  
  n <- length(surv_active)
  
  if (n < start_cycle) return(surv_active)
  
  # Check inputs are valid
  stopifnot(start_cycle <= finish_cycle)
  if(length(surv_active) != length(surv_ref)) {
    if (length(surv_active) > length(surv_ref)) {
      new_index <- (length(surv_ref)+1):length(surv_active)
      surv_ref[new_index] <- surv_active[new_index]
    } else {
      new_index <- (length(surv_active)+1):length(surv_ref)
      surv_active[new_index] <- surv_ref[new_index]
    }
  }
  
  res <- numeric(n)
  
  # Copy survival curve exactly before `start_cycle`
  res[1:(start_cycle-1)] <- surv_active[1:(start_cycle-1)]
  
  # If the method is a then the operation is very simple - simply interpolate
  # absolute survival from the active curve to the reference curve using the
  # start and end values
  
  if (method == "a") {
    
    # Instant waning - people instantly die when the treatment effect starts to
    # disappears all at once...
    if (start_cycle == finish_cycle) {
      res[start_cycle:length(res)] <- pmin(
        surv_ref[start_cycle:length(res)],
        surv_active[start_cycle:length(res)]
      )
      return(res)
    } else {
      # Gradual waning but on absolute survival - linearly interpolate absolute
      # survival between the two curves and return the result. Perform interpolation
      # take the values we need, use them to populate the rest of the vector.
      stopifnot(start_cycle < finish_cycle)
      t <- seq_along(surv_active) - 1 + wcc
      w_t <- approx(
        x    = c(0, start_cycle, finish_cycle, Inf),
        y    = c(1, 1, 0, 0),
        xout = t[1:(n-1)]
      )$y[start_cycle:n]
      
      if (if_cross_use_worst) {
        surv_ref <- pmin(surv_active,surv_ref)
      }
      
      res[start_cycle:n] <- 
        (surv_active[start_cycle:n] * (w_t)) + 
        (surv_ref[start_cycle:n] * (1-w_t))
      return(res)
    }
    
  } else if (method == "h") {
    
    
    # Calculate per cycle hazard
    h_active <- log(surv_active[1:(n-1)]) - log(surv_active[2:n])
    h_ref    <- log(surv_ref[1:(n-1)]) - log(surv_ref[2:n])
    
    if (any(h_ref < h_active)) {
      total_warnings <- sum(h_ref < h_active)
      cycle_first_warning <- which.max(h_ref < h_active)
      warning_msg <- paste(
        "Hazard rate in the active treatment is more than the hazard rate in the",
        "reference treatment in",
        total_warnings,
        "cycles, the first of which is",
        (cycle_first_warning - 1)
      )
      if (if_cross_use_worst) {
        warning_msg <- paste0(warning_msg," . Using the higher of the two hazards where they cross!")
        h_ref <- pmax(h_ref, h_active)
      }
      warning(warning_msg, immediate. = TRUE)
    }
    
    # Calculate waning
    t <- seq_along(surv_active) - 1 + wcc
    r <- if (start_cycle < finish_cycle) {
      approx(
        x    = c(0, start_cycle, finish_cycle, Inf),
        y    = c(1, 1, 0, 0),
        xout = t[1:(n-1)]
      )$y
    } else {
      rep(c(1, 0), times = c(start_cycle, n - start_cycle - 1))
    }
    
    # Calculate adjusted hazards
    h_adj <- r * h_active + (1 - r) * h_ref
    p_adj <- exp(-h_adj)
    
    # Produce adjusted survivor function
    s_adj <- cumprod(c(1,p_adj))
    
    # QC NOTE for TS: the above is identical to the below commented out lines 
    # s_adj <- Reduce(
    #   f          = `*`,
    #   x          = p_adj,
    #   init       = 1,
    #   accumulate = TRUE
    # )
    
    res[start_cycle:n] <- s_adj[start_cycle:n]
    
    return(res)
  } else {
    stop("Method should be h or a - this error shouldn't be possible!")
  }
}


# convenience functions ---------------------------------------------------

#' Translator function to return the argument needed for the waning function on method
#' using the dropdown option from excel
#' 
#' @param m either `ref trt hazard` or `ref trt abs surv` according to named range `apply_waning_to` in Excel
#' 
f_misc_twaning_methodTranslate <- Vectorize(function(m) {
  switch (m,
          "ref trt hazard" = return("h"),
          "ref trt abs surv" = return("a"),
          "0" = NA
  )
},USE.NAMES = FALSE)



# Function to apply in the cost-effectiveness model -----------------------

#' Function to apply treatment effect waning to all PLMTEs according to the table
#' provided in the excel inputs, under named range "R_table_TE_waning_settings".
#' 
#' @param st_list the list of all extrapolations AFTER propagation of relative efficacy, i.e. `p$surv$st`
#' @param tab_waning named range `R_table_TE_waning_settings` from excel
#' @param tab_eff_set named range `R_table_eff_data_settings` from excel
#'  
#' 
f_surv_twaning_apply <- function(st_list, tab_waning, tab_eff_set, verbose = FALSE) {
  
  # abbreviate the waning table and translate the method:
  t_w <- data.table(tab_waning)
  t_w$row_orig <- 1:nrow(t_w)
  t_e <- data.table(tab_eff_set)
  
  # Tidy up the excel table, filtering it down to just the yes:
  t_w$method <- f_misc_twaning_methodTranslate(t_w$apply.to)
  t_w <- t_w[!is.na(apply.to) & apply.waning == "Yes",]
  t_w$row_func <- 1:nrow(t_w)
  
  # filter the efficacy table down to just destinations and origins, as that's
  # all we need it for (d for destination o for origin):
  d_o <- t_e[,list(
    Include.in.this.analysis.,
    Population,
    Treatment.line,
    Molecule,
    End.point,
    Origin.population,
    Origin.line,
    Origin.treatment,
    Origin.trial,
    Origin.endpoint
  )]
  
  # Cycle down the rows of t_w, updating st_list each time using d_o to give us
  # the location of the reference curve and t_w to give us the location of the
  # active curve to change in st_list
  
  return(Reduce(
    x = 1:nrow(t_w),
    init = st_list,
    accumulate = FALSE,
    f = function(prev, tw_row) {
      
      tw <- t_w[tw_row,]
      
      dN <- list(
        pop      = tw$Population,
        line     = tw$Treatment.line,
        mol      = tw$Treatment,
        endpoint = tw$End.point
      )
      
      # Use the numeric destination id's to filter down d_o to the destination
      # (Note that there's only one for each dest, and it's not sensitive to trial)
      omatch <- d_o[Population == dN$pop & Treatment.line == dN$line & Molecule == dN$mol & End.point == dN$endpoint,]
      
      # if this curve isn't included in the analysis then nope
      if (omatch$Include.in.this.analysis. == "No") return(prev)
      
      # Use the reference table to populate the destination trial (for completeness and QC)
      dN$trial = omatch$Origin.trial
      
      # Use this information to derive the location within st_list for the
      # active survival to be "twaned", and then the reference curve to use:
      d <- list(
        pop = paste0("pop_",dN$pop),
        line = paste0("line_",dN$line),
        mol = paste0("mol_",dN$mol),
        trial = paste0("trial_",dN$trial),
        endpoint = paste0("endpoint_",dN$endpoint)
      )
      o <- list(
        pop = paste0("pop_",omatch$Origin.population),
        line = paste0("line_",omatch$Origin.line),
        mol = paste0("mol_",omatch$Origin.treatment),
        trial = paste0("trial_",omatch$Origin.trial),
        endpoint = paste0("endpoint_",omatch$Origin.endpoint)
      )
      
      # if verbose, then return a nice message:
      if (verbose) {
        f_misc_colcat(paste0(
          "TxWaning #",tw$row_orig," (",tw_row,"). d=.$",
          paste(unlist(d),collapse = "$"),
          " | o=.$",
          paste(unlist(o),collapse = "$"),
          " | m=",
          tw$apply.to,
          " | s=",
          tw$start,
          " | e=",
          tw$end
        ))
      }
      
      # Get our plmte to update:
      plmte <- f_misc_get_plmte(st_list,d)
      
      # Get the active curve and the reference curve
      ref_curve <- f_misc_get_plmte(st_list,o)$st
      
      stopifnot(!is.null(ref_curve))
      
      # perform treatment effect waning adjustment on the active curve, replacing
      # the element "st" in plmte with the result
      plmte$st <- treatment_effect_waning_with_absolute(
        surv_active  = plmte$st,
        surv_ref     = ref_curve,
        start_cycle  = tw$start,
        finish_cycle = tw$end,
        method       = tw$method)
      
      # slot the plmte back into prev, now that it has been udpated, and return
      prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
      return(prev)
    }
  ))
}


# Testing ground: ---------------------------------------------------------

if (FALSE) {
  
  # Run the model until you populate p$surv$st (i.e. the propagator). good idea
  # to then output that to a file so you can play with it without having to
  # rerun everything each time.
  p <- list()
  p$surv <- list()
  p$surv$st <- readRDS("./1_Data/example_st.rds")
  
  # Example inputs:
  st_act <- p$surv$st$pop_0$line_1$mol_1$trial_2$endpoint_0$st
  st_ref <- p$surv$st$pop_0$line_1$mol_7$trial_2$endpoint_0$st
  wan_st <- i$R_table_TE_waning_settings[Population == 0 & Treatment.line == 1 & Treatment == 1 & End.point == 0,]$start
  wan_end <- i$R_table_TE_waning_settings[Population == 0 & Treatment.line == 1 & Treatment == 1 & End.point == 0,]$end
  method<- "h"
  
  st_wan <- treatment_effect_waning_with_absolute(
    surv_active = st_act,
    surv_ref    = st_ref,
    start_cycle = wan_st,
    finish_cycle = wan_end,
    apply_waning = TRUE,
    method = method)
  
  plot(st_wan, type="l")
  lines(st_act, col="green")
  lines(st_ref, col="red")
  
}

