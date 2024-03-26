# ~~ Reference curves -----------------------------------------------------

# reference curves require generating sets of parameter estimates using
# the parameters and their variance-covariance matrices. To that end, we define 
# several functions:
# 
#   - Efficiently drawing from multivariate normal using eigenvalue appraoch (per MASS package)
#   - Generating one set of reference curves across all PLMTEs
# 

#' Function which takes uniform random draws and produces mutlinorminv results
#' 
#' @param rands uniform random draws to be entered into this function. Ensures reproducibility
#' @param mu named vector of values. if no names then the matrix output will have none
#' @param sigma variance-covariance matrix
#' 
#' @details credit to `MASS::mvrnorm` which has the code but doesn't let you insert
#'          uniform draws.
#' 
#' 
f_qmultinorminv <- function(rands,mu,sigma) {
  p    <- length(mu)
  n    <- length(rands) / p
  X    <- matrix(qnorm(rands), n)
  eS   <- eigen(sigma)
  ev   <- eS$values
  out <- t(mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X))
  colnames(out) <- names(mu)
  return(out)
}

#' Function to quickly translate from description to number using the excel 
#' lookup tables (provided in the function)
f_psa_trans_desc2id <- function(id_list, lookups) list(
  pop      = paste0("pop_",lookups$pop[Description == id_list$pop,]$Number),
  line     = paste0("line_",lookups$line[Description == id_list$line,]$Number),
  mol      = paste0("mol_",lookups$mol[Description == id_list$mol,]$Number),
  trial    = paste0("trial_",lookups$trial[Description == id_list$tr,]$Number),
  endpoint = paste0("endpoint_",lookups$endpoint[Description == id_list$endpoint,]$Number)
)
f_psa_trans_id2numb <- function(id_list) list(
  pop      = as.numeric(gsub("pop_","",id_list$pop)),
  line     = as.numeric(gsub("line_","",id_list$line)),
  mol      = as.numeric(gsub("mol_","",id_list$mol)),
  trial    = as.numeric(gsub("trial_","",id_list$trial)),
  endpoint = as.numeric(gsub("endpoint_","",id_list$endpoint))
)


#' Function to draw parametric model parameters for all reference curves for all
#' distributions
#' 
#' 
#' 
f_PSA_drawFSParams <- function(surv_regs, n_psa, lookups, return_rands = FALSE, verbose = FALSE) {
  lapply(surv_regs, function(risk_pop) {
    lapply(risk_pop, function(tr_line) {
      lapply(tr_line, function(mol) {
        lapply(mol, function(trial) {
          lapply(trial, function(endpoint) {
            if (is.null(endpoint$fs_fits)) {
              return(NULL)
            } else {
              fits <- .subset2(endpoint,"fs_fits")
              # Since we are only going to use one of these distributions, 
              # we can apply one set of random numbers 
              max_length <- max(lengths(fits))
              rands      <- runif(max_length * n_psa)
              plmte_id   <- f_psa_trans_desc2id(
                id_list = endpoint[c("pop", "line", "mol", "tr", "endpoint")],
                lookups = lookups
              )
              if(verbose) f_misc_colcat(paste0(
                "RNG for PSA - reference curve parameters | ",
                paste(unlist(plmte_id), collapse = " | ")
              ))
              
              fs_nam <- names(fits)
              names(fs_nam) <- fs_nam
              
              out <- lapply(fs_nam, function(fs_finam) {
                
                fs_fit <- .subset2(.subset2(endpoint,"fs_fits"),fs_finam)
                
                if("logical" %in% class(.subset2(fs_fit,"vcov"))) {
                  warning(paste0(
                    "There is no variance-covariance matrix available for ",
                    paste(unlist(plmte_id), collapse = " | "),
                    " for ",fs_finam,".",
                    " I cannot compute probabilsitic parameters and therefore have to use the mean!")
                  )
                  return(matrix(rep(.subset2(fs_fit,"coefs"),n_psa),nrow=n_psa,byrow = TRUE))
                }
                f_qmultinorminv(
                  rands = rands[1:(length(.subset2(fs_fit,"coefs")) * n_psa)],
                  mu    = .subset2(fs_fit,"coefs"),
                  sigma = .subset2(fs_fit,"vcov")
                )
              })
              # If the user wants the uniform draws to replicate, then it can be done.
              if (return_rands) {
                return(list(draws=out, id=plmte_id, rands = rands))
              } else {
                return(list(draws = out, id = plmte_id))
              }
            }
          })
        })
      })
    })
  })
}



#' Function to filter down the PSA PSM parameters to the distributions selected
#' in the excel book (and therefore `i`)
f_psa_surv_params_filter <- function(psa_psm_param, excel_eff_table, lookups) {
  
  # If it's not a data.table, make it one:
  if (!"data.table" %in% class(excel_eff_table)) excel_eff_table <- data.table(excel_eff_table)
  
  # Filter down to curve selections:
  curve_select <- excel_eff_table[Include.in.this.analysis. == "Yes" & Effectiveness.data.source == "Trial survival analysis", list(
    Population,
    Treatment.line,
    Molecule,
    End.point,
    Curve.fit..for.survival.analysis.
  )]
  
  # Now, we can use each id for each set of draws to filter down this table
  # to the relevant row. With that row we can select the appropriate 
  # parametric model and ONLY compute extrapolations for that one, then take
  # the colsums or trapz after extrapolating only the relevant model
  lapply(psa_psm_param, function(risk_pop) {
    lapply(risk_pop, function(tr_line) {
      lapply(tr_line, function(mol) {
        lapply(mol, function(trial) {
          lapply(trial, function(endpoint) {
            if (is.null(endpoint)) {
              return(NULL)
            } else {
              
              # Use the function we made above to get numbers we can use to filter
              # the table:
              idn <- f_psa_trans_id2numb(endpoint$id)
              tabmatch <- curve_select[Population == idn$pop & Treatment.line == idn$line & Molecule == idn$mol & End.point == idn$endpoint,]
              
              if (nrow(tabmatch) == 0) return (NULL)
              
              curve <- lookups$dist[match(tabmatch$Curve.fit..for.survival.analysis., lookups$dist$Description),]$RCC_input_desc
              endpoint$draws <- endpoint$draws[[curve]]
              endpoint$id$dist <- curve
              return(endpoint)
            }
          })
        })
      })
    })
  })
}



#' Computes the approximate lambda rate (1/AUC) using either fast or accurate
#' method. Respects the time unit of the CE model in cycles, so conversion comes
#' after.
#' 
#' @param psa_params the result of function f_PSA_drawFSParams, which draws from a multivariate normal 
#' @param method either sum or trap. sum sums the columns, trap uses the trapezoidal rule
#' 
#' @details Cycle down all levels of PLMTE. if there are parameters there, convert
#'          those parameters into approximate exponental rates by extrapolating the lines
#'          and then calculating AUC using sum or fast method. This is computed in
#'          model cycles, meaning the lambda for the exponential is in model cycles as
#'          the time unit of the analysis.
#' 
f_psa_approx_lambda <- function(psa_params, method = "sum", th, disc=NULL) {
  # Checking that the method is one of those that are allowed:
  stopifnot(method %in% c("sum", "trap"))
  
  # Set up the time vector to avoid repeatedly making it thousands of times:
  t <- 0:th
  
  if (is.null(disc)) {
    disc <- rep(1,th+1)
  } else {
    stopifnot(length(disc)== th+1)
  }
  
  # Cycle through the levels of PLMTE, if nothing to do do nothing, otherwise
  lapply(psa_params, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        lapply(mol, function(tr) {
          lapply(tr, function(plmte) {
            
            if (is.null(plmte)) return(NULL)
            
            # So if there's something in this plmte it will be a list of
            # different parametric distributions. each element is a set of
            # parameters per row, with each row being one PSA iteration.
            
            
            dr <- .subset2(plmte,"draws")
            par_nam <- colnames(dr)
            drt <- t(dr)
            dis <- plmte$id$dist
            npsa <- ncol(drt)
            
            # Two methods - one for speed one for accuracy:
            if (method == "sum") {
              vapply(
                X = 1:npsa,
                FUN = function(cyc) 1 / sum(f_extrapolate(t, drt[, cyc], dis) * disc),
                FUN.VALUE = numeric(1)
              )
            } else {
              # More precise method using trapezoidal integration 
              vapply(
                X = 1:npsa,
                FUN = function(cyc) 1 / trapz(t,f_extrapolate(t, drt[, cyc], dis) * disc),
                FUN.VALUE = numeric(1)
              )
            }
          })
        })
      })
    })
  })  
}


#' Function to apply `lambda` and `t` to extrapolate an exponential line
f_psa_exp <- function(t,lamda) exp(-lamda*t)



#' function to compute the first TP from a lambda for an exponential and therefore
#' all TPs as they are time invariant.
#' 
#' @param PSA_est_Lambdas result of function `f_psa_approx_lambda`
#' 
#' @details Lambda can be translated to transition probability TP like so:
#'  
#'  TP_t = 1-s(t)/s(t-1)
#'  Given that the rate is constant with exponential, all TP for a given endpoint
#'  are equal to TP_t, that is TP = 1-s(t)/s(t-1). 
#' 
#'  At cycle 1 (given cycle 0 is model start), TP_t = 1-(s(t) / s(t-1)), but s(t-1)
#'  is known to be 1. Therefore TP = 1-s(1).
#' 
#'  Thus, the entire set of lambdas can be converted to TPs by cycling through them
#'  and applying 1-f_psa_exp(1,lambda)
#'  
f_psa_lambda2TP <- function(PSA_est_Lambdas) {
  lapply(PSA_est_Lambdas, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        lapply(mol, function(tr) {
          lapply(tr, function(plmte) {
            if(is.null(plmte)) {
              return(NULL)
            } else {
              1-f_psa_exp(1,plmte)
            }
          })
        })
      })
    })
  })  
}
f_psa_lambda2St <- function(PSA_est_Lambdas,t) {
  TH <- length(t)
  lapply(PSA_est_Lambdas, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        lapply(mol, function(tr) {
          lapply(tr, function(plmte) {
            if(is.null(plmte)) {
              return(NULL)
            } else {
              vapply(X = 1:length(plmte), FUN.VALUE = numeric(TH), FUN = function(psa_it) f_psa_exp(t,plmte[psa_it]))
            }
          })
        })
      })
    })
  })  
}




#' Function taking extrapolated survival in different treatment lines and 
#' either computing TP out of that state or lambda for exponential approximation.
f_psa_collapse_st_lambda2lplus <- function(st, th, disc, dfacQ = NULL) {
  
  # If the user wants discounted lambda rates they have to provide the discount factor!
  if (disc) stopifnot(!is.null(dfacQ))
  
  # Drill down from the top level to the line level. if 1L compute TP otherwise
  # compute lambda:
  lapply(st, function(popu) {
    
    na_li <- names(popu)
    names(na_li) <- na_li
    
    lapply(na_li, function(li_lab) {
      
      li <- popu[[li_lab]]
      
      # li <- st$pop_0$line_1
      
      # Now, if we are in first-line, we do something differently as we 
      # simply calculate TP
      # 
      # If we are in 2L+ we compute lambda (1/int(s(t)))
      if (li_lab == "line_1") {
        # This is 1st line, so keep the shape and compute the TP for it
        lapply(li, function(mol) {
          lapply(mol, function(tr) {
            lapply(tr, function(endp) {
              # If there's no s(t) to manipulate, do nothing
              if (!"st" %in% names(endp)) {
                return(endp)
              } else {
                endp$tp <- 1 - (endp$st / shift(endp$st,fill = 1))
                endp$st <- NULL
                return(endp)
              }
            })
          })
        })
      } else {
        # this is the 2L+, so compute lambda
        lapply(li, function(mol) {
          lapply(mol, function(tr) {
            lapply(tr, function(endp) {
              # cat(paste0(paste(endp$dest,collapse = "$"),"\n"))
              # If there's no s(t) to manipulate, do nothing
              if(!"include" %in% names(endp)) {
                return(endp)
              } else if (endp$include == FALSE) {
                return(endp)
              } else if (!"st" %in% names(endp)) {
                if(length(endp$fp$HR)==th+1) endp$fp$HR <- NULL
                return(endp)
              } else if (all(length(endp$st)==1,is.na(endp$st[1]))) {
                if(length(endp$fp$HR)==th+1) endp$fp$HR <- NULL
                return(endp)
              } else {
                # If discounted lambda, use the discount factor, if not, don't!
                if (disc) {
                  endp$lambda <- 1 / trapz(0:th,endp$st * dfacQ)
                } else {
                  endp$lambda <- 1 / trapz(0:th,endp$st)
                }
                if(length(endp$fp$HR)==th+1) endp$fp$HR <- NULL
                endp$st <- NULL
                return(endp)
              }
            })
          })
        })
      }
    })
  })
}




#' Function to generate all PSA iterations for all hazard ratios in one go within
#' the PLMTE structure.
#' 
#' @details cycles down the excel table and for each row makes npsa normal draws
#' for the hazard ratio to apply in that position. These HRs are then entered
#' into the efficacy networks during the assumption step. 
#' 
f_psa_assumptionsTab_genHRs <- function(excel_efficacy_table,npsa) {
  excel_efficacy_table$r <- 1:nrow(excel_efficacy_table)
  tab <- data.table(excel_efficacy_table)[Effectiveness.data.source == "Apply HR to",]
  
  # Cycle down the table one row at a time generating HRs for all PSA iterations
  # for all of those which are informed by those HRs. This assumes that all HRs
  # are independent from each other
  # 
  lapply(1:nrow(tab), function(HR_row) {
    list2env(as.list(tab[HR_row,]),envir = environment())
    list(
      dest = list(
        pop      = paste0("pop_",      Population),
        line     = paste0("line_",     Treatment.line),
        mol      = paste0("mol_",      Molecule),
        trial    = paste0("trial_",    Origin.trial),
        endpoint = paste0("endpoint_", End.point)
      ),
      orig = list(
        pop      = paste0("pop_",      Origin.population),
        line     = paste0("line_",     Origin.line),
        mol      = paste0("mol_",      Origin.treatment),
        trial    = paste0("trial_",    Origin.trial),
        endpoint = paste0("endpoint_", Origin.endpoint)
      ),
      # hr = rnorm(npsa,HR.to.apply,(HR.95..CI..UCL. - HR.to.apply)/1.96)
      hr = rlnorm(npsa, log(HR.to.apply),  estSDlog(HR.to.apply, HR.95..CI..LCL., HR.95..CI..UCL.))
      
    )
    
  })
}


f_misc_plmte2StringID <- function(plmte_id_list) paste(unlist(plmte_id_list),collapse = "$")


#' function to go into demographics and pull out just the PSA iteration
#' requested
f_psa_get_it_demo <- function(demo,it) {
  lapply(demo, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(cata) {
        if (class(cata) == "numeric") {
          return(cata[it])
        } else {
          cata
        }
      })
    })
  })
}




#' Function for the lambda approximation which provides vector of costs for 1L
#' and per-cycle costs for 2L+
#' 
#' @details This function should be used to collapse the cost object for ONE 
#' PSA RUN. Each PSA version of the cost object p$costs$mk should contain the
#' exact same structure as the deterministic model. This function collapses
#' the data to just what's needed for one PSA iteration.
#' 
f_psa_lambda_cost <- function(cost) {
  lapply(cost, function(category) {
    li_lab <- names(category)
    names(li_lab) <- li_lab
    lapply(li_lab, function(li_nam) {
      li <- category[[li_nam]]
      
      # if first-line we need the full vector, if not we need per cycle:
      if (li_nam == "line_1") {
        return(li)
      } else {
        return(lapply(li,mean))
      }
    })
  })
}


#' Lambda approximation version of TP calculator
f_psa_pf_mk_ComputeTPs <- function(tp_list) {
  out <- lapply(1:length(tp_list), function(trt_li) {
    
    # Figure out the time horizon:
    tl    <- tp_list[[trt_li]]
    
    # In this version of the function, if treatment line is 1, then we do what
    # we did in the deterministic model, otherwise the rates are all fixed
    # from lambda
    if (trt_li < length(tp_list)) {
      
      if (trt_li == 1) {
        # Pull out the endpoints for ease of reading. These are time-varying
        tp_TTD <- tl$TTD$tp
        tp_TTP <- tl$TTP$tp
        tp_PFS <- tl$PFS$tp
      } else {
        # lambda-based - TP=1-(s(t) / s(t-1)); s(t-1) = 1, so 1-s(t)
        lambda_list <- tl$OS$lambda
        tp_TTD <- 1-f_psa_exp(1,tl$TTD$lambda)
        tp_TTP <- 1-f_psa_exp(1,tl$TTP$lambda)
        tp_PFS <- 1-f_psa_exp(1,tl$PFS$lambda)
      }
      
      # Apply assumptions to figure out TP from each state to each other state:
      disc <- tp_TTD
      
      # ON-treatment transition probabilities
      
      # probability of death directly from on treatment state:
      death_on <- 1 - tp_TTP - (1-tp_PFS)
      
      # Going directly onto next line from the on treatment state
      
      # The probability of discontinuation removing the probability of any other event
      disc_only <- disc - tp_TTP - death_on
      disc_only_adj <- pmax(disc_only, 0)
      
      # The probability of going to next therapy but not death directly from on-treatment
      next_on  <- tp_TTP - (disc_only_adj - disc_only)
      
      # the probability of staying on treatment is simply the inverse of the probability
      # of going off treatment for any reason
      stay_on  <- 1 - disc
      
      # The probability of discontinuation when on treatment
      disc_on  <- disc_only_adj
      
      # Just to note some equivalence here:
      # 1 - next_on - disc - stay_on == 1 - tp_TTP - (1-tp_PFS)
      
      next_off  <- next_on
      death_off <- death_on
      stay_off  <- 1 - death_on - next_on
      
    } else {
      # Lambda approximate rates in the last line
      tp_OS  <- 1-f_psa_exp(1,tl$OS$lambda)
      
      stay_on  <- 1-tp_OS
      disc_on  <- 0
      next_on  <- 0
      death_on <- tp_OS
      
      next_off  <- 0
      death_off <- 0
      stay_off  <- 0
    }
    
    
    
    # return a list of TPs for each of the 7 required to compute the Markov model
    return(list(
      disc_on   = disc_on,
      stay_on   = stay_on,
      next_on   = next_on,
      death_on  = death_on,
      stay_off  = stay_off,
      next_off  = next_off,
      death_off = death_off
    ))
    
  })
  names(out) <- c(paste0("L",1:(length(tp_list)-1)),"NT")
  return(out)
}



f_psa_lambda_M_compile <- function(TPM, max_active, TPs) {
  
  TH <- length(TPs$L1$disc_on)
  nline <- length(TPs)
  
  # Strip down TPM to the size required
  included_index <- c(1:((nline-1) * 2), ((max_active+1)*2)-1,(max_active+1)*2)
  ATL <- nline-1
  
  # First, we strip down TPM to suit the number of ATLs we have
  m <- TPM[included_index,included_index]
  dm <- dim(m)[[1]]
  
  # Now that we have the empty version of the matrix to fill in, let's populate
  # it with some transitions.
  # 
  # We leave first-line transitions blank because they get populated at the end
  # of the function (when we expand it out to per-cycle). 2L+ are all time-invariant
  
  # Any line can be NT, so we use the names within TP to capture that, and create
  # a matrix of co-ordinates and values to put in those co-ordinates, then put
  # them there in one command
  coord <- do.call(
    rbind,
    lapply(2:nline, function(tx_line) {
      # top left element of this set is linex2 - 1
      tl <- (tx_line * 2) - 1
      tr <- TPs[[tx_line]]
      if (tx_line == nline) {
        # This line is BSC. There is ONLY the probability of death or stay
        # and there is no "next". In this situation only top left and top left+1
        # are the co-ordinates:
        matrix(
          c(
            tl  , tl  , .subset2(tr,"stay_on"),
            tl  , tl+1, .subset2(tr,"death_on"),
            tl+1, tl+1, 1
          ),
          ncol=3,
          byrow = TRUE
        )
      } else {
        # this is not the last line of therapy, so we need the full suite:
        matrix(
          c(
            tl  , tl  , .subset2(tr,"stay_on"),
            tl  , tl+1, .subset2(tr,"disc_on"),
            tl  , tl+2, .subset2(tr,"next_on"),
            tl  , dm  , .subset2(tr,"death_on"),
            tl+1, tl+1, .subset2(tr,"stay_off"),
            tl+1, tl+2, .subset2(tr,"next_off"),
            tl+1, dm  , .subset2(tr,"death_off")
          ),
          ncol=3,
          byrow = TRUE
        )
      }
      
    })
  )
  m[coord[,1:2]] <- coord[,3]
  
  L1_trans <- .subset2(TPs,"L1")
  list2env(L1_trans, envir = environment())
  # m is now the time-invariant part of the matrix which can be replicated TH
  # times:
  return(lapply(1:TH, function(cyc) {
    
    # Like with the time invariant, but with the time-varying bit, make a 
    # co-ordinate matrix:
    coord <- matrix(
      c(
        1, 1 , .subset2(stay_on,cyc),
        1, 2 , .subset2(disc_on,cyc),
        1, 3 , .subset2(next_on,cyc),
        1, dm, .subset2(death_on,cyc),
        2, 2 , .subset2(stay_off,cyc),
        2, 3 , .subset2(next_off,cyc),
        2, dm, .subset2(death_off,cyc)
      ),
      ncol=3,
      byrow = TRUE
    )
    m[coord[,1:2]] <- coord[,3]
    return(m)
  }))
}




#' Patient flow calculator for probabilistic sensitivity analysis using lambda approximation
#' to remove the need for tunnel states
#' 
#' @details this function has less comments in it, see f_pf_computePF_mk for 
#' fully explained. Where this function differs from f_pf_computePF_mk there are
#' comments in here. 
#' 
f_psa_pf_computePF_mkLambda <- function(pops,
                                        basic,
                                        demo,
                                        sequences,
                                        survival,
                                        costs,
                                        util,
                                        ae,
                                        eff_table,
                                        verbose,
                                        include_plots,
                                        just_nlines,
                                        just_seq) {
  
  # Some basic inputs for easy reference and avoiding repitition in the code:
  dfac_q <- .subset2(basic,"discFacQ")
  dfac_c <- .subset2(basic,"discFacC")
  
  # abbreviate cost per cycle (cpc) and costs on initation (coi). we only need the
  # mols which we are using for this sequence. one off are for all.
  cpc  <- .subset2(costs,"per_cycle")
  coo <- .subset2(costs,"one_off")
  coi  <- .subset2(coo,"Prog")
  cod  <- .subset2(coo,"Death")
  
  # Getting more specific:
  cpc_d     <-.subset2(cpc,"drug")
  cpc_a     <-.subset2(cpc,"admin")
  cpc_m_on  <-.subset2(cpc,"mru_on")
  cpc_m_off <-.subset2(cpc,"mru_off")
  
  
  # Now simulate the treatment sequences:
  max_active <- max(sapply(sequences$qc, ncol)) - 1
  lookup <- basic$lookup
  pop_map <- lookup$pop_map
  th <- basic$th+1
  
  # Full size empty TPM, then expand out to list of TH length
  TPM <- matrix(
    0,
    nrow = ((max_active+1)*2), 
    ncol = ((max_active+1)*2),
    dimnames = lapply(1:2, function(x) {
      c(
        unlist(lapply(1:max_active, function(trt_line) {
          paste0("L",trt_line,c("_on", "_off"))
        })),
        "BSC", "dead"
      )
    })
  )
  
  # population level inputs, these change at each population:
  lapply(pops, function(popu) {
    # population ids
    popu_n <- as.integer(gsub("pop_","",popu))
    overall_pop <- as.list(pop_map[Overall.population.number == popu_n,])
    rpop_n <- overall_pop$Risk.population.number
    rpop <- paste0("pop_",rpop_n)
    spop_n <- overall_pop$Sequencing.population.number
    spop <- paste0("pop_",spop_n)
    
    # demographics
    demog <- demo[[rpop]]
    
    # Sequence list derivation
    L1_inc_rpop <- sort(unique(eff_table[Treatment.line == 1 &
                                           Population == rpop_n &
                                           Include.in.this.analysis. == "Yes", ]$Molecule))
    L2p_inc_rpop <- lapply(2:max_active, function(line_n) {
      sort(unique(eff_table[Treatment.line == line_n &
                              Population == 0 & Include.in.this.analysis. == "Yes", ]$Molecule))
    })
    trt_filter <- c(list(L1_inc_rpop),L2p_inc_rpop)
    rm(L1_inc_rpop)
    rm(L2p_inc_rpop)
    ateachline <- sequences$n[[spop]]
    ateachline$nlines <- ncol(ateachline)-rowSums(is.na(ateachline))
    molecule_list_full <- lapply(2:(ncol(ateachline)-1), function(n_lines) {
      ateachline[nlines==n_lines,1:n_lines]
    })
    sequence_list <- lapply(molecule_list_full, function(lines_amount) {
      Reduce(
        x = 1:(ncol(lines_amount)-1),
        init = lines_amount,
        accumulate = FALSE,
        f = function(prev, treatment_line) {
          
          # Find out which molecules are in that line column and filter the table
          molecules_allowed <- trt_filter[[treatment_line]]
          which_to_allow <- which(prev[[treatment_line]] %in% molecules_allowed)
          return(prev[which_to_allow,])
        }
      )
    })
    
    # account for running subset
    if (!is.null(just_nlines)) {
      sequence_list <- sequence_list[just_nlines]
      names(sequence_list) <- paste0("active_lines_",just_nlines)
    } else {
      names(sequence_list) <- paste0("active_lines_",1:length(sequence_list))
    }
    
    # gpop matrix to multiply HRQoL by over time
    util_gpop_mat <- matrix(
      util$gpop[[rpop]], 
      nrow = basic$th+1,
      ncol = (max_active*2) +1 + 1
    )
    
    # Sequence "block" populations - treatments with specific amounts of ATLs:
    return(lapply(sequence_list, function(amount_of_lines) {
      
      seq_id <- 1:nrow(amount_of_lines)
      names(seq_id) <- paste0("seq_",seq_id)
      if (!is.null(just_seq)) {seq_id <- seq_id[just_seq]}
      
      
      # Treatment sequence level: specific population, ATLs, cycling down:
      lapply(seq_id, function(seq_n) {
        
        number_of_lines <- ncol(amount_of_lines)
        
        # Pull out this specific treatment sequence:
        trt_seq <- structure(paste0("mol_",amount_of_lines[seq_n,]),.Names=names(amount_of_lines[seq_n,]))
        trts_n <- as.integer(gsub("mol_","",trt_seq))
        
        trt_seq_plotLab <- trt_seq
        names(trt_seq_plotLab) <- gsub("line_","L",names(trt_seq))
        names(trt_seq_plotLab)[length(names(trt_seq_plotLab))] <- "NT"
        
        treatment_names <- basic$lookup$ipd$mol[match(trts_n,basic$lookup$ipd$mol$Number),Description]
        treatment_abbr  <- basic$lookup$ipd$mol[match(trts_n,basic$lookup$ipd$mol$Number),RCC_input_desc]
        
        if (verbose) cat(paste0(
          "PSA #", util$hsuv$iteration[1],
          " | ", popu,
          " | ", paste(treatment_names,collapse = "->"),"\n"
        ))
        
        
        # In this case we provide survival$tp instead of survival$st. This returns
        # the PLMTE format extraps, but instead of st for each of them there is
        # tp at first-line and lambda at 2L+. 
        extraps <- f_seq_extrapCollector(
          treatment_sequence = trt_seq,
          st                 = survival$tp,
          lookups            = lookup,
          pop_n              = rpop_n,
          pop_0_2Lp          = TRUE
        )
        
        # For the lambda approximation, the TP calculations are the same
        # methodologically but the way to GET those TPs is different
        TPs      <- f_psa_pf_mk_ComputeTPs(tp_list = extraps$st)
        
        # Correct for non-finite TPs:
        TPs <- lapply(TPs, function(li) {
          lapply(li, function(tp) {
            tp[!is.finite(tp)] <- 0
            tp
          })
        })
        
        # For computing L1 off trt entrants
        L1_onoff_adjustment <- TPs$L1$death_on + TPs$L1$next_on
        
        # Compile TH Ms
        M <- f_psa_lambda_M_compile(
          TPM = TPM,
          max_active = max_active,
          TPs = TPs
        )
        
        # Set up the baseline population (all on 1L on trt)
        p_0 <- numeric(dim(M[[1]])[[1]])
        p_0[1] <- 1
        
        # Undiscounted trace:
        TRACE <- matrix(
          data = unlist(
            Reduce(
              x = 1:p$basic$th,
              init = p_0,
              accumulate = TRUE,
              f = function(prev_pop, cyc) {
                prev_pop %*% M[[cyc]]
              }
            ),
            use.names = FALSE
          ),
          byrow = TRUE,
          ncol = dim(M[[1]])[[1]],
          dimnames = list(NULL,dimnames(M[[1]])[[1]])
        )
        
        tTRACE <- t(TRACE)
        
        M_nam <- dimnames(M[[1]])[[1]]
        
        # compute the % of the initial cohort transitioning into each state
        # at each cycle beyond the 1st state. This is used to assign one-off
        # costs/QALY losses.
        ENTRANTS <- matrix(
          unlist(lapply(1:ncol(tTRACE), function(cyc) {
            dm <- dim(M[[cyc]])[1]
            coord <- rbind(
              cbind(1:(dm-1),2:dm),
              cbind(seq(1,dm-2,2),seq(1,dm-2,2)+2)
            )
            coord <- coord[order(coord[,1]),]
            coord <- cbind(coord,M[[cyc]][coord])
            up_one_trans <- unlist(lapply(unique(coord[,2]), function(x) sum(coord[coord[,2]==x,3])))
            p_tm1 <- tTRACE[1:(nrow(tTRACE)-1),cyc]
            return(as.numeric(p_tm1 * up_one_trans))
          }),use.names = FALSE),
          nrow = ncol(tTRACE),
          byrow = TRUE,
          dimnames = list(NULL,M_nam[2:length(M_nam)])
        )
        
        # Discounted versions:
        DTRACEQ    <- TRACE * dfac_q
        DTRACEC    <- TRACE * dfac_c
        tDTRACEQ   <- t(DTRACEQ)
        tDTRACEC   <- t(DTRACEC)
        DENTRANTSQ <- ENTRANTS * dfac_q
        DENTRANTSC <- ENTRANTS * dfac_c
        
        # drop some values to save memory
        rm(extraps, M, TPs, tTRACE, M_nam)
        
        # Generate something similar to the consolidated_trace object
        # we had for the deterministic model:
        CT <- list(
          OS = rowSums(TRACE[,1:(length(p_0)-1)]),
          full_lines = TRACE,
          entrants = ENTRANTS
        )
        CT_disc <- list(
          OS_q = rowSums(DTRACEQ[,1:(length(p_0)-1)]),
          OS_c = rowSums(DTRACEC[,1:(length(p_0)-1)]),
          full_lines_q = DTRACEQ,
          full_lines_c = DTRACEC,
          entrants_q = DENTRANTSQ,
          entrants_c = DENTRANTSC
        )
        
        # We already made the cost elements once so now we use them
        
        # cycle through the treatment sequence, computing the different cost components:
        
        pf_costs <- lapply(1:length(trt_seq), function(line_n) {
          mol <- trt_seq[line_n]
          
          if (line_n == 1) {
            # first-line - vector multiplied by cost
            nam_mol <- names(mol)
            
            tr_1l_on <- TRACE[,1]
            tr_1l_off <- TRACE[,2]
            
            drug_cost    <- cpc_d[[nam_mol]][[mol]] * tr_1l_on
            admin_cost   <- cpc_a[[nam_mol]][[mol]] * tr_1l_on
            ruc_on       <- cpc_m_on[[nam_mol]][[mol]] * tr_1l_on
            ruc_off      <- cpc_m_off[[nam_mol]][[mol]] * tr_1l_off
            
            # Discounted version:
            dtr_1l_on <- TRACE[,1] * dfac_c
            dtr_1l_off <- TRACE[,2] * dfac_c
            
            ddrug_cost    <- cpc_d[[nam_mol]][[mol]] * dtr_1l_on
            dadmin_cost   <- cpc_a[[nam_mol]][[mol]] * dtr_1l_on
            druc_on       <- cpc_m_on[[nam_mol]][[mol]] * dtr_1l_on
            druc_off      <- cpc_m_off[[nam_mol]][[mol]] * dtr_1l_off
            
            rm(tr_1l_on,tr_1l_off,dtr_1l_on, dtr_1l_off)
            
          } else if (line_n < length(trt_seq)) {
            
            # In the lambda approximated world, the costs are time invariant in 2L+
            nam_mol <- names(mol)
            tr_on  <- TRACE[,(line_n*2)-1]
            tr_off <- TRACE[,line_n*2]
            e_on   <- ENTRANTS[,(line_n*2)-2]
            
            dc     <- cpc_d[[nam_mol]][[mol]]
            ac     <- cpc_a[[nam_mol]][[mol]]
            ru_on  <- cpc_m_on[[nam_mol]][[mol]]
            ru_off <- cpc_m_off[[nam_mol]][[mol]]
            
            drug_cost  <- tr_on * dc
            admin_cost <- tr_on * ac
            ruc_on     <- (tr_on * ru_on) + (e_on * coi)
            ruc_off    <- tr_off * ru_off
            
            dtr_on  <- DTRACEC[,(line_n*2)-1]
            dtr_off <- DTRACEC[,line_n*2]
            de_on   <- DENTRANTSC[,(line_n*2)-2]
            
            ddrug_cost  <- dtr_on * dc
            dadmin_cost <- dtr_on * ac
            druc_on     <- (dtr_on * ru_on) + (de_on * coi)
            druc_off    <- dtr_off * ru_off
            
            rm(tr_on, tr_off,e_on, dtr_on,dc,ac,ru_on,ru_off)
            
          } else {
            
            # Future switch? Add in one-off costs upon entering this treatment line:
            # In future, link this to the excel switch for applying this cost
            # upon starting BSC.
            if (TRUE) {
              # Sometimes BSC is coming after any active treatment line.
              # The cost should come from the last active treatment line's BSC
              # values:
              if (is.null(cpc$mru_on[[names(mol)]])) {
                # reduce line id by 1 just for pulling out the right data
                # from cpc:
                namline_n <- as.numeric(gsub("line_","",names(mol)))
                line_nam_temp <- paste0("line_",namline_n-1)
                
                # Now use this adjusted name to pull out the right values
                mru_on_tx <- cpc$mru_on[[line_nam_temp]][[mol]]
              } else {
                mru_on_tx <- cpc$mru_on[[names(mol)]][[mol]]
              }
            }
            
            ruc_on  <- (mru_on_tx * TRACE[,"BSC"]) + (ENTRANTS[,"BSC"] * coi)
            druc_on  <- (mru_on_tx * DTRACEC[,"BSC"]) + (DENTRANTSC[,"BSC"] * coi)
            
            drug_cost <- rep(0,th)
            admin_cost <- rep(0,th)
            ddrug_cost <- rep(0,th)
            dadmin_cost <- rep(0,th)
            ruc_off  <- rep(0,th)
            druc_off <- ruc_off
            
          }
          
          # Ok now no matter how many lines and their order we've done drug cost
          # with full memory for our treatment sequence -_-
          
          return(list(
            molecule = mol,
            undisc = list(
              drug = drug_cost,
              admin = admin_cost,
              mru_on = ruc_on,
              mru_off = ruc_off
            ),
            disc = list(
              drug = ddrug_cost,
              admin = dadmin_cost,
              mru_on = druc_on,
              mru_off = druc_off
            )
          ))
        })
        names(pf_costs) <- names(trt_seq)
        
        # We now no longer need TRACE (or discounted trace) and can drop it to free up (a lot of) memory
        rm(TRACE, DTRACEC, DTRACEQ)
        
        # Finally, add in end of life costs by multiplying the entrants to death from consolidated
        # trace by the cost upon death
        pf_costs$eol <- list(
          undisc = CT$entrants[,ncol(CT$entrants)] * cod,
          disc   = CT_disc$entrants_c[,ncol(CT_disc$entrants_c)] * cod
        )
        
        # QALYs and AEs are calculated using the consolidated traces: pf_qalys
        
        # Get the HSUVs for the active treatment line by treatment status:
        hsuv_active <- unlist(lapply(1:(number_of_lines-1), function(line_n) {
          
          # ASSUMEs POPULATION 0 IF LATER LINES!!!!
          if (line_n == 1) {
            line_hsuv <- as.list(util$hsuv[Population == rpop_n & Treatment.line == line_n & Molecule == trts_n[line_n],list(OnTxt,OffTxt)])
          } else {
            line_hsuv <- as.list(util$hsuv[Population == 0 & Treatment.line == line_n & Molecule == trts_n[line_n],list(OnTxt,OffTxt)])
          }
          
          lab_line <- paste0("L",line_n)
          names(line_hsuv) <- paste0(lab_line,c("_on","_off"))
          unlist(line_hsuv)
        }))
        # Get the HSUV for BSC:
        hsuv_bsc <- structure(
          util$hsuv[Population == 0 & Treatment.line == number_of_lines & Molecule == trts_n[number_of_lines],]$OnTxt,
          .Names = paste0("L",number_of_lines,"_on")
        )
        
        # Stick them together to form the vector of HSUVs per possible state for
        # this treatment sequence:
        hsuv_unadj <- c(hsuv_active,hsuv_bsc,"dead"=0)
        
        # Expand this to the time horizon by multiplying column-wise by the gpop line.
        # The result is U, which we can element-wise multiply by our consolidated trace
        # to get our undiscounted QALYs :)
        U <- util_gpop_mat[,1:length(hsuv_unadj)] %*% diag(hsuv_unadj)
        
        rm(hsuv_active)
        rm(hsuv_bsc)
        rm(hsuv_unadj)
        
        # Undiscounted QALYs (note AEs are separate as they're a separate thing):
        pf_qalys <- list(
          undisc = (CT$full_lines * .subset2(basic,"cl_y")) * U,
          disc   = (CT_disc$full_lines_q * .subset2(basic,"cl_y")) * U
        )
        
        rm(U)
        
        # AEs undiscounted
        
        # pf_ae
        # Same as above, compute your AE stuff here
        
        if(ae$approach[1] == "one-off") {
          
          ae_impact_atl <- do.call(
            rbind,
            lapply(1:(number_of_lines-1), function(line_n) {
              c(cost = ae$per_cycle[ae$per_cycle$trt == treatment_abbr[line_n] & ae$per_cycle$line == line_n]$cost,
                qaly = ae$per_cycle[ae$per_cycle$trt == treatment_abbr[line_n] & ae$per_cycle$line == line_n]$QALYs,
                dur  = mean(
                  ae$one_off$duration_weeks[
                    ae$one_off$Molecule == trts_n[line_n]
                    & ae$one_off$Treatment.line == min(line_n,2)
                  ]))
            })
          )
          
          ae_impact_atl <- data.table(ae_impact_atl)
          ae_oo_cost <- ae_impact_atl$cost * ae_impact_atl$dur
          ae_oo_qaly <- ae_impact_atl$qaly * ae_impact_atl$dur
          
          ent <- data.table(as.matrix(CT$entrants))
          ent <- rbindlist(list(data.frame(t(structure(rep(0,ncol(ent)),.Names=colnames(ent)))),ent))
          ent <- cbind(L1_on = c(1,rep(0,nrow(ent)-1)),ent)
          
          dent_c <- data.table(as.matrix(CT_disc$entrants_c))
          dent_c <- rbindlist(list(data.frame(t(structure(rep(0,ncol(dent_c)),.Names=colnames(dent_c)))),dent_c))
          dent_c <- cbind(L1_on = c(1,rep(0,nrow(dent_c)-1)),dent_c)
          
          dent_q <- data.table(as.matrix(CT_disc$entrants_q))
          dent_q <- rbindlist(list(data.frame(t(structure(rep(0,ncol(dent_q)),.Names=colnames(dent_q)))),dent_q))
          dent_q <- cbind(L1_on = c(1,rep(0,nrow(dent_q)-1)),dent_q)
          
          nonzero_ind <- 1:number_of_lines + 0:(number_of_lines-1)
          
          nstate <- ncol(CT$full_lines)
          n_cyc  <- nrow(CT$full_lines)
          nam_st <- colnames(CT$full_lines)
          
          # clever trick to expand a vector
          ae_oo_cost <- "[<-"(numeric(nstate), nonzero_ind, c(ae_oo_cost,0))
          ae_oo_qaly <- "[<-"(numeric(nstate), nonzero_ind, c(ae_oo_qaly,0))
          
          # multiply payoff by entrants
          pf_ae <- list(
            undisc = list(
              costs = matrix(
                as.matrix(ent)[1:th,] %*% diag(ae_oo_cost),
                ncol = nstate,
                nrow = n_cyc,
                dimnames = list(NULL,nam_st)
              ),
              qalys = matrix(
                as.matrix(ent)[1:th,] %*% diag(ae_oo_qaly),
                ncol = nstate,
                nrow = n_cyc,
                dimnames = list(NULL,nam_st)
              )
            ),
            disc = list(
              costs = matrix(
                as.matrix(dent_c)[1:th,] %*% diag(ae_oo_cost),
                ncol = nstate,
                nrow = n_cyc,
                dimnames = list(NULL,nam_st)
              ),
              qalys = matrix(
                as.matrix(dent_q)[1:th,] %*% diag(ae_oo_qaly),
                ncol = nstate,
                nrow = n_cyc,
                dimnames = list(NULL,nam_st)
              )
            )
          )
          
          
        } else {
          # compute cost per cycle on treatment and qalys lost per cycle on treatment:
          ae_cpc_ontrt <- unlist(lapply(1:number_of_lines, function(act_line_n) {
            as.list(ae$per_cycle[line == act_line_n & molecule == trts_n[act_line_n],])$cost
          }))
          ae_qlpc_ontrt <- unlist(lapply(1:number_of_lines, function(act_line_n) {
            as.list(ae$per_cycle[line == act_line_n & molecule == trts_n[act_line_n],])$QALYs
          }))
          # Make a correction for if we have 4 atl's:
          if (number_of_lines == max_active+1) {
            ae_cpc_ontrt[5] <- ae$per_cycle[line == max_active & molecule == trts_n[number_of_lines],]$cost
            ae_qlpc_ontrt[5] <- ae$per_cycle[line == max_active & molecule == trts_n[number_of_lines],]$QALYs
          }
          nonzero_ind <- 1:number_of_lines + 0:(number_of_lines-1)
          
          # pad the cpc and qlpc with 0s:
          ae_cpc <- "[<-"(numeric(ncol(CT$full_lines)), nonzero_ind, ae_cpc_ontrt)
          ae_qlpc <- "[<-"(numeric(ncol(CT$full_lines)), nonzero_ind, ae_qlpc_ontrt)
          
          names(ae_cpc) <- colnames(CT$full_lines)
          names(ae_qlpc) <- colnames(CT$full_lines)
          
          # Implement adverse events (yay)
          pf_ae <- list(
            undisc = list(
              costs = matrix(
                CT$full_lines %*% diag(ae_cpc),
                ncol = ncol(CT$full_lines),
                nrow = nrow(CT$full_lines),
                dimnames = list(NULL,colnames(CT$full_lines))
              ),
              qalys = matrix(
                CT$full_lines %*% diag(ae_qlpc),
                ncol = ncol(CT$full_lines),
                nrow = nrow(CT$full_lines),
                dimnames = list(NULL,colnames(CT$full_lines))
              )
            ),
            disc = list(
              costs = matrix(
                (CT$full_lines * dfac_c) %*% diag(ae_cpc),
                ncol = ncol(CT$full_lines),
                nrow = nrow(CT$full_lines),
                dimnames = list(NULL,colnames(CT$full_lines))
              ),
              qalys = matrix(
                CT_disc$full_lines %*% diag(ae_qlpc),
                ncol = ncol(CT$full_lines),
                nrow = nrow(CT$full_lines),
                dimnames = list(NULL,colnames(CT$full_lines))
              )
            )
          )
        }
        
        if(include_plots) {
          f_plot_mk_draw_consol_trace(consol_trace = CT$full_lines,
                                      treatment_names = treatment_names,
                                      tmax = 15)
        }
        
        
        return(list(
          population   = popu,
          trt_nam      = treatment_names,
          trace_consol = CT$full_lines,
          costs        = pf_costs,
          qalys        = pf_qalys,
          aes          = pf_ae
        ))
      })
    }))
  })
}





# Post run analysis -------------------------------------------------------




#' Function to generate probabilistic draws for the subsequent treatment proportions
#' 
#' @param subsTx_table The named range `R_table_sub_txts_prop_n_costs` from excel
#' @param sims default 10000 - the amount of samples to draw.
#' @param lookups `p$basic$lookup` containing the lookup tables to make the columns to match model results
#' 
sub_tx_PSA_samples <- function(subsTx_table, sims = 10000, lookups, PSA = TRUE) {
  
  # strip the table down to what gets used in the CE model to reduce unecessary
  # stuff
  subsTx_table <- data.table(subsTx_table)
  subsTx_table <- subsTx_table[!is.na(Population),list(Population,Line.1,Line.2,Line.3,Line.4,Adj.proportion.given.line.1,n)]
  
  # Make the columns that match with the model results:
  
  lu_mol <- lookups$ipd$mol
  
  # Weighting table:
  # The subsTx table only needs to be computed once, so just get it done:
  subsTx_table$L1 <- lu_mol[match(subsTx_table$Line.1,RCC_input_desc,nomatch = NA),]$Number
  subsTx_table$L2 <- lu_mol[match(subsTx_table$Line.2,RCC_input_desc,nomatch = NA),]$Number
  subsTx_table$L3 <- lu_mol[match(subsTx_table$Line.3,RCC_input_desc,nomatch = NA),]$Number
  subsTx_table$L4 <- lu_mol[match(subsTx_table$Line.4,RCC_input_desc,nomatch = NA),]$Number
  subsTx_table$L5 <- 999
  subs    <-
    subsTx_table[!is.na(Population), list(
      Population,
      L1,
      L2,
      L3,
      L4,
      L5,
      Adj.proportion.given.line.1,
      n
    )]
  subs$trt_n <- do.call(paste, c(subs[,paste0("L",1:5),with=FALSE], sep="→"))
  subs$trt_n <- gsub("→NA","",subs$trt_n)
  subs <- subs[,list(Population,L1,trt_n,Adj.proportion.given.line.1,n)]
  
  # split the table by population and first-line therapy to give us something
  # to lapply through:
  split_table <- split(subs, by=c("Population","L1"))
  
  # Each element in the lapply is then one pairing of population and first-line
  # therapy. We can then generate our dirichlet draws for it as we know 
  # that the column Adj.proportion.given.line.1 sums to 1. Therefore we can
  # generate diriclet draws from it by multiplying it by N
  rbindlist(lapply(split_table, function(popL1) {
    
    # sum(popL1$Adj.proportion.given.line.1) will always be 1. Therefore, draw
    # the n column then back calculate the proportion
    
    if (PSA == TRUE) {
    draws <- gtools::rdirichlet(sims,alpha = popL1$n)
    } else {
    draws <- popL1$Adj.proportion.given.line.1
    }
    
    # QC proof: the below is true:
    # c(t(draws))[1:ncol(draws)] == draws[1,]
    
    id <- as.list(popL1[,list(Population,L1,trt_n)])
    
    # the PSA results use oo_pop as a numeric
    
    id$oo_pop <- as.numeric(gsub("pop","",id$Population))
    id$Population <- NULL
    
    id <- as.data.table(lapply(id, function(col) rep(col,sims)))
    id$iteration <- rep(1:sims, each=nrow(popL1))
    id$Adj.proportion.given.line.1 <- c(t(draws))
    
    # QC test 2: all grouped sums are equal to 1 (or extremely close due to floating point errors)
    # all(round(id[,.(tst = sum(Adj.proportion.given.line.1)),by=list(Population,L1,iteration)]$tst,13) == 1)
    id
    
  }))
}


f_psa_computeWAModelRes <- function(R_table_sub_txts_prop_n_costs, sims, lookups, psa_results, PSA = TRUE) {
  
  
  weighting_table <- sub_tx_PSA_samples(
    subsTx_table = data.table(R_table_sub_txts_prop_n_costs),
    sims = sims,
    lookups = lookups, 
    PSA = PSA
  )
  
  lu_pop <- lookups$pop_map
  
  weighting_table$oo_pop <- lu_pop[match(weighting_table$oo_pop, Sequencing.population.number,nomatch = NA),]$Overall.population.number
  
  add_pop_2 <- weighting_table[oo_pop==1]
  add_pop_2$oo_pop <- 2
  
  add_pop_5 <- weighting_table[oo_pop==4]
  add_pop_5$oo_pop <- 5
  
  weighting_table <- rbind(weighting_table, add_pop_2, add_pop_5)
  
  # merge in the weightings - note that some rows disappear because there's not a row
  # in the weighting table for it!
  # 
  # To show these rows up, use this:
  
  
    res_tab <-
      merge.data.table(
        weighting_table,
        psa_results,
        by = c("oo_pop", "L1", "trt_n", "iteration"),
        all.y = TRUE
      )
    
  
  res_tab$Adj.proportion.given.line.1[is.na(res_tab$Adj.proportion.given.line.1)] <- 0
  
  # QC test: tst should be 1 in the below
  # res_tab[,.(tst = sum(Adj.proportion.given.line.1)),by=list(oo_pop,L1,iteration)]
  
  # EW edit: testing taking into account cPAS vs List price
  # res_tab[,.(tst = sum(Adj.proportion.given.line.1)),by=list(oo_pop,L1,iteration,dd_drug_price_options)]
  
  # Apply the weightings for all iterations all in one go to all the outputs!
  res_w <- res_tab[,`:=`(
    qaly        = qaly * Adj.proportion.given.line.1,
    mol_0       = mol_0 * Adj.proportion.given.line.1,
    mol_1       = mol_1 * Adj.proportion.given.line.1,
    mol_2       = mol_2 * Adj.proportion.given.line.1,
    mol_3       = mol_3 * Adj.proportion.given.line.1,
    mol_4       = mol_4 * Adj.proportion.given.line.1,
    mol_5       = mol_5 * Adj.proportion.given.line.1,
    mol_6       = mol_6 * Adj.proportion.given.line.1,
    mol_7       = mol_7 * Adj.proportion.given.line.1,
    mol_8       = mol_8 * Adj.proportion.given.line.1,
    mol_9       = mol_9 * Adj.proportion.given.line.1,
    mol_10      = mol_10 * Adj.proportion.given.line.1,
    mol_11      = mol_11 * Adj.proportion.given.line.1,
    mol_12      = mol_12 * Adj.proportion.given.line.1,
    mol_999     = mol_999 * Adj.proportion.given.line.1,
    other_costs = other_costs * Adj.proportion.given.line.1,
    LY          = LY * Adj.proportion.given.line.1
  )]
  
  # now sum them up:
    res_w_sum <- res_w[,.(
      qaly        = sum(qaly),
      mol_0       = sum(mol_0),
      mol_1       = sum(mol_1),
      mol_2       = sum(mol_2),
      mol_3       = sum(mol_3),
      mol_4       = sum(mol_4),
      mol_5       = sum(mol_5),
      mol_6       = sum(mol_6),
      mol_7       = sum(mol_7),
      mol_8       = sum(mol_8),
      mol_9       = sum(mol_9),
      mol_10      = sum(mol_10),
      mol_11      = sum(mol_11),
      mol_12      = sum(mol_12),
      mol_999     = sum(mol_999),
      other_costs = sum(other_costs),
      LY          = sum(LY),
      Adj.proportion.given.line.1 = sum(Adj.proportion.given.line.1)
    ), by = list(oo_pop,L1,iteration, dd_drug_price_options)]
    
    #NOTE from EW to DL/DB 24-8-23:  I've added in dd_drug_price_options here as the
    #output from the HPC for the PSA includes this.  Let me know if this creates
    #problems elsewhere (recommend dd_drug_price_options added in as a column 
    #in res_tab if so)
    
  avg_res_w <- res_w_sum[,.(
    mean_qaly        = mean(qaly),
    lb_qaly = quantile(qaly,prob=0.025),
    ub_qaly = quantile(qaly,prob=0.975),
    mean_mol_0       = mean(mol_0),
    lb_mol_0 = quantile(mol_0,prob=0.025),
    ub_mol_0 = quantile(mol_0,prob=0.975),
    mean_mol_1       = mean(mol_1),
    lb_mol_1 = quantile(mol_1,prob=0.025),
    ub_mol_1 = quantile(mol_1,prob=0.975),
    mean_mol_2       = mean(mol_2),
    lb_mol_2 = quantile(mol_2,prob=0.025),
    ub_mol_2 = quantile(mol_2,prob=0.975),
    mean_mol_3       = mean(mol_3),
    lb_mol_3 = quantile(mol_3,prob=0.025),
    ub_mol_3 = quantile(mol_3,prob=0.975),
    mean_mol_4       = mean(mol_4),
    lb_mol_4 = quantile(mol_4,prob=0.025),
    ub_mol_4 = quantile(mol_4,prob=0.975),
    mean_mol_5       = mean(mol_5),
    lb_mol_5 = quantile(mol_5,prob=0.025),
    ub_mol_5 = quantile(mol_5,prob=0.975),
    mean_mol_6       = mean(mol_6),
    lb_mol_6 = quantile(mol_6,prob=0.025),
    ub_mol_6 = quantile(mol_6,prob=0.975),
    mean_mol_7       = mean(mol_7),
    lb_mol_7 = quantile(mol_7,prob=0.025),
    ub_mol_7 = quantile(mol_7,prob=0.975),
    mean_mol_8       = mean(mol_8),
    lb_mol_8 = quantile(mol_8,prob=0.025),
    ub_mol_8 = quantile(mol_8,prob=0.975),
    mean_mol_9       = mean(mol_9),
    lb_mol_9 = quantile(mol_9,prob=0.025),
    ub_mol_9 = quantile(mol_9,prob=0.975),
    mean_mol_10      = mean(mol_10),
    lb_mol_10 = quantile(mol_10,prob=0.025),
    ub_mol_10 = quantile(mol_10,prob=0.975),
    mean_mol_11      = mean(mol_11),
    lb_mol_11 = quantile(mol_11,prob=0.025),
    ub_mol_11 = quantile(mol_11,prob=0.975),
    mean_mol_12      = mean(mol_12),
    lb_mol_12 = quantile(mol_12,prob=0.025),
    ub_mol_12 = quantile(mol_12,prob=0.975),
    mean_mol_999     = mean(mol_999),
    lb_mol_999 = quantile(mol_999,prob=0.025),
    ub_mol_999 = quantile(mol_999,prob=0.975),
    mean_other_costs = mean(other_costs),
    lb_other_costs = quantile(other_costs,prob=0.025),
    ub_other_costs = quantile(other_costs,prob=0.975),
    mean_LY          = mean(LY),
    lb_LY = quantile(LY,prob=0.025),
    ub_LY = quantile(LY,prob=0.975),
    Adj.proportion.given.line.1 = mean(Adj.proportion.given.line.1)
  ), by = list(oo_pop,L1, dd_drug_price_options)]
  
  
  avg_res_w <- avg_res_w[Adj.proportion.given.line.1 != 0,]
  
  return(list(
    weightings = weighting_table,
    weighted = res_w,
    weighted_average = res_w_sum,
    mean_weighted_average = avg_res_w
  ))
  
}


