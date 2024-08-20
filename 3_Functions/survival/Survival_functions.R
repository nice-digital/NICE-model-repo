

# Cleaning functions ------------------------------------------------------

# This function gets the KM data into the proper format to use with the survminer package
# This packages requires that column names are the same (you have to specify in advance)
# t_multiplier can be used to change the time unit of the analysis.
# time unit has been specified in weeks

f_ce_km_MakeDatSurvFriendly <- function(Data_required, time_column, event_column, t_multiplier = 1) {
  dat           <- Data_required[,c(time_column,event_column),with =FALSE]
  colnames(dat) <- c("t", "ec")
  dat[,"t"] <- dat[,"t"] * t_multiplier
  return(dat)
}



# Extrapolation functions -------------------------------------------------


# ~ Routing function ------------------------------------------------------


#' router function for extrapolations. Extrapolation functions called are defined
#' after this function so you can read things in order.
#' 
#' @param v_cycles vector of cycles
#' @param coefficients coefficients for this distribution
#' @param distribution distribution identifier
#' 
f_extrapolate <- function(v_cycles, coefs, distribution) {
  
  if(is.null(coefs)) return(NULL)
  
  if (!distribution %in% c("gengamma","exp","weibull","lnorm","gamma","gompertz","llogis","lognormal")) stop("distribution not listed")
  
  # Check that the coefs match the distribution stated
  
  if(distribution == "gengamma") {if (sum(names(coefs) != c("mu","sigma","Q")) > 0) {stop("incorrect coefficients")} 
  } else if(distribution == "exp") {if (length(coefs) > 1) {stop("incorrect coefficients")}
  } else if(distribution == "weibull") {if (sum(names(coefs) != c("shape", "scale")) > 0) {stop("incorrect coefficients")}
  } else if(distribution == "lnorm") {if (sum(names(coefs) != c("meanlog", "sdlog")) > 0) {stop("incorrect coefficients")}
  } else if(distribution == "gamma") {if (sum(names(coefs) != c("shape", "rate")) > 0) {stop("incorrect coefficients")}
  } else if(distribution == "gompertz") {if (sum(names(coefs) != c("shape", "rate")) > 0) {stop("incorrect coefficients")}
  } else if(distribution == "llogis") {if (sum(names(coefs) != c("shape", "scale")) > 0) {stop("incorrect coefficients")}
  }
  
  # select the correct distribution to apply and produce table vector of S(t) per cycle
  
  if (distribution == "exp") {return(function_apply_exp(v_cycles, coefs))}
  if (distribution == "weibull") {return(function_apply_weibull(v_cycles, coefs))}
  if (distribution == "gengamma") {return(function_apply_gengamma(v_cycles, coefs))}
  if (distribution == "lnorm") {return(function_apply_lnorm(v_cycles, coefs))}
  if (distribution == "gamma") {return(function_apply_gamma(v_cycles, coefs))}
  if (distribution == "gompertz") {return(function_apply_gompertz(v_cycles, coefs))}
  if (distribution == "llogis") {return(function_apply_llogis(v_cycles, coefs))}
  
  
}


# ~ extrapolation functions -----------------------------------------------


# these functions each produce a vector of S(t) per cycle for the selected distriubtion

function_apply_exp <- function(v_cycles, rate) {
  output<- exp(-1* c(exp(rate)) * v_cycles)
  output
}

function_apply_weibull <- function(v_cycles, coefs) {
  output<- 1 - pweibull(v_cycles, exp(coefs[1]), exp(coefs[2]))
  output
}

function_apply_lnorm <- function(v_cycles, coefs) {
  output<- 1 - plnorm(v_cycles, coefs[1], 1/exp(-1*coefs[2]))
  output
}

function_apply_llogis <- function(v_cycles, coefs) {
  output<- 1/(1+(v_cycles*exp(-1*coefs[2]))^(1/exp(-1*coefs[1])))
  output
}

function_apply_gengamma <- function(v_cycles, coefs) {
  pgamma <- pgamma(((-1*coefs[3])^-2)*exp(-1*coefs[3]*-((log(v_cycles)-(coefs[1]))/exp(coefs[2]))), ((-1*coefs[3])^-2), scale = 1)
  output<- if (coefs[3]<0) {pgamma} else {1-pgamma}
  output
}

function_apply_gompertz <- function(v_cycles, coefs) {
  output<- exp((-1/coefs[1])*exp(coefs[2])*(exp(coefs[1]*v_cycles)-1))
  output
}

function_apply_gamma <- function(v_cycles, coefs) {
  output<- 1-((exp(gammaln(exp(coefs[1])))*pgamma(exp(coefs[2])*v_cycles,exp(coefs[1]),1))/gamma(exp(coefs[1])))
  output
}


# plotting functions ------------------------------------------------------



# function to produce survival analysis output plots

f_surv_plot <- function(v_years, Plot_data, ipd, xlim, km_time_mult) {
  
  km_line_data <- f_ce_km_gen_KMPlotData(
    DATA     = ipd,
    timevar  = "timew",
    eventvar = "event_censor"
  )
  km_line_data$dist <- "KM"
  km_line_data$t <- km_line_data$t * km_time_mult
  
  dat <- data.table(Plot_data)
  dists <- dimnames(dat)[[2]]
  dat$t <- v_years
  dat <- melt.data.table(
    data = dat,
    measure.vars = dists,
    variable.name = "dist",
    value.name = "s_t"
  )
  
  dat <- dat[t <= xlim,]
  dat <- rbind(km_line_data,dat)
  
  # colours
  colour <- c(
    KM       = "black",
    gengamma = "blue",
    exp      = "red",
    weibull  = "green",
    lnorm    = "purple",
    gamma    = "yellow",
    gompertz = "orange",
    llogis   = "brown"
  )
  size <- c(
    KM       = 1.5,
    gengamma = 0.5,
    exp      = 0.5,
    weibull  = 0.5,
    lnorm    = 0.5,
    gamma    = 0.5,
    gompertz = 0.5,
    llogis   = 0.5
  )
  
  
  ggplot(data = dat, aes(x = t, y = s_t, colour = dist)) + 
    geom_line() +    
    ylim(0,1) +
    xlim(0,xlim) +
    labs(x ="Years", y = "Proportion Surviving") +
    theme_bw() + theme(legend.position = "bottom") + 
    scale_colour_manual(name="Curve Fit",values = colour) +
    scale_size_manual(name="Curve Fit",values = size) +
    scale_x_continuous(expand = expansion(mult = c(0,0))) +
    scale_y_continuous(expand = expansion(mult = c(0,0))) +
    guides(
      color = guide_legend(title = NULL),
      linetype = "none",
      size = "none"
    )
}



# Takes the survival models and draws Kaplan Meier with a risk table under it.
#   Then also takes the survival extrapolations and overlays them

f_extrap_plot <- function(SurvEstimate, Data_required, curvefits_data, time_vector, xlim, break_by) {
  
  # bit of reformatting
  dat   <- data.table(curvefits_data)
  dists <- dimnames(dat)[[2]]
  dat$t <- time_vector
  dat   <- melt.data.table(
    data = dat,
    measure.vars = dists,
    variable.name = "dist",
    value.name = "s_t"
  )
  
  cols  <- c("KM"       = "#030303",
             "gengamma" = "#FF4040",
             "exp"      = "#0000FF",
             "weibull"  = "#D2691E",
             "lnorm"    = "#EE6A50",
             "gamma"    = "#556B2F",
             "gompertz" = "#FF7F00",
             "llogis"   = "#00CDCD")
  thick  <- c("KM"       = 1.5,
              "gengamma" = 0.5,
              "exp"      = 0.5,
              "weibull"  = 0.5,
              "lnorm"    = 0.5,
              "gamma"    = 0.5,
              "gompertz" = 0.5,
              "llogis"   = 0.5)
  legend  <- c("KM",
             "Generalised Gamma",
             "Exponential",
             "Weibull",
             "Lognormal",
             "Gamma",
             "Gompertz",
             "Loglogistic")
  
  legend  <- c("Exponential",       # sorted by alphabetical order
               "Gamma",
               "Generalised Gamma",
               "Gompertz",
               "KM",
               "Lognormal",
               "Loglogistic",
               "Weibull")
  
  
  # plot <- ggsurvplot(
  plot <- ggsurvplot(
    fit                   = SurvEstimate,
    data                  = Data_required,
    combine               = TRUE,
    censor                = TRUE,
    risk.table            = TRUE,
    conf.int              = TRUE,
    break.x.by            = break_by,
    break.y.by            = 0.1,
    xlim                  = c(0, xlim),
    xlab                  = "Years",
    size                  = 0.72,
    legend.title          = '',
    legend.labs           = c("KM"),
    risk.table.y.text.col = FALSE
  )
  
  plot$plot <- plot$plot +
    geom_line(aes(
      x      = t,
      y      = s_t,
      colour = dist
    ),
    data = dat,
    alpha = 0.8,
    linewidth = 0.3) +
    scale_color_manual(values = cols, labels = legend) +    
    scale_size_manual(values = thick) +
    scale_y_continuous(expand = expansion(c(0,0.05))) +
    scale_x_continuous(expand = expansion(c(0.05,0.01))) +
    guides(
      color = guide_legend(title = NULL),
      linetype = "none",
      size = "none"
    ) 
    suppressWarnings(suppressMessages(plot))
}




# Makes KM plotting data to put onto a plot without having to us survminer.
# 
# Could use this instead of function_extrap_plot
# 
f_ce_km_gen_KMPlotData <- function(DATA, timevar, eventvar) {
  #make KM  testing ground
  #load example data and do survfit
  x <- as.data.table(DATA)
  y <- Surv(time = unlist(x[,timevar, with =FALSE],use.names = F), unlist(x[,eventvar,with =FALSE],use.names = F))
  z <- survfit(y ~ 1)
  
  #make adjustments to data required for graph (i.e. double both datasets)
  t   <- rep(z$time,2)
  s_t <- rep(z$surv,2)
  
  #sort parameters as required (by time ascending, by survival descending), adding extra data points
  t   <- append(0,t[order(t)])
  s_t <- append(1,append(1,s_t[order(s_t,decreasing = TRUE)]))[1:length(t)]
  
  #put into a dataframe and return as output of function
  df <- data.table(t = t, s_t = s_t)
  return(df)
}



# Run all models for all data  ----------------------------------------------------------------

#' Function to run TSD14 survival analysis for all available data by going through
#' all possible populations, lines, regimen, trials and endpoints, returning NULL
#' if the data is empty and a set of outputs if there is data.
#' 
#' Contains options for producing plots (considerable impact on speed, but useful for reporting)
#' 
#' @param r_pld All of the survival analyses in the PLMTE format. `i$surv$pld` in the model structure
#' @param id Identifiers for ipd - `i$id$ipd` in the model structure
#' @param lookups lookup tables to translate between numbers and labels: `i$lookup$ipd` in the model structure
#' @param distnames `i$distnames` in the model structure
#' @param cl_y cycle length in years `p$basic$cl_y` in the model structure
#' @param xlim_survplots_yr x-limit for the plots: `p$misc$plot$xlim_survplots_yr` in the model structure
#' @param t_yr vector of time in years: `p$basic$t_yr` in the model structure
#' @param draw_plots logical for drawing plots
#' @param verbose logical for extra console output that can help with debugging
#' @param min_obs cutoff for running regressions. `N<min_obs` won't run.
#' 
#' 
f_surv_runAllTSD14 <-
  function(r_pld,
           id,
           lookups,
           distnames,
           t_cyc,
           cl_y              = NULL,
           xlim_survplots_yr = NULL,
           t_yr              = NULL,
           draw_plots        = FALSE,
           verbose           = FALSE,
           min_obs           = 30) {
    
    if (draw_plots) {
      if(is.null(cl_y)) stop("If draw_plots is TRUE, then cl_y (cycle length in years) must be provided")
      if(is.null(xlim_survplots_yr)) stop("If draw_plots is TRUE, then xlim_survplots_yr (x axis limit for plot in years) must be provided")
      if(is.null(t_yr)) stop("If draw_plots is TRUE, then t_yr (vector of time in years to TH) must be provided")
    }
    
    
  lapply(id$pop, function(trial_population) {
    lapply(id$line, function(l) {
      lapply(id$mol, function(m) {
        lapply(id$trial, function(tr) {
          lapply(id$endpoint, function(en) {
            
            
            # Filter down to the parameters avove associated with this combination:
            ipd <- r_pld[population == trial_population & line==l & molecule==m & trial==tr & endpoint==en,list(timew,event_censor)]
            
            names(ipd) <- c("t","e")
            
            if(all(verbose, nrow(ipd) > 0)) {
              cat(paste0(
                "Survival analysis - population: ", lookups$pop[Number      == trial_population, Description],
                "\t line: "                       , lookups$line[Number     == l               , Description],
                "\t molecule: "                   , lookups$mol[Number      == m               , Description],
                "\t trial: "                      , tr,
                "\t endpoint: "                   , lookups$endpoint[Number == en              , Description], "\n"
              ))
            }
            
            
            # If that dataset is empty (i.e., has no rows), then there's no data for it. return
            # nothing
            if (nrow(ipd) == 0) {
              return(list(
                pop      = lookups$pop[     Number == trial_population,Description],
                line     = lookups$line[    Number == l,Description],
                mol      = lookups$mol[     Number == m,Description],
                tr       = lookups$trial[     Number == tr,Description],
                endpoint = lookups$endpoint[Number == en,Description],
                ipd      = NULL,
                fs_fits  = NULL,
                gof      = NULL,
                st       = NULL,
                plot     = NULL
              ))
            } else if (nrow(ipd) < min_obs) {
              
              warning(paste0(
                lookups$pop[Number      == trial_population, Description],
                " population "     , lookups$line[Number     == l               , Description],
                " " , lookups$mol[Number      == m               , Description],
                " " , lookups$endpoint[Number == en              , Description],
                " from "    , lookups$trial[Number == tr              , Description],
                " has ",nrow(ipd), "(<", min_obs, ") observations (i.e. < the minimum set by you!). Skipping this PLMTE."
              ),immediate. = TRUE)
              
              return(list(
                pop      = lookups$pop[     Number == trial_population,Description],
                line     = lookups$line[    Number == l,Description],
                mol      = lookups$mol[     Number == m,Description],
                tr       = lookups$trial[     Number == tr,Description],
                endpoint = lookups$endpoint[Number == en,Description],
                ipd      = NULL,
                fs_fits  = NULL,
                gof      = NULL,
                st       = NULL,
                plot     = NULL
              ))
              
              
            } else {
              
              # If there IS data, continue on to use that data to produce full survival analysis results.
              # 
              # Remember that the IPD is sensitive so should not ever be saved to disk. Being stored in i
              # and never saved is the way to do this.
              
              # Even within the multiple nesting by population, line, molecule, trial and endpoint,
              # we need ANOTHER (6th) layer - the parametric model:
              
              
              
              
              
              
              fs_fits <- lapply(distnames, function(dist) {  # applying all parametric survival curves in the list of distNames
                fs_fit <- flexsurvreg(
                  formula = Surv(t, e) ~ 1,
                  data = ipd,
                  dist = dist
                )
                return(list(
                  coefs = coefficients(fs_fit),                                         # coefficients for the fitted model
                  vcov  = vcov(fs_fit),                                                 # variance covariance matrix for the fitted model
                  fit   = c(AIC= AIC(fs_fit), BIC=BIC(fs_fit), logLik = logLik(fs_fit)) # goodness of fit statistics for the fitted model
                ))
              })
              
              gof <- do.call(rbind, lapply(distnames, function(dist) fs_fits[[dist]]$fit))
              
             
              st <- matrix(
                unlist(lapply(distnames, function(dist) {
                  f_extrapolate(t_cyc, fs_fits[[dist]]$coefs, dist)
                })),
                ncol = length(distnames),
                dimnames = list(NULL, distnames),
                byrow = FALSE
              )
              
              if (draw_plots) {
                
                plot <- {
                  # First the IPD is produced in a format that survminer will accept. Data must all be 
                  # the same format with the same column names. 
                  # this assumes no covariate adjustment.
                  
                  sm_ipd <- f_ce_km_MakeDatSurvFriendly(
                    Data_required = ipd,
                    time_column   = "t",                 # note that this is taking IPD in weeks 
                    event_column  = "e",
                    t_multiplier  = cl_y             # data in weeks, cycle length in plot years
                  )
                  
                  # get the survival analysis in the form we need for survminer
                  # and make the extrapolations we need for survminer
                  
                  form          <- Surv(t, ec) ~ 1
                  sm_surv_est   <- surv_fit(formula = form, data = sm_ipd)
                  
                  # make the plot with the input data:
                  survival_plot <- suppressMessages(f_extrap_plot(
                    SurvEstimate   = sm_surv_est,
                    Data_required  = sm_ipd,
                    curvefits_data = st,
                    time_vector    = t_yr,
                    xlim           = xlim_survplots_yr,   
                    break_by       = round(xlim_survplots_yr/8,0) 
                  ))
                  list(
                    ipd     = sm_ipd,
                    formula = form,
                    plot    = survival_plot
                  )
                }
              } else {
                plot <- NULL
              }
              
              # Now that we've done everything for this dataset, return a list of the stuff
              # we need for it:
              return(list(
                pop      = lookups$pop[     Number == trial_population,Description],
                line     = lookups$line[    Number == l,Description],
                mol      = lookups$mol[     Number == m,Description],
                tr       = lookups$trial[     Number == tr,Description],
                endpoint = lookups$endpoint[Number == en,Description],
                ipd      = ipd,
                fs_fits  = fs_fits,
                gof      = gof,
                st       = st,
                plot     = plot
              ))
            }
          })
        })
      })
    })
  })
}


# Reporting ---------------------------------------------------------------

#' Function to automatically bring together graphs and goodness of fit for all TSD14 models
#' 
#' @param fs_res the results from f_surv_runAllTSD14
#' @param id the identifiers for all the numerical indices to separate out populations, lines etc
#' 
f_surv_makeTSD14Report <- function(fs_res, id, lookup) {
  
  doc_surv <- read_docx() %>%
    body_add_toc(level = 4)
  
  for(subgroup in names(id$pop)) {
    doc_surv <- doc_surv %>%
      body_add_par(paste0(lookup$pop[Number == id$pop[subgroup],Description], " population"),style = "heading 1")
    for (l in names(id$line)) {
      doc_surv <- doc_surv %>%
        body_add_break() %>%
        body_add_par(lookup$line[Number == id$line[l],Description],style = "heading 2")
      for (m in names(id$mol)) {
        for (tr in names(id$trial)) {
          for (en in names(id$endpoint)) {
            # cat(paste0(paste(subgroup, l, m, tr, en, collapse = "\t"),"\n"))
            if (!is.null(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$ipd)) {
              
              row_to_bold <- which.min(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$gof[,"AIC"])
              
              gof_tab <- round(as.data.frame(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$gof), 2) %>%
                rownames_to_column(var = "Distribution") %>%
                mutate(Distribution = case_when(
                  Distribution == "gengamma" ~ "Generalised gamma",
                  Distribution == "exp"      ~ "Exponential",
                  Distribution == "weibull"  ~ "Weibull",
                  Distribution == "lnorm"    ~ "Log-normal",
                  Distribution == "gamma"    ~ "Gamma",
                  Distribution == "gompertz" ~ "Gompertz",
                  Distribution == "llogis"   ~ "Log-logistic"
                )) %>%
                flextable() %>%
                bold(i = row_to_bold, bold = TRUE, part = "body") %>%
                flextable::set_header_labels(
                  AIC = "A.I.C.",
                  BIC = "B.I.C.",
                  logLik = "Log likelihood"
                ) %>%
                autofit()
                
              
              doc_surv <- doc_surv %>%
                body_add_par("",style = "Normal") %>%
                body_add_par("",style = "Normal") %>%
                body_add_par(paste0(
                  fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$mol,
                  ", ",
                  fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$endpoint
                ), style = "heading 3") %>%
                body_add_par("",style = "Normal") %>%
                body_add_normal(
                  paste0(
                    "The Figure below shows ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$endpoint,
                    " in the ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$pop),
                    " population - ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$line),
                    " ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$mol,
                    " (trial ",
                    tr,
                    "). A plot is provided for visual fit, alongside goodness-of-fit statistics for decision making"
                  )
                ) %>%
                body_add_par("",style = "Normal") %>%
                body_add_figure_legend(
                  legend = paste0(
                    "Extrapolation - ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$endpoint,
                    " in the ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$pop),
                    " population - ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$line),
                    " ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$mol,
                    " (trial ",
                    tr,
                    ")"
                  ),
                  bookmark = gsub("_", "", paste0("fig_", subgroup, l, m, tr, en)),
                  
                ) %>%
                body_add_plot(print(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$plot$plot), width = 6) %>%
                body_add_par("", style = "Normal") %>%
                body_add_table_legend(
                  paste0(
                    "Goodness of fit - ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$endpoint,
                    " in the ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$pop),
                    " population - ",
                    tolower(fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$line),
                    " ",
                    fs_res[[subgroup]][[l]][[m]][[tr]][[en]]$mol,
                    " (trial ",
                    tr,
                    ")"
                  ),
                  bookmark = gsub("_", "", paste0("table_", subgroup, l, m, tr, en))
                ) %>%
                body_add_flextable(gof_tab,
                  topcaption = TRUE,
                  split = TRUE
                ) %>%
                body_add_par("", style = "Normal") %>%
                body_add_par("", style = "Normal")
            }
          }
        }
      }
    }
  }
  return(doc_surv)
}


# Comparative efficacy ----------------------------------------------------

#' Function to generate an empty list representing an evidence network with space for 
#' application of hazard ratios to populate modeled comparisons between different
#' drugs across all populations, lines, molecules, trials and endpoints
#' 
#' 
#' @param id a list containing identifiers for the different PLMTE
#' @param lookup a translator between numerical ids and descriptions of them 
#'
#'
f_NMA_generateNetwork <- function(id,lookup) {
  lapply(id$pop, function(trial_population) {
    lapply(id$line, function(l) {
      lapply(id$mol, function(m) {
        lapply(id$trial, function(tr) {
          lapply(id$endpoint, function(en) {
            list(
              dest = list(
                pop      = paste0("pop_"      ,trial_population),
                line     = paste0("line_"     ,l),
                mol      = paste0("mol_"      ,m),
                trial    = paste0("trial_"    ,tr),
                endpoint = paste0("endpoint_" ,en)
              ),
              orig = list(
                pop      = NULL,
                line     = NULL,
                mol      = NULL,
                trial    = NULL,
                endpoint = NULL,
                dist     = NULL,
                source   = NULL
              ),
              hr = 1,
              fp = list()
            )
          })
        })
      })
    })
  })
}


#' Function to pull out one sample from CODA in the RCC format. Notice that the
#' naming of the columns MUST match exactly. Spelling errors will break the function.
#' 
#' @param dat data (CODA sample), run must have column called "Run"
#' @param n which iteration to pull for each grouping
#' 
f_PHNMA_getSample <- function(dat, n) dat[Run == n, .(HR = mean(HR)), by = list(Population,
                                                                                Line,
                                                                                Molecule,
                                                                                Endpoint,
                                                                                Reference.treatment,
                                                                                Reference.trial)]

# e.g. Get the 2nd sample HRs
# hr_table <- f_PHNMA_getSample(i$PHNMA$data,2)




#' Function to "link" the result of f_NMA_generateNetwork by entering hazard ratios
#' and origin identifiers to the foundations it lays
#' 
#' @param network the result of f_NMA_generateNetwork
#' @param hr_table a table of HRs in line with the formatting of the CODA sample for the RCC project. Numbering MUST align with id for f_NMA_generateNetwork
#' @param rel_eff_table the named range R_table_eff_data_settings from the excel front end
#' 
f_NMA_linkPHNMA <- function(network, hr_table) {
  
  # Simply put, this function takes each row of the NMA HR table and puts that information
  # along with origin information into network.
  
  Reduce(
    x = 1:nrow(hr_table),
    init = network,
    accumulate = FALSE,
    f = function(prev, hr_tab_row) {
      
      dat <- as.list(hr_table[hr_tab_row,])
      
      # Establish the destination, d as the labeling that is used to identify
      # list elements. trial_0 is asserted because application of comp eff isn't really
      # associated with any trial, so internal numbering within group will always
      # be 0
      
      d <- list(
        pop      = paste0("pop_"     , .subset2(dat,"Population")),
        line     = paste0("line_"    , .subset2(dat,"Line")),
        mol      = paste0("mol_"     , .subset2(dat,"Molecule")),
        trial    = paste0("trial_"   , .subset2(dat,"Reference.trial")), # trial_0 assumed for HR applied - could be a flag like NMA but that creates complication later
        endpoint = paste0("endpoint_", .subset2(dat,"Endpoint"))
      )
      
      # Extract the PLMTE for this destination, so that we have a list to manipulate
      # directly without fear of overwriting anything :)
      
      plmte <- f_misc_get_plmte(prev,d)
      
      # Derive the origin from dat, these values can then go in plmte$orig
      
      o <- list(
        pop      = paste0("pop_"     , .subset2(dat,"Population")),
        line     = paste0("line_"    , .subset2(dat,"Line")),
        mol      = paste0("mol_"     , .subset2(dat,"Reference.treatment")),
        trial    = paste0("trial_"   , .subset2(dat,"Reference.trial")),
        endpoint = paste0("endpoint_", .subset2(dat,"Endpoint"))
      )
      
      # Place the relevant information into the destination PLMTE. Whether or not
      # this information is used in the end is decided in the function which follows
      # this one, called f_NMA_AddAssumptionsToNetwork
      
      if (!is.null(o$pop))      plmte$orig$pop      <- o$pop
      if (!is.null(o$line))     plmte$orig$line     <- o$line
      if (!is.null(o$mol))      plmte$orig$mol      <- o$mol
      if (!is.null(o$trial))    plmte$orig$trial    <- o$trial
      if (!is.null(o$endpoint)) plmte$orig$endpoint <- o$endpoint
      if (!is.null(dat$HR))     plmte$hr            <- dat$HR
      
      # Take the updated plmte, being careful not to overwrite anything and 
      # slot it back into its destination, returning the updated network.
      
      prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
      return(prev)
    }
  )
}


#' Function to "link" the result of f_NMA_generateNetwork by entering fractional
#' polynomial network meta analysis (FPNMA) parameters table (which contains numerical
#' identifiers for population line molecule trial and endpoint, PLMTE)
#' 
#' @param network the result of f_NMA_generateNetwork
#' @param hr_table a table of HRs parameters in line with the formatting of the CODA sample for the RCC project. Numbering ids MUST align with id for f_NMA_generateNetwork!
#' @param rel_eff_table the named range R_table_eff_data_settings from the excel front end for any other bits/assumptions that're required
#' 
f_NMA_linkFPNMA <- function(network, destinations, hr_table, time_horizon) {
  
  # This function takes each item of hr_list and puts that information
  # along with origin information into network.
  
  Reduce(
    x = 1:nrow(destinations),
    init = network,
    accumulate = FALSE,
    f = function(prev, dest_row_number) {
      
      dat <- as.list(destinations[dest_row_number,])
      
      # Establish the destination, d as the labeling that is used to identify
      # list elements. trial_0 is asserted because application of comp eff isn't really
      # associated with any trial, so internal numbering within group will always
      # be 0
      
      d <- list(
        pop      = paste0("pop_"     , .subset2(dat,"Population")),
        line     = paste0("line_"    , .subset2(dat,"Line")),
        mol      = paste0("mol_"     , .subset2(dat,"Molecule")),
        trial    = paste0("trial_"   , .subset2(dat,"Reference.trial")), # trial_0 assumed for HR applied - could be a flag like NMA but that creates complication later
        endpoint = paste0("endpoint_", .subset2(dat,"Endpoint"))
      )
      dN <- list(
        pop      = .subset2(dat,"Population"),
        line     = .subset2(dat,"Line"),
        mol      = .subset2(dat,"Molecule"),
        trial    = .subset2(dat,"Reference.trial"),
        endpoint = .subset2(dat,"Endpoint")
      )
      
      # Extract the PLMTE for this destination, so that we have a list to manipulate
      # directly without fear of overwriting anything :)
      plmte <- prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]
      
      # Derive the origin from dat, these values can then go in plmte$orig
      o <- list(
        pop      = paste0("pop_"     , .subset2(dat,"Population")),
        line     = paste0("line_"    , .subset2(dat,"Line")),
        mol      = paste0("mol_"     , .subset2(dat,"Reference.treatment")),
        trial    = paste0("trial_"   , .subset2(dat,"Reference.trial")),
        endpoint = paste0("endpoint_", .subset2(dat,"Endpoint"))
      )
      
      # Filter down hr_table to matching destination given origin reference treatment
      # and keep the time-varying HRs.
      fp <-
        list(HR = hr_table[Population == dN$pop &
                             Line == dN$line &
                             Molecule == dN$mol &
                             Reference.treatment == .subset2(dat,"Reference.treatment") &
                             Reference.trial == dN$trial & 
                             Endpoint == dN$endpoint, HR])
      
      if(length(fp$HR) < time_horizon+1) {
        fp_extend <- rep(fp$HR[length(fp$HR)],((time_horizon+1)-length(fp$HR)))
        fp$HR <- c(fp$HR,fp_extend)[1:(time_horizon+1)]
      }
      
      # Place the relevant information into the destination PLMTE. Whether or not
      # this information is used in the end is decided in the function which follows
      # this one, called f_NMA_AddAssumptionsToNetwork
      if (!is.null(o$pop))      plmte$orig$pop      <- o$pop
      if (!is.null(o$line))     plmte$orig$line     <- o$line
      if (!is.null(o$mol))      plmte$orig$mol      <- o$mol
      if (!is.null(o$trial))    plmte$orig$trial    <- o$trial
      if (!is.null(o$endpoint)) plmte$orig$endpoint <- o$endpoint
      plmte$fp <- fp
      
      # Take the updated plmte, being careful not to overwrite anything and 
      # slot it back into its destination, returning the updated network.
      prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
      return(prev)
    }
  )
}


#' Function to take results after running f_NMA_generateNetwork and then f_NMA_linkPHNMA
#' and then apply assumption-based HRs (e.g. same as another treatment/endpoint etc)
#' 
#' WARNING: I expect the table in excel to be complete and to provide the HRs 
#'          required, rather than working them out based on the drop down selection
#'          The function below simply adds the HR in question to the efficacy network
#'          list. It won't work it out for you!!
#'          
#'          
#' NOTE: This function is ONLY for those rows that are NOT trial survival analysis,
#'       FP NMA or PH NMA. Everything else is assumption based, and that is what
#'       this function applies!
#'          
#' @param network `p$releff$network`
#' @param phnma_table `p$releff$CODA$PH`
#' @param fpnma_table `p$releff$CODA$FP`
#' @param fpnma_destinations `p$releff$fp_dest`
#' @param excel_table `data.table(i$R_table_eff_data_settings)`
#' @param trial_flag `i$List_eff_datasources[1]`
#' @param fpnma_flag `i$List_eff_datasources[3]`
#' @param phnma_flag `i$List_eff_datasources[2]`
#' @param et_flag `i$List_eff_datasources[4]`
#' @param ahr_flag `i$List_eff_datasources[5]`
#' @param verbose `qc_mode`
#' 
f_NMA_AddAssumptionsToNetwork <-
  function(network,
           phnma_table,
           fpnma_table,
           fpnma_destinations,
           excel_table,
           trial_flag = "Trial survival anlaysis",
           fpnma_flag = "FP_NMA",
           phnma_flag = "PH_NMA",
           et_flag = "Assume equal to",
           ahr_flag = "Apply HR to",
           verbose = FALSE,
           psa_flag = FALSE) {

  Reduce(
    x = 1:nrow(excel_table),
    init = network,
    accumulate = FALSE,
    f = function(prev, excel_tab_row) {
      
      dat <- as.list(excel_table[excel_tab_row,])
      
      # add a message showing progress and what's happening:
      if (all(verbose,dat$Include.in.this.analysis. == "Yes")) {
        cat(paste0(
          paste0(
            "R_table_eff_data_settings row ",excel_tab_row," - dest: pop_",
            dat$Population,
            "$line_",
            dat$Treatment.line,
            "$mol_",
            dat$Molecule,
            "$trial_",
            dat$Origin.trial,
            "$endpoint_",
            dat$End.point,
            " | method: ",
            dat$Effectiveness.data.source,
            " | orig: pop_",
            dat$Origin.population,
            "$line_",
            dat$Origin.line,
            "$mol_",
            dat$Origin.treatment,
            "$trial_",
            dat$Origin.trial,
            "$endpoint_",
            dat$Origin.endpoint,
            "\n"
          )
        ))
      }
      
      # Cycle through the different possibilities, handling the case for different
      # dropdown options one at a time. The cases are:
      # 
      # - trial_flag: directly based on PLD extrapolation. sampling is done using parameters outside of this function.
      # - phnma_flag: PHNMA. HRs are linked corresponding to the CODA sample, but overrides are to be enacted. CODA sampling is done externally
      # - fpnma_flag: PHNMA. time-varying HRs are linked corresponding to the input data, but overrides are to be enacted. sampling done within.
      # - Assume equal to: put in origin and set HR to 1, simple :). doesn't vary in PSA
      # - apply HR to: put in origin and set HR to the "HR to apply" column in the table. also simple. for PSA use bounds externally, replacing the values in the table
      
      d <- list(
        pop      = paste0("pop_",dat$Population),
        line     = paste0("line_",dat$Treatment.line),
        mol      = paste0("mol_",dat$Molecule),
        trial    = paste0("trial_",dat$Origin.trial), 
        endpoint = paste0("endpoint_",dat$End.point)
      )
      dN <- list(
        pop      = dat$Population,
        line     = dat$Treatment.line,
        mol      = dat$Molecule,
        trial    = dat$Origin.trial, 
        endpoint = dat$End.point
      )
      # If the user in Excel has excluded by selecting no in the inclusion dropdown,
      # exclude it!
    if (dat$Include.in.this.analysis. == "No") {
        prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$include <- FALSE
        if (all(verbose,dat$Include.in.this.analysis. == "No")) {
          cat(paste0(
            paste0(
              "R_table_eff_data_settings row ",excel_tab_row," - dest: pop ",
              dat$Population,
              " line ",
              dat$Treatment.line,
              " mol ",
              dat$Molecule,
              " trial ",
              dat$Origin.trial,
              " endp ",
              dat$End.point,
              " | method: EXCLUDED FROM ANALYSIS BY EXCEL. INCLUDE IN THIS ANALYSIS SET TO NO FOR THIS PLMTE",
              "\n"
            )
          ))
        }
        return(prev)
      } else {
        plmte <- f_misc_get_plmte(prev,d)
        plmte$include <- TRUE
      }
      
      # Do a quick check on whether effectiveness data source has gone blank for
      # some reason!
      if (all(
        is.null(plmte$orig$source),
        !is.null(dat$Effectiveness.data.source)
      )) {
        plmte$orig$source <- dat$Effectiveness.data.source
      }
      
      
      if (plmte$orig$source == trial_flag) {
        
        # Simple. extrapolation based directly on PLD. The user in excel has indicated
        # that this is one of the reference curves, overriding everything else.
        # Origin should be the same as destination.
        
        d <- list(
          pop      = paste0("pop_",dat$Population),
          line     = paste0("line_",dat$Treatment.line),
          mol      = paste0("mol_",dat$Molecule),
          trial    = paste0("trial_",dat$Origin.trial), 
          endpoint = paste0("endpoint_",dat$End.point)
        )
        dN <- list(
          pop      = dat$Population,
          line     = dat$Treatment.line,
          mol      = dat$Molecule,
          trial    = dat$Origin.trial, 
          endpoint = dat$End.point
        )
        
        # populate the origin
        o <- list(
          pop      = paste0("pop_",dat$Origin.population),
          line     = paste0("line_",dat$Origin.line),
          mol      = paste0("mol_",dat$Origin.treatment),
          trial    = paste0("trial_",dat$Origin.trial), 
          endpoint = paste0("endpoint_",dat$Origin.endpoint)
        )
        o$dist   <- dat$Curve.fit..for.survival.analysis.
        o$source <- dat$Effectiveness.data.source
        plmte$orig <- o
        plmte$dest <- d
        plmte$fp <- list()
        
        # slot the plmte back into the list and return the result:
        prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
        return(prev)
        
        
      } else if (plmte$orig$source %in% c(phnma_flag, fpnma_flag)) {
        
        d <- list(
          pop      = paste0("pop_",dat$Population),
          line     = paste0("line_",dat$Treatment.line),
          mol      = paste0("mol_",dat$Molecule),
          trial    = paste0("trial_",dat$Origin.trial), 
          endpoint = paste0("endpoint_",dat$End.point)
        )
        dN <- list(
          pop      = dat$Population,
          line     = dat$Treatment.line,
          mol      = dat$Molecule,
          trial    = dat$Origin.trial, 
          endpoint = dat$End.point
        )
        oN <- list(
          pop      = dat$Origin.population,
          line     = dat$Origin.line,
          mol      = dat$Origin.treatment,
          trial    = dat$Origin.trial, 
          endpoint = dat$Origin.endpoint
        )
        
        
        
        # Now, a check for an error in the excel file:
        if (dN$mol == dat$Origin.treatment) {
          warning(paste0(
            "In row ",
            excel_tab_row,
            ", origin molecule (",
            dN$mol,
            ") is set to the SAME as destination molecule (",
            dat$Origin.treatment,
            ") in Excel. Please fix your human error. I cannot generate s(t)",
            " for ", paste(d,collapse = "$"), "."
          ))
          return(prev)
        }
        
        
        if (dat$Effectiveness.data.source == phnma_flag) {
          # Match on the phnma table as there's only 1 row per destination
          tab_nma_match <- phnma_table[Population == dN$pop &
                                         Line     == dN$line &
                                         Molecule == dN$mol &           # destination molecule to populate
                                         Reference.treatment== oN$mol & # origin molecule to apply HR to
                                         Endpoint == dN$endpoint,]
        } else {
          # Match on the "destinations" for the FPNMA to populate across trials:
          tab_nma_match <- fpnma_destinations[ Population           == dN$pop &
                                                Line                == dN$line &
                                                Molecule            == dat$Molecule &
                                                Reference.treatment == dat$Origin.treatment &
                                                Reference.trial     == oN$trial &
                                                Endpoint            == dN$endpoint, ]
        }
        
        # If there is one hazard ratio, put it in:
        if (nrow(tab_nma_match) == 1) {
          # Take the ORIGIN from the excel table, superseding the linkage from
          # the original PHNMA.
          # 
          # EXPLANATION: If the NMA is linking A to B, but the user wants to 
          #              apply the HR linking A to B to extrapolations derived
          #              from C, then the A vs B HR should be applied to dataset
          #              C. If it should be applied to reference curve B, then
          #              the data in the Excel table should have origin columns
          #              which match with the CODA sample for the NMA. If they
          #              are different, (e.g. apply HR from trials to RWE-based
          #              extrapolations per our base case!) then we should
          #              use the origin from the excel file and apply it to the
          #              HR that's given from the CODA sample!!!
          
          plmte$orig$pop      <- paste0("pop_",dat$Origin.population)
          plmte$orig$line     <- paste0("line_",dat$Origin.line)
          plmte$orig$mol      <- paste0("mol_",dat$Origin.treatment)
          plmte$orig$trial    <- paste0("trial_",dat$Origin.trial)
          plmte$orig$endpoint <- paste0("endpoint_",dat$Origin.endpoint)
          plmte$orig$source   <- dat$Effectiveness.data.source
          
          # Take the HR from the PH NMA, irrespective of what its reference trial
          # is. Note that the distributional selection DOES NOT MATTER, as 
          # during propagation the original reference curve will propagate via
          # those rows in excel that are set to "Trial survival analysis". See above
          # 
          # e.g. prev$pop_1$line_1$mol_7$trial_0$endpoint_0 should have no HR but should
          # have PLD-based extraps.
          # 
          
          if(dat$Effectiveness.data.source == phnma_flag) {
            plmte$hr <- tab_nma_match$HR
          }
          
          # Now that we've updated the origin for this plmte and have added the HR
          # that we want to apply, we can slot the plmte back into prev, and return
          # that as we're done with this row.
          
          prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
          return(prev)
          
        } else if (nrow(tab_nma_match) == 0) {
          warning(paste0(
            "Row: ",
            excel_tab_row,
            " Dest: ",
            paste(unlist(d), collapse = "$"),
            " is using ",
            dat$Effectiveness.data.source,
            ", but there is no NMA efficacy data which corresponds to it in the CODA sample!"
          ))
          return(prev)
        }
        
      } else if (dat$Effectiveness.data.source %in% c(et_flag, ahr_flag)) {
        
        # Either we're applying a HR from the excel table, or we're assuming equal.
        # Either way, we need the destination, the plmte and the origin:
        d <- list(
          pop      = paste0("pop_",dat$Population),
          line     = paste0("line_",dat$Treatment.line),
          mol      = paste0("mol_",dat$Molecule),
          trial    = paste0("trial_",dat$Origin.trial), # NOTE THE ASSUMPTION HERE ON trial==origin.trial
          endpoint = paste0("endpoint_",dat$End.point)
        )
        o <- list(
          pop      = paste0("pop_",dat$Origin.population),
          line     = paste0("line_",dat$Origin.line),
          mol      = paste0("mol_",dat$Origin.treatment),
          trial    = paste0("trial_",dat$Origin.trial), 
          endpoint = paste0("endpoint_",dat$Origin.endpoint),
          dist     = NULL,
          source   = dat$Effectiveness.data.source
        )
        
        # override the info that's been put in there for origin and destination
        # from the NMAs by the linkages being put in from the excel file. this
        # then allows the user to apply relative efficacy estimated by the NMA
        # to data not used in the NMA (e.g. apply NMA-based HRs to RWE rather than
        # trial PLD)
        plmte$orig <- o
        plmte$dest <- d
        
        # If it's a hazard ratio and not assume the same, then put in the HR
        if (dat$Effectiveness.data.source != et_flag) {
          # Note for PSA we could use rlnorm or such here, including deriving meanlog and sdlog from bounds.
          if (psa_flag) {
            plmte$hr <- rnorm(1,dat$HR.to.apply,(dat$HR.95..CI..UCL.-dat$HR.to.apply)/1.96)
          } else {
            plmte$hr <- dat$HR.to.apply 
          }
        } else {
          plmte$hr <- 1
        }
        
        plmte$fp <- list()
        
        # either way slot in the origin information directly (we don't need anything more)
        prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
        
        return(prev)
      } else {
        stop(paste0(
          "In the effectivness settings in excel, every row with 'Yes'",
          " for inclusion MUST have an effectiveness data source that is ",
          trial_flag, ", ",phnma_flag,", ",fpnma_flag,", ",et_flag,", or ", ahr_flag,
          ". However, pop ", dat$Population.name, " line ", dat$Treatment.line, " molecule ",
          dat$Molecule, " endpoint ", dat$End.point.name, " has a value of ",dat$Effectiveness.data.source,
          " there! Please find and fix in the table"
          ))
      }

    }
  )
}



#' Function to implement a fractional polynomial network meta analytic result (time-varying
#' hazard ratio) by multiplying the baseline hazard of the reference curve by the HR_t
#' and then recalculating the survival curve afterwards.
#' 
#' @param ref_curve the reference curve to apply the time-varying HRs to
#' @param HR the time-varying hazard ratios to apply
#' 
f_FPNMA_implement <- function(ref_curve, HR) {
  
  # Making sure the two vectors have the same length:
  stopifnot(length(ref_curve) == length(HR)) 
  
  # Calculate cumulative h(t), reference baseline h(t), apply HR, cumsum, calculate s(t)
  ref_cum_haz      <- -log(ref_curve) 
  ref_baseline_haz <- c(0,diff(ref_cum_haz))
  int_baseline_haz <- ref_baseline_haz * HR 
  int_cum_haz      <- cumsum(int_baseline_haz)
  return(exp(-int_cum_haz))
}


#' Function to cycle through all possible pop line mol trial endpoint getting the
#' selected distribution for data-based, using that to apply hazard ratios to populate
#' "destination" entries
#' 
#' As the comparative efficacy network provided and the TSD14 survival analysis all
#' follow identical structures (nesting of population line molecule trial endpoint)
#' and naming within those list structures, we can leverage this to essentially
#' replicate network but with an entry (wherever possible) for survival at time t
#' or "st" for short.
#' 
#' This uses the "dist" entry in "origin" and "destination" entries to populate 
#' those extrapolations that are not informed by either the NMA or assumption. 
#' 
#' The origin and destination lists for each pop line mol trial endpoint are then used
#' to figure out the distributional selection for "origin" extrapolations (i.e. where
#' the data is coming from). The function then takes the correct
#' distribution from the correct place (the "origin" location in extraps) and
#' applies the correct hazard ratio to it (the HR within the "destination" location) 
#' to get the resulting extrapolated survival, to then go into the destination
#' under entry st.
#' 
#' In other words, this function gets the right information from the right place
#' and uses it to generate all possible extrapolations given the data that is 
#' provided. 
#' 
#' This function DOES NOT check which treatment pathways can and cannot be
#' simulated given the evidence provided. This is left to another function and
#' this one simply proliferates what can be proliferated!
#' 
#' 
#' @param network list resulting from f_NMA_generateNetwork, f_NMA_linkPHNMA, and then f_NMA_AddAssumptionsToNetwork
#' @param extraps result of applying f_surv_getExtrapolations to the TSD14 output from f_surv_runAllTSD14
#' @param dos degrees of separation between original data and extrapolation (e.g. PFS for A (data) to PFS for B (HR1) to PFS for C (HR2) requires 2)
#' @param excel_table named range `R_table_eff_data_settings` from excel, as a `data.table` object
#' @param verbose extra console output if true
#' @param dist_lookups lookup table for distributions
#' @param psa_lambda_flag flag for using lambda approximation to improve efficiency in the PSA. default FALSE
#' @param psa_iteration optional - if `psa_lambda_flag` is TRUE or PSA is true (to be built) it is required.
#' @param th `NULL` except for PSA with lambda, which needs TH to extrapolate the curve should be th not the number of cycles in the model (i.e. that minus one)
#' @param psa_params optional - if `psa_lambda_flag` is TRUE or PSA is true (to be built) it is required.
#' 
#' 
f_releff_PropNetwork <-
  function(network,
           extraps,
           excel_table,
           dos = 5,
           verbose = FALSE,
           dist_lookups,
           psa_lambda_flag = FALSE,
           psa_iteration = NULL,
           psa_params = NULL,
           th = NULL) {
    
  
  dl <- data.table(dist_lookups)
  excel_table$xl_rn <- 1:nrow(excel_table)
  ref_curves <- excel_table[Include.in.this.analysis. == "Yes" & Effectiveness.data.source == "Trial survival analysis",]
  included   <- excel_table[Include.in.this.analysis. == "Yes",]
  releff     <- excel_table[Include.in.this.analysis. == "Yes" & Effectiveness.data.source != "Trial survival analysis",]
  
  # Add some row numbers for ease of cross referencing things
  ref_curves$ref_rn <- 1:nrow(ref_curves)
  included$inc_rn   <- 1:nrow(included)
  releff$rel_rn     <- 1:nrow(releff)
  
  # The reason there are 2 nested reduce statements here is that we need to
  # cycle through degrees of separation, for each degree of separation we then
  # need to go down the list of those PLMTEs that are included (according to
  # excel_table) testing whether or not something can be done at that dos, then
  # doing it if possible, leaving it if there's nothing at the origin (i.e. 1 dos away)

  
  
  # Reduce 1: degree of separation - essentially just repeating the inner reduce
  #           dos times, but keeping a cumulative result
  Reduce(
    x = 1:dos,
    init = network,
    accumulate = FALSE,
    f = function(prev_network, deg_of_sep) {
      
      # The first time, do
      # prev_network <- network 
      # 
      # to feed in the original network and 
      # deg_of_sep <- 1
      # 
      # To do the inner loop, run the inside, then do
      # 
      # prev_network <- prev_network_this
      # deg_of_sep <- deg_of_sep + 1
      # 
      # then run again. 
      
      # It hugely simplifies the degree of separation issue if all of the "starting
      # points" are entered into the final network first. this allows relative
      # efficacy to fan out from all INITIAL ORIGIN s(t) right from the start
      # without having to wait. This hugely simplifies the process at higher
      # dos, which is then simply is there a ref curve, if no do nothing, if so
      # apply releff.
      
      if (deg_of_sep == 1) {
        prev_network_this <- Reduce(
          x = 1:nrow(ref_curves),
          init = prev_network,
          accumulate = FALSE,
          f = function(prev, ref_curve_row) {
            
            dat <- ref_curves[ref_curve_row,]
            
            # get the destination so we know where to put the extrapolation
            d <- list(
              pop = paste0("pop_",dat$Population),
              line = paste0("line_",dat$Treatment.line),
              mol = paste0("mol_",dat$Molecule),
              trial = paste0("trial_",dat$Origin.trial),
              endpoint = paste0("endpoint_",dat$End.point)
            )
            o <- list(
              pop = paste0("pop_",dat$Origin.population),
              line = paste0("line_",dat$Origin.line),
              mol = paste0("mol_",dat$Origin.treatment),
              trial = paste0("trial_",dat$Origin.trial),
              endpoint = paste0("endpoint_",dat$Origin.endpoint)
            )
            
            # pull out that plmte from prev
            plmte <- f_misc_get_plmte(prev, d)
            
            # Make sure excel is determinining the origin, superseding
            # if it came from the NMA
            plmte$orig[c("pop","line","mol","trial","endpoint")] <- o
            
            if(verbose) {
              cat(
                paste0(
                  "ref curve #",
                  ref_curve_row,
                  ", dest: ",
                  paste(d, collapse = "$"),
                  " orig: ",
                  paste(o, collapse = "$"),
                  "\n"
                )
              )
            }
            
            # pull out the distributional selection from it:
            if (is.null(plmte$orig$dist)) plmte$orig$dist <- dat$Curve.fit..for.survival.analysis.
            dst <- dl[Description == plmte$orig$dist,]$RCC_input_desc
            
            if (is.null(plmte$orig$source)) plmte$orig$source <- dat$Effectiveness.data.source
            
            # check extrapolations available
            # If it's the PSA and we're using lambda approximation, then we
            # have a THxnPSA matrix here (all exponential), 
            # otherwise we have a THx7 matrix by dist
            if (psa_lambda_flag) {
              # When it's the PSA and the lambda approximation method is being used
              # then extraps is a set of lambdas, each a fixed rate for that PLMTE.
              # Extrapolate that to get curves.
              
              # in first-line therapy, we preserve the distributional shape
              if (d$line == "line_1") {
                par <- f_misc_get_plmte(psa_params,d)
                dr  <- .subset2(par,"draws")
                dis <- par$id$dist
                if (dis == "exp") {
                  plmte$st <- f_extrapolate(0:th,dr[psa_iteration],dis)
                } else {
                  plmte$st <- f_extrapolate(0:th,dr[psa_iteration,],dis)
                }
              } else {
                # Otherwise we use the lambda we calculated earlier to remove the
                # need for tunnel states for the PSA only.
                plmte$st <- f_psa_exp(0:th,f_misc_get_plmte(extraps,o)[psa_iteration])
              }
              plmte$populated <- TRUE
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
              return(prev)
            } else if (all(class(f_misc_get_plmte(extraps,o)) %in% c("matrix", "array"))) {
              plmte$st <- f_misc_get_plmte(extraps,o)[,dst]
              plmte$populated <- TRUE
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
              return(prev)
            } else {
              warning(paste0("ref curve #",ref_curve_row, ", destination: ", paste(d,collapse="$"), " is set to from trial analysis but there's no extrapolations at that location"))
              return(prev)
            }
          }
        )
        return(prev_network_this)
      } else {
        # We've already inserted all the reference curves, so we can fan out
        # from each of them one dos at a time.
        
        # Reduce 2: inclusion table - cycling through all the PLMTEs which are
        #           included with a "Yes" in the column in excel efficacy settings
        #           sheet. Each one gets tested for being populated already, then
        #           method, then the origin is extracted. if no data in the origin
        #           nothing can be done. if data in the origin then apply the
        #           correct method.
        
        # To run line by line:
        # 
        # prev_network <- prev_network_this
        # deg_of_sep <- deg_of_sep + 1
        # prev <- prev_network
        # 
        # set the releff_row to be the row in releff you want to apply
        # 
        # Whilst inside, to look up a row (row numbers were added at the beginning):
        # 
        # releff[Population == 0 & Treatment.line == 4 & Molecule == 5 & End.point == 2, ]
        # 
        # Just change the numbers to fit. put that number in releff_row to go 
        # through that line of the table
        # 
        # releff_row <- releff[Population == 0 & Treatment.line == 2 & Molecule == 0 & End.point == 2, ]$rel_rn
        # 
        
        prev_network_this <- Reduce(
          # x = 1:63,
          x = 1:nrow(releff),
          init = prev_network,
          accumulate = FALSE,
          f = function(prev, releff_row) {
            dat <- releff[releff_row, ]
            
            # Get the destination from the table. Note the assumption on trial
            d <- list(
              pop     = paste0("pop_"      , dat$Population),
              line    = paste0("line_"     , dat$Treatment.line),
              mol     = paste0("mol_"      , dat$Molecule),
              trial   = paste0("trial_"    , dat$Origin.trial),
              endpoint = paste0("endpoint_", dat$End.point)
            )
            
            # pull out the PLMTE
            plmte <- f_misc_get_plmte(prev, d)
            
            # Check if done already, if so return prev to avoid repetition
            if (all("populated" %in% names(plmte), plmte$populated == TRUE)) return(prev)
            
            # get information on the origin from plmte, then pull out that data
            o <- plmte$orig
            orig_data <- f_misc_get_plmte(prev, o)
            
            # Error check - if the plmte and origin plmte are identical then 
            # something is wrong (shouldn't ever be the case for non trial data analysis):
            if (identical(orig_data,plmte)) {
              if (plmte$orig$source == "Assume equal to") {
                stop(paste0(
                  "Circular reference at: ",
                  paste(unlist(d), collapse = "$"),
                  " (entry ", dat$xl_rn, ". See Column BQ in Excel effectiveness settings page)",
                  ". Method is 'Assume equal to' and origin=destination."
                ))
              }
            }
            
            
            # Check if the origin is empty for this destination. shouldn't happen
            # much but worth a failsafe:
            if (any(
              is.null(o$pop),
              is.null(o$line),
              is.null(o$mol),
              is.null(o$trial),
              is.null(o$endpoint),
              is.null(o$source)
            )) {
              warning(paste0(
                "Entry at ",
                paste(d$pop, d$line, d$mol, d$trial, d$endpoint, sep = "$"),
                " has missing origin data. This should've been entered via f_NMA_AddAssumptionsToNetwork",
                " earlier. Attempting to recover..."
              ),immediate. = TRUE)
              
              # Attempt to recover by getting the settings from the excel table
              # by matching up the entries
              dN <- list(
                pop     = dat$Population,
                line    = dat$Treatment.line,
                mol     = dat$Molecule,
                trial   = dat$Origin.trial,
                endpoint = dat$End.point
              )
              o <- as.list(included[Population     == dat$Population & 
                                 Treatment.line == dat$Treatment.line &
                                 Molecule       == dat$Molecule &
                                 Origin.trial   == dat$Origin.trial &
                                 End.point    == dat$End.point 
                               , list(Origin.population,Origin.line,Origin.treatment,Origin.trial, Origin.endpoint,Effectiveness.data.source)])
              names(o)   <- c("pop","line","mol","trial","endpoint","source")
              o$pop      <- paste0("pop_",o$pop)
              o$line     <- paste0("line_",o$line)
              o$mol      <- paste0("mol_",o$mol)
              o$trial    <- paste0("trial_",o$trial)
              o$endpoint <- paste0("endpoint_",o$endpoint)
              
              
              # slot them into plmte
              plmte$orig$pop <- o$pop
              plmte$orig$line <- o$line
              plmte$orig$mol <- o$mol
              plmte$orig$trial <- o$trial
              plmte$orig$endpoint <- o$endpoint
              plmte$orig$source <- o$source
              
              # update prev to keep a permanent record of what we've done,
              # making sure to only replace those items which exist
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$pop <- o$pop
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$line <- o$line
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$mol <- o$mol
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$trial <- o$trial
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$endpoint <- o$endpoint
              
              # get information on the origin from plmte, then pull out that data
              o <- plmte$orig
              orig_data <- f_misc_get_plmte(prev, o)
            }
            
            # if that data has s(t) in it, we can do something, otherwise we
            # can't.
            
            if (!"populated" %in% names(orig_data)) {
              # We can't do anything here yet. If orig is busted we already tried
              # to fix it above.
              return(prev)
            } else {
              # there is data in the orig_data to apply relative efficacy to.
              # we can apply it
              
              # limit data to the model time horizon
              
              orig_data$st <- orig_data$st[1:(p$basic$th+1)]
              plmte$fp$HR <- plmte$fp$HR[1:(p$basic$th+1)]
              
              if (verbose) {
                cat(
                  paste0(
                    "dos ",
                    deg_of_sep,
                    " | releff row: ",
                    releff_row,
                    " | excel row: ",
                    dat$xl_rn,
                    " | D: ",
                    paste(d$pop, d$line, d$mol, d$trial, d$endpoint, sep = "$"),
                    "\t | O: ",
                    paste(o$pop, o$line, o$mol, o$trial, o$endpoint, sep = "$"),
                    ": Method = ",
                    o$source,
                    "\n"
                  )
                )
              }
              method <- o$source
              
              # note that if source was missing, we already caught that above
              # and tried to fix it. if it's still blank it's blank in excel.
              if (is.null(method)) {
                warning(paste0(
                  "dos ",
                  deg_of_sep,
                  " | row: ",
                  releff_row,
                  " | Destination: ",
                  paste(d$pop, d$line, d$mol, d$trial, d$endpoint, sep = "$"),
                  " Origin: ",
                  paste(o$pop, o$line, o$mol, o$trial, o$endpoint, sep = "$"),
                  ": o has no method in it, trying to recover..."
                ),immediate. = TRUE)
                
                o$source <- dat$Effectiveness.data.source
                method   <- dat$Effectiveness.data.source
                prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]]$orig$source <- dat$Effectiveness.data.source
                
                if (is.null(o$source)) stop("Failed...row in excel has no method, network has no method. shouldn't happen!")
                
                return(prev)
              }
              
              # Now we can get on with applying relative efficacy
              if (method == "Assume equal to") {
                # if assuming equal to, then set equal to, update populated, slot in, done.
                plmte$st <- orig_data$st
                if (length(plmte$st) > 0) {
                  plmte$populated <- TRUE
                } else {
                  plmte$populated <- NULL
                }
              } else if (method %in% c("Apply HR to", "PH_NMA")) {
                # simply apply hazard ratio, update, slot in, return.
                plmte$st        <- orig_data$st ^ plmte$hr
                if (length(plmte$st) > 0) {
                  plmte$populated <- TRUE
                } else {
                  plmte$populated <- NULL
                }
              } else if (method == "FP_NMA") {
                # Apply FP NMA using the function, update, slot in, return.
                # To plot:
                # plot(orig_data$st,type="l")
                # lines(f_FPNMA_implement(orig_data$st, plmte$fp$HR),col="red")
                plmte$st <- f_FPNMA_implement(orig_data$st, plmte$fp$HR)
                if (length(plmte$st) > 0) {
                  plmte$populated <- TRUE
                } else {
                  plmte$populated <- NULL
                }
              } else{
                stop(
                  paste0(
                    "Destination ",
                    paste(unlist(d), collapse = "$"),
                    " origin ",
                    paste(unlist(o), collapse = "$"),
                    " has an invalid dropdown menu selection in the 'Effectiveness data source' column in Excel 'Effectiveness settings' sheet"
                  )
                )
              }
              # slot in the updated PLMTE and return the result.
              prev[[d$pop]][[d$line]][[d$mol]][[d$trial]][[d$endpoint]] <- plmte
              return(prev)
            }
            
          }
        )
        return(prev_network_this)
      }
    })
}





# Misc functions ----------------------------------------------------------

#' Function to add the estimated hazard of two curves. Useful for getting
#' TTD not censoring for death for the partitioned survival model.
f_surv_hazardadd <- function(s_1t, s_2t) {
  s_1t * s_2t
}

#' Calculate probability of dying in each cycle from the survival curve
f_surv_get_q_t <- function(s_t) {
  data.table::nafill(1 - (data.table::shift(s_t, type = "lead") / s_t), "locf")
}

#' Estimate the hazard function from the survival curve
f_surv_get_h_t <- function(s_t, cl_yr = 1/52) {
  data.table::nafill((log(s_t) - log(data.table::shift(s_t, type = "lead"))), "locf") / cl_yr
}


# Computing extrapolations in CE model ----------------------------------------------------------


#' Function to go into the regression structure generated by f_surv_runAllTSD14 and
#' pull out just the extrapolations (st for survival at time t) for all PLMTEs
#' 
#' @param regs regression list resulting from f_surv_runAllTSD14
#' 
f_surv_getExtrapolations <- function(regs) {
  lapply(regs, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        lapply(mol, function(tr) {
          lapply(tr, function(plmte) {
            plmte$st
          })
        })
      })
    })
  })
}

#' Function which checks whether a treatment sequence can be estimated given the
#' extrapolations which are currently available
#' 
#' @param treatment_sequence a string vector containing the molecule identifiers (e.g. mol_1 mol_4 mol_7 mol_999)
#' @param st calcualted extrapolations with comparative efficacy applied and before transformation into TPs
#' @param lookups vlookup table list from excel (i.e. `i$lookup`)
#' @param pop_n population as a number (0, 1, 2 or 3)
#' 
f_seq_extrapCollector <- function(treatment_sequence, st, lookups, pop_n = 0, pop_0_2Lp = TRUE, required_endpoints = paste0("endpoint_",0:4)) {
  
  len_ts <- length(treatment_sequence)
  
  # Make the name for the population
  risk_pop <- paste0("pop_",pop_n)
  
  # labels to cycle through for lapply operations for plmte stuff
  line_labs <- structure(paste0("line_",1:len_ts),names=paste0("line_",1:len_ts))
  endp_labs <- structure(paste0("endpoint_",0:4),names=paste0("endpoint_",0:4))
  
  # bodge error for now, genericisation must allow subgroups at later lines.
  if (pop_0_2Lp == FALSE) stop("pop_0_2Lp must be TRUE. Function forces assumption of 2L+ risk population 0. contact developers if you want something else!")
  
  # Derive which are populated and which aren't, simplifying into tables:
  
  which_populated <- lapply(line_labs, function(li) {
    
    if (li == "line_1") {
      i_st <- st[[risk_pop]][[li]]
    } else {
      i_st <- st[["pop_0"]][[li]]
    }
    
    lapply(i_st, function(mol) {
      do.call(
        rbind,
        lapply(mol, function(tr) {
          tst <- unlist(lapply(tr, function(plmte) {
            "populated" %in% names(plmte)
          }))
          # Now, every so often not all endpoints make it this far
          # so those that aren't there need to get added back in with a FALSE
          if (all(required_endpoints %in% names(tr))) {
            return(tst)
          } else {
            # put a FALSE in the place that was missing and then ensure it's
            # in the right order
            which_names_missing <- required_endpoints[which(!required_endpoints %in% names(tr))]
            tst[which_names_missing] <- FALSE
            return(tst[required_endpoints])
          }
        })
      )
    })
  })
  
  # Narrow it down to only those relevant to this pathway:
  # shorten the treatment sequence name and derive a numerical version, then trt names!
  ts     <- treatment_sequence
  ts_n   <- as.integer(gsub(pattern = "mol_",replacement = "",ts))
  ts_lab <- names(lookups$trt[which(lookups$trt %in% ts_n)])
  
  availability_tables <- lapply(1:len_ts, function(trt_line) {
      which_populated[[attr(treatment_sequence[trt_line], "names")]][[treatment_sequence[trt_line]]]
    })
  names(availability_tables) <- ts_lab
  
  # Now identify the trial that's recorded for each PLMTE for this TS.
  which_trials <- lapply(availability_tables, function(trt_line) {
    rn <- rownames(trt_line)
    
    if (length(rn) == 1) {
      temp_vec <- ifelse(trt_line,as.numeric(gsub("trial_","",rn))+1,NA)[rn,]
      temp_vec[!is.na(temp_vec)]
    } else {
      temp_tab <- do.call(
        rbind,
        lapply(1:nrow(trt_line), function(the_row) {
          ifelse(trt_line[the_row,],as.numeric(gsub("trial_","",rownames(trt_line)[the_row]))+1,NA)
        })
      )
      unlist(lapply(1:ncol(trt_line), function(the_column) {
        temp_tab[which(!is.na(temp_tab[,the_column])),the_column]
      }))
    }
  })
  
  # error check - there should only ever be one trial basis for PLMTE. if more than
  # one then we need a load of inputs adding to excel and a much more complicated function.
  for (tr in 1:len_ts) {
    if (any(lapply(which_trials[[tr]],length) > 1)) {
      stop(
        paste0(
          "In treatment sequence ",
          paste(treatment_sequence, collapse = " "),
          ", Multiple trials available to inform ",
          tr,
          "L extrapolations! these are ",
          which_trials[[tr]],
          "...FUNCTION AND EXCEL NEED TO BE UPDATED TO COPE WITH THIS!!!"
        )
      )
    }
  }
  
  # pull out the s(t) for each plmte which is to be used for this treatment
  # sequence:
  
  st_ts <- lapply(1:len_ts, function(trt_line) {
    pop_txt  <- ifelse(trt_line == 1, risk_pop, "pop_0") # later lines assumed risk pop 0 (all patients)
    line_txt <- paste0("line_",trt_line)
    tr_ids   <- names(which_trials[[trt_line]])
    mol_txt  <- treatment_sequence[line_txt]
    
    st_line <- lapply(1:length(tr_ids), function(endp_n) {
      tr_txt   <- paste0("trial_",which_trials[[trt_line]][endp_n]-1)
      endp_txt <- paste0("endpoint_",endp_n-1)
      st[[ifelse(trt_line == 1, risk_pop, "pop_0")]][[line_txt]][[mol_txt]][[tr_txt]][[endp_txt]]
    })
    names(st_line) <- lookups$ipd$endpoint[match(1:length(tr_ids)-1,Number),]$Description
    
    return(st_line)
    
  })
  names(st_ts) <- paste0("line_",1:len_ts)
  
  # So now we've collected all of the relevant extrapolations for this TS. 
  # The next step is to figure out which method/assumptions can/should be applied to 
  # derive transition probabilities between all of the model health states
  # for this treatment pathway. However, this will be handled in another function
  # as this one is already long.
  
  return(list(
    availability = availability_tables,
    trials       = which_trials,
    st           = st_ts
  ))
}


#' Function to maximize estimated hazards between two curves.
#'
#' takes s_t and the reference s_t to maximize hazards against.
#' 
#' As the function with the greatest mortality probability will also have the
#' greatest hazard, uses mortality probabilities as they should be quicker to
#' calculate
f_surv_hazardmax <- function(s_t, reference_s_t) {
  th_c    <- length(s_t)
  h_t_fixed <- pmax(f_surv_get_q_t(s_t)[1:th_c],f_surv_get_q_t(reference_s_t)[1:th_c])
  return(cumprod(1-h_t_fixed))
}

#' Hazard minimizer instead of maximizer
f_surv_hazardmin <- function(s_t, reference_s_t) {
  th_c      <- length(s_t)
  h_t_fixed <- pmin(f_surv_get_q_t(s_t)[1:th_c],f_surv_get_q_t(reference_s_t)[1:th_c])
  return(cumprod(1-h_t_fixed))
}



#' Function to apply general population mortality adjustment to all OS lines
f_surv_gpopadjust <- function(st,gpop, method = "hazardmax", verbose=FALSE) {
  
  # Cycle through all PLMT's only for endpoint_0, applying hazard max with
  # the corresponding gpop OS line
  npopu        <- names(st)
  names(npopu) <- npopu
  lapply(npopu, function(this_npopu) {
    
    nline        <- names(st[[this_npopu]])
    names(nline) <- nline
    
    lapply(nline, function(this_nline) {
      
      nmol        <- names(st[[this_npopu]][[this_nline]])
      names(nmol) <- nmol
      
      lapply(nmol, function(this_nmol) {
        
        ntrial <- names(st[[this_npopu]][[this_nline]][[this_nmol]])
        names(ntrial) <- ntrial
        
        lapply(ntrial, function(this_ntrial) {
          
          nendp <- names(st[[this_npopu]][[this_nline]][[this_nmol]][[this_ntrial]])
          names(nendp) <- nendp
          
          lapply(nendp, function(this_nendp) {
            
            plmte <- st[[this_npopu]][[this_nline]][[this_nmol]][[this_ntrial]][[this_nendp]]
            
            # If there's no survival extrapolation here, then there's nothing to do
            if (!"st" %in% names(plmte)) {
              return(plmte)
            } else {
              # Basically, if it's empty already, then there's nothing to do!
              if (is.null(plmte$dest)) {
                if (verbose) cat(paste0(
                  "This PLMTE has no destination, doing nothing. (plmte = ",
                  paste(c(this_npopu,this_nline,this_nmol,this_ntrial,this_nendp), collapse = "$"),
                  ")\n"
                ))
                return(plmte)
              }
              
              d  <- plmte$dest
              if (verbose) cat(paste0("Gpop adj: ",paste(unlist(d),collapse = " | "),"\n"))
              gp <- gpop[[d$pop]][[d$line]]$os
              if (method == "hazardmax") {
                plmte$st <- f_surv_hazardmax(plmte$st, gp)
              } else {
                plmte$st <- pmin(plmte$st, gp)
              }
              
              if (length(plmte$st) == 1) {
                plmte$st <- NULL
                return(plmte)
              } else {
                plmte$gpop_adjusted <- TRUE
                plmte$gpop_method   <- method
                return(plmte)
              }
            }
          })
        })
      })
    })
  })
}

#' hazard max/abs min for PFS vs OS to prevent curve crossing
f_surv_PFSxOS <- function(st,method = "hazardmax") {
  
  # Cycle through all PLMT's only for endpoint_0, applying hazard max with
  # the corresponding gpop OS line
  
  Reduce(
    x = 1,
    init = st,
    accumulate = FALSE,
    f = function(prev, dos) {
      
      lapply(st, function(popu) {
        lapply(popu, function(li) {
          lapply(li, function(mol) {
            lapply(mol, function(plmt) {
              
              OS  <- plmt$endpoint_0
              PFS <- plmt$endpoint_1
              
              # If there's no extrapolations to adjust then don't.
              if (!"st" %in% names(OS)) {
                return(plmt)
              }
              
              if (method == "hazardmax") {
                PFS$st <- f_surv_hazardmax(PFS$st, OS$st)
              } else {
                PFS$st <- pmin(PFS$st, OS$st)
              }
              PFS$curve_cross        <- TRUE
              PFS$curve_cross_method <- method
              plmt$endpoint_1        <- PFS
              
              # Now that we updated the OS for this PLMT, we can just return it
              return(plmt)
            })
          })
        })
      })
    }
  )
}


#' hazard min/abs max for TTD vs PFS to prevent curve crossing. TTP should be 
#' the same as or above PFS at all times, since it censors for death and PFS
#' includes death as an event. It should be impossible for TTP to be below PFS
#' 
#' Also note that as TTP (and TTD) ignore death events, they can be above OS
#' 
f_surv_PFSxTTP <- function(st, method = "hazardmax") {
  # Cycle through all PLMT's only for endpoint_0, applying hazard max with
  # the corresponding gpop OS line
  
  lapply(st, function(popu) {
    lapply(popu, function(li) {
      lapply(li, function(mol) {
        lapply(mol, function(plmt) {
          PFS <- .subset2(plmt, "endpoint_1")
          TTP <- .subset2(plmt, "endpoint_3")
          if (!"st" %in% names(PFS)) {
            return(plmt)
          }
          if (method == "hazardmax") {
            TTP$st <- f_surv_hazardmin(PFS$st, TTP$st)
          } else if (method == "abs") {
            TTP$st <- pmax(PFS$st, TTP$st)
          } else {
            stop("the argument 'method' should be 'hazards' or 'abs'")
          }
          PFS$curve_cross        <- TRUE
          PFS$curve_cross_method <- method
          plmt$endpoint_3        <- TTP
          
          # Now that we updated the OS for this PLMT, we can just return it
          return(plmt)
        })
      })
    })
  })
}

#' Although time to discontinuation can have lower hazard than PFS, it can't
#' have lower hazard than OS, since you can't be dead but on treatment.
f_surv_TTDxOS <- function(st, method = "hazardmax") {
  
  # Cycle through all PLMT's only for endpoint_0, applying hazard max with
  # the corresponding gpop OS line
  
  return(Reduce(
    x = 1,
    init = st,
    accumulate = FALSE,
    f = function(prev, dos) {
      
      lapply(st, function(popu) {
        lapply(popu, function(li) {
          lapply(li, function(mol) {
            lapply(mol, function(plmt) {
              
              OS <- .subset2(plmt,"endpoint_0")
              TTD <- .subset2(plmt,"endpoint_2")
              
              if (!"st" %in% names(OS)) {
                return(plmt)
              }
              
              if (method == "hazardmax") {
                TTD$st <- f_surv_hazardmax(TTD$st,OS$st)
              } else if (method == "abs") {
                TTD$st <- pmin(TTD$st,OS$st)
              } else {
               stop("the argument 'method' should be 'hazards' or 'abs'") 
              }
              TTD$curve_cross        <- TRUE
              TTD$curve_cross_method <- method
              plmt$endpoint_2        <- TTD
              
              # Now that we updated the OS for this PLMT, we can just return it
              return(plmt)
            })
          })
        })
      })
    }
  ))
}

f_surv_adjuvant_HR <- function(st,adjuvant_impact,demo_table, lookup, verbose = FALSE) {
  
  # Cycle through all PLMT's only for endpoint_0, applying hazard ratio associated with
  # adjuvant treatment
  
  adjuvant_impact <- as.data.table(adjuvant_impact)
  
  # Make a multi-id, which can id entries inside of st in two different ways:
  # 
  #  - by number for risk population (0, 1, 2)
  #  - by "label" for risk population (pop_0, pop_1, pop_2)
  #  
  # The reason for doing this is that the table adjuvant_impact has the numbers
  # but the output list needs to have the pop_0 type labelling to match the rest
  # of the PLMTE relational database type system we've set up.
  # 
  pop_names <- names(st)
  pop_id <- as.numeric(gsub("pop_","",names(st)))
  names(pop_id) <- names(st)
  
  
  # Cycle through pop_id, using the numbers or the labels to pull the right
  # set of data:
 
  lapply(pop_id, function(popu) {
    
    # simplify things by pulling out the label, position and number as separate
    # things to keep everything simple in the code that follows.
    p_pos <- which(pop_id %in% popu)
    p_lab <- names(pop_id)[p_pos]
    p_num <- as.numeric(popu)
    
    # Use the lookup table to also pull up the name of the risk population
    # as used within Excel:
    p_nam <- lookup$ipd$pop[Number == p_num,]$Description
    
    # now, do the same thing by treatment line, we only need numbers and 
    # labels this time though:
    line_names <- names(st[[p_lab]])
    line_id    <- as.numeric(gsub("line_", "", line_names))
    names(line_id) <- line_names
    
    # Cycle through the line ids. We still need numbers to filter the demo table and
    # ids to name the list:
    lapply(line_id, function(li) {
      
      # Number for this line. we need the label for getting the names in the next 
      # level:
      l_num <- as.numeric(li)
      l_pos <- which(line_id %in% l_num)
      l_lab <- names(line_id)[l_pos]
      
      # Note that the demo table in excel uses p_nam and NOT the number! 
      prop_adj <- demo_table[Treatment.line == l_num & Population == p_nam]$Prior.IO...in.12.months.Mean
      
      mol_names <- names(st[[p_lab]][[l_lab]])
      mol_id    <- as.numeric(gsub("mol_", "", mol_names))
      names(mol_id) <- mol_names
      
      lapply(mol_id, function(mol) {
        
        # now number position label and name for molecule: 
        m_num <- as.numeric(gsub("mol_", "", mol))
        m_pos <- which(mol_id %in% m_num)
        m_lab <- names(mol_id)[m_pos]
        
        # The table in excel named range R_table_prior_IO_impact_eff uses 
        # the treatment name rather than molecule number, so we need to lookup
        # the name
        m_nam <- lookup$ipd$mol[Number == m_num,]$RCC_input_desc
        
        t_labs <- names(st[[p_lab]][[l_lab]][[m_lab]])
        t_nums <- as.numeric(gsub("trial_", "", t_labs))
        
        # Make the cycling index for trials. the function structure
        # lets you make the values and the names in one go.
        trial_id <- structure(t_nums,.Names=t_labs)
        
        lapply(trial_id, function(tr) {
          
          # Make an id set for the trial being considered:
          t_num <- as.numeric(gsub("trial_", "", tr))
          t_pos <- which(trial_id %in% t_num)
          t_lab <- names(trial_id)[t_pos]
          
          # calculate the impact of prior adjuvant therapy (HR assumed 1 for patients with no prior adjuvant)
          impact_adj <- prop_adj * adjuvant_impact[Treatments == m_nam,]$Prior.adj.impact_HR_mean + (1 - prop_adj)
          
          # The impact is applied like a hazard ratio to the element "st" inside
          # of each PLMTE. Therefore, we need to cycle through each endpoint
          # one at a time, checking whether element st even exists, then applying
          # the hazard ratio to it if it does.
          
          # cycle through this plmt, one endpoint at a time, applying HR if
          # it's not 1 and s(t) exists and is finite.
          plmt <- lapply(st[[p_lab]][[l_lab]][[m_lab]][[t_lab]], function(plmte) {
            
            if ("st" %in% names(plmte)) {
              if (all(is.finite(plmte$st), is.finite(impact_adj))) { #, impact_adj != 1
                
                if (verbose) {
                  f_misc_colcat(
                    paste0(
                      "Applying adjuvant HRs | ",
                      paste(p_lab,l_lab,m_lab,t_lab),
                      " | HR=",
                      impact_adj
                    )
                  )
                }
                
                # If the entry st (survival at time t) exists in this PLMTE and
                # it is finite (i.e. not NA or something) then apply the HR to it
                # otherwise don't. Obviously if the HR is 1 then computation is 
                # unecessary.
                plmte$st <- plmte$st ^ impact_adj
                
                # Now that we're done editing this PLMTE, we return the updated
                # version:
                return(plmte)
              }
            } else {
              # there's nothing to do because there's no entry in this PLMTE
              # for st (survival at time t), so there's nothing to apply a HR to
              return(plmte)
            }
          })
          # Now that we've edited this plmt, we can return the updated result
          return(plmt)
        })
      })
    })
  })
  
}


# Stopping rules ----------------------------------------------------------

#' Function to apply stopping rules by directly modifying TTD to make it 0 at 
#' the start of the given cycle (e.g. if 104 weeks given, it should be 0 from the
#' 105th week onwards, i.e. patients get treated in week 103 starting from 0, but
#' then not from 104 onwards)
#' 
f_surv_apply_stopping_rules <- function(st, tot_table, lookups) {
  
  # filter down the table to just those with a max TxDur in cycles that have Yes for inclusion:
  txDur <- tot_table[Treatment.given.for.fixed.duration...Yes.No. == "Yes" & If.yes..number.of.cycles != 0 ,]
  
  # Cycle down the table using the Reduce method, only changing TTD in corresponding
  # DESTINATION PLMTs. We don't have TOT lines, so TTD is the only thing we can 
  # manipulate here. However, it may be that we have to create TOT for the PS model
  # applying these. 
  
  Reduce(
    x = 1:nrow(txDur),
    init = st,
    accumulate = FALSE,
    f = function(prev, txDur_row) {
      
      dat <- as.list(txDur[txDur_row,])
      
      # Round the max treatment duration up to the next cycle. In reality these
      # should all be integer values for. The +1 is for cycle 0 as it will stop
      # treatment at the start of the cycle given.
      max_txDur <- ceiling(dat$If.yes..number.of.cycles) + 1
      
      # the destination is by population line and molecule only across all trials
      # and endpoints. Therefore, PLM rather than PLMTE
      d <- list(
        pop  = paste0("pop_",lookups$ipd$pop[Description == dat$Population,"Number"]),
        line = paste0("line_",dat$Treatment.line),
        mol  = paste0("mol_",dat$Molecule)
      )
      
      # Cycle through all relevant trials for the plm, applying the stopping rule to the 
      # TTD line in each by setting all TTD values at or beyond the stopping rule
      # time to 0
      
      plm_updated <- lapply(prev[[d$pop]][[d$line]][[d$mol]], function(plmt) {
        
        # For this plmt, check if TTD is populated or not. if not, then just return it
        # unchanged, if so, apply the stopping rule and return the result.
        
        if ("populated" %in% names(plmt$endpoint_2)) {
          
          # Apply the stopping rule at the nth cycle (taking 0 into account)
          plmt$endpoint_2$st[max_txDur:length(plmt$endpoint_2$st)] <- 0
          plmt$endpoint_2$stopping_rule    <- TRUE
          plmt$endpoint_2$stopping_rule_cl <- max_txDur
          
          return(plmt)
        } else {
          return(plmt)
        }
      })
      
      # Now simply inject the updated plm into prev, and return the updated list,
      # allowing us to then apply the next row of the Excel table in the next 
      # go round.
      prev[[d$pop]][[d$line]][[d$mol]] <- plm_updated
      return(prev)
    }
  )
  
}



# QC functions ------------------------------------------------------------

# plots to look at abs survival and estimated hazard

f_qc_surv_ExtrapPlot <- function(st,popu,li,mo,tr,t_yr,th) {
  
  endpoints <- structure(
    .Data = c("endpoint_0", "endpoint_1", "endpoint_2", "endpoint_3"),
    .Names  = c("OS", "PFS", "TTD", "TTP")
  )
  
  p_dat <- rbindlist(lapply(1:length(endpoints), function(endp) {
    
    if (length(st[[popu]][[li]][[mo]][[tr]][[endpoints[endp]]]$st[1]) > 0) {
      data.table(t_yr = t_yr,
                 s_t  = st[[popu]][[li]][[mo]][[tr]][[endpoints[endp]]]$st,
                 endp = names(endpoints[endp]))
    }
  }))
 
  if(is.null(p_dat)) {
    return(NULL)
  } else {
    ggplot(p_dat, aes(x = t_yr, y = s_t, colour = endp)) + 
      geom_line() + 
      theme_classic() +
      theme(legend.position = "bottom", legend.title=element_blank()) + 
      labs(title = NULL, x = "Time (years)", y = "% Survival") + 
      scale_x_continuous(expand = expansion(mult = c(0,0.05))) + 
      scale_y_continuous(labels = scales::percent)
  }
}



f_qc_surv_gethead <- function(x, p,l,m,t,len = 10) {
  do.call(cbind,lapply(x[[paste0("pop_",p)]][[paste0("line_",l)]][[paste0("mol_",m)]][[paste0("trial_",t)]], function(endp) {
    head(endp$st,len)
  })
  )
}

f_qc_lookupplm <- function(tab, p,l,m) {
  tab[Population == p & Treatment.line == l & Molecule == m, ]
}
f_qc_lookupplm_e <- function(tab, p,l,m,e) {
  tab[Population == p & Treatment.line == l & Molecule == m & End.point == e, ]
}


f_qc_surv_EstHazPlot <- function(st,gpop,popu,li,mo,tr,t_yr,th) {
  endpoints <- structure(
    .Data = c("endpoint_0", "endpoint_1", "endpoint_2", "endpoint_3"),
    .Names  = c("OS", "PFS", "TTD", "TTP")
  )
  
  p_dat <- do.call(rbind,lapply(1:length(endpoints), function(endp) {
    if (!is.na(st[[popu]][[li]][[mo]][[tr]][[endpoints[endp]]]$st[1])) {
      data.table(t_yr = t_yr,
                 est_h_t  = f_surv_get_h_t(st[[popu]][[li]][[mo]][[tr]][[endpoints[endp]]]$st, t_yr[2]),
                 endp = names(endpoints[endp]))
    }
  }))
  
  p_dat <- rbind(p_dat,data.table(t_yr = t_yr,est_h_t = f_surv_get_h_t(gpop[[popu]][[li]]$os, t_yr[2]),endp="gpop"))
  if(nrow(p_dat) == length(t_yr)) {
    return(NULL)
  } else {
    ggplot(p_dat, aes(x = t_yr, y = est_h_t, colour = endp)) + 
      geom_line() + 
      theme_classic() +
      theme(legend.position = "bottom", legend.title=element_blank()) + 
      labs(title = NULL, x = "Time (years)", y = "Transition probability") + 
      scale_x_continuous(expand = expansion(mult = c(0,0.05)))
  }
}



