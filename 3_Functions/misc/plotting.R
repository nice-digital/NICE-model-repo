#' Function to make a state residency plot - useful for the PS model and checking that
#' states sum to 1. The integral of each line is lifetime expected time in state `:)`
#' 
#' @param sr_list a list of s_t per endpoint to be included in the plt
#' @param t_yr time in years per cycle
#' 
f_plot_srPlot <- function(sr_list, t_yr, plot_x_lim) {
  p_dat <- data.table::rbindlist(lapply(structure(names(sr_list)[1:(length(names(sr_list))-1)],.Names=names(sr_list)[1:(length(names(sr_list))-1)]), function(endpoint) {
    data.table(
      t_yr = t_yr,
      s_t  = .subset2(sr_list,endpoint),
      endp = endpoint
    )
  }))[t_yr <= plot_x_lim,]
  ggplot(p_dat, aes(x = t_yr, y = s_t, colour=endp)) + 
    geom_line() + 
    theme(legend.position = "bottom") +
    theme_classic() + 
    scale_x_continuous(expand = expansion(mult = c(0,0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05)))
}



# Markov trace plots ------------------------------------------------------

#' Function to draw a consolidated trace as a plot
#' 
#' @param consol_trace consolidated trace produced during ST model (added up columns)
#' @param treatment_names names to allocate to the legend for each line
#' @tmax time max. Function assumes t is in years with weekly cycle length.
#' 
#' 
f_plot_mk_draw_consol_trace <- function(consol_trace,treatment_names , tmax = 15) {
  plot_dt   <- data.table(consol_trace)
  th <- nrow(plot_dt)
  plot_dt$w <- 0:(th-1)
  plot_dt$y <- plot_dt$w / 52.17857
  plot_dt$w <- NULL
  plot_dt   <- melt.data.table(plot_dt, id.vars = "y")
  
  trt_seq_plotLab <- c(paste0("L",1:(length(treatment_names)-1)),"NT")
  names(trt_seq_plotLab) <- trt_seq_plotLab
  
  plot_dt   <- Reduce(
    x = 1:length(treatment_names),
    init = plot_dt,
    accumulate = FALSE,
    f = function(prev, trt_n) {
      
      if (trt_n < length(treatment_names)) {
        which_on  <- which(prev$variable == paste0(names(trt_seq_plotLab)[trt_n],"_on"))
        which_off <- which(prev$variable == paste0(names(trt_seq_plotLab)[trt_n],"_off"))
        prev[which_on,]$variable  <- paste0(treatment_names[trt_n]," (on treatment)")
        prev[which_off,]$variable <- paste0(treatment_names[trt_n]," (off treatment)")
      } else {
        # this is BSC
        which_ones  <- which(prev$variable == "BSC")
        prev[which_ones,]$variable  <- treatment_names[trt_n]
      }
      
      return(prev)
    }
  )
  
  # We really did it! A fully working trace calculator for a sequencing model in R :)
  return(ggplot(plot_dt[y < tmax,],aes(x = y, y = value, colour = variable)) + 
           geom_line() + 
           theme_classic() + 
           theme(legend.position = "bottom") + 
           scale_x_continuous(expand = expansion(mult = c(0,0.05))) +
           scale_y_continuous(labels = scales::percent,expand = expansion(mult = c(0,0.05))) + 
           labs(x = "Years from baseline", y = "% of baseline cohort") +
           guides(
             colour = guide_legend(title = NULL)
           ))
}
