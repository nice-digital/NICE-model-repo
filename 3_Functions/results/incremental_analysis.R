#' Function to produce incremental analysis results
#' 
#' @param res_tab the results from the markov model (e.g. res$mk$disc$pop_3$res)
#' @param produce_plot whether or not to produce an efficiency frontier plot
#' 
#' 
#' @details The plot has a lot of dots on it so sequences are numbered. the numbers
#' correspond to the column `r` in the expanded_results table. one can take the
#' plot and change the labelling or produce a separate plot from the various tables
#' returned.
#' 
#' 
f_res_mk_incremental <- function(res_d, res_ud,lu_mol, produce_plot = TRUE, no_active_lines = NULL, output_weighted = "No") {
  
  out <- list()
  
  # Sort by cost in descending order:
  if ("trt" %in% colnames(res_ud)) {
    ru <- res_ud[,-c("costs","qalys","trt")]
  } else {
    ru <- res_ud[,-c("costs","qalys")]
  }
  r <- merge.data.table(res_d,ru)
  r <- r[order(costs),]
  
  # Translate molecule number to full name to make the output more readable.
  # Only do this when we're looking at 1st line treatments (weighted models):
  if ("L1" %in% colnames(r)) {
    r$L1 <- lu_mol[match(r$L1, lu_mol$Number)]$Description
  } else {
    r$L1 <- gsub("→","|",r$trt_n)
  }
  
  # Filter according to the number of active lines
  if(output_weighted == "No") {
    if (!is.null(no_active_lines)) {
      r <- r[str_count(r$trt_n, pattern = "→")<(no_active_lines+1)]
    }
  }
  
  # Identify strictly dominated strategies, i.e. for each row are there any
  # other rows which have both lower costs and higher QALYs
  r$str_dom <- FALSE
  
  # If for any treatment a worse treatment is more expensive, then it's dominated strictly
  
  str_dom_list <- lapply(nrow(r):1, function(tx) {
    
    c <- r$costs[tx]
    q <- r$qalys[tx]
    
    if(tx == nrow(r)) {
      return(q < max(r$qalys))
    } else {
      dat <- r[-tx,]
      
      dat$t_c <- dat$costs <= c
      dat$t_q <- dat$qalys >= q
      
      dat$tt <- dat$t_c & dat$t_q
      
      return(any(dat$tt))
    }
  })
  
  r$str_dom <- unlist(str_dom_list)[nrow(r):1]
  
  
  if(sum(!r$str_dom) == 1) {
    r <- r[str_dom == FALSE,]
    warning(paste0("Molecule ", r$L1, " as first line therapy strictly dominates all other strategies..."))
    return(r)
  }
  
  # not strictly dominated strategies
  stdom  <- r[str_dom == TRUE,]
  nsddom <- r[str_dom == FALSE,]
  
  
  # Extended dominance is iterative by nature, requiring pairwise ICERs to be 
  # calcualted repeatedly after eliminating strategies one at a time:
  
  nsddom$extdom  <- FALSE
  nsddom$ic      <- 0
  nsddom$iq      <- 0
  nsddom$il      <- 0
  nsddom$ICER    <- 0
  nsddom$ICER[length(nsddom$ICER)]   <- NA
  nsddom <- nsddom[order(qalys),]
  nsddom$r <- 1:nrow(nsddom)
  
  out$not_strictly_dominated <- nsddom
  
  # Order checks for QC - they all go in increasing order in cost and qalys (otherwise
  # some would be strdom by a cheaper and better strat)
  
  extdom <- Reduce(
    x = 1:10,
    accumulate = FALSE,
    init = nsddom,
    function(not_strictly_dominated, dos) {
      
      
      after_this_pass <- Reduce(
        x = 1:(nrow(not_strictly_dominated)-1),
        init = not_strictly_dominated,
        accumulate = FALSE,
        f = function(prev, comp_tx) {
          
          
          
          # For each strategy, calculate the lowest ICER compared to any other
          # which is more effective than it is, which should be the rows below it.
          # this can be done in reverse order cycling through interventions from the
          # bottom of the table upwards, calculating the lowest ICER to any other strategy
          # that at the time isn't dominated.
          
          # If this strategy is extendedly dominated, just return the table as we 
          # don't need to do anything further
          if (prev$extdom[comp_tx] == TRUE) return(prev)
          
          ned <- prev[extdom == FALSE & r %in% comp_tx:nrow(prev),]
          
          # Pull out the comparator:
          comp <- prev[r == comp_tx,]
          
          # Calculate pairwise ICERs against all cheaper and less effective treatments
          # compared to this one, that are not themselves dominated:
          ic    <- ned[r > comp_tx]$costs - comp$costs
          iq    <- ned[r > comp_tx]$qalys - comp$qalys
          il    <- ned[r > comp_tx]$ly - comp$ly
          
          ICERs <- ic/iq
          
          ned$ic <- c(NA,ic)
          ned$iq <- c(NA,iq)
          ned$il <- c(NA,il)
          ned$ICER <- c(NA,ICERs)
          which_best_ICER <- ned[which.min(ICER),]$r
          
          # Locate the best ICER (efficiency frontier) for this iteration
          ned$best_ICER <- FALSE 
          ned[r == which_best_ICER,]$best_ICER <- TRUE
          
          # Figure out if there is any extended dominance to add to the table
          # prev$testboi <- between(prev$r,which_best_ICER,int_tx-1,incbounds = FALSE)
          if (which_best_ICER > (comp_tx+1)) {
            prev$extdom[data.table::between(prev$r,comp_tx,which_best_ICER,incbounds = FALSE)] <- TRUE
            prev$ic[data.table::between(prev$r,comp_tx,which_best_ICER,incbounds = FALSE)] <- NA
            prev$iq[data.table::between(prev$r,comp_tx,which_best_ICER,incbounds = FALSE)] <- NA
            prev$il[data.table::between(prev$r,comp_tx,which_best_ICER,incbounds = FALSE)] <- NA
            prev$ICER[data.table::between(prev$r,comp_tx,which_best_ICER,incbounds = FALSE)] <- NA
            if (which_best_ICER + 1 < nrow(prev)) {
              prev$ICER[(which_best_ICER+1):nrow(prev)] <- NA
            }
            prev[r == which_best_ICER,]$ic <- ned[best_ICER == TRUE,]$ic
            prev[r == which_best_ICER,]$iq <- ned[best_ICER == TRUE,]$iq
            prev[r == which_best_ICER,]$il <- ned[best_ICER == TRUE,]$il
            prev[r == which_best_ICER,]$ICER <- ned[best_ICER == TRUE,]$ICER
          } else {
            # it's just the next treatment along in the reduced table. there's no
            # extended dominance to add
            prev[r == which_best_ICER,]$ic   <- ned[r == which_best_ICER,]$ic
            prev[r == which_best_ICER,]$iq   <- ned[r == which_best_ICER,]$iq
            prev[r == which_best_ICER,]$il   <- ned[r == which_best_ICER,]$il
            prev[r == which_best_ICER,]$ICER <- ned[r == which_best_ICER,]$ICER
          }
          return(prev)
        })
      
      return(after_this_pass)
    }
  )
  
  out$expanded_results <- extdom
  
  # Non-dominated strategies and ICERs between them
  reduced_table <- extdom[extdom == FALSE,]
  
  out$non_dominated <- reduced_table[,-c("str_dom","extdom","r")]
  
  if (produce_plot) {
    out$p <- ggplot(nsddom, aes(x = qalys, y = costs, colour = as.factor(L1), label = r)) + 
      geom_point() +
      theme_classic() +
      # ggrepel::geom_text_repel(max.overlaps = 100, alpha = 0.2) +
      geom_line(data = reduced_table, aes(x=qalys,y=costs,colour=NULL)) +
      # ggrepel::geom_label_repel(
      #   data = reduced_table,
      #   # arrow = arrow(ends = "last",type = "closed"),
      #   aes(
      #     x = qalys,
      #     y = costs,
      #     colour = NULL,
      #     label = as.factor(L1),
      #   ))+
      theme(legend.position = "bottom") + 
      scale_x_continuous(limits = c(0,max(nsddom$qalys)), expand = expansion(mult = c(0,0.05))) +
      scale_y_continuous(
        limits = c(0, max(nsddom$costs)),
        expand = expansion(mult = c(0, 0.05)),
        labels = label_dollar(prefix = "£")
      ) +
      labs(x= "QALYs", y = "Costs") +
      theme(legend.title = element_blank())
  }
  
  return(out)
}





# testing ground ----------------------------------------------------------

if (FALSE) {
  f_res_mk_incremental(res_tab = res$mk$disc$pop_1$res)
  f_res_mk_incremental(res$mk$disc$pop_2$res)
  f_res_mk_incremental(res$mk$disc$pop_3$res)
  f_res_mk_incremental(res$mk$disc$pop_4$res)
  f_res_mk_incremental(res$mk$disc$pop_5$res)
  f_res_mk_incremental(res$mk$disc$pop_6$res)
}

