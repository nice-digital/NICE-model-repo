### This code produces CEAFs with the cost axes correctly outputting following HPC run issues with encoding

# Set population

popu <- "pop_3"

# Set nsdom and reduced_table for the population

nsddom <- res$weighted_incremental[[popu]]$not_strictly_dominated
reduced_table <- res$weighted_incremental[[popu]]$expanded_results[extdom == FALSE,]

# run plot

ggplot(nsddom, aes(x = qalys, y = costs, colour = as.factor(L1), label = r)) + 
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

