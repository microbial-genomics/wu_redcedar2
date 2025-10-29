### combined parameter plots
heatmap_pcc
pcc_plot_maintext
parameter_evolution_plot

stacked_sensitivity_posteriors_plot <- plot_grid(
                                      pcc_plot_maintext, 
                                      parameter_evolution_plot, 
                                      ncol = 1, labels = c("A", "B"))

#library(cowplot)
library(patchwork)
stacked_sensitivity_posteriors_plot <- heatmap_pcc | (pcc_plot_maintext / parameter_evolution_plot)

stacked_sensitivity_posteriors_plot <- stacked_sensitivity_posteriors_plot + 
                                  plot_annotation(tag_levels = "A") +
                                  plot_layout(widths = c(1, 4))
stacked_sensitivity_posteriors_plot


png(paste(file.path(graphics_dir, "stacked_sensitivity_posteriors_plot_maintext.png")), 
    width=10, height=9, units="in", res=300)
  stacked_sensitivity_posteriors_plot
dev.off()


### ribbon and scatter plots
stacked_daily_ribbon_plot_subset
combined_rain_flow_conc_obs_sims_plot

combined_ribbon_scatter_plots <- plot_grid(stacked_daily_ribbon_plot_subset,
                                                   combined_rain_flow_conc_obs_sims_plot,
                                                   rel_heights = c(2,1),
                                                   ncol = 1)
combined_ribbon_scatter_plots

ggsave(paste(file.path(graphics_dir, "combined_ribbon_scatter_plots.png")), 
       plot = combined_ribbon_scatter_plots, 
       width=12, height=12, units="in", dpi=300)



# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")