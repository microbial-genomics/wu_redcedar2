pcc_plot_maintext
parameter_evolution_plot

stacked_sensitivity_posteriors_plot <- plot_grid(
                                      pcc_plot_maintext, 
                                      parameter_evolution_plot, 
                                      ncol = 1, labels = c("A", "B"))

stacked_sensitivity_posteriors_plot 

png(paste(file.path(graphics_dir, "stacked_sensitivity_posteriors_plot_maintext.png")), 
    width=10, height=9, units="in", res=300)
  stacked_sensitivity_posteriors_plot
dev.off()





stacked_daily_ribbon_plot

# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")