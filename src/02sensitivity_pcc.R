# source directory for hpc pcc output
hpc_data_sensitivity

# daily pcc results for 45 parameters (hpc output does not have colnames)
load(file.path(hpc_data_sensitivity, "bac_pcc.RData"))
dim(bac_pcc)
# View(bac_pcc)
load(file.path(hpc_data_sensitivity, "flows_pcc.RData"))
dim(flows_pcc)
load(file.path(hpc_data_sensitivity, "flux_pcc.RData"))
dim(flux_pcc)

# sensitivity analysis inputs and evaluation
load(file.path(hpc_data_sensitivity, "bac_sensitivity.RData"))
str(bac_cal1)
sim_parameters <- bac_cal1$parameter$value
dim(sim_parameters)
#View(sim_parameters)

#fix problematic colnames
colnames(sim_parameters)[2] <- "SOL_K"
colnames(sim_parameters)[3] <- "SOL_AWC"
full_parameter_list <- colnames(sim_parameters)

### resolving parameter names
# sensitivity output list
### resolving parameter names
colnames(sim_parameters)
# par_bound 
colnames(par_bound)
#View(cbind(colnames(sim_parameters), colnames(par_bound)))
# the parameters that are carried forward from the sensitivity analysis
simulated_parameter_list
sensitive_parameter_indices <- match(simulated_parameter_list, colnames(sim_parameters))
cbind(simulated_parameter_list, colnames(sim_parameters)[sensitive_parameter_indices])

# 
sim_bac_concs <- bac_cal1$simulation$bac_out
dim(sim_bac_concs)
#View(sim_bac_concs)
colnames(sim_bac_concs)

# 
sim_flows <- bac_cal1$simulation$q_out 
dim(sim_flows)
#View(sim_flows)
colnames(sim_flows)

#simulation is for the 365 days of 2013
#drop date column for sim_bac_concs
sim_bac_concs2 <- sim_bac_concs[,-c(1)]
sim_dates <- sim_bac_concs[,1]
#View(sim_bac_concs2)
dim(sim_bac_concs2)
#transpose bacteria concentration data frame
sim_bac_concs3 <- t(sim_bac_concs2)
dim(sim_bac_concs3)
# pcc on row medians
sim_bac_medians <- rowMedians(sim_bac_concs3)
pcc(sim_parameters, sim_bac_medians)
# pcc on row maxs
sim_bac_maxs <- rowMaxs(sim_bac_concs3)
pcc(sim_parameters, sim_bac_maxs)


# drop dates for sim_flows
sim_flows2 <- sim_flows[,-c(1)]
dim(sim_flows2)
#transpose flows data frame
sim_flows3 <- t(sim_flows2)
dim(sim_flows3)
# pcc on row medians
sim_flows_medians <- rowMedians(sim_flows3)
pcc(sim_parameters, sim_flows_medians)
# pcc on row maxs
sim_flows_maxs <- rowMaxs(sim_flows3)
pcc(sim_parameters, sim_flows_maxs)

# flux
sim_flux3 <- sim_bac_concs3*sim_flows3*10^4
dim(sim_flux3)
#pcc on row medians
sim_flux_medians <- rowMedians(sim_flux3)
pcc(sim_parameters, sim_flux_medians)
#pcc on row maxs
sim_flux_maxs <- rowMaxs(sim_flux3)
pcc(sim_parameters, sim_flux_maxs)



dim(bac_pcc)
colnames(bac_pcc) <- colnames(sim_parameters)
bac_pcc2 <- as.data.frame(bac_pcc)
colnames(flows_pcc) <- colnames(sim_parameters)
flows_pcc2 <- as.data.frame(flows_pcc)
colnames(flux_pcc) <- colnames(sim_parameters)
flux_pcc2 <- as.data.frame(flux_pcc)

# Add a dataset label and bind them together
bac_pcc2$Dataset <- "Bacteria PCC"
flows_pcc2$Dataset <- "Flow PCC"
flux_pcc2$Dataset <- "Flux PCC"
combined_pcc <- bind_rows(bac_pcc2, flows_pcc2, flux_pcc2)

# Convert to long format
long_combined_pcc <- pivot_longer(combined_pcc, 
                              cols = -Dataset, 
                              names_to = "Parameter", 
                              values_to = "Value")

### calculate mean absolute value pccs for each of the 45 initial parameters
#View(combined_pcc)
colnames(combined_pcc)
dim(combined_pcc) #11595 46
#11595/3 = 3865
unique(combined_pcc$Dataset)
bacteria_combined_pcc <- combined_pcc[1:3865,1:45]
dim(bacteria_combined_pcc)
bacteria_mav_pcc <- sapply(bacteria_combined_pcc, function(x) mean(abs(x)))
bacteria_mav_pcc

flow_combined_pcc <- combined_pcc[3866:7730,1:45]
dim(flow_combined_pcc)
flow_mav_pcc <- sapply(flow_combined_pcc, function(x) mean(abs(x)))

flux_combined_pcc <- combined_pcc[7731:11595,1:45]
dim(flux_combined_pcc)
flux_mav_pcc <- sapply(flux_combined_pcc, function(x) mean(abs(x)))

combined_mav_pcc <- as.data.frame(bind_rows(bacteria_mav_pcc, flow_mav_pcc, flux_mav_pcc))
rownames(combined_mav_pcc) <- c("Bacteria_MAV_PCC", "Flow_MAV_PCC", "Flux_MAV_PCC")
combined_mav_pcc


library(DT)
library(htmlwidgets)
library(webshot2)
# Calculate quantile breaks from all cell values
combined_mav_pcc <- combined_mav_pcc %>% mutate_if(is.numeric, signif, digits = 3)
breaks <- quantile(unlist(combined_mav_pcc), probs = c(0.33, 0.66), na.rm = TRUE)

dt_mav_pcc <- datatable(combined_mav_pcc) %>%
  formatStyle(
    columns = names(combined_mav_pcc),
    backgroundColor = styleInterval(
      breaks,
      c("lightblue", "orange", "red")
    )
  )

saveWidget(dt_mav_pcc, file.path(graphics_dir, "dt_mav_pcc.html"))
webshot(file.path(graphics_dir, "dt_mav_pcc.html"), file.path(graphics_dir, "dt_mav_pcc.png"), vwidth = 4000, vheight = 40)



# which parameters from the sensitivity analysis are kept?
simulated_parameter_list #18
full_parameter_list #45
alphabetical_parameter_list <- sort(full_parameter_list)
alphabetical_parameter_list
alphabetical_sensitive_parameter_positions <- sort(match(simulated_parameter_list, alphabetical_parameter_list))
print(alphabetical_sensitive_parameter_positions)

# color them differently (red) than the ones that are dropped
my_labels <- alphabetical_parameter_list
highlight_indices <- alphabetical_sensitive_parameter_positions
label_colors <- ifelse(seq_along(alphabetical_parameter_list) %in% highlight_indices, "red", "black")
# Create labels with <span> for coloring
my_labels_colored <- paste0(
  "<span style='color:", label_colors, "'>", my_labels, "</span>"
)

# Make the stacked violin plot of all 45 parameters for the appendix
pcc_plot_appendix1 <- ggplot(long_combined_pcc, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_grid(rows = vars(Dataset), scales = "free_y", switch = "y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
pcc_plot_appendix1

png(paste(file.path(graphics_dir, "pcc_violin_plot_appendix_45parameters_dont_use.png")), 
    width=11, height=8, units="in", res=300)
   print(pcc_plot_appendix1)
dev.off()

my_palette <- colorRampPalette(brewer.pal(9, "GnBu"))
colors45 <- my_palette(90)[11:55]
pcc_plot_appendix2 <- ggplot(long_combined_pcc, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_grid(rows = vars(Dataset), scales = "free_y", switch = "y") +
  scale_fill_manual(values = colors45) +
  stat_summary(fun = median, geom = "point", color = "red", size = 1) +
  scale_x_discrete(labels = my_labels_colored) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
pcc_plot_appendix2


png(paste(file.path(graphics_dir, "pcc_plot_appendix_45parameters.png")), 
    width=11, height=8, units="in", res=300)
   print(pcc_plot_appendix2)
dev.off()

#####
# filter long_combined_pcc to sensitive parameters for main text figure
colnames(long_combined_pcc)
unique(long_combined_pcc$Parameter)

simulated_parameter_list #18

short_combined_pcc <- long_combined_pcc %>%
  filter(Parameter %in% simulated_parameter_list)
unique(short_combined_pcc$Parameter)

# Make the stacked violin plot of the 18 selected sensitive parameters for the main text
my_palette <- colorRampPalette(brewer.pal(9, "GnBu"))
colors18 <- my_palette(40)[11:28]
pcc_plot_maintext <- ggplot(short_combined_pcc, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_grid(rows = vars(Dataset), scales = "free_y", switch = "y") +
  scale_fill_manual(values = colors18) +
  stat_summary(fun = median, geom = "point", color = "red", size = 1) +
  theme_bw() +
  labs(
    title = "Initial Sensitivity Analyses",
    x = "Parameter",
    y = "Sensitivity"
  ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
pcc_plot_maintext

png(paste(file.path(graphics_dir, "pcc_plot_maintext_18parameters.png")), 
    width=11, height=8, units="in", res=300)
  print(pcc_plot_maintext)
dev.off()





######################################################################
colnames(flux_pcc)
flux_pcc_sensitive <- flux_pcc[,c(1,4,6,7,8,9,10,11,14,17,19,20,22,23,26,27,28,30,40,41)]
colnames(flux_pcc_sensitive)

#vioplot(flux_pcc_sensitive, col="lightblue")

# long for ggplot
long_flux_pcc_sensitive <- as.data.frame(flux_pcc_sensitive) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "pcc")

# sort by median of absolute values of pcc
long_flux_pcc_sensitive <- long_flux_pcc_sensitive %>%
  mutate(parameter = fct_reorder(parameter, abs(pcc), .fun = median))

ggplot(long_flux_pcc_sensitive, aes(x = parameter, y = pcc)) +
  geom_violin(trim = FALSE, color = "black", fill="lightblue", alpha = 0.7, linewidth = 0.8, width=1.1) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  labs(title = "Violin Plot of Sensitivity for Selected Parameters",
       x = "Parameter",
       y = "PCC Sensitivity") +
  theme_bw()



#create time series
bac_dates <- ts(sim_dates, start=c(2004, 1), end=c(2014, 3865), frequency=3865)


##simple ggplot
#dim(bac_pcc)
#pdf("pcc_ts_5000_45-0.2.pdf",width=11,height=8, onefile=TRUE)
#  for(i in 1:45){
#    data <- data.frame(
#      day = as.Date("2014-01-01") - 0:3864,
#      bac_value = bac_pcc[,i],
#      flows_value = flows_pcc[,i],
#      flux_value = flux_pcc[,i]
#    )
#    p <- ggplot(data, aes(x=day)) +
#      geom_line(aes(y=bac_value), color = "darkred") + 
#      geom_line(aes(y=flows_value), color = "steelblue", linetype="twodash") +
#      geom_line(aes(y=flux_value), color = "orange") +
#      xlab("") +
#      ylab("pcc sensitivity") +
#      ylim(-0.2,0.2)
#    print(p + ggtitle(colnames(bac_pcc)[i]))
#  }
#dev.off()

# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")
gc()
