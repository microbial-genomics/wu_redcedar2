length(sim_date) # 3865

dim(concentration_daily_quantile_results) # [1] 3865    8
dim(flow_daily_quantile_results2) # [1] 3865    8

dim(pcp_obs)

sim_date[1]
pcp_obs$date[1]
concentration_daily_quantile_results[1,]

sim_date[3865]
pcp_obs$date[3865]

# 3 day moving average of precip data
#length()
pcp_obs$precipitation_ma3 <- zoo::rollmean(pcp_obs$precipitation, k = 3, fill = NA, align = "right")

rainfall_response_data <- as.data.frame(cbind(sim_date, pcp_obs$precipitation[1:3865], concentration_daily_quantile_results$four_q, flow_daily_quantile_results2$four_q))
dim(rainfall_response_data)
colnames(rainfall_response_data) <- c("date", "rainfall", "bacteria_conc", "flow")
min(rainfall_response_data$bacteria_conc, na.rm=TRUE)
max(rainfall_response_data$bacteria_conc, na.rm=TRUE)
rainfall_response_data$bacteria_conc[rainfall_response_data$bacteria_conc < 1] <- 1
min(rainfall_response_data$bacteria_conc, na.rm=TRUE)

# create new columns for rainfall_response_data corresponding to 3-day moving average of precip and flow
rainfall_response_data$rainfall_ma3 <- zoo::rollmean(rainfall_response_data$rainfall, k = 3, fill = NA, align = "right")
rainfall_response_data$rainfall_ma3
rainfall_response_data$flow_ma3 <- zoo::rollmean(rainfall_response_data$flow, k = 3, fill = NA, align = "right")
rainfall_response_data$flow_ma3

daily_rain_flow_sims_plot <- ggplot(rainfall_response_data, aes(x = rainfall, y = flow)) +
  geom_point(color = "darkgray")+
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Flow (cfs))") +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "red", size = 1.2)+
  theme_classic()

daily_rain_conc_sims_plot <- ggplot(rainfall_response_data, aes(x = rainfall, y = bacteria_conc)) +
  geom_point(color = "darkgray")+
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "red", size = 1.2) +
  theme_classic()

daily_flow_conc_sims_plot <- ggplot(rainfall_response_data, aes(x = flow, y = bacteria_conc)) +
  geom_point(color = "darkgray")+
  xlab("log10(Flow (cfs))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_x_log10() +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "red", size = 1.2)+
  theme_classic()

daily_rain_flow_sims_plot
daily_rain_conc_sims_plot
daily_flow_conc_sims_plot
rain_flow_conc_plot <- plot_grid(daily_rain_flow_sims_plot, 
                                       daily_rain_conc_sims_plot, 
                                       daily_flow_conc_sims_plot,
                                       ncol = 3, labels = c("A", "B", "C"))

rain_flow_conc_plot
ggsave(paste(file.path(graphics_dir, "daily_sims_rain_flow_conc_plot.png")), 
       plot = rain_flow_conc_plot, 
       width=11, height=8, units="in", dpi=300)


# now with observations
dim(filtered_obs_merged)
head(filtered_obs_merged)
colnames(filtered_obs_merged)

daily_rain_flow_obs_plot <- ggplot(filtered_obs_merged, aes(x = precipitation, y = discharge)) +
  geom_point(color = "black")+
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Flow (cfs))") +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "blue", size = 1.2)+
  theme_classic()

daily_rain_conc_obs_plot <- ggplot(filtered_obs_merged, aes(x = precipitation, y = bacteria)) +
  geom_point(color = "black")+
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "blue", size = 1.2) +
  theme_classic()

daily_flow_conc_obs_plot <- ggplot(filtered_obs_merged, aes(x = discharge, y = bacteria)) +
  geom_point(color = "black")+
  xlab("log10(Flow (cfs))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_x_log10() +
  scale_y_log10()+
  geom_smooth(method = "loess", se = FALSE, color = "blue", size = 1.2)+
  theme_classic()

daily_rain_flow_obs_plot
daily_rain_conc_obs_plot
daily_flow_conc_obs_plot
rain_flow_conc_obs_plot <- plot_grid(daily_rain_flow_obs_plot, 
                                 daily_rain_conc_obs_plot, 
                                 daily_flow_conc_obs_plot,
                                 ncol = 3, labels = c("A", "B", "C"))
rain_flow_conc_obs_plot
ggsave(paste(file.path(graphics_dir, "daily_obs_rain_flow_conc_plot.png")), 
       plot = rain_flow_conc_obs_plot, 
       width=11, height=8, units="in", dpi=300)


combined_rain_flow_conc_obs_plot <- plot_grid(daily_rain_flow_sims_plot, 
                                              daily_rain_conc_sims_plot, 
                                              daily_flow_conc_sims_plot,
                                              daily_rain_flow_obs_plot, 
                                               daily_rain_conc_obs_plot, 
                                               daily_flow_conc_obs_plot,
                                               ncol = 3, labels = c("A", "B", "C",
                                                                    "D", "E", "F"))
combined_rain_flow_conc_obs_plot
ggsave(paste(file.path(graphics_dir, "combined_daily_rain_flow_conc_plot.png")), 
       plot = combined_rain_flow_conc_obs_plot, 
       width=11, height=8, units="in", dpi=300)


# merge simulation and observation data on date
class(rainfall_response_data$date)
colnames(rainfall_response_data) # sims
class(filtered_obs_merged$date) <- "numeric"
colnames(filtered_obs_merged) # obs
# merge observations on date
merge_list <- list(rainfall_response_data, filtered_obs_merged)
sims_obs_merged <- merge_list %>% reduce(full_join, by = "date")
head(sims_obs_merged)
min(sims_obs_merged$discharge, na.rm=TRUE)

# rainfall vs flow
colnames(sims_obs_merged)
rainfall_flow_plot <- ggplot(sims_obs_merged, aes(x = rainfall)) +
  geom_point(aes(y = flow), color = "lightblue") +
  geom_point(aes(y = discharge), color = "#FF6666") +
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Flow (cfs))") +
  scale_y_log10()+
  geom_smooth(aes(y = flow), method = "loess", se = FALSE, color = "blue", size = 1.2)+
  geom_smooth(aes(y = discharge), method = "loess", se = FALSE, color = "red", size = 1.2)+
  theme_classic()
rainfall_flow_plot

# rainfall vs conc
rainfall_conc_plot <- ggplot(sims_obs_merged, aes(x = rainfall)) +
  geom_point(aes(y = bacteria_conc), color = "lightblue") +
  geom_point(aes(y = bacteria), color = "#FF6666") +
  xlab("log10(Rainfall (mm))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_y_log10()+
  geom_smooth(aes(y = bacteria_conc), method = "loess", se = FALSE, color = "blue", size = 1.2)+
  geom_smooth(aes(y = bacteria), method = "loess", se = FALSE, color = "red", size = 1.2)+
  theme_classic()

rainfall_conc_plot

# flow vs conc
flow_conc_plot <- ggplot(sims_obs_merged, aes(x = flow)) +
  # Points: use 'Data Type Point' and 'Data Type Line'
  geom_point(aes(y = bacteria_conc, color = "Simulated"), size = 1.2) +
  geom_point(aes(y = bacteria, color = "Observed"), size = 1.2) +
  geom_smooth(aes(y = bacteria_conc, color = "Simulated loess"), method = "loess", se = FALSE, size = 1.2) +
  geom_smooth(aes(y = bacteria, color = "Observed loess"), method = "loess", se = FALSE, size = 1.2) +
  xlab("log10(Flow (cfs))") +
  ylab("log10(Concentration (MPN/100 ml))") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    name = "Legend",
    values = c(
      "Simulated" = "lightblue",
      "Observed" = "#FF6666",
      "Simulated loess" = "blue",
      "Observed loess" = "red"
    )
  ) +
  theme_classic()
flow_conc_plot

combined_rain_flow_conc_obs_sims_plot <- plot_grid(rainfall_flow_plot, 
                                              rainfall_conc_plot, 
                                              flow_conc_plot,
                                              ncol = 3, 
                                              rel_widths = c(1, 1, 1.5),
                                              labels = c("D", "E", "F"))
combined_rain_flow_conc_obs_sims_plot

ggsave(paste(file.path(graphics_dir, "combined_daily_rain_flow_conc_obs_sims_plot.png")), 
       plot = combined_rain_flow_conc_obs_sims_plot, 
       width=11, height=6, units="in", dpi=300)




###################################
length(na.omit(sims_obs_merged$rainfall))
length(na.omit(sims_obs_merged$precipitation))

precip_df <- data.frame(
  value = c(na.omit(sims_obs_merged$rainfall), na.omit(sims_obs_merged$precipitation)),
  group = c(rep("Sim",length(na.omit(sims_obs_merged$rainfall))), rep("Obs",length(na.omit(sims_obs_merged$precipitation)))) )

colnames(precip_df)

# Overlapping kernel density plot
ggplot(precip_df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4) +           # alpha controls transparency for overlap
  labs(x = "Value", y = "Density", fill = "Group") +
  theme_minimal()

ggplot(precip_df, aes(x = value, fill = group, color = group)) +
  geom_histogram(aes(y=after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Scaled Dual Histogram",
       x = "Value",
       y = "Density",
       fill = "Group") +
  theme_minimal()
# Boxplot
ggplot(precip_df, aes(x = group, y = log(value), fill = group)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Boxplot Comparison",
       x = "Group",
       y = "Value") +
  theme_minimal()
# Violin plot
ggplot(precip_df, aes(x = group, y = log(value), fill = group)) +
  geom_violin(alpha = 0.7) +
  labs(title = "Violin Plot Comparison",
       x = "Group",
       y = "Value") +
  theme_minimal()


# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")
