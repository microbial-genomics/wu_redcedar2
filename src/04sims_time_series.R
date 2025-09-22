#####################################
#####time variation quantile plots####
####################################
# simulated concentrations for final generation
daily_conc_sims <- bac_cal_output11$simulation$bac_out
dim(daily_conc_sims) #3865
class(daily_conc_sims)


df_unique <- df[!(df$col %in% df$col[duplicated(df$col)]), ]


#### daily bacteria concentration #### 
## load bac_cal_output data
#load("E:/boxplot/daily_conc.RData")
#load("E:/boxplot/bac_obs.RData")

#conc<-bac_cal_output$simulation$bac_out
#class(bac_cal_output)
#class(conc)
sim_date <- daily_conc_sims$date
class(sim_date)
length(sim_date)
sim_date[duplicated(sim_date)] # no duplicates

# find the keepers ...
generation_stats
# nses_parameters files has all the keepers and the rejected simulations
nses_parameters11
dim(nses_parameters11)
str(bac_cal_output11)
str(bac_cal_output11$parameter$values)
str(bac_cal_output11$parameter$definition)
str(bac_cal_output11$simulation$bac_out)
str(bac_cal_output11$simulation$q_out)

# strip out numbers of keepers
keepers_text <- rownames(nses_parameters11)
keepers_indices <- as.integer(sub("^0+", "", substr(keepers_text, 6, 9)))
keepers_indices # use this to subset the full generation of sims

#### daily concentrations  #### 
# extract the concentration data for the keepers only
dim(daily_conc_sims)
conc2 <- daily_conc_sims[,c(2:2139)]
#dim(conc2) #3865 days 2138 sims
conc3 <- t(conc2)
conc3
dim(conc3) # 2138 sims 3865 days
colnames(conc3)<-as.Date(sim_date,"%Y%m%d")
class(colnames(conc3))
class(colnames(conc3))<-"Date"
conc3
head(conc3)
colnames(conc3)
class(conc3)
### re-assign and subset conc3 to be the reduced set of sims
dim(conc3)
# View(conc3)
conc4 <- conc3[keepers_indices,]
dim(conc4) # 200 (sims) 3865 (days)

# calculate concentration quantiles for each day
concentration_daily_quantile_results <- data.frame(matrix(nrow = 3865, ncol = 8))
colnames(concentration_daily_quantile_results) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:3865) {
  first_q <- quantile(conc4[,i], c(0.001))
  sec_q <- quantile(conc4[,i], c(0.023))
  third_q <- quantile(conc4[,i], c(0.159))
  four_q <- quantile(conc4[,i], c(0.5))
  fif_q <- quantile(conc4[,i], c(0.841))
  six_q <- quantile(conc4[,i], c(0.977))
  sev_q <- quantile(conc4[,i], c(0.99))
  concentration_daily_quantile_results[i, 1] <- i
  concentration_daily_quantile_results[i, 2] <- first_q
  concentration_daily_quantile_results[i, 3] <- sec_q
  concentration_daily_quantile_results[i, 4] <- third_q
  concentration_daily_quantile_results[i, 5] <- four_q
  concentration_daily_quantile_results[i, 6] <- fif_q
  concentration_daily_quantile_results[i, 7] <- six_q
  concentration_daily_quantile_results[i, 8] <- sev_q
}
head(concentration_daily_quantile_results)
dim(concentration_daily_quantile_results)
#View(concentration_daily_quantile_results)

concentration_daily_quantile_results2 <- cbind(sim_date, concentration_daily_quantile_results[,c(2:8)])
head(concentration_daily_quantile_results2)
dim(concentration_daily_quantile_results2)

# truncate to min = 1
concentration_daily_quantile_results2[concentration_daily_quantile_results2<1] <- 1
dim(concentration_daily_quantile_results2)
# we subset the dates to this period for reasons I do not recall
conc_daily_quantile_results <- concentration_daily_quantile_results2[c(2193:3865),]
range(conc_daily_quantile_results$first_q)
# check on bacteria observations
dim(bac_obs)

### narrow in on desired time range
dim(bac_obs)
head(bac_obs)
#View(bac_obs)
# 184 = April 7, 2013
# 335 last observation on May 1, 2014
dim(conc_daily_quantile_results)
conc_daily_quantile_results$sim_date[1193]
# 1193 - April 7, 2013
conc_daily_quantile_results$sim_date[1582]
# 1582 = May 1, 2014

# create ribbon plots
# https://www.geeksforgeeks.org/r-language/combine-and-modify-ggplot2-legends-with-ribbons-and-lines/
# daily concentrations subset for April 7, 2013 - May 1, 2014
daily_concs_quantiles_ribbon_plot_subset <- ggplot(conc_daily_quantile_results[1193:1582,], aes(x = sim_date)) +
  geom_ribbon(aes(ymin = first_q, ymax = sev_q, fill = "Median \u00B1 3 SD"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sec_q, ymax = six_q, fill = "Median \u00B1 2 SD"), alpha = 0.4) +
  geom_ribbon(aes(ymin = third_q, ymax = fif_q, fill = "Median \u00B1 1 SD"), alpha = 0.6) +
  geom_line(aes(y = four_q), color = "slateblue4", linewidth = 0.8, alpha=0.8) +
  geom_point(data = bac_obs[c(184:335),], aes(x=date, y=bacteria, color = "Bacteria"))+ #c(155:334)
  labs(color = "Observations") +
  scale_fill_manual(name = "Simulations", 
                    values = c("Median \u00B1 3 SD" = "skyblue",
                               "Median \u00B1 2 SD" = "skyblue3",
                               "Median \u00B1 1 SD" = "skyblue4")) +
  xlab(" ") +
  ylab("Concentration (MPN/100 ml)") +
  scale_x_date(
    date_breaks = "2 months",      # 3 month ticks
    date_labels = "%b %Y"          # month abbreviation, year
  ) +
  scale_y_log10() +
  labs(title="Daily Concentration")+
  theme_bw()

daily_concs_quantiles_ribbon_plot_subset
ggsave(paste(file.path(graphics_dir, "daily_concs_quantiles_ribbon_plot_subset.png")), 
       plot = daily_concs_quantiles_ribbon_plot_subset, 
       width=8, height=5, units="in", dpi=300)

# daily concentrations all dates
daily_concs_quantiles_ribbon_plot_all <- ggplot(concentration_daily_quantile_results2, aes(x = sim_date)) +
  geom_ribbon(aes(ymin = first_q, ymax = sev_q, fill = "Median \u00B1 3 SD"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sec_q, ymax = six_q, fill = "Median \u00B1 2 SD"), alpha = 0.4) +
  geom_ribbon(aes(ymin = third_q, ymax = fif_q, fill = "Median \u00B1 1 SD"), alpha = 0.6) +
  geom_line(aes(y = four_q), color = "slateblue4", linewidth = 0.8, alpha=0.8) +
  geom_point(data = bac_obs, aes(x=date, y=bacteria, color = "Bacteria"))+ #c(155:334)
  labs(color = "Observations") +
  scale_fill_manual(name = "Simulations", 
                    values = c("Median \u00B1 3 SD" = "skyblue",
                               "Median \u00B1 2 SD" = "skyblue3",
                               "Median \u00B1 1 SD" = "skyblue4")) +
  xlab(" ") +
  ylab("Concentration (MPN/100 ml)") +
  scale_x_date(
    date_breaks = "12 months",      # 12 month ticks
    date_labels = "%b %Y"          # month abbreviation, year
  ) +
  scale_y_log10() +
  labs(title="Daily Concentration")+
  theme_bw()

daily_concs_quantiles_ribbon_plot_all
ggsave(paste(file.path(graphics_dir, "daily_concs_quantiles_ribbon_plot_all.png")), 
       plot = daily_concs_quantiles_ribbon_plot_all, 
       width=11, height=3, units="in", dpi=300)


#### daily flow  #### 
# simulated flows of final generation
daily_flow_sims <- bac_cal_output11$simulation$q_out
class(daily_flow_sims)

flow_sim_date <- daily_flow_sims$date
daily_flow_sims2 <- daily_flow_sims[,c(2:2139)]
dim(daily_flow_sims2) # 3865 days 2138 sims
#dim(conc2) #3865*1000
daily_flow_sims3 <- t(daily_flow_sims2)
dim(daily_flow_sims3) # 2138 sims 3865 days

colnames(daily_flow_sims3)<-as.Date(flow_sim_date,"%Y%m%d")
class(colnames(daily_flow_sims3))
class(colnames(daily_flow_sims3))<-"Date"
daily_flow_sims3
dim(daily_flow_sims3)
head(daily_flow_sims3)
daily_flow_sims4 <- daily_flow_sims3[keepers_indices,]
dim(daily_flow_sims4) # 200 (sims) 3865 (days)

# calculate flow quantiles for each day
flow_daily_quantile_results <- data.frame(matrix(nrow = 3865, ncol = 8))
colnames(flow_daily_quantile_results) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:3865) {
  first_q <- quantile(daily_flow_sims4[,i], c(0.001))
  sec_q <- quantile(daily_flow_sims4[,i], c(0.023))
  third_q <- quantile(daily_flow_sims4[,i], c(0.159))
  four_q <- quantile(daily_flow_sims4[,i], c(0.5))
  fif_q <- quantile(daily_flow_sims4[,i], c(0.841))
  six_q <- quantile(daily_flow_sims4[,i], c(0.977))
  sev_q <- quantile(daily_flow_sims4[,i], c(0.99))
  flow_daily_quantile_results[i, 1] <- i
  flow_daily_quantile_results[i, 2] <- first_q
  flow_daily_quantile_results[i, 3] <- sec_q
  flow_daily_quantile_results[i, 4] <- third_q
  flow_daily_quantile_results[i, 5] <- four_q
  flow_daily_quantile_results[i, 6] <- fif_q
  flow_daily_quantile_results[i, 7] <- six_q
  flow_daily_quantile_results[i, 8] <- six_q
  i=i+1  
}
head(flow_daily_quantile_results)

flow_daily_quantile_results2 <- cbind(sim_date, flow_daily_quantile_results[,c(2:8)])
head(flow_daily_quantile_results2)
dim(flow_daily_quantile_results2) #3865 8

# don't truncate
min(flow_daily_quantile_results2$first_q)
dim(flow_daily_quantile_results2)

### narrow in on desired time range
dim(bac_obs)
#View(bac_obs)
# 184 = April 7, 2013
# 335 last observation on May 1, 2014



# we subset the dates to this period for reasons I do not recall

flow_daily_quantile_results <- flow_daily_quantile_results2[c(2193:2793),]
#View(flow_daily_quantile_results)
range(flow_daily_quantile_results$first_q)

# check on flow observations
flow_obs <- q_obs # [4018,2]
flow_obs_daily <- flow_obs$discharge 
obs_flow_xts <- as.xts(flow_obs$discharge,order.by=as.Date(flow_obs$date))
dim(obs_flow_xts) #4018 1
head(obs_flow_xts)
#View(obs_flow_xts)
#colnames(obs_flow_xts) <- c("date", "flow")
obs_flow_xts[3415,]
# 3415 - April 7, 2013
obs_flow_xts[3774,]
# 3774 = May 1, 2014
colnames(obs_flow_xts) <- c("flow")

# add empirical flows observations to the flow quantile data frame
dim(flow_daily_quantile_results2)
colnames(flow_daily_quantile_results2)
dim(obs_merged)
colnames(obs_merged)
obs_merged$date[3865]
flow_daily_quantile_results2$sim_date[3865]
#View(cbind(obs_merged$date[1:3865], flow_daily_quantile_results2$sim_date[1:3865]))
#View(obs_merged)
colnames(flow_daily_quantile_results2)
flow_daily_quantile_results2$sim_date[1]
flow_daily_quantile_results2 <- cbind(sim_date, flow_daily_quantile_results2[,c(2:8)])
flow_daily_quantile_results2 <- cbind(flow_daily_quantile_results2, obs_merged[1:3865,])
#View(flow_daily_quantile_results2)

# create ribbon plots
# https://www.geeksforgeeks.org/r-language/combine-and-modify-ggplot2-legends-with-ribbons-and-lines/
colnames(flow_daily_quantile_results2)
dim(flow_daily_quantile_results2)

daily_flow_quantiles_ribbon_plot <- ggplot(flow_daily_quantile_results2[3415:3774,], aes(x = sim_date)) +
  geom_ribbon(aes(ymin = first_q, ymax = sev_q, fill = "Median \u00B1 3 SD"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sec_q, ymax = six_q, fill = "Median \u00B1 2 SD"), alpha = 0.4) +
  geom_ribbon(aes(ymin = third_q, ymax = fif_q, fill = "Median \u00B1 1 SD"), alpha = 0.6) +
  geom_line(aes(y = four_q, color = "Simulated median"), linewidth = 0.8, alpha = 0.8) +
  geom_line(aes(y = discharge, color = "Observed"), linewidth = 0.6) +
  scale_fill_manual(name = "Simulations", 
                    values = c("Median \u00B1 3 SD" = "skyblue", "Median \u00B1 2 SD" = "skyblue3", "Median \u00B1 1 SD" = "skyblue4")) +
  scale_color_manual(
    name = "Line Type",
    values = c("Simulated median" = "slateblue4", "Observed" = "red2"),
    labels = c("Observed", "Simulated median")
  ) +
  xlab(" ") +
  ylab("Flow (cfs)") +
  scale_x_date(
    date_breaks = "2 months",      
    date_labels = "%b %Y"        
  ) +
  scale_y_log10() +
  labs(title = "Daily Flow") +
  theme_bw()

daily_flow_quantiles_ribbon_plot



stacked_daily_ribbon_plot <- plot_grid(daily_concs_quantiles_ribbon_plot, 
          daily_flow_quantiles_ribbon_plot, 
          ncol = 1, labels = c("A", "B"))

stacked_daily_ribbon_plot

png(paste(file.path(graphics_dir, "daily_stacked_simulation_ribbon_plot_appendix.png")), 
    width=8, height=11, units="in", res=300)
  stacked_daily_ribbon_plot
dev.off()


#### weekly flow ####
flow_obs <- q_obs # [4018,2]
flow_obs_daily <- flow_obs$discharge 
obs_flow_xts <- as.xts(flow_obs$discharge,order.by=as.Date(flow_obs$date))
dim(obs_flow_xts)
flow_obs_weekly <- as.data.frame(apply.weekly(obs_flow_xts, mean)) #575

flow_all_days <- bac_cal_output$simulation$q_out # [3865,1451]
# dim(flow_all_days)
# adds date and flow fields also
flow_daily_temp <- right_join(flow_all_days, flow_obs, by="date") #[4018,3540]
# dim(flow_daily_temp) #
flow_daily <- flow_daily_temp[,-which((colnames(flow_daily_temp)=="bacteria" | 
                                         colnames(flow_daily_temp)=="discharge"))] #[4018,3539]
nsim_cols <- ncol(flow_daily)
flow_daily_data <- as.xts(flow_daily[2:nsim_cols],order.by=as.Date(flow_daily$date))  #[4018,3538]
flow_daily <- flow_daily[,-1] #[4018,3538]
flow_weekly <- as.data.frame(apply.weekly(flow_daily_data,mean)) #[575,3538]
flow_weekly<-flow_weekly[,c(1:1000)] #[575,1000]
dim(flow_weekly)

flow3<-t(flow_weekly)
dim(flow3)
head(flow3)
date<-colnames(flow3)[1:550]
date<- as.Date(date, format="%Y-%m-%d")


result <- data.frame(matrix(nrow = 550, ncol = 8))
colnames(result) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:550) {
  first_q <- quantile(flow3[,i], c(0.001))
  sec_q <- quantile(flow3[,i], c(0.023))
  third_q <- quantile(flow3[,i], c(0.159))
  four_q <- quantile(flow3[,i], c(0.5))
  fif_q <- quantile(flow3[,i], c(0.841))
  six_q <- quantile(flow3[,i], c(0.977))
  sev_q <- quantile(flow3[,i], c(0.99))
  result[i, 1] <- i
  result[i, 2] <- first_q
  result[i, 3] <- sec_q
  result[i, 4] <- third_q
  result[i, 5] <- four_q
  result[i, 6] <- fif_q
  result[i, 7] <- six_q
  result[i, 8] <- six_q
  i=i+1  
}
head(result)

summary(result)

result2<-cbind(date,result[,c(2:8)])
head(result2)
dim(result2)


result2[result2<1] <-1
result2

class(result2$first_q)
range(result2$first_q)
date2<-row.names(flow_obs_weekly)
date2<-as.Date(date2)

flow_obs_weekly<-cbind(date2, flow_obs_weekly)
flow_obs_weekly
colnames(flow_obs_weekly)<-c("date","discharge")
dim(flow_obs_weekly)

ggplot() +
  geom_line(data = result2[c(250:500),], aes(x = date, y = log10(first_q), colour = "first_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date, y = sec_q, colour = "sec_q")) +
  # geom_line(data = result2, aes(x = date, y = third_q, colour = "third_q")) +
  geom_line(data = result2[c(250:500),], aes(x = date, y = log10(four_q), colour = "four_q")) +
  #geom_line(data = result2, aes(x = date,y = fif_q, colour = "fif_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date,y = six_q, colour = "six_q")) +
  geom_line(data = result2[c(250:500),], aes(x = date,y = log10(sev_q), colour = "sev_q")) +
  geom_point(data =flow_obs_weekly[c(250:500),], aes(x=date, y=log10(discharge), colour = "discharge"))+
  scale_colour_manual("", 
                      breaks = c("first_q",   "four_q",  "sev_q",  "discharge"),
                      values =c("grey70", "#CC79A7", "grey30",  "darkblue")) +
  xlab(" ") +
  scale_y_continuous("log10(discharge) (cms)") + 
  labs(title="Weekly_flow")+
  theme_bw()



#### weekly bacteria concentration #### 
load("E:/boxplot/weekly_conc.RData")
bac_obs_daily <- bac_obs$bacteria #335#removed the highest number on 2006/9/18
obs_data_xts <- as.xts(bac_obs$bacteria,order.by=as.Date(bac_obs$date))
bac_obs_weekly <- as.data.frame(apply.weekly(obs_data_xts, mean)) #204
bac_obs_weekly
dim(bac_obs_weekly)



bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,1571]
# dim(bac_sims_all_days)
# adds date and bacteria fields also
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
#dim(bac_sims_daily_temp)
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
#dim(bac_sims_daily)
#colnames(bac_sims_daily)[1]
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
#head(colnames(bac_sims_daily)) #  date 
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,1571]
dim(bac_sims_daily) #[336,1570]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[203,1570]

conc_weekly <-t(bac_sims_weekly)
conc_weekly<-conc_weekly[c(1:1000),]
# dim(conc_weekly) #1000*203
# head(flow3)
# colnames(conc3)
# class(conc3)

# day1=conc3[,1]
# q1<-quantile(conc3[,1], c(0.001, 0.023, 0.159, 0.5, 0.841, 0.977, 0.99))
# q1

result <- data.frame(matrix(nrow = 203, ncol = 8))
colnames(result) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:203) {
  first_q <- quantile(conc_weekly[,i], c(0.001))
  sec_q <- quantile(conc_weekly[,i], c(0.023))
  third_q <- quantile(conc_weekly[,i], c(0.159))
  four_q <- quantile(conc_weekly[,i], c(0.5))
  fif_q <- quantile(conc_weekly[,i], c(0.841))
  six_q <- quantile(conc_weekly[,i], c(0.977))
  sev_q <- quantile(conc_weekly[,i], c(0.99))
  result[i, 1] <- i
  result[i, 2] <- first_q
  result[i, 3] <- sec_q
  result[i, 4] <- third_q
  result[i, 5] <- four_q
  result[i, 6] <- fif_q
  result[i, 7] <- six_q
  result[i, 8] <- six_q
  i=i+1  
}
head(result)

date<-row.names(bac_sims_weekly)
date
date<- as.Date(date, format="%Y-%m-%d")
date
result2<-cbind(date,result[,c(2:8)])
head(result2)
dim(result2)


result2[result2<1] <-1
result2

range(result2$first_q)

bac_obs_weekly<-cbind(date, bac_obs_weekly)
bac_obs_weekly
dim(bac_obs_weekly)
colnames(bac_obs_weekly)<-c("date","bacteria")
head(bac_obs_weekly)

ggplot() +
  geom_line(data = result2, aes(x = date, y = log10(first_q), colour = "first_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date, y = sec_q, colour = "sec_q")) +
  # geom_line(data = result2, aes(x = date, y = third_q, colour = "third_q")) +
  geom_line(data = result2, aes(x = date, y = log10(four_q), colour = "four_q")) +
  #geom_line(data = result2, aes(x = date,y = fif_q, colour = "fif_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date,y = six_q, colour = "six_q")) +
  geom_line(data = result2, aes(x = date,y = log10(sev_q), colour = "sev_q")) +
  geom_point(data =bac_obs_weekly, aes(x=date, y=log10(bacteria), colour = "bacteria"))+
  scale_colour_manual("", 
                      breaks = c("first_q",   "four_q",  "sev_q",  "bacteria"),
                      values =c("grey70", "#CC79A7", "grey30",  "red")) +
  xlab(" ") +
  scale_y_continuous("log10(E. coli) (MPN/100ml)") + 
  labs(title="Weekly_concentration")+
  theme_bw()

# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")
gc()
