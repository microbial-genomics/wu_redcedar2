#####################################
#####time variation quantile plots####
####################################
# simulated concentrations for final generation
daily_conc_sims <- bac_cal_output11$simulation$bac_out
dim(daily_conc_sims) #3865
class(daily_conc_sims)

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

# find the keepers ...
# THIS IS IMPORTANT, CRAP SIMS ARE STILL INCLUDED
generation_stats

nses_parameters11
dim(nses_parameters11)

str(bac_cal_output11)
str(bac_cal_output11$parameter$values)
str(bac_cal_output11$parameter$definition)
str(bac_cal_output11$simulation$bac_out)
str(bac_cal_output11$simulation$q_out)

# date<-as.vector(conc[,1])
# length(date)
# nrow(date)
# dim(date)
# class(date)
conc2 <- daily_conc_sims[,c(2:2139)]
#conc2
#dim(conc2) #3865*1000
conc3<-t(conc2)
conc3
dim(conc3)
colnames(conc3)<-as.Date(sim_date,"%Y%m%d")
class(colnames(conc3))
class(colnames(conc3))<-"Date"
conc3
head(conc3)
colnames(conc3)
class(conc3)

dim(conc3)
quantile_results <- data.frame(matrix(nrow = 3865, ncol = 8))
colnames(quantile_results) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:3865) {
  first_q <- quantile(conc3[,i], c(0.001))
  sec_q <- quantile(conc3[,i], c(0.023))
  third_q <- quantile(conc3[,i], c(0.159))
  four_q <- quantile(conc3[,i], c(0.5))
  fif_q <- quantile(conc3[,i], c(0.841))
  six_q <- quantile(conc3[,i], c(0.977))
  sev_q <- quantile(conc3[,i], c(0.99))
  quantile_results[i, 1] <- i
  quantile_results[i, 2] <- first_q
  quantile_results[i, 3] <- sec_q
  quantile_results[i, 4] <- third_q
  quantile_results[i, 5] <- four_q
  quantile_results[i, 6] <- fif_q
  quantile_results[i, 7] <- six_q
  quantile_results[i, 8] <- sev_q
}
head(quantile_results)
#View(quantile_results)

quantile_results2 <- cbind(date, quantile_results[,c(2:8)])
head(quantile_results2)
dim(quantile_results2)
# plot(x=result2$date, y=result2$first_q)
# plot(x=result2$date, y=result2$sec_q)
# result3<-merge(bac_obs, result2, by="date", incomparables = NA)
# dim(result3)




# The summary data frame ds is used to plot larger red points on top
# of the raw data. Note that we don't need to supply `data` or `mapping`
# # in each layer because the defaults from ggplot() are used.
# p <-ggplot() +
#   # geom_point(data=result2, aes(date,first_q),colour = "#999999",size =1)+ 
#   # geom_point(data=result2,aes(date,sec_q),colour = "#E69F00",size =1)+ 
#   # geom_point(data=result2,aes(date,third_q),colour ="#56B4E9",size =1)+ 
#   # geom_point(data=result2,aes(date,four_q),colour =  "#009E73",size =1)+ 
#   # geom_point(data=result2,aes(date,fif_q),colour ="#F0E442",size =1)+ 
#   # geom_point(data=result2,aes(date,six_q),colour =  "#0072B2",size =1)+ 
#   geom_line(data=result2,aes(date,sev_q),colour ="#D55E00",size =1)+
#   theme()+
#   theme(panel.spacing = unit(0.2, "lines"),
#         legend.position = "bottom")+
#   labs( x = "Date (yyyy)", 
#         y = "Concentration") 
# p  



quantile_results2[quantile_results2<1] <- 1
conc_daily_quantile_results <- quantile_results2[c(2193:2793),]
range(conc_daily_quantile_results$first_q)

dim(bac_obs)

# still working this out
# https://www.geeksforgeeks.org/r-language/combine-and-modify-ggplot2-legends-with-ribbons-and-lines/
daily_concs_quantiles_ribbon_plot <- ggplot(conc_daily_quantile_results, aes(x = date)) +
  geom_ribbon(aes(ymin = first_q, ymax = sev_q), fill = "Median+-3s", alpha = 0.2) +
  geom_ribbon(aes(ymin = sec_q, ymax = six_q), fill = "Median+-2s", alpha = 0.4) +
  geom_ribbon(aes(ymin = third_q, ymax = fif_q), fill = "Median+-1s", alpha = 0.6) +
  geom_line(aes(y = four_q), color = "slateblue4", linewidth = 0.8, alpha=0.8) +
  geom_point(data = bac_obs[c(155:180),], aes(x=date, y=bacteria, colour = "bacteria"))+
  scale_fill_manual(name = "Ribbon", 
                    values = c("Median+-3s" = "skyblue",
                               "Median+-2s" = "skyblue3",
                               "Median+-1s" = "skyblue4")) +
                    #) +       # manual legend
  #scale_fill_manual(name = "Ribbon", values = c()) +       # manual legend
  #scale_fill_manual(name = "Ribbon", values = c() +       # manual legend
  #scale_color_manual(name = "Lines", values = c("Lower" = "blue", "Upper" = "red")) +    # manual legend
  xlab(" ") +
  # scale_y_continuous("Concentration (log10)") + 
  scale_x_date(
    date_breaks = "2 months",      # 3 month ticks
    date_labels = "%b %Y"          # month abbreviation, year
  ) +
  scale_y_log10() +
  labs(title="Daily Concentration")+
  theme_bw()

daily_concs_quantiles_ribbon_plot


#### daily flow  #### 
# simulated flows of final generation
daily_flow_sims <- bac_cal_output11$simulation$q_out
dim(daily_flow_sims)
class(daily_flow_sims)
#load("E:/boxplot/daily_flow.RData")
#flow<-bac_cal_output$simulation$q_out

date <- daily_flow_sims$date
daily_flow_sims2 <- daily_flow_sims[,c(2:2139)]
#conc2
#dim(conc2) #3865*1000
daily_flow_sims3 <- t(daily_flow_sims2)
# conc3
# dim(conc3)
colnames(daily_flow_sims3)<-as.Date(date,"%Y%m%d")
class(colnames(daily_flow_sims3))
class(colnames(daily_flow_sims3))<-"Date"
daily_flow_sims3
dim(daily_flow_sims3)
head(daily_flow_sims3)
# colnames(conc3)
# class(conc3)

# day1=conc3[,1]
# q1<-quantile(conc3[,1], c(0.001, 0.023, 0.159, 0.5, 0.841, 0.977, 0.99))
# q1

result <- data.frame(matrix(nrow = 3865, ncol = 8))
colnames(result) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:3865) {
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

result2<-cbind(date,result[,c(2:8)])
head(result2)
dim(result2)
# plot(x=result2$date, y=result2$first_q)
# plot(x=result2$date, y=result2$sec_q)
# result3<-merge(bac_obs, result2, by="date", incomparables = NA)
# dim(result3)

result2[result2<1] <-1
result2

range(result2$first_q)

ggplot() +
  geom_line(data = result2[c(2193:2793),], aes(x = date, y = log10(first_q), colour = "first_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date, y = sec_q, colour = "sec_q")) +
  # geom_line(data = result2, aes(x = date, y = third_q, colour = "third_q")) +
  geom_line(data = result2[c(2193:2793),], aes(x = date, y = log10(four_q), colour = "four_q")) +
  #geom_line(data = result2, aes(x = date,y = fif_q, colour = "fif_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date,y = six_q, colour = "six_q")) +
  geom_line(data = result2[c(2193:2793),], aes(x = date,y = log10(sev_q), colour = "sev_q")) +
  geom_point(data =q_obs[c(2193:2793),], aes(x=date, y=log10(discharge), colour = "discharge"))+
  scale_colour_manual("", 
                      breaks = c("first_q",   "four_q",  "sev_q",  "discharge"),
                      values =c("grey70", "#CC79A7", "grey30",  "darkblue")) +
  xlab(" ") +
  scale_y_continuous("log10(discharge) (cms)") + 
  labs(title="Daily_flow")+
  theme_bw()



#### weekly flow ####
load("E:/boxplot/weekly_flow.RData")
load(file="E:/boxplot/q_obs.RData")
flow_obs <- q_obs # [4018,2]
flow_obs_daily <- flow_obs$discharge 
obs_flow_xts <- as.xts(flow_obs$discharge,order.by=as.Date(flow_obs$date))
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


###########################
#### boxplot ####
#####################
####mnse_mean_daily####
##run this one first, otherwise, it will overwrite other result###
load("E:/boxplot/daily_mean.RData")
load("E:/boxplot/flux_obs.RData")
load("E:/boxplot/bac_obs.RData")
load("E:/boxplot/q_obs.RData")

sim_bac <- bac_cal_output$simulation$bac_out
# merge simulated and observed bacteria concentrations, calculate nses for all sims
mnse_bac_daily<- right_join(sim_bac,bac_obs,by="date")%>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~mNSE(.x, bac_obs$bacteria))
sim_q <- bac_cal_output$simulation$q_out
# merge simulated and observed flows, calculate nses for all sims
mnse_flow_daily <- right_join(sim_q,q_obs,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-discharge) %>%
  map_dbl(., ~mNSE(.x, q_obs$discharge))
date <-bac_cal_output$simulation$bac_out$date
flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
sim_flux<- cbind(date, flux_sim)
#merge simulated and observed fluxes, calculate nses for all sims
mnse_flux_daily <-  right_join(sim_flux, flux_obs, by = "date") %>%
  dplyr::select(-date) %>%dplyr::select(-flux) %>%
  map_dbl(., ~mNSE(.x, flux_obs$flux))

mnse_mean_daily <- rowMeans(cbind(mnse_bac_daily, mnse_flow_daily, mnse_flux_daily))

####mnse_mean_weekly####
##run mean second, otherwise, it will rewrite bac_nse and q_nse
load("E:/boxplot/weekly_mean.RData")
load("E:/boxplot/flux_obs_w.RData")
load("E:/boxplot/bac_obs_w.RData")
load("E:/boxplot/q_obs_w.RData")
load("E:/boxplot/flux_obs.RData")
load("E:/boxplot/bac_obs.RData")
load("E:/boxplot/q_obs.RData")
flow_obs <- right_join(q_obs, bac_obs, by="date")
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,1503]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
#dim(bac_flows_daily_temp)
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[203,2214]
# dim(bac_flows_weekly)
bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
#dim(bac_fluxes_weekly) #[204,2214]
mnse_flux_weekly <- mapply(mNSE, bac_fluxes_weekly, flux_obs_weekly)
mnse_bac_weekly <- mapply(mNSE, bac_sims_weekly, bac_obs_weekly)
mnse_flow_weekly <- mapply(mNSE, bac_flows_weekly, flow_obs_weekly)
mnse_mean_weekly<-rowMeans(cbind(mnse_bac_weekly, mnse_flow_weekly, mnse_flux_weekly))


####mnse_mean_monthly####
##run mean second, otherwise, it will rewrite bac_nse and q_nse
load("E:/boxplot/monthly_mean.RData")
load("E:/boxplot/flux_obs_m.RData")
load("E:/boxplot/bac_obs_m.RData")
load("E:/boxplot/q_obs_m.RData")
load("E:/boxplot/flux_obs.RData")
load("E:/boxplot/bac_obs.RData")
load("E:/boxplot/q_obs.RData")
flow_obs <- right_join(q_obs, bac_obs, by="date")
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,1503]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_monthly <- as.data.frame(apply.monthly(bac_sims_daily_data,mean)) #[204,2214]
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
#dim(bac_flows_daily_temp)
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_monthly <- as.data.frame(apply.monthly(bac_flows_daily_data,mean)) #[203,2214]
# dim(bac_flows_monthly)
bac_fluxes_monthly <- bac_sims_monthly * bac_flows_monthly * 10^4
#dim(bac_fluxes_monthly) #[204,2214]
mnse_flux_monthly <- mapply(mNSE, bac_fluxes_monthly, flux_obs_monthly)
mnse_bac_monthly <- mapply(mNSE, bac_sims_monthly, bac_obs_monthly)
mnse_flow_monthly <- mapply(mNSE, bac_flows_monthly, flow_obs_monthly)

mnse_mean_monthly<-rowMeans(cbind(mnse_bac_monthly, mnse_flow_monthly, mnse_flux_monthly))

#####mnse_flux_daily####
##run flux second, otherwise, it will rewrite bac_nse and q_nse
load("E:/boxplot/daily_flux.RData")
load("E:/boxplot/flux_obs.RData")
sim_bac <- bac_cal_output$simulation$bac_out
sim_q <- bac_cal_output$simulation$q_out
date <-bac_cal_output$simulation$bac_out$date
flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
sim_flux<- cbind(date, flux_sim)
#merge simulated and observed fluxes, calculate nses for all sims

mnse_flux_daily <-  right_join(sim_flux, flux_obs, by = "date") %>%
  dplyr::select(-date) %>%dplyr::select(-flux) %>%
  map_dbl(., ~mNSE(.x, flux_obs$flux))

#####mnse_flux_weekly####
##run flux second, otherwise, it will rewrite bac_nse and q_nse
load("E:/boxplot/weekly_flux.RData")
load("E:/boxplot/flux_obs_w.RData")
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
#dim(bac_flows_all_days)
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
#dim(bac_flows_daily_temp)
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[203,2214]
# dim(bac_flows_weekly)
bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
#dim(bac_fluxes_weekly) #[204,2214]
mnse_flux_weekly <- mapply(mNSE, bac_fluxes_weekly, flux_obs_weekly)


####mnse_flux_monthly####
##run flux second, otherwise, it will rewrite bac_nse and q_nse
load("E:/boxplot/monthly_flux.RData")
load("E:/boxplot/flux_obs_m.RData")
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_monthly <- as.data.frame(apply.monthly(bac_sims_daily_data,mean)) #[204,2214]
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
#dim(bac_flows_all_days)
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
#dim(bac_flows_daily_temp)
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_monthly <- as.data.frame(apply.monthly(bac_flows_daily_data,mean)) #[203,2214]
# dim(bac_flows_monthly)
bac_fluxes_monthly <- bac_sims_monthly * bac_flows_monthly * 10^4
#dim(bac_fluxes_monthly) #[204,2214]
mnse_flux_monthly <- mapply(mNSE, bac_fluxes_monthly, flux_obs_monthly)


#mnse_bac_daily
load("E:/boxplot/bac_obs.RData")
load("E:/boxplot/daily_conc.RData")
sim_bac <- bac_cal_output$simulation$bac_out
# merge simulated and observed bacteria concentrations, calculate nses for all sims
mnse_bac_daily<- right_join(sim_bac,bac_obs,by="date")%>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~mNSE(.x, bac_obs$bacteria))

####mnse_bac_weekly####
load("E:/boxplot/bac_obs_w.RData")
load("E:/boxplot/weekly_conc.RData")

bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
# adds date and bacteria fields also
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
#dim(bac_sims_daily_temp)
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]

mnse_bac_weekly <- mapply(mNSE, bac_sims_weekly, bac_obs_weekly)

####mnse_bac_monthly####
load("E:/boxplot/bac_obs_m.RData")
load("E:/boxplot/monthly_conc.RData")

bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
# adds date and bacteria fields also
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
# dim(bac_sims_daily_temp)
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_monthly <- as.data.frame(apply.monthly(bac_sims_daily_data,mean)) 

mnse_bac_monthly <- mapply(mNSE, bac_sims_monthly, bac_obs_monthly)
range(mnse_bac_monthly)


#####mnse_flow_daily####
load("E:/boxplot/q_obs.RData")
load("E:/boxplot/daily_flow.RData")
sim_q <- bac_cal_output$simulation$q_out
# merge simulated and observed flows, calculate nses for all sims

mnse_flow_daily <- right_join(sim_q,q_obs,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-discharge) %>%
  map_dbl(., ~mNSE(.x, q_obs$discharge))

####mnse_flow_weekly####
load("E:/boxplot/q_obs_w.RData")
load("E:/boxplot/weekly_flow.RData")
flow_obs <- right_join(q_obs, bac_obs, by="date")
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] 
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
#head(colnames(bac_flows_daily)) #  date 
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) 
bac_flows_daily <- bac_flows_daily[,-1] #
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[203,2214]
#dim(bac_flows_weekly)
# bac_flows_monthly <- as.data.frame(apply.monthly(bac_flows_daily_data,mean)) #51

mnse_flow_weekly<- mapply(mNSE, bac_flows_weekly, flow_obs_weekly)


####mnse_flow_monthly####
load("E:/boxplot/q_obs_m.RData")
load("E:/boxplot/monthly_flow.RData")
flow_obs <- right_join(q_obs, bac_obs, by="date")
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] 
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
#head(colnames(bac_flows_daily)) #  date 
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) 
bac_flows_daily <- bac_flows_daily[,-1] #
bac_flows_monthly <- as.data.frame(apply.monthly(bac_flows_daily_data,mean)) #51

mnse_flow_monthly<- mapply(mNSE, bac_flows_monthly, flow_obs_monthly)

####boxplot####
daily<- replicate(1000, "Daily")
weekly<-replicate(1000,"Weekly")
monthly<-replicate(1000,"Monthly")
bacteria<- replicate(1000, "bacteria")
flow<- replicate(1000, "flow")
flux<- replicate(1000, "flux")
mean<- replicate(1000, "mean")

data1<-cbind(mnse_bac_daily[c(1:1000)], daily, bacteria)
data2<-cbind(mnse_flow_daily[c(1:1000)], daily, flow)
data3<-cbind(mnse_flux_daily[c(1:1000)], daily, flux)
data4<-cbind(mnse_mean_daily[c(1:1000)], daily, mean)

data5<-cbind(mnse_bac_weekly[c(1:1000)], weekly, bacteria)
data6<-cbind(mnse_flow_weekly[c(1:1000)], weekly, flow)
data7<-cbind(mnse_flux_weekly[c(1:1000)], weekly, flux)
data8<-cbind(mnse_mean_weekly[c(1:1000)], weekly, mean)

data9<-cbind(mnse_bac_monthly[c(1:1000)], monthly, bacteria)
data10<-cbind(mnse_flow_monthly[c(1:1000)], monthly, flow)
data11<-cbind(mnse_flux_monthly[c(1:1000)], monthly, flux)
data12<-cbind(mnse_mean_monthly[c(1:1000)], monthly, mean)

data <-as.data.frame(rbind(data1,data2,data3,data4,data5,data6, data7,data8,data9,data10,data11,data12))
dim(data)
head(data)
colnames(data)<-c("mnse", "interval","target")
# is.numeric(data$mnse)
data$mnse<-as.numeric(levels(data$mnse))[data$mnse]
data$interval<-as.factor(data$interval)


ggplot(data, aes(x=target, y=mnse,color=target)) + 
  geom_boxplot()+
  theme_classic()


data %>%
  arrange(interval)%>%
  mutate(interval = factor(interval,levels = c("Daily","Weekly","Monthly"))) %>% 
  ggplot(aes(x=target, y=mnse, fill=interval)) +
  geom_boxplot()+
  labs(title="Modified nse summary",x="Targets", y = "mNSE")+
  # scale_y_discrete(limits=c("0.2","0.4","0.6","0.8","1"))+
  theme_bw()


data %>%
  arrange(interval)%>%
  mutate(interval = factor(interval,levels = c("Daily","Weekly","Monthly"))) %>% 
  ggplot(aes(x=interval, y=mnse, fill=target)) +
  geom_boxplot()+
  labs(title="Modified nse summary, first 1000 simulations",x="Time Interval", y = "mNSE")+
  # scale_y_discrete(limits=c("0.2","0.4","0.6","0.8","1"))+
  theme_bw()
