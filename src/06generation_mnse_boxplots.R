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

# For memory used by R
cat("Current R memory used: ", mem_used()/1e+9, "\n")
# Cross platform
m <- ps::ps_system_memory()
cat("Total RAM:", m$total/1E+9, "\n")
cat("Available RAM:", m$avail/1E+9, "\n")
