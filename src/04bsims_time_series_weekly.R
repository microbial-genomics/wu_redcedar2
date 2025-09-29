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
