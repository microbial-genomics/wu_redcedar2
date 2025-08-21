library(xts)
library(tibble)
#glimpse(df)


load(file="/work/OVERFLOW/RCR/sim53/bac_obs.RData")
data <- as.xts(bac_obs$bacteria,order.by=as.Date(bac_obs$date))
bac_weekly<- apply.weekly(data,mean)
bacteria<- coredata(bac_weekly)
class(bacteria)
bacteria <-bacteria[,1]
date <- index(bac_weekly)
bac_obs_weekly<-tibble(date,bacteria)
save(bac_obs_weekly, file="/work/OVERFLOW/RCR/sim53/bac_obs_w.RData")
load(file="/work/OVERFLOW/RCR/sim53/bac_obs_w.RData")

load(file="/work/OVERFLOW/RCR/sim53/q_obs.RData")
data2 <- as.xts(q_obs$discharge,order.by=as.Date(q_obs$date))
q_weekly<- apply.weekly(data2,mean)
flow<- coredata(q_weekly)
class(flow)
discharge <-flow[,1]
date <- index(q_weekly)
q_obs_weekly<-tibble(date,discharge)
save(q_obs_weekly, file="/work/OVERFLOW/RCR/sim53/q_obs_w.RData")
load(file="/work/OVERFLOW/RCR/sim53/q_obs_w.RData")



obs<-q_obs
obs$bacteria<-"NA"
df<-left_join(obs,bac_obs,by="date")
df <-subset(df,select=-c(3))
colnames(df)[3]<-"bacteria"
obs_daily <-df
save(obs_daily, file="/work/OVERFLOW/RCR/sim53/obs_daily.RData")
load(file="/work/OVERFLOW/RCR/sim53/obs_daily.RData")



bac_daily <-subset(obs_daily,select=c(1,3))
bac_daily
data <- as.xts(bac_daily$bacteria,order.by=as.Date(bac_daily$date))
bac_weekly<- apply.weekly(data,mean,na.rm=TRUE)
bacteria<- coredata(bac_weekly)
bacteria <-bacteria[,1]
date <- index(bac_weekly)
bac_obs_weekly<-tibble(date,bacteria)
save(bac_obs_weekly, file="/work/OVERFLOW/RCR/sim53/bac_obs_w.RData")
load(file="/work/OVERFLOW/RCR/sim53/bac_obs_w.RData")



################################
#####load simulation data########
################################
#convert daily sims to weekly
load("/work/OVERFLOW/RCR/sim55/bac_cal5.RData")
data3 <- as.xts(bac_cal_output$simulation$bac_out,order.by=as.Date(bac_cal_output$simulation$bac_out$date))
bac_cal<- apply.weekly(data3,mean, na.rm=TRUE)
bac_cal_w<- coredata(bac_cal) #xts function
bac_cal_w <-as_tibble(bac_cal_w)
bac_cal_w$date<-index(bac_cal)
save(bac_cal_w, file="/work/OVERFLOW/RCR/sim55/bac_cal5_w.RData")
load(file="/work/OVERFLOW/RCR/sim55/bac_cal5_w.RData")


nse_bac_w <- right_join(bac_cal_w,bac_obs_weekly,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs_weekly$bacteria))

sort(nse_bac_w, decreasing = T) %>% enframe()




