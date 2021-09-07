
module load intel/19.0.5
module load R/3.6.2
module load geos/3.8.0
module load gdal-2.4.3/intel-19.0
module load proj-5.2.0/intel-19.0
module load udunits-2.2.26/intel-19.0




library(devtools)
library(dplyr)
library(dygraphs)
#library(fast) #error
library(forcats)
library(ggplot2)
library(hydroGOF)
library(lhs)
library(lubridate)
library(mapview)
library(plotly)
library(purrr)
library(sensitivity)
library(sf)
library(tibble)
library(tidyr)
library(fitdistrplus)
library(truncnorm)
library(xts)


base_dir <- file.path("/work", "OVERFLOW", "RCR", "sim55-weekly","sim55-nse-flow")
  data_in_dir <- base_dir
  graphics_dir <- base_dir
  src_dir <- base_dir
  project_path <- base_dir
  swat_path <- base_dir

source(file.path(src_dir, "abc_functions.R"))

load_observations()

load ("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/bac_cal14.RData")
iter <- 14


bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]


###### subset simulated flow to observed days and average by week
# extract only the observed days from the simulated daily output
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
#dim(bac_flows_all_days)
# adds date and flow fields also
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
#dim(bac_flows_daily_temp)
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
#dim(bac_flows_daily)
# then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
#head(colnames(bac_flows_daily)) #  date 
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[204,2214]
t(bac_flows_weekly)
#dim(bac_flows_weekly)
###### calculate simulated flux data for observed days and average by week
#dim(bac_sims_weekly)
#dim(bac_flows_weekly)
bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
#dim(bac_fluxes_weekly) #[204,2214]


###Flow
###### calculate various nses for daily data 
nse_flow_daily <- calculate_nse_flow_daily(iter, bac_cal_output, q_obs)
# calculate various modified nses for daily data
mnse_flow_daily <- calculate_mnse_flow_daily(iter, bac_cal_output, q_obs)
# calculate various nses for weekly data
nse_flow_weekly <- calculate_nse_flow_weekly(iter, bac_flows_weekly, flow_obs_weekly)
# calculate various modified nses for weekly data
mnse_flow_weekly <- calculate_mnse_flow_weekly(iter, bac_flows_weekly, flow_obs_weekly)


######################################################################
############### bacteria simulation and observation plot##############
#####################################################################
##load bac_obs.RData
##load bac_cal*.RData(bac_cal_output) this is the file generated from swatplusr package
##calculate nse_bac, use abc_functions.R
##identify iter
bac_sim <- bac_cal_output$simulation$bac_out

#get the run_*** number
plot <-right_join(bac_sim,bac_obs,by="date")%>%
  dplyr::select(date, run_0808)%>%
left_join(., bac_obs, by ="date")%>%
  rename (bac_obs=bacteria)%>%
  gather(., key= "variable", value="bacteria",-date)

ggplot(data = plot)+
  geom_line(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  geom_point(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  scale_color_manual(values = c("black", "tomato3")) +
  scale_y_log10()+
  theme_bw()+
  ggtitle("mnse=0.305,opt(conc. weekly)")
#ggsave("**********.pdf")


sapply(my.data, class)
str(my.data)



####count the non zero point in simulation output #####
bac <- bac_cal_output$simulation$bac_out$run_0001
length(which(bac!=0))


######################################################################
###############WEEKLY bacteria simulation and observation plot##############
#####################################################################
##load bac_obs_w.RData
##load bac_cal*.RData(bac_cal_output) this is the file generated from swatplusr package
##calculate nse_bac, use abc_functions.R
##identify iter
load("/work/OVERFLOW/RCR/sim55/bac_cal13.RData")
data3 <- as.xts(bac_cal_output$simulation$bac_out,order.by=as.Date(bac_cal_output$simulation$bac_out$date))
bac_cal<- apply.weekly(data3,mean, na.rm=TRUE)
bac_cal_w<- coredata(bac_cal)
bac_cal_w <-as_tibble(bac_cal_w)
bac_cal_w$date<-index(bac_cal)
save(bac_cal_w, file="/work/OVERFLOW/RCR/sim55/bac_cal13_w.RData")
load(file="/work/OVERFLOW/RCR/sim55/bac_cal13_w.RData")
load(file="/work/OVERFLOW/RCR/sim53/bac_obs_w.RData")

sort(nse_bac, decreasing = T) %>% enframe()

#get the run_*** number
#run_0808
bac_sim <- bac_cal_output$simulation$bac_out
bac_plot <-right_join(bac_sim,bac_obs,by="date")%>%
  dplyr::select(date, run_0808)%>%
  left_join(.,bac_obs, by ="date")%>%
  rename (bac=bacteria)%>%
  gather(., key= "variable", value="bacteria",-date)


ggplot(data = bac_plot)+
  geom_line(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  geom_point(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  scale_color_manual(values = c("black", "tomato3")) +
  theme_bw()
ggsave("bac_sim_obs_gen*.pdf")

#################################################################
#######WEEKLY plot log bac simulation and log bac observation ##############
#################################################################
sort(nse_bac_w, decreasing = T) %>% enframe()

bac_sim0 <-bac_cal_w[,c(1,1065+1)]
l_run_1065<-log10(bac_sim0$run_1065)
l_bac_sim_w <-tibble(bac_sim0[,1],l_run_1065)
###need to refine this code#####
load(file="/work/OVERFLOW/RCR/sim52.3/l_bac_obs.RData")
####generate log bacteria observation RData
date<-bac_obs_weekly$date
bacteria<-log10(bac_obs_weekly$bacteria)
l_bac_obs_w <-tibble(date, bacteria)
save(l_bac_obs_w,file="/work/OVERFLOW/RCR/sim53/l_bac_obs_w.RData")
load(file="/work/OVERFLOW/RCR/sim53/l_bac_obs_w.RData")

l_bac_plot <-right_join(l_bac_sim_w,l_bac_obs_w,by="date")%>%
  dplyr::select(date, l_run_1065)%>%
  left_join(.,l_bac_obs_w, by ="date")%>%
  rename (l_bac_obs=bacteria)%>%
  gather(., key= "variable", value="bacteria",-date)


ggplot(data =l_bac_plot)+
  geom_line(aes(x = date, y =bacteria, col = variable, lty = variable)) +
  geom_point(aes(x = date, y =bacteria, col = variable, lty = variable)) +
  scale_color_manual(values = c("black", "tomato3")) +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
 # ggtitle("weekly_gene5_nse=0.015")+
  theme_bw()

########################################################
#########Observation bacteria with discharge plot##################
#######################################################
###method1
ggplot() +
  geom_line(mapping = aes(x = q_obs$date, y = q_obs$discharge), size = 0.7, color = "grey50") +
  geom_point(mapping = aes(x = bac_obs$date, y = bac_obs$bacteria*0.01), color = "tomato3") +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  scale_y_continuous(name = "discharge (cms)",
                     sec.axis = sec_axis(~./0.01, name = "bacteria(MPN/100ml)")) +
  theme(
    axis.title.y = element_text(color = "grey70"),
    axis.title.y.right = element_text(color = "tomato3"))+
  ggtitle("observation bac vs. flow")

##method2

ggplot() +
  geom_line(mapping = aes(x = q_obs$date, y = q_obs$discharge, colour="discharge"), size = 0.7) +
  geom_point(mapping = aes(x = bac_obs$date, y = bac_obs$bacteria*0.02, colour= "bacteria")) +
  scale_y_continuous(sec.axis = sec_axis(~./0.02, name = "bacteria(MPN/100ml)")) +
  scale_color_manual(values = c("tomato3","blue"))+
  labs(y = "discharge (cms)",
       x= "Date",
       colour="Parameter") +
  theme(legend.position = c(0.8,0.9))


######################################################################
############### discharge simulation and observation plot##############
####################################################################

sort(nse_q, decreasing = T) %>% enframe()
#run_1261 nse=0.537
#run_0763 nse=0.464

q_plot <-bac_cal_output$simulation$q_out%>%dplyr::select(date, run_0763)%>%
  left_join(., q_obs, by ="date")%>% 
  rename (q_obs=discharge)%>%gather(., key= "variable", value="discharge",-date)

ggplot(data = q_plot) +
  geom_line(aes(x = date, y = discharge, col = variable, lty = variable)) +
  scale_x_date(name = "date", date_breaks = "1 year",date_labels = "%Y") +
  scale_color_manual(values = c("black", "tomato3")) +
  labs(y = "discharge (cms)",
       x= "Date")+
  ggtitle("opt(weekly,flow); mnse=0.464")+
  ylim(0,100)+
  theme_bw()


ggsave("q_obs_sim_gene11.pdf")
ggsave("/home/hwu/wu_redcedar2/graphics/sim53/q_obs_sim_gene11.pdf")
######################################################################
############### flux simulation and observation plot##############
####################################################################
flux_sim <- sim_bac[c(97:167),c(-1)]*
  sim_q[c(97:167), c(-1)]*10^4

flux_plot <-cbind(flux_sim[,11805], flux_obs[,1])

ggplot() +
  geom_line(mapping = aes(x = bac_obs$date, y = flux_plot[,1]), size = 0.7, colour = "black") +
  geom_point(mapping = aes(x = bac_obs$date, y = flux_plot[,2]), colour = "green") +
  #scale_y_continuous(sec.axis = sec_axis(~., name = "flux")) +
  #scale_color_manual(values = c("tomato3","blue"))+
  labs(y = "flux",
       x= "Date",
       colour="Parameter") +
  theme(legend.position = c(0.8,0.9))
#########################################################
#########plot bacteria with discharge (simulation data)####################
###########################################################
###Method1
bac <- bac_cal_output$simulation$bac_out%>%dplyr::select(date, run_1065 )
names(bac)[2]<- paste("bacteria")

q  <- bac_cal_output$simulation$q_out%>%dplyr::select(date, run_1065)
names(q)[2] <- paste("discharge")


bac_q_plot <- bac %>% left_join(., q, by ="date" )

ggplot(bac_q_plot,aes(date,discharge)) +
  geom_line(aes(y = discharge), size =0.5, color = "blue") +
  geom_point(aes(y = bacteria/100), size = 0.3, color = "tomato3") +
  scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
  scale_y_continuous(name = "discharge (cms)", limits=c(0,30),
                     sec.axis = sec_axis(~(.*100), name = "bacteria (MPN/100ml)")) +
  theme(
    axis.title.y = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "tomato3"))

ggsave("q_bac_sim_gen10_02779.pdf")
ggsave("/home/hwu/wu_redcedar2/graphics/sim53/q_bac_sim_gen10_02779.pdf")

#

  
#################################################################
############precipitation with discharge plot###########################
######################################################################
precipitation <- read.table("pcp_obs.txt")[,2]
date0414 <- read.table("pcp_obs.txt")[,1]
date = as.Date(date0414, format="%m/%d/%Y")
pcp_obs <- tibble(date, precipitation)
save(pcp_obs, file ='/work/OVERFLOW/RCR/sim53/pcp_obs.RData')

###method 1
ggplot() +
  geom_line(mapping = aes(x = q_obs$date, y = q_obs$discharge), size = 0.5, 
            color = "cyan3") +
  geom_bar(mapping = aes(x = pcp_obs$date, y = pcp_obs$precipitation*10),
           stat = "identity", fill = "black",width=4) +
  scale_x_date(name = "date", date_breaks = "1 year",date_labels = "%y") +
  scale_y_continuous(name = "discharge",
                     sec.axis = sec_axis(~./10, name = "precipitation")) +
  theme(
    axis.title.y = element_text(color = "cyan3"),
    axis.title.y.right = element_text(color = "black"))

###method 2

obs<-cbind(q_obs,pcp_obs[2])

p <- ggplot(obs, aes(x=date))+
  geom_line(aes(y= discharge, colour = "discharge"))+
  geom_col(aes(y=precipitation*10, colour = "precipitation")) +
  scale_y_continuous(sec.axis =sec_axis(~./10, name="precipitation(mm)")) +
  scale_colour_manual(values = c("blue","black"))+
  labs(y = "discharge (cms)",
       x= "Date",
       colour="Parameter") +
  theme(legend.position = c(0.8,0.9))

p

##############################################
####################Boxplot nse###############
##############################################

####opt(weekly,mean)-daily data###############################
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/bac_cal14.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/bac_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/q_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/flux_obs.RData")

sim_bac <- bac_cal_output$simulation$bac_out
nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~mNSE(.x, bac_obs$bacteria))

sim_q <- bac_cal_output$simulation$q_out
nse_q <- right_join(sim_q,q_obs,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-discharge) %>%
  map_dbl(., ~mNSE(.x, q_obs$discharge))

date <-bac_cal_output$simulation$bac_out$date
flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
sim_flux<- cbind(date, flux_sim)
nse_flux <-  right_join(sim_flux, flux_obs, by = "date") %>%
  dplyr::select(-date) %>%dplyr::select(-flux) %>%
  map_dbl(., ~mNSE(.x, flux_obs$flux))

nse_mean <- rowMeans(cbind(nse_bac, nse_q, nse_flux))##select first 2000
nse_mean<-nse_mean[c(1:1500)]
sort(nse_mean, decreasing = T) %>% enframe()

####opt(conc,weekly)-daily_data###############################
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-conc/bac_cal24.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-conc/bac_obs.RData")
sim_bac <- bac_cal_output$simulation$bac_out
# merge simulated and observed bacteria concentrations, calculate nses for all sims
nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~mNSE(.x, bac_obs$bacteria))##select first 2000 #the same with log??
nse_bac<-nse_bac[c(1:1500)]
sort(nse_bac, decreasing = T) %>% enframe() 


######opt(flow,weekly)-daily-data################
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/bac_cal15.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/q_obs.RData")
sim_q <- bac_cal_output$simulation$q_out
# merge simulated and observed flows, calculate nses for all sims
nse_q <- right_join(sim_q,q_obs,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-discharge) %>%
  map_dbl(., ~mNSE(.x, q_obs$discharge))##select first 2000
nse_q<-nse_q[c(1:1500)]
sort(nse_q, decreasing = T) %>% enframe()

####opt(flux,weekly)-daily-data###############################
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flux/bac_cal12.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flux/flux_obs.RData")
sim_bac <- bac_cal_output$simulation$bac_out
sim_q <- bac_cal_output$simulation$q_out
# calculate the simulated fluxes from concs and flows for all sims
# flux_sim <- sim_bac[,c(-1)]*
#   sim_q[, c(-1)]*10^4
date <-bac_cal_output$simulation$bac_out$date
flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
sim_flux<- cbind(date, flux_sim)
#merge simulated and observed fluxes, calculate nses for all sims
nse_flux <-  right_join(sim_flux, flux_obs, by = "date") %>%
  dplyr::select(-date) %>%dplyr::select(-flux) %>%
  map_dbl(., ~mNSE(.x, flux_obs$flux))##select first 2000
nse_flux<-nse_flux[c(1:1500)]
sort(nse_flux, decreasing = T) %>% enframe() 

#####Boxplot###
nse_comp <- tibble(flow = nse_q,
                   bacteria = nse_bac,
                   flux =nse_flux,
                   mean= nse_mean)
nse_comp%>% 
  gather(key = "weekly", value = "mnse" ) %>% 
    ggplot(data = .) +
  geom_boxplot(aes(x = weekly, y = mnse), fill = "grey") +
  theme_bw()+
  ylim(0,0.5)+
  ggtitle("opt (weekly) first 1500 simulation daily data") 


###########################################
#############weekly_weekly data################
####opt(mean,weekly)-weekly data###########
###################################
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/bac_cal14.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/bac_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/q_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-mean/flux_obs.RData")
#bacteria
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
#dim(bac_sims_daily) #[336,2214]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]
#flow
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" |colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[204,2214]
bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
nse_flux <- mapply(mNSE, bac_fluxes_weekly, flux_obs_weekly)
nse_q <- mapply(mNSE, bac_flows_weekly, flow_obs_weekly)
nse_bac <- mapply(mNSE, bac_sims_weekly, bac_obs_weekly)
nse_mean <- rowMeans(cbind(nse_bac, nse_q, nse_flux))
nse_mean<-nse_mean[1:1500]

####opt(conc,weekly)-weekly data###########
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-conc/bac_cal24.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-conc/bac_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-conc/bac_obs_w.RData")

bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]

nse_bac <- mapply(mNSE, bac_sims_weekly, bac_obs_weekly)
#sort(nse_bac, decreasing = T) %>% enframe() 
nse_bac<-nse_bac[1:1500]

####opt(flow,weekly)-weekly data###########
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/bac_cal15.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/q_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flow/q_obs_w.RData")
flow_obs <- right_join(q_obs, bac_obs, by="date") #needed to reduce the number of flow observations
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                   colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
#dim(bac_flows_daily) #[336,2214]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[204,2214]
nse_q <- mapply(mNSE, bac_flows_weekly, flow_obs_weekly)##NSE has higher score than mNSE
#sort(nse_q, decreasing = T) %>% enframe() 
nse_q<-nse_q[1:1500]


####opt(flux,weekly)-weekly data###########
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flux/bac_cal12.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flux/flux_obs.RData")
load("/work/OVERFLOW/RCR/sim55-weekly/sim55-nse-flux/flux_obs_w.RData")
#bacteria
bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
#dim(bac_sims_daily) #[336,2214]
bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]
#flow
bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" |colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[204,2214]
bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
nse_flux <- mapply(mNSE, bac_fluxes_weekly, flux_obs_weekly)
#sort(nse_flux, decreasing = T) %>% enframe() 
nse_flux<-nse_flux[1:1500]


#############
###box-plot###
##############
nse_comp <- tibble(flow = nse_q,
                   bacteria = nse_bac,
                   flux =nse_flux,
                   mean= nse_mean)
nse_comp%>% 
  gather(key = "weekly", value = "mnse" ) %>% 
  ggplot(data = .) +
  geom_boxplot(aes(x = weekly, y = mnse), fill = "grey") +
  theme_bw()+
  ylim(0,0.5)+
  ggtitle("opt (weekly) first 1500 simulation weekly data") 





