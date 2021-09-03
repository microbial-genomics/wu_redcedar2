
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

load ("/work/OVERFLOW/RCR/sim55/bac_cal7.RData")
iter <- 7

nse_bac <- calculate_nse_bac(iter, bac_cal_output, bac_obs)
nse_q <- calculate_nse_q(iter, bac_cal_output, q_obs)
nse_flux <- calculate_nse_flux(iter, bac_cal_output, flux_obs)


######################################################################
############### bacteria simulation and observation plot##############
#####################################################################
##load bac_obs.RData
##load bac_cal*.RData(bac_cal_output) this is the file generated from swatplusr package
##calculate nse_bac, use abc_functions.R
##identify iter
nse_bac <- calculate_nse_bac(iter, bac_cal_output, bac_obs)
sort(nse_bac, decreasing = T) %>% enframe()
#get the run_*** number
bac_plot <-right_join(bac_cal_output$simulation$bac_out,bac_obs,by="date")%>%
  dplyr::select(date, run_1065)%>%
left_join(., bac_obs, by ="date")%>%
  rename (bac_obs=bacteria)%>%
  gather(., key= "variable", value="bacteria",-date)


ggplot(data = bac_plot)+
  geom_line(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  geom_point(aes(x = date, y = bacteria, col = variable, lty = variable)) +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  scale_color_manual(values = c("black", "tomato3")) +
  theme_bw()
ggsave("bac_sim_obs_gen*.pdf")




#################################################################
#######plot log bac simulation and log bac observation ##############
#################################################################
bac_sim0 <-bac_cal_output$simulation$bac_out[,c(1,1065+1)]
l_run_1065<-log10(bac_sim0$run_1065)
bac_sim <-cbind(bac_sim0[,1],l_run_1065)
l_bac_sim <-bac_sim
###need to refine this code#####
load(file="/work/OVERFLOW/RCR/sim52.3/l_bac_obs.RData")

l_bac_plot <-right_join(l_bac_sim,l_bac_obs,by="date")%>%
  dplyr::select(date, l_run_1065)%>%
left_join(.,l_bac_obs, by ="date")%>%
  rename (l_bac=l_bacteria)%>%
  gather(., key= "variable", value="l_bacteria",-date)


ggplot(data =l_bac_plot)+
  geom_line(aes(x = date, y =l_bacteria, col = variable, lty = variable)) +
  geom_point(aes(x = date, y =l_bacteria, col = variable, lty = variable)) +
  scale_color_manual(values = c("black", "tomato3")) +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  #ggtitle("sim55_generation5_bacnse=-0.0012")+
  theme_bw()


ggsave("bac_sim_obs_gen7.pdf")
ggsave("/home/hwu/wu_redcedar2/graphics/sim53/bac_sim_obs_gen7.pdf")


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


nse_bac_w <- right_join(bac_cal_w,bac_obs_weekly,by="date") %>%
  dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs_weekly$bacteria))

sort(nse_bac_w, decreasing = T) %>% enframe()
#get the run_*** number
bac_plot <-right_join(bac_cal_w,bac_obs_weekly,by="date")%>%
  dplyr::select(date, run_1065)%>%
  left_join(., bac_obs_weekly, by ="date")%>%
  rename (bac_obs_weekly=bacteria)%>%
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
  geom_line(mapping = aes(x = q_obs$date, y = q_obs$discharge), size = 0.7, color = "blue") +
  geom_point(mapping = aes(x = bac_obs$date, y = bac_obs$bacteria*0.002), color = "tomato3") +
  scale_x_date(name = "date",date_breaks = "1 year",date_labels = "%Y") +
  scale_y_continuous(name = "discharge (cms)",
                     sec.axis = sec_axis(~./0.002, name = "bacteria(MPN/100ml)")) +
  theme(
    axis.title.y = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "tomato3"))

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

q_plot <-bac_cal_output$simulation$q_out%>%dplyr::select(date, run_1065)%>%
  left_join(., q_obs, by ="date")%>% 
  rename (q_obs=discharge)%>%gather(., key= "variable", value="discharge",-date)

ggplot(data = q_plot) +
  geom_line(aes(x = date, y = discharge, col = variable, lty = variable)) +
  scale_x_date(name = "date", date_breaks = "1 year",date_labels = "%Y") +
  scale_color_manual(values = c("black", "tomato3")) +
  ggtitle("sim53.2")+
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

  

