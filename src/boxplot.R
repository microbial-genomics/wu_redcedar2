######load modules in HPC terminal Bash
###############################
module load intel/19.0.5
module load R/3.6.2
module load geos/3.8.0
module load gdal-2.4.3/intel-19.0
module load proj-5.2.0/intel-19.0
module load udunits-2.2.26/intel-19.0
module load rstudio




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
require(data.table) # load it






###daily bacteria concentration
conc<-bac_cal_output$simulation$bac_out
class(bac_cal_output)
conc<-as.data.table(conc)
class(conc)
date<-conc$date
class(date)
length(date)

# date<-as.vector(conc[,1])
# length(date)
# nrow(date)
# dim(date)
# class(date)
conc2<-conc[,c(2:1001)]
#conc2
#dim(conc2) #3865*1000
conc3<-t(conc2)
conc3
dim(conc3)
colnames(conc3)<-as.Date(date,"%Y%m%d")
class(colnames(conc3))
class(colnames(conc3))<-"Date"
conc3
head(conc3)
colnames(conc3)
class(conc3)

day1=conc3[,1]
q1<-quantile(conc3[,1], c(0.001, 0.023, 0.159, 0.5, 0.841, 0.977, 0.99))
q1

result <- data.frame(matrix(nrow = 3865, ncol = 8))
colnames(result) <- c("i", "first_q","sec_q","third_q","four_q","fif_q","six_q","sev_q")
i=1
for (i in 1:3865) {
  first_q <- quantile(conc3[,i], c(0.001))
  sec_q <- quantile(conc3[,i], c(0.023))
  third_q <- quantile(conc3[,i], c(0.159))
  four_q <- quantile(conc3[,i], c(0.5))
  fif_q <- quantile(conc3[,i], c(0.841))
  six_q <- quantile(conc3[,i], c(0.977))
  sev_q <- quantile(conc3[,i], c(0.99))
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



result2[result2<1] <-1
result2

range(result2$first_q)

ggplot() +
  geom_line(data = result2[c(2193:2793),], aes(x = date, y = log10(first_q), colour = "first_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date, y = log10(sec_q), colour = "sec_q")) +
# geom_line(data = result2, aes(x = date, y = log10(third_q), colour = "third_q")) +
  geom_line(data = result2[c(2193:2793),], aes(x = date, y = log10(four_q), colour = "four_q")) +
  #geom_line(data = result2, aes(x = date,y = log10(fif_q), colour = "fif_q")) +
  # geom_line(data = result2[c(2193:2793),], aes(x = date,y = log10(six_q), colour = "six_q")) +
  geom_line(data = result2[c(2193:2793),], aes(x = date,y = log10(sev_q), colour = "sev_q")) +
  geom_point(data =bac_obs[c(155:180),], aes(x=date, y=log10(bacteria), colour = "bacteria"))+
  scale_colour_manual("", 
                      breaks = c("first_q",   "four_q",  "sev_q",  "bacteria"),
                      values =c("grey70", "#CC79A7", "grey30",  "red")) +
  xlab(" ") +
  scale_y_continuous("Concentration (log10)") + 
  labs(title="Daily_concentration")+
  theme_bw()

