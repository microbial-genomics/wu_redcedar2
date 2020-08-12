library(SWATplusR)
library(dplyr)
library(dygraphs)
library(ggplot2)
library(lubridate)
library(mapview)
library(plotly)
library(sf)
library(tibble)
library(tidyr)
library(purrr)
library(lhs)
library(hydroGOF)
library(fast)
library(forcats)
library(lubridate)
library(sensitivity)
library(Hmisc)




load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/q_obs2.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/flux_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs2.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal12.RData')

nse_bac <- right_join(bac_cal1$simulation$bac_out,bac_obs,by="date") %>%
  select(-date) %>% select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs$bacteria))
  sort(nse_cal1, decreasing = T) %>% enframe()

nse_q <- right_join(bac_cal1$simulation$q_out,q_obs,by="date") %>%
  select(-date) %>% select(-discharge) %>%
  map_dbl(., ~NSE(.x, q_obs$discharge))
  sort(nse_cal1_q, decreasing = T) %>% enframe()



flux_sim <- bac_cal1$simulation$bac_out[c(97:167),c(-1)]*
            bac_cal1$simulation$q_out [c(97:167), c(-1)]*10^4
nse_flux <- flux_sim %>%
     map_dbl(., ~NSE(.x, flux_obs[,1]))   
  sort(nse_cal1_flux, decreasing = T) %>% enframe()


###combine parameters and three nses
par.sim.nse <- cbind(par_sim,nse_q,nse_bac,nse_flux)
dim(par.sim.nse)
#[1] 50000    21




######################
### i-- parameter index
#### sim_par: 16288*18
### wt =NSIi/sum(nse)*1
range(nse_q)
length(nse_q[nse_q>0]) ###16288 simulation out of 50000 nse_q great than zero
par.sim.nse.2 <- par.sim.nse[order(par.sim.nse$nse_q)]
par.sim.nse.2 <- par.sim.nse[order(-nse_q),]
head(par.sim.nse.2,10)

par.sim.nse.3 <- par.sim.nse.2[c(1:16288),]
df <- par.sim.nse.3 
wt_q <- df$nse_q/(sum(df$nse_q))
quantiles.probs <- c(0.001,0.023,0.025,0.05,0.159,0.25,0.5,0.75,0.841,0.95,0.975,0.977,0.999)
nquantiles <- length(quantiles.probs)
par_sim_q <- array(data=NA, dim=c(18,nquantiles))

  for(i in 1:18){
    print(i)
   par_sim_q[i,] <- wtd.quantile(df[,i],weights=wt_q, probs=quantiles.probs,normwt=TRUE)  
  }
###
###
###assign colnames and rownames for the dataframe
colnames(par_sim_q) <-c("1%", "2.3%","2.5%","5%", "15.9%", "25%","50%","75%","84.1%","95%","97.5%","97.7%","99.9%")
rownames(par_sim_q) <- colnames(par_sim)

write.csv(par_sim_q, file ='/home/hwu/wu_redcedar2/data_out/weighted.pars.q.csv')



##########################################################################
range(nse_bac)
head(par.sim.nse[order(-nse_bac),],10)
par.sim.bac.nse <- par.sim.nse[order(-nse_bac),]
df.bac <- par.sim.bac.nse[c(1:16288),]
range(df.bac$nse_bac) ## -77.38714327  -0.07664198
colnames(df.bac)
nse_bac2<- df.bac[,20]+78 ## adjust the value to above zero
wt_bac<- nse_bac2/sum(nse_bac2)

par_sim_bac <- array(data=NA, dim=c(18,nquantiles))

  for(i in 1:18){
    print(i)
   par_sim_bac[i,] <- wtd.quantile(df.bac[,i],weights=wt_bac, probs=quantiles.probs,normwt=TRUE)  
  }


colnames(par_sim_bac) <-c("1%", "2.3%","2.5%","5%", "15.9%", "25%","50%","75%","84.1%","95%","97.5%","97.7%","99.9%")
rownames(par_sim_bac) <- colnames(par_sim)

write.csv(par_sim_bac, file ='/home/hwu/wu_redcedar2/data_out/weighted.pars.bac.csv')


############################################################################

range(nse_flux)
head(par.sim.nse[order(-nse_flux),],10)
par.sim.flux.nse <- par.sim.nse[order(-nse_flux),]
df.flux <- par.sim.bac.nse[c(1:16288),]
range(df.flux$nse_flux) ##  -638.2132142    0.3854417
colnames(df.flux)
nse_flux2<- df.flux[,21]+639 ## adjust the value to above zero
wt_flux<- nse_flux2/sum(nse_flux2)

par_sim_flux <- array(data=NA, dim=c(18,nquantiles))

  for(i in 1:18){
    print(i)
   par_sim_flux[i,] <- wtd.quantile(df.flux[,i],weights=wt_flux, probs=quantiles.probs,normwt=TRUE)  
  }


colnames(par_sim_flux) <-c("1%", "2.3%","2.5%","5%", "15.9%", "25%","50%","75%","84.1%","95%","97.5%","97.7%","99.9%")
rownames(par_sim_flux) <- colnames(par_sim)

write.csv(par_sim_flux, file ='/home/hwu/wu_redcedar2/data_out/weighted.pars.flux.csv')












