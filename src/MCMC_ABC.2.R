
load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal13.RData')


#6.Use the truncated normal distributions from 4) to set up the next round of simulations. Simulate with these inputs. For each simulation, calculate the average_nse. We will only keep individual simulations if the average_nse is higher that the first_quartile_average_nse (calculated in step 5) from the last set of simulations. 


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



load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/q_obs2.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/flux_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs2.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal13.RData')


sim_bac<- bac_cal1$simulation$bac_out
sim_q <- bac_cal1$simulation$q_out
sim_pars <- bac_cal1$parameter$values



nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
  select(-date) %>% select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs$bacteria))

sort(nse_bac, decreasing = T) %>% enframe()

nse_q <- right_join(sim_q,q_obs,by="date") %>%
  select(-date) %>% select(-discharge) %>%
  map_dbl(., ~NSE(.x, q_obs$discharge))

sort(nse_q, decreasing = T) %>% enframe()


flux_sim <- sim_bac[c(97:167),c(-1)]*
  sim_q  [c(97:167), c(-1)]*10^4
nse_flux <- flux_sim %>%
  map_dbl(., ~NSE(.x, flux_obs[,1]))

sort(nse_flux, decreasing = T) %>% enframe()


nse_average <- matrix(data=NA, nrow=50000, ncol=1)
for(i in 1:50000){
  print(i)
  nse_average[i] <- mean(c(nse_bac[i], nse_q[i], nse_flux[i]))
}
sort(nse_average, decreasing = T) %>% enframe()



run<-c(1:50000)

nse_mean <-cbind(run,nse_average)
colnames(nse_mean)<-c("run","nse_average")


valid.n <-head(which(nse_mean[,2]>-110.4801),n=10000)
sim_pars_2 <-sim_pars[valid.n,]

nse_ave.2 <-nse_average[valid.n]

#8.	Calculate the updated unweighted kernel densities based on these new 10k simulations and fit to the normal distribution, truncate at the original range limits for each parameter.
kde_mcabc <- sim_pars.2 %>% 
  gather(key = "par", value = "parameter_range")

ggplot(data = kde_mcabc) +
  geom_density(aes(x = parameter_range)) +
  facet_wrap(.~par, nrow=5, scales = "free") +
  theme_bw()
#ggsave("/home/hwu/wu_redcedar2/graphics/kde_mcabc.2.pdf")
# ###
#9.	Now use these new 10k simulations to calculate the updated first_quartile_average_nse, this will be the average of the 2500 and 2501st highest average_nse.
mean(nse_ave.2[2500:2501])
#[1] -51.68906


##################
################
#############
library(fitdistrplus)
library(truncnorm)

fln_CN2 <- fitdist(sim_pars_2$CN2, "norm")
fln_GWQMN <- fitdist(sim_pars_2$GWQMN, "norm")
fln_ALPHA_BNK <- fitdist(sim_pars_2$ALPHA_BNK, "norm")
fln_CH_K2 <- fitdist(sim_pars_2$CH_K2, "norm")
fln_CH_N2 <- fitdist(sim_pars_2$CH_N2, "norm")
fln_TRNSRCH <- fitdist(sim_pars_2$TRNSRCH, "norm")
fln_CH_N1 <- fitdist(sim_pars_2$CH_N1, "norm")
fln_CH_K1 <- fitdist(sim_pars_2$CH_K1, "norm")
fln_RCHRG_DP <- fitdist(sim_pars_2$RCHRG_DP, "norm")
fln_SFTMP <- fitdist(sim_pars_2$SFTMP, "norm")
fln_SMTMP <- fitdist(sim_pars_2$SMTMP, "norm")
fln_DEP_IMP <- fitdist(sim_pars_2$DEP_IMP, "norm")
fln_DDRAIN <- fitdist(sim_pars_2$DDRAIN, "norm")
fln_GDRAIN <- fitdist(sim_pars_2$GDRAIN, "norm")
fln_BACTKDQ <- fitdist(sim_pars_2$BACTKDQ, "norm")
fln_BACT_SWF<- fitdist(sim_pars_2$BACT_SWF, "norm")
fln_THBACT <- fitdist(sim_pars_2$THBACT, "norm")
fln_WDPRCH <- fitdist(sim_pars_2$WDPRCH, "norm")




 # fln_CN2
# SCS runoff curve number f
 # estimate   Std. Error
 # mean -0.02112694 0.0013701323
 # sd    0.13701323 0.0009685976
CN2<-rtruncnorm(50000, min(sim_pars_2$CN2), max(sim_pars_2$CN2), +
                  mean = -0.02112694, sd =  0.13701323)

# fln_GWQMN
# estimate  Std. Error
# mean 0.5249057 0.005570374
# sd   0.5570374 0.003938792
GWQMN<-rtruncnorm(50000, min(sim_pars_2$GWQMN), max(sim_pars_2$GWQMN), mean =0.5249057, sd = 0.5570374)


# fln_CH_N1 
# estimate   Std. Error
# mean 0.10006564 0.0002335225
# sd   0.02335225 0.0001637648
CH_N1<-rtruncnorm(50000, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = 0.10006564, sd = 0.02335225)


# fln_ALPHA_BNK
# Baseflow alpha factor for bank storage
# estimate  Std. Error
# mean 0.5003083 0.002359214
# sd   0.2359214 0.001668081
ALPHA_BNK<-rtruncnorm(50000, min(sim_pars_2$ALPHA_BNK), max(sim_pars_2$ALPHA_BNK), mean = 0.5003083, sd = 0.2359214)


# fln_CH_K2
# Effective hydraulic conductivity in main channel alluvium
# estimate Std. Error
# mean 26.54609 0.11625638
# sd   11.62564 0.08220568
CH_K2<-rtruncnorm(50000, min(sim_pars_2$CH_K2), max(sim_pars_2$CH_K2), mean = 26.54609, sd = 11.62564)


# fln_CH_N2
# estimate   Std. Error
# mean 0.10150102 0.0002367507
# sd   0.02367507 0.0001660660
CH_N2<-rtruncnorm(50000, min(sim_pars_2$CH_N2), max(sim_pars_2$CH_N2), mean = 0.10150102, sd = 0.02367507)


# fln_TRNSRCH
# estimate   Std. Error
# mean 0.16563411 0.0006859082
# sd   0.06859082 0.0004845465
TRNSRCH<-rtruncnorm(50000, min(sim_pars_2$TRNSRCH), max(sim_pars_2$TRNSRCH), mean = 0.16563411, sd = 0.06859082)

# fln_CH_N1
# estimate   Std. Error
# mean 0.10006564 0.0002335225
# sd   0.02335225 0.0001637648
CH_N1<-rtruncnorm(50000, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = 0.10006564, sd = 0.02335225)

# fln_CH_K1
# estimate Std. Error
# mean 159.28917  0.6879786
# sd    68.79785  0.4864743
CH_K1<-rtruncnorm(50000, min(sim_pars_2$CH_K1), max(sim_pars_2$CH_K1), mean = 159.28917, sd = 68.79785)


# fln_RCHRG_DP
# estimate  Std. Error
# mean 0.5129635 0.002378551
# sd   0.2378551 0.001681756
RCHRG_DP<-rtruncnorm(50000, min(sim_pars_2$RCHRG_DP), max(sim_pars_2$RCHRG_DP), mean = 0.5129635, sd = 0.2378551)


# fln_SFTMP
# estimate Std. Error
# mean 0.06362268 0.02349272
# sd   2.34927153 0.01661184
SFTMP<-rtruncnorm(50000, min(sim_pars_2$SFTMP), max(sim_pars_2$SFTMP), mean = 0.06362268, sd = 2.34927153)

# fln_SMTMP
# estimate Std. Error
# mean 0.1168178 0.02352965
# sd   2.3529649 0.01663796
SMTMP<-rtruncnorm(50000, min(sim_pars_2$SMTMP), max(sim_pars_2$SMTMP), mean = 0.1168178, sd = 2.3529649)


# fln_DEP_IMP
# estimate Std. Error
# mean 3505.098  13.291228
# sd   1329.247   9.398318
DEP_IMP<-rtruncnorm(50000, min(sim_pars_2$DEP_IMP), max(sim_pars_2$DEP_IMP), mean = 3505.098, sd = 1329.247)


# fln_DDRAIN
# estimate Std. Error
# mean 925.1803   4.687875
# sd   468.7940   3.314828
DDRAIN<-rtruncnorm(50000, min(sim_pars_2$DDRAIN), max(sim_pars_2$DDRAIN), mean = 925.1803, sd = 468.7940)


# fln_GDRAIN
# estimate Std. Error
# mean 54.10903  0.2317882
# sd   23.17882  0.1638990
GDRAIN<-rtruncnorm(50000, min(sim_pars_2$GDRAIN), max(sim_pars_2$GDRAIN), mean = 54.10903, sd = 23.17882)

# fln_BACTKDQ
# estimate Std. Error
# mean 280.9932  1.0867222
# sd   108.6722  0.7684278
BACTKDQ<-rtruncnorm(50000, min(sim_pars_2$BACTKDQ), max(sim_pars_2$BACTKDQ), mean = 280.9932, sd = 108.6722)


# fln_BACT_SWF
# estimate  Std. Error
# mean 0.4526996 0.002327963
# sd   0.2327963 0.001645982
BACT_SWF<-rtruncnorm(50000, min(sim_pars_2$BACT_SWF), max(sim_pars_2$BACT_SWF), mean = 0.4526996, sd = 0.2327963)

# fln_THBACT
# estimate Std. Error
# mean 3.818817 0.02322536
# sd   2.322536 0.01642280
THBACT<-rtruncnorm(50000, min(sim_pars_2$THBACT), max(sim_pars_2$THBACT), mean = 3.818817, sd = 2.322536)

# fln_WDPRCH
# estimate  Std. Error
# mean 0.5255286 0.002302976
# sd   0.2302976 0.001628312
WDPRCH<-rtruncnorm(50000, min(sim_pars_2$WDPRCH), max(sim_pars_2$WDPRCH), mean = 0.5255286, sd = 0.2302976)



load(file='/work/OVERFLOW/RCR/MSU/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')


path <- "/work/OVERFLOW/RCR/calibration/MSU"


pars <- tibble(
  "CN2.mgt|change = relchg"= CN2,
  "GWQMN.gw|change = relchg" = GWQMN,
  "ALPHA_BNK.rte|change = absval" = ALPHA_BNK,
  "CH_K2.rte|change = absval" = CH_K2,
  "CH_N2.rte|change = absval" = CH_N2,
  "TRNSRCH.bsn|change = absval" = TRNSRCH,
  "CH_N1.sub|change = absval" = CH_N1,
  "CH_K1.sub|change = absval" = CH_K1,
  "RCHRG_DP.gw|change = absval" = RCHRG_DP,
  "SFTMP.bsn|change = absval"= SFTMP,
  "SMTMP.bsn|change = absval"= SMTMP,
  "DEP_IMP.hru|change = absval"= DEP_IMP,
  "DDRAIN.mgt|change = absval"= DDRAIN,
  "GDRAIN.mgt|change = absval"= GDRAIN,
  "BACTKDQ.bsn|change = absval" = BACTKDQ,
  "BACT_SWF.bsn|change = absval" = BACT_SWF,
  "THBACT.bsn|change = absval"= THBACT,
  "WDPRCH.bsn|change = absval"= WDPRCH)

bac_cal1 <- run_swat2012(project_path = path,
                         output = list(q_out = define_output(file = "rch",
                                                             variable = "FLOW_OUT",
                                                             unit = 4),
                                       bac_out = define_output(file = "rch",
                                                               variable = "BACTP_OUT",
                                                               unit = 4)),
                         parameter = pars,
                         start_date = "2011-01-01",
                         end_date = "2013-12-31",
                         years_skip = 2,
                         n_thread = 32)


save(bac_cal1, file='/work/OVERFLOW/RCR/calibration/MSU/bac_cal14.RData')



