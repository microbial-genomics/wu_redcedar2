
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

median_score = -17.21075
valid.n <-head(which(nse_mean[,2]>median_score),n=5000)
sim_pars_2 <-sim_pars[valid.n,]

nse_ave.2 <-nse_average[valid.n]

#8.	Calculate the updated unweighted kernel densities based on these new 10k simulations and fit to the normal distribution, truncate at the original range limits for each parameter.
kde_mcabc <- sim_pars_2 %>% 
  gather(key = "par", value = "parameter_range")

ggplot(data = kde_mcabc) +
  geom_density(aes(x = parameter_range)) +
  facet_wrap(.~par, nrow=5, scales = "free") +
  theme_bw()
#ggsave("/home/hwu/wu_redcedar2/graphics/kde_mcabc.2.pdf")
# ###
#9.	Now use these new 10k simulations to calculate the updated first_quartile_average_nse, this will be the average of the 2500 and 2501st highest average_nse.
median_score <- mean(nse_ave.2[1250:1251])
#[1] -5.864447


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
# estimate  Std. Error
# mean -0.04183794 0.001904808
# sd    0.13469026 0.001346568

n_sims =100000
CN2_mean <- fln_CN2$estimate[1]
CN2_sd <- fln_CN2$estimate[2]
CN2<-rtruncnorm(n_sims, min(sim_pars_2$CN2), max(sim_pars_2$CN2),mean = CN2_mean, sd =  CN2_sd)

# fln_GWQMN
# #     estimate  Std. Error
# mean 0.4499712 0.007905126
# sd   0.5589768 0.005589688
GWQMN_mean <- fln_GWQMN$estimate[1]
GWQMN_sd <- fln_GWQMN$estimate[2]
GWQMN<-rtruncnorm(n_sims, min(sim_pars_2$GWQMN), max(sim_pars_2$GWQMN), mean =GWQMN_mean, sd = GWQMN_sd)


# fln_CH_N1 
# estimate   Std. Error
# mean 0.09991472 0.0003314811
# sd   0.02343925 0.0002324755
CH_N1_mean <- fln_CH_N1$estimate[1]
CH_N1_sd <- fln_CH_N1$estimate[2]
CH_N1<-rtruncnorm(n_sims, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = CH_N1_mean, sd = CH_N1_sd)


# fln_ALPHA_BNK
# Baseflow alpha factor for bank storage
# estimate  Std. Error
# mean 0.4911617 0.003369632
# sd   0.2382690 0.002382501
ALPHA_BNK_mean <- fln_ALPHA_BNK$estimate[1]
ALPHA_BNK_sd <- fln_ALPHA_BNK$estimate[2]
ALPHA_BNK<-rtruncnorm(n_sims, min(sim_pars_2$ALPHA_BNK), max(sim_pars_2$ALPHA_BNK), mean = ALPHA_BNK_mean, sd = ALPHA_BNK_sd)


# fln_CH_K2
# Effective hydraulic conductivity in main channel alluvium
# estimate Std. Error
# mean 27.90991  0.1635768
# sd   11.56662  0.1156662
CH_K2_mean <- fln_CH_K2$estimate[1]
CH_K2_sd <- fln_CH_K2$estimate[2]
CH_K2<-rtruncnorm(n_sims, min(sim_pars_2$CH_K2), max(sim_pars_2$CH_K2), mean = CH_K2_mean, sd = CH_K2_sd)


# fln_CH_N2
# estimate   Std. Error
# mean 0.10415042 0.0003312391
# sd   0.02342214 0.0002323030
CH_N2_mean <- fln_CH_N2$estimate[1]
CH_N2_sd <- fln_CH_N2$estimate[2]
CH_N2<-rtruncnorm(n_sims, min(sim_pars_2$CH_N2), max(sim_pars_2$CH_N2), mean = CH_N2_mean, sd = CH_N2_sd)


# fln_TRNSRCH
# estimate   Std. Error
# mean 0.17602745 0.0009794462
# sd   0.06925731 0.0006919234
TRNSRCH_mean <- fln_TRNSRCH$estimate[1]
TRNSRCH_sd <- fln_TRNSRCH$estimate[2]
TRNSRCH<-rtruncnorm(n_sims, min(sim_pars_2$TRNSRCH), max(sim_pars_2$TRNSRCH), mean = TRNSRCH_mean, sd = TRNSRCH_sd)

# fln_CH_N1
# estimate   Std. Error
# mean 0.09991472 0.0003314811
# sd   0.02343925 0.0002324755
CH_N1_mean <- fln_CH_N1$estimate[1]
CH_N1_sd <- fln_CH_N1$estimate[2]
CH_N1<-rtruncnorm(n_sims, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = CH_N1_mean, sd = CH_N1_sd)

# fln_CH_K1
# estimate Std. Error
# mean 158.93894  0.9599238
# sd    67.87685  0.6787685
CH_K1_mean <- fln_CH_K1$estimate[1]
CH_K1_sd <- fln_CH_K1$estimate[2]
CH_K1<-rtruncnorm(n_sims, min(sim_pars_2$CH_K1), max(sim_pars_2$CH_K1), mean = CH_K1_mean, sd = CH_K1_sd)


# fln_RCHRG_DP
# estimate  Std. Error
# mean 0.4808134 0.003439768
# sd   0.2432284 0.002432099
RCHRG_DP_mean <- fln_RCHRG_DP$estimate[1]
RCHRG_DP_sd <- fln_RCHRG_DP$estimate[2]
RCHRG_DP<-rtruncnorm(n_sims, min(sim_pars_2$RCHRG_DP), max(sim_pars_2$RCHRG_DP), mean = RCHRG_DP_mean, sd = RCHRG_DP_sd)


# fln_SFTMP
# estimate Std. Error
# mean 0.1230878 0.03301559
# sd   2.3345547 0.02334553
SFTMP_mean <- fln_SFTMP$estimate[1]
SFTMP_sd <- fln_SFTMP$estimate[2]
SFTMP<-rtruncnorm(n_sims, min(sim_pars_2$SFTMP), max(sim_pars_2$SFTMP), mean = SFTMP_mean, sd = SFTMP_sd)

# fln_SMTMP
# estimate Std. Error
# mean 0.1653056 0.03339191
# sd   2.3611645 0.02361163
SMTMP_mean <- fln_SMTMP$estimate[1]
SMTMP_sd <- fln_SMTMP$estimate[2]
SMTMP<-rtruncnorm(n_sims, min(sim_pars_2$SMTMP), max(sim_pars_2$SMTMP), mean = SMTMP_mean, sd = SMTMP_sd)


# fln_DEP_IMP
# estimate Std. Error
# mean 3667.422   18.83298
# sd   1331.808   13.31908
DEP_IMP_mean <- fln_DEP_IMP$estimate[1]
DEP_IMP_sd <- fln_DEP_IMP$estimate[2]
DEP_IMP<-rtruncnorm(n_sims, min(sim_pars_2$DEP_IMP), max(sim_pars_2$DEP_IMP), mean = DEP_IMP_mean, sd = DEP_IMP_sd)


# fln_DDRAIN
# estimate Std. Error
# mean 902.7242   6.710620
# sd   474.4889   4.744931
DDRAIN_mean <- fln_DDRAIN$estimate[1]
DDRAIN_sd <- fln_DDRAIN$estimate[2]
DDRAIN<-rtruncnorm(n_sims, min(sim_pars_2$DDRAIN), max(sim_pars_2$DDRAIN), mean = DDRAIN_mean, sd = DDRAIN_sd)


# fln_GDRAIN
# estimate Std. Error
# mean 53.36850  0.3311553
# sd   23.41622  0.2341622
GDRAIn_mean <- fln_GDRAIN$estimate[1]
GDRAIn_sd <- fln_GDRAIN$estimate[2]
GDRAIN<-rtruncnorm(n_sims, min(sim_pars_2$GDRAIN), max(sim_pars_2$GDRAIN), mean = GDRAIn_mean, sd = GDRAIn_sd)

# fln_BACTKDQ
# estimate Std. Error
# mean 284.6590   1.550853
# sd   109.6619   1.096619
BACTKDQ_mean <- fln_BACTKDQ$estimate[1]
BACTKDQ_sd <- fln_BACTKDQ$estimate[2]
BACTKDQ<-rtruncnorm(n_sims, min(sim_pars_2$BACTKDQ), max(sim_pars_2$BACTKDQ), mean = BACTKDQ_mean, sd = BACTKDQ_sd)


# fln_BACT_SWF
# estimate  Std. Error
# mean 0.4287223 0.003284534
# sd   0.2322516 0.002322323
BACT_SWF_mean <- fln_BACT_SWF$estimate[1]
BACT_SWF_sd <- fln_BACT_SWF$estimate[2]
BACT_SWF<-rtruncnorm(n_sims, min(sim_pars_2$BACT_SWF), max(sim_pars_2$BACT_SWF), mean = BACT_SWF_mean, sd = BACT_SWF_sd)

# fln_THBACT
# estimate Std. Error
# mean 1.384695 0.02367388
# sd   1.673996 0.01673994
THBACT_mean <- fln_THBACT$estimate[1]
THBACT_sd <- fln_THBACT$estimate[2]
THBACT<-rtruncnorm(n_sims, min(sim_pars_2$THBACT), max(sim_pars_2$THBACT), mean = THBACT_mean, sd = THBACT_sd)

# fln_WDPRCH
# estimate  Std. Error
# mean 0.5643253 0.003093901
# sd   0.2187718 0.002187513
WDPRCH_mean <- fln_WDPRCH$estimate[1]
WDPRCH_sd <- fln_WDPRCH$estimate[2]
WDPRCH<-rtruncnorm(n_sims, min(sim_pars_2$WDPRCH), max(sim_pars_2$WDPRCH), mean = WDPRCH_mean, sd = WDPRCH_sd)



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



load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal14.RData')


