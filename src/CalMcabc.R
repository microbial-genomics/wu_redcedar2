

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

load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal12.RData')




# steps for the Approximate Bayesian Computation with sequential Monte Carlo simulation (ABC-MC):##############

#Step 1 You have already done 50k simulations from the uniform priors. We are going to take the first 10k (not the best 10k, do not sort) of these as the first generation of simulations. We are simply throwing away the last 40k and not using them for the official simulations for the manuscript. This is necessary to honor the assumptions behind ABC-MC.

sim_bac<- bac_cal1$simulation$bac_out[,1:10001]
sim_q <- bac_cal1$simulation$q_out[,1:10001]
sim_pars <- bac_cal1$parameter$values[1:10000,]

#Step 2 Calculate the average of the 3 nses for each of the 10k simulations: nse_average = mean(nse_conc, nse_flow, nse_ flux)

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


nse_average <- matrix(data=NA, nrow=10000, ncol=1)
	for(i in 1:10000){
		print(i)
  nse_average[i] <- mean(c(nse_bac[i], nse_q[i], nse_flux[i]))
					}
sort(nse_average, decreasing = T) %>% enframe()


#Step3.	Take the top 5k simulations based on the nse_average, discard the other 5K simulations.
#run<-colnames(sim_bac)

#run <- run[2:10001]
#nse_average <- nse_average
#nse_mean <-cbind(run,nse_average)
#top5k_sim <- nse_mean[order(-nse_average,run),][1:5000,]
#valid_sims <-c(top5k[,1])
#sim_bac_2 <-sim_bac[valid_sims]
#sim_q_2 <-sim_q[valid_sims]


run<-c(1:10000)

nse_average <- nse_average
nse_mean <-cbind(run,nse_average)
top2.5k_par <- nse_mean[order(-nse_average,run),][1:2500,]
median(top2.5k_par[,2])


valid_pars <-c(top2.5k_par[,1])
sim_pars_2 <- sim_pars[valid_pars,]

##Step 4 Use the 5k set of inputs associated with these top 5k simulations. Calculate the unweighted kernel densities using the kde package and fit to a normal distribution, truncate at the range limits for each parameter.


# kde_mcabc <- sim_pars_2 %>% 
#   #filter(nse > -10) %>%
#   gather(key = "par", value = "parameter_range")

#      ggplot(data = kde_mcabc) +
#    geom_density(aes(x = parameter_range)) +
#    facet_wrap(.~par, nrow=5, scales = "free") +
#     theme_bw()
# ggsave("/home/hwu/wu_redcedar2/graphics/kde_mcabc.pdf")
# ###


 ####.	

	
#5.	Find the median_average_nse of these 2.5k nse_averages, this will be the average of the 1225 and 1226th highest average_nse.
median(top2.5k_par[,2])
#[1] -17.21075



#6.	Use the truncated normal distributions from 4) to set up the next round of simulations. Simulate with these inputs. For each simulation, calculate the average_nse. Only keep the simulation if the average_nse is higher that the median_average_nse from the last set of simulations.
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




fln_CN2
# SCS runoff curve number f
#        estimate  Std. Error
#mean -0.03060313 0.003394195
#sd    0.16970975 0.002399683

CN2<-rtruncnorm(50000, min(sim_pars_2$CN2), max(sim_pars_2$CN2), mean = -0.03060313, sd = 0.16970975)


 fln_GWQMN
      # estimate  Std. Error
# mean 0.4915682 0.013951768
# sd   0.6975884 0.009865299

GWQMN<-rtruncnorm(50000, min(sim_pars_2$GWQMN), max(sim_pars_2$GWQMN), mean =0.4915682, sd = 0.6975884)


fln_CH_N1 
#       estimate   Std. Error
#mean 0.09968964 0.0005745893
#sd   0.02872947 0.0004040831

CH_N1<-rtruncnorm(50000, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = 0.09968964, sd = 0.02872947)

fln_ALPHA_BNK
# Baseflow alpha factor for bank storage
#       estimate  Std. Error
# mean 0.4906096 0.005878505
# sd   0.2939252 0.004156514
ALPHA_BNK<-rtruncnorm(50000, min(sim_pars_2$ALPHA_BNK), max(sim_pars_2$ALPHA_BNK), mean = 0.4906096, sd = 0.2939252)


fln_CH_K2
# Effective hydraulic conductivity in main channel alluvium
#      estimate Std. Error
# mean 27.43193  0.2815924
# sd   14.07962  0.1991159
CH_K2<-rtruncnorm(50000, min(sim_pars_2$CH_K2), max(sim_pars_2$CH_K2), mean = 27.43193, sd = 14.07962)

 fln_CH_N2
#        estimate   Std. Error
# mean 0.10236356 0.0005781796
# sd   0.02890898 0.0004066355
CH_N2<-rtruncnorm(50000, min(sim_pars_2$CH_N2), max(sim_pars_2$CH_N2), mean = 0.10236356, sd = 0.02890898)


 fln_TRNSRCH
#        estimate  Std. Error
# mean 0.16805179 0.001698463
# sd   0.08492317 0.001200246
TRNSRCH<-rtruncnorm(50000, min(sim_pars_2$TRNSRCH), max(sim_pars_2$TRNSRCH), mean = 0.16805179, sd = 0.08492317)

 fln_CH_N1
       # estimate   Std. Error
# mean 0.09968964 0.0005745893
# sd   0.02872947 0.0004040831
CH_N1<-rtruncnorm(50000, min(sim_pars_2$CH_N1), max(sim_pars_2$CH_N1), mean = 0.09968964, sd = 0.02872947)

 fln_CH_K1
      # estimate Std. Error
# mean 159.53393   1.693993
# sd    84.69961   1.197834
CH_K1<-rtruncnorm(50000, min(sim_pars_2$CH_K1), max(sim_pars_2$CH_K1), mean = 159.53393, sd = 84.69961)


fln_RCHRG_DP
      # estimate  Std. Error
# mean 0.5015587 0.005983016
# sd   0.2991508 0.004230418
RCHRG_DP<-rtruncnorm(50000, min(sim_pars_2$RCHRG_DP), max(sim_pars_2$RCHRG_DP), mean = 0.5015587, sd = 0.2991508)


 fln_SFTMP
      estimate Std. Error
# mean 0.1093467 0.05705058
# sd   2.8525289 0.04034083
SFTMP<-rtruncnorm(50000, min(sim_pars_2$SFTMP), max(sim_pars_2$SFTMP), mean = 0.1093467, sd = 2.8525289)


 fln_SMTMP
      # estimate Std. Error
# mean 0.1483035 0.05855348
# sd   2.9276738 0.04140354
SMTMP<-rtruncnorm(50000, min(sim_pars_2$SMTMP), max(sim_pars_2$SMTMP), mean = 0.1483035, sd = 2.9276738)


 fln_DEP_IMP
     # estimate Std. Error
# mean 3556.082   33.39346
# sd   1669.921   23.61274
DEP_IMP<-rtruncnorm(50000, min(sim_pars_2$DEP_IMP), max(sim_pars_2$DEP_IMP), mean = 3556.082, sd = 1669.921)


 fln_DDRAIN
     # estimate Std. Error
# mean 919.3622  11.782497
# sd   589.0822   8.330957
DDRAIN<-rtruncnorm(50000, min(sim_pars_2$DDRAIN), max(sim_pars_2$DDRAIN), mean = 919.3622, sd = 589.0822)


 fln_GDRAIN
     # estimate Std. Error
# mean 53.26543  0.5793765
# sd   28.96882  0.4096810
GDRAIN<-rtruncnorm(50000, min(sim_pars_2$GDRAIN), max(sim_pars_2$GDRAIN), mean = 53.26543, sd = 28.96882)


 fln_BACTKDQ
     # estimate Std. Error
# mean 284.7103   2.681929
# sd   134.0963   1.896409
BACTKDQ<-rtruncnorm(50000, min(sim_pars_2$BACTKDQ), max(sim_pars_2$BACTKDQ), mean = 284.7103, sd = 134.0963)


 fln_BACT_SWF
      # estimate  Std. Error
# mean 0.4383358 0.005708378
# sd   0.2854189 0.004036210
BACT_SWF<-rtruncnorm(50000, min(sim_pars_2$BACT_SWF), max(sim_pars_2$BACT_SWF), mean = 0.4383358, sd = 0.2854189)

 fln_THBACT
     # estimate Std. Error
# mean 3.234899 0.06148666
# sd   3.074333 0.04347761
THBACT<-rtruncnorm(50000, min(sim_pars_2$THBACT), max(sim_pars_2$THBACT), mean = 3.234899, sd = 3.074333)

 fln_WDPRCH
      # estimate  Std. Error
# mean 0.5312934 0.005626855
# sd   0.2813427 0.003978561
WDPRCH<-rtruncnorm(50000, min(sim_pars_2$WDPRCH), max(sim_pars_2$WDPRCH), mean = 0.5312934, sd = 0.2813427)



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


  save(bac_cal1, file='/work/OVERFLOW/RCR/calibration/MSU/bac_cal13.RData')













#7.	Keep simulating until we get 5k new simulations with an average_nse greater than the median_average_nse calculated for the previous set of 5k simulation results.
#8.	Repeat step 4 to calculate the updated unweighted kernel densities and fit to the normal distribution.
#9.	Now repeat step 5 for the new set of 5k simulations. Calculate the updated median_average_nse.
#10.	Repeat steps 6)-9) over and over until the median_average_nse fails to improve by X% versus the median_average_nse from the last generation. We have not explicitly defined what this percentage is just yet. It will probably be something like 1% or less.



