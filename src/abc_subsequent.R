#check to make sure required packages are installed
list.of.packages <- c("devtools","dplyr","dygraphs","fast","forcats","ggplot2",
                      "hydroGOF","lhs","lubridate","mapview","plotly",
                      "purrr","sensitivity","sf","SWATplusR","tibble",
                      "tidyr","fitdistrplus","truncnorm")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)>0) {install.packages(new.packages)}
 
#installed.packages()
#load relevant packages
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

#SWATPlusR needs to be installed via devtools
# can require manual installation of tidy, etc packages with reboots
# devtools::install_github("chrisschuerz/SWATplusR")
library(SWATplusR)

# the swat executable swat2012_rev670 needs to be copied to the run directory
# and you must have exe privileges on it

if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  base_dir <- file.path("c:", "git", "wu_redcedar2")
  data_in_dir <- file.path(base_dir, "data_in")
  graphics_dir <- file.path(base_dir, "graphics")
}else{
  base_dir <- file.path("/work", "OVERFLOW", "stp", "MSU")
  data_in_dir <- base_dir
  graphics_dir <- base_dir
}

#load outside data
load(file= file.path(data_in_dir,'bac_obs.RData'))
load(file = file.path(data_in_dir,'flux_obs.RData'))
load(file = file.path(data_in_dir,'pcp_obs.RData'))
load(file = file.path(data_in_dir,'pcp_obs2.RData'))
load(file= file.path(data_in_dir, 'q_obs.RData'))
load(file= file.path(data_in_dir, 'q_obs2.RData'))


# ### initial run
# pars <- tibble("CN2.mgt|change = relchg"= runif(nsims,-0.25,0.1),
#                "SOL_K(1).sol|change = relchg" = runif(nsims,-0.8,0.8),
#                "SOL_AWC(1).sol|change = relchg" = runif(nsims,-0.8,2),
#                "OV_N.hru|change = relchg" = runif(nsims,-0.8,2),
#                "ALPHA_BF.gw|change = relchg" = runif(nsims,-0.3,0.3),
#                "GW_DELAY.gw|change = relchg" = runif(nsims,-15,10),
#                "GWQMN.gw|change = relchg" = runif(nsims,-0.5, 2),
#                "HRU_SLP.hru|change = relchg" = runif(nsims,-15,10),
#                "SLSUBBSN.hru|change = relchg" = runif(nsims,-0.5, 3),	
#                "ALPHA_BNK.rte|change = absval" = runif(nsims, 0, 1),
#                "CH_K2.rte|change = absval" = runif(nsims, 0, 500),
#                "CH_N2.rte|change = absval" = runif(nsims, 0, 0.3),
#                "ESCO.bsn |change = absval" = runif(nsims, 0, 1),
#                "EPCO.bsn|change = absval" = runif(nsims, 0, 1),
#                "TRNSRCH.bsn|change = absval" = runif(nsims, 0, 1),
#                "SURLAG.bsn|change = absval" = runif(nsims, 1, 24),
#                "CH_N1.sub|change = absval" = runif(nsims, 0.01, 30),
#                "CH_K1.sub|change = absval" = runif(nsims, 0, 300),
#                "REVAPMN.gw |change = absval" = runif(nsims, 0, 1000),
#                "GW_REVAP.gw|change = absval" = runif(nsims, 0.02, 0.2),
#                "RCHRG_DP.gw|change = absval" = runif(nsims, 0, 1),
#                "GW_SPYLD.gw|change = absval" = runif(nsims, 0, 0.4))
# 
# path <- "/work/OVERFLOW/stp/MSU"
# 
# swat_output0 <- run_swat2012(project_path = path,
#                              output = list(q_out = define_output(file = "rch",
#                                                                  variable = "FLOW_OUT",
#                                                                  unit = 4),
#                                            bac_out = define_output(file = "rch",
#                                                                    variable = "BACTP_OUT",
#                                                                    unit = 4)),
#                              parameter = pars,
#                              start_date = "2011-01-01",
#                              end_date = "2013-12-31",
#                              years_skip = 2,
#                              n_thread = 32)
# 
# save_file <- paste(path,"/swat_output0.RData",sep="")
# save(swat_output0, file=save_file)


## start the loop here
for(iter in 23:26){

  ### subsequent runs
  #load in last set of simulations
  # list with parameter and simulation elements
  bac_cal_filename <- paste('bac_cal', iter, '.RData', sep="")
  rdata_file_in <- file.path(data_in_dir, bac_cal_filename)
  print(paste("loading data file: ", rdata_file_in))
  load(file = rdata_file_in)
  
  #determine number of simulations last time
  dim(bac_cal1$simulation$bac_out)
  
  #set variables
  previous_nsims <- ncol(bac_cal1$simulation$bac_out) - 1
  print(paste("the last generation ", iter, " had ", previous_nsims, " sims"))
  n_to_keep <- 5000 #number to keep each generation
  
  #load the simulated concentrations, flows, inputs for last simulations
  sim_bac<- bac_cal1$simulation$bac_out
  sim_q <- bac_cal1$simulation$q_out
  sim_pars <- bac_cal1$parameter$values
  
  # merge simulated and observed bacteria concentrations, calculate nses for all sims
  nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
    dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
    map_dbl(., ~NSE(.x, bac_obs$bacteria))
  sort(nse_bac, decreasing = T) %>% enframe()
  
  # merge simulated and observed flows, calculate nses for all sims
  nse_q <- right_join(sim_q,q_obs,by="date") %>%
    dplyr::select(-date) %>% dplyr::select(-discharge) %>%
    map_dbl(., ~NSE(.x, q_obs$discharge))
  sort(nse_q, decreasing = T) %>% enframe()
  
  # calculate the simulated fluxes from concs and flows for all sims
  flux_sim <- sim_bac[c(97:167),c(-1)]*
    sim_q[c(97:167), c(-1)]*10^4
  
  #merge simulated and observed fluxes, calculate nses for all sims
  nse_flux <- flux_sim %>%
    map_dbl(., ~NSE(.x, flux_obs[,1]))
  sort(nse_flux, decreasing = T) %>% enframe()
  
  #calculate average nse of conc, flow and flux for all sims
  nse_mean_calc <- matrix(data=NA, nrow=previous_nsims, ncol=1)
  for(i in 1:previous_nsims){
    #print(i)
    nse_mean_calc[i] <- mean(c(nse_bac[i], nse_q[i], nse_flux[i]))
  }
  sort(nse_mean_calc, decreasing = T) %>% enframe()
  
  #create new df with average of the 3 nses for all sims
  run <- c(1:previous_nsims)
  nse_mean <-cbind(run, nse_mean_calc)
  colnames(nse_mean)<-c("run", "nse_mean")
  
  #should make this a 3D dataframe across the generations

  
  #use the median score from the last generation to sort the keepers
  #kluge until we write a tracker csv
  if(iter==15){previous_median_score = -1.679666} #from previous generation -- 14}
  if(iter==16){previous_median_score = -0.84941146160688} #from previous generation -- 15}
  if(iter==23){previous_median_score = -0.0161426045727871}
  all_keepers <- which(nse_mean[,2] > previous_median_score)
  n_all_keepers <- length(all_keepers)
  valid_keepers <- head(all_keepers, n = n_to_keep)
  sim_pars_keepers <- sim_pars[valid_keepers,] # this dimension is problematic for saving with everything else
  nse_mean_keepers <- nse_mean[valid_keepers,2]
  
  #8.	Calculate the updated unweighted kernel densities based on these new 10k simulations 
  #and fit to the normal distribution, truncate at the original range limits for each parameter.
  kde_next_gen <- sim_pars_keepers %>% 
    gather(key = "par", value = "parameter_range")
  
  # print the distributions
  ggplot(data = kde_next_gen) +
    geom_density(aes(x = parameter_range)) +
    facet_wrap(.~par, nrow=5, scales = "free") +
    theme_bw()
  density_plot_filename <- paste("kde_mcabc_gen", iter, ".pdf", sep="")
  ggsave(file.path(graphics_dir, density_plot_filename))
  
  # ###
  #9.	Now use these new 10k simulations to calculate the updated first_quartile_average_nse, 
  # this will be the average of the 2500 and 2501st highest average_nse.
  new_median_score <- median(nse_mean_keepers)
  proportion_kept <- n_all_keepers/previous_nsims
  #[1] -0.9092294
  
  # print results
  print(paste("Generation ", iter))
  print(paste("median score for the last generation was:", previous_median_score))
  print(paste("generation x:",n_all_keepers, "of", format(previous_nsims,scientific=F), " simulations kept; proportion kept =", round(proportion_kept,4)))
  print(paste("best kept mean nse for this generation is:", max(round(nse_mean_keepers,4))))
  print(paste("best bacteria nse for this generation is:", max(round(nse_bac,4))))
  print(paste("best flow nse for this generation is:", max(round(nse_q,4))))
  print(paste("best flux nse for this generation is:", max(round(nse_flux,4))))
  print(paste("median mean nse score for this generation is:", round(new_median_score,4)))
  

  # fit the new distributions assuming normality from the keepers
  fitted_CN2 <- fitdist(sim_pars_keepers$CN2, "norm")
  fitted_GWQMN <- fitdist(sim_pars_keepers$GWQMN, "norm")
  fitted_ALPHA_BNK <- fitdist(sim_pars_keepers$ALPHA_BNK, "norm")
  fitted_CH_K2 <- fitdist(sim_pars_keepers$CH_K2, "norm")
  fitted_CH_N2 <- fitdist(sim_pars_keepers$CH_N2, "norm")
  fitted_TRNSRCH <- fitdist(sim_pars_keepers$TRNSRCH, "norm")
  fitted_CH_N1 <- fitdist(sim_pars_keepers$CH_N1, "norm")
  fitted_CH_K1 <- fitdist(sim_pars_keepers$CH_K1, "norm")
  fitted_RCHRG_DP <- fitdist(sim_pars_keepers$RCHRG_DP, "norm")
  fitted_SFTMP <- fitdist(sim_pars_keepers$SFTMP, "norm")
  fitted_SMTMP <- fitdist(sim_pars_keepers$SMTMP, "norm")
  fitted_DEP_IMP <- fitdist(sim_pars_keepers$DEP_IMP, "norm")
  fitted_DDRAIN <- fitdist(sim_pars_keepers$DDRAIN, "norm")
  fitted_GDRAIN <- fitdist(sim_pars_keepers$GDRAIN, "norm")
  fitted_BACTKDQ <- fitdist(sim_pars_keepers$BACTKDQ, "norm")
  fitted_BACT_SWF<- fitdist(sim_pars_keepers$BACT_SWF, "norm")
  fitted_THBACT <- fitdist(sim_pars_keepers$THBACT, "norm")
  fitted_WDPRCH <- fitdist(sim_pars_keepers$WDPRCH, "norm")
  
  #reset nsims based on acceptance frequency of last generation
  new_nsims <- max(10000, round(n_to_keep/proportion_kept)*2)
  print(paste("next round we will do", new_nsims, "simulations"))
  
  #hard coding truncated parameters based on values from Sensitivity.R
  CN2_mean <- fitted_CN2$estimate[1]
  CN2_sd <- fitted_CN2$estimate[2]
  print(paste("generation", iter, "CN2", round(CN2_mean,3), round(CN2_sd,3)))
  
  CN2 <- rtruncnorm(new_nsims, -0.25, 0.1, mean = CN2_mean, sd =  CN2_sd)
  
  # fitted_GWQMN
  # #     estimate  Std. Error
  # mean 0.4499712 0.007905126
  # sd   0.5589768 0.005589688
  GWQMN_mean <- fitted_GWQMN$estimate[1]
  GWQMN_sd <- fitted_GWQMN$estimate[2]
  print(paste("generation", iter, "GWQMN", round(GWQMN_mean,3), round(GWQMN_sd,3)))
  GWQMN <- rtruncnorm(new_nsims, -0.5, 2, mean =GWQMN_mean, sd = GWQMN_sd)

  # fitted_ALPHA_BNK
  # Baseflow alpha factor for bank storage
  # estimate  Std. Error
  # mean 0.4911617 0.003369632
  # sd   0.2382690 0.002382501
  ALPHA_BNK_mean <- fitted_ALPHA_BNK$estimate[1]
  ALPHA_BNK_sd <- fitted_ALPHA_BNK$estimate[2]
  print(paste("generation", iter, "ALPHA_BNK", round(ALPHA_BNK_mean,3), round(ALPHA_BNK_sd,3)))
  ALPHA_BNK <- rtruncnorm(new_nsims, 0, 1, mean = ALPHA_BNK_mean, sd = ALPHA_BNK_sd)

  # fitted_CH_K2
  # Effective hydraulic conductivity in main channel alluvium
  # estimate Std. Error
  # mean 27.90991  0.1635768
  # sd   11.56662  0.1156662
  CH_K2_mean <- fitted_CH_K2$estimate[1]
  CH_K2_sd <- fitted_CH_K2$estimate[2]
  print(paste("generation", iter, "CH_K2", round(CH_K2_mean,3), round(CH_K2_sd,3)))
  CH_K2 <- rtruncnorm(new_nsims, 0, 500, mean = CH_K2_mean, sd = CH_K2_sd)
  
  # fitted_CH_N2
  # estimate   Std. Error
  # mean 0.10415042 0.0003312391
  # sd   0.02342214 0.0002323030
  CH_N2_mean <- fitted_CH_N2$estimate[1]
  CH_N2_sd <- fitted_CH_N2$estimate[2]
  print(paste("generation", iter, "CH_N2", round(CH_N2_mean,3), round(CH_N2_sd,3)))
  CH_N2 <- rtruncnorm(new_nsims, 0, 0.3, mean = CH_N2_mean, sd = CH_N2_sd)
  
  # fitted_TRNSRCH
  # estimate   Std. Error
  # mean 0.17602745 0.0009794462
  # sd   0.06925731 0.0006919234
  TRNSRCH_mean <- fitted_TRNSRCH$estimate[1]
  TRNSRCH_sd <- fitted_TRNSRCH$estimate[2]
  print(paste("generation", iter, "TRNSRCH", round(TRNSRCH_mean,3), round(TRNSRCH_sd,3)))
  TRNSRCH <- rtruncnorm(new_nsims, 0, 1, mean = TRNSRCH_mean, sd = TRNSRCH_sd)

  # fitted_CH_N1 
  # estimate   Std. Error
  # mean 0.09991472 0.0003314811
  # sd   0.02343925 0.0002324755
  CH_N1_mean <- fitted_CH_N1$estimate[1]
  CH_N1_sd <- fitted_CH_N1$estimate[2]
  print(paste("generation", iter, "CH_N1", round(CH_N1_mean,3), round(CH_N1_sd,3)))
  CH_N1 <- rtruncnorm(new_nsims, 0.01, 30, mean = CH_N1_mean, sd = CH_N1_sd)
  
  # fitted_CH_K1
  # estimate Std. Error
  # mean 158.93894  0.9599238
  # sd    67.87685  0.6787685
  CH_K1_mean <- fitted_CH_K1$estimate[1]
  CH_K1_sd <- fitted_CH_K1$estimate[2]
  print(paste("generation", iter, "CH_K1", round(CH_K1_mean,3), round(CH_K1_sd,3)))
  CH_K1 <- rtruncnorm(new_nsims, 0, 300, mean = CH_K1_mean, sd = CH_K1_sd)
  
  # fitted_RCHRG_DP
  # estimate  Std. Error
  # mean 0.4808134 0.003439768
  # sd   0.2432284 0.002432099
  RCHRG_DP_mean <- fitted_RCHRG_DP$estimate[1]
  RCHRG_DP_sd <- fitted_RCHRG_DP$estimate[2]
  print(paste("generation", iter, "RCHRG_DP", round(RCHRG_DP_mean,3), round(RCHRG_DP_sd,3)))
  RCHRG_DP <- rtruncnorm(new_nsims, 0, 1, mean = RCHRG_DP_mean, sd = RCHRG_DP_sd)
  
  # fitted_SFTMP
  # estimate Std. Error
  # mean 0.1230878 0.03301559
  # sd   2.3345547 0.02334553
  SFTMP_mean <- fitted_SFTMP$estimate[1]
  SFTMP_sd <- fitted_SFTMP$estimate[2]
  print(paste("generation", iter, "SFTMP", round(SFTMP_mean,3), round(SFTMP_sd,3)))
  SFTMP <- rtruncnorm(new_nsims, -5, 5, mean = SFTMP_mean, sd = SFTMP_sd)
  
  # fitted_SMTMP
  # estimate Std. Error
  # mean 0.1653056 0.03339191
  # sd   2.3611645 0.02361163
  SMTMP_mean <- fitted_SMTMP$estimate[1]
  SMTMP_sd <- fitted_SMTMP$estimate[2]
  print(paste("generation", iter, "SMTMP", round(SMTMP_mean,3), round(SMTMP_sd,3)))
  SMTMP <- rtruncnorm(new_nsims, -5, 5, mean = SMTMP_mean, sd = SMTMP_sd)
  
  # fitted_DEP_IMP
  # estimate Std. Error
  # mean 3667.422   18.83298
  # sd   1331.808   13.31908
  DEP_IMP_mean <- fitted_DEP_IMP$estimate[1]
  DEP_IMP_sd <- fitted_DEP_IMP$estimate[2]
  print(paste("generation", iter, "DEP_IMP", round(DEP_IMP_mean,3), round(DEP_IMP_sd,3)))
  DEP_IMP <- rtruncnorm(new_nsims, 0, 6000, mean = DEP_IMP_mean, sd = DEP_IMP_sd)
  
  # fitted_DDRAIN
  # estimate Std. Error
  # mean 902.7242   6.710620
  # sd   474.4889   4.744931
  DDRAIN_mean <- fitted_DDRAIN$estimate[1]
  DDRAIN_sd <- fitted_DDRAIN$estimate[2]
  print(paste("generation", iter, "DDRAIN", round(DDRAIN_mean,3), round(DDRAIN_sd,3)))
  DDRAIN <- rtruncnorm(new_nsims, 0, 2000, mean = DDRAIN_mean, sd = DDRAIN_sd)
  
  # fitted_GDRAIN
  # estimate Std. Error
  # mean 53.36850  0.3311553
  # sd   23.41622  0.2341622
  GDRAIN_mean <- fitted_GDRAIN$estimate[1]
  GDRAIN_sd <- fitted_GDRAIN$estimate[2]
  print(paste("generation", iter, "GDRAIN", round(GDRAIN_mean,3), round(GDRAIN_sd,3)))
  GDRAIN <- rtruncnorm(new_nsims, 0, 100, mean = GDRAIN_mean, sd = GDRAIN_sd)
  
  # fitted_BACTKDQ
  # estimate Std. Error
  # mean 284.6590   1.550853
  # sd   109.6619   1.096619
  BACTKDQ_mean <- fitted_BACTKDQ$estimate[1]
  BACTKDQ_sd <- fitted_BACTKDQ$estimate[2]
  print(paste("generation", iter, "BACTKDQ", round(BACTKDQ_mean,3), round(BACTKDQ_sd,3)))
  BACTKDQ <- rtruncnorm(new_nsims, 0, 500, mean = BACTKDQ_mean, sd = BACTKDQ_sd)
  
  # fitted_BACT_SWF
  # estimate  Std. Error
  # mean 0.4287223 0.003284534
  # sd   0.2322516 0.002322323
  BACT_SWF_mean <- fitted_BACT_SWF$estimate[1]
  BACT_SWF_sd <- fitted_BACT_SWF$estimate[2]
  print(paste("generation", iter, "BACT_SWF", round(BACT_SWF_mean,3), round(BACT_SWF_sd,3)))
  BACT_SWF <- rtruncnorm(new_nsims, 0, 1, mean = BACT_SWF_mean, sd = BACT_SWF_sd)
  
  # fitted_THBACT
  # estimate Std. Error
  # mean 1.384695 0.02367388
  # sd   1.673996 0.01673994
  THBACT_mean <- fitted_THBACT$estimate[1]
  THBACT_sd <- fitted_THBACT$estimate[2]
  print(paste("generation", iter, "THBACT", round(THBACT_mean,3), round(THBACT_sd,3)))
  THBACT <- rtruncnorm(new_nsims, 0, 10, mean = THBACT_mean, sd = THBACT_sd)
  
  # fitted_WDPRCH
  # estimate  Std. Error
  # mean 0.5643253 0.003093901
  # sd   0.2187718 0.002187513
  WDPRCH_mean <- fitted_WDPRCH$estimate[1]
  WDPRCH_sd <- fitted_WDPRCH$estimate[2]
  print(paste("generation", iter, "WDPRCH", round(WDPRCH_mean,3), round(WDPRCH_sd,3)))
  WDPRCH <- rtruncnorm(new_nsims, 0, 1, mean = WDPRCH_mean, sd = WDPRCH_sd)
  
  
  

  # create tibble for next generation of monte carlo
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
  
  # run swat
  bac_cal1 <- run_swat2012(project_path = base_dir,
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
  
  print(paste("swat runs finished for generation ", iter+1))
  rdata_out_filename <- paste('bac_cal', iter+1, '.RData', sep = "")
  rdata_file_out <- file.path(data_in_dir, rdata_out_filename)
  save(bac_cal1, file=rdata_file_out)
  print(paste("radata file for generation ", iter+1, " saved to ", rdata_file_out))
  previous_median_score = new_median_score
  print(paste("the median score cutoff used to create proposals for generation ", iter+1, " was ", previous_median_score))
}
