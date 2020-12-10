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

#setup directory structure
#set paths for local machine or hpc
# we are dumping everything in root directory on hpc
huiyun <- FALSE #Huiyun set to true when you are running this code
print("load support functions")
if(huiyun){
  base_dir <- file.path("/work", "OVERFLOW", "RCR", "stp", "MSU")
  data_in_dir <- base_dir
  graphics_dir <- base_dir
  src_dir <- base_dir
  project_path <- base_dir
  swat_path <- base_dir
}else{
  if(Sys.info()[4]=="LZ2626UTPURUCKE"){
    base_dir <- file.path("c:", "git", "wu_redcedar2")
    data_in_dir <- file.path(base_dir, "data_in")
    graphics_dir <- file.path(base_dir, "graphics")
    src_dir <- file.path(base_dir, "src")
    project_path <- base_dir
    swat_path <- base_dir
  }else{
    base_dir <- file.path("/work", "OVERFLOW", "stp", "MSU")
    data_in_dir <- base_dir
    graphics_dir <- base_dir
    src_dir <- base_dir
    project_path <- base_dir
    swat_path <- base_dir
  }
}

# source support functions
source(file.path(src_dir, "abc_functions.R"))

#load outside data
load_observations()

# preset the generations to be simulated
# stargen = 0 means starting from scratch
startgen <- 0
ngens <- 40

# start the loop here
for(iter in startgen:ngens){
  #will be NA for generation 0 because not needed, this routine will clear out the vector though
  cutoff_median_score <- get_cutoff_median_score(iter)
  #  first (zeroeth) generation
  if(iter==0){
    # every generation will have 5000 accepted particles
    # the median score of these 5000 particles will be used as the cutoff for the next generation
    # the first generation will be 5000 without rejection
    nsims <- 5000
    pars_initial <- create_tibble_initial(nsims)
    #simulate_generation_zero(nsims, swat_path, base_dir, pars_initial)
    # run the initial set of swat simulations
    bac_cal_output <- simulate_generation_zero(nsims, swat_path, base_dir, pars_initial)
    # save output to disk
    save_bac_cal_output(iter, bac_cal_output)
    #get parameter names and values
    sim_pars <- bac_cal_output$parameter$values
    # calculate various nses
    nse_bac <- calculate_nse_bac(iter, bac_cal_output, bac_obs)
    nse_q <- calculate_nse_q(iter, bac_cal_output, q_obs)
    nse_flux <- calculate_nse_flux(iter, bac_cal_output, bac_obs, q_obs)
    # calculate nse means
    nse_mean <- calculate_nse_mean(iter, nse_bac, nse_q, nse_flux)
    #combine nses w parameters
    nses_w_parameters <- cbind(nse_bac, nse_flux, nse_q, nse_mean, sim_pars)
    # calculate median score
    next_median_score <- median(nses_w_parameters$nse_mean)
    # save next median score for future use
    save_cutoff_median_score(iter, data_in_dir, next_median_score)
    # log results
    log_results(iter, cutoff_median_score, nsims, nsims, nse_mean,
                nse_bac, nse_q, nse_flux, next_median_score)
    # update and save parameter inputs
    fitted_parameter_list <- fit_normal_parameters(sim_pars)
    save_fitted_parameter_list(iter, data_in, fitted_parameter_list)
    # calculate nsims for next generation
    #TODO
    # sample from input distributions for next generation
    parameter_input_sims <- sample_truncated_normals(iter, new_nsims, fitted_parameter_list)
    save_parameter_input_sims(iter, data_in, parameter_input_sims) #TODO
  }
}else{
  #
  ### subsequent runs
  #load in last set of simulations
  # list with parameter and simulation elements
  load_previous_swat_simulations(iter, data_in_dir)

  
  
  #determine number of simulations last time
  length(bac_cal_output$simulation$bac_out)
  #set variables
  previous_nsims <- ncol(bac_cal_output$simulation$bac_out) - 1
  print(paste("the last generation ", iter-1, " had ", previous_nsims, " sims"))
  n_to_keep <- 5000 #number to keep each generation
  
  #load the simulated concentrations, flows, inputs for last simulations
  sim_bac<- bac_cal1$simulation$bac_out
  sim_q <- bac_cal1$simulation$q_out
  sim_pars <- bac_cal1$parameter$values
  
  # merge simulated data with above observations, calculate nses
  nse_simulations_v_observations()
  nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
    dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
    map_dbl(., ~NSE(.x, bac_obs$bacteria))
  sort(nse_bac, decreasing = T) %>% enframe()
  
  nse_q <- right_join(sim_q,q_obs,by="date") %>%
    dplyr::select(-date) %>% dplyr::select(-discharge) %>%
    map_dbl(., ~NSE(.x, q_obs$discharge))
  sort(nse_q, decreasing = T) %>% enframe()
  
  flux_sim <- sim_bac[c(97:167),c(-1)]*
    sim_q[c(97:167), c(-1)]*10^4
  
  nse_flux <- flux_sim %>%
    map_dbl(., ~NSE(.x, flux_obs[,1]))
  sort(nse_flux, decreasing = T) %>% enframe()
    
  #calculate average nse of conc, flow and flux for all sims
  calculate_mean_nses()
  nse_mean_calc <- matrix(data=NA, nrow=previous_nsims, ncol=1)
  for(i in 1:previous_nsims){
    nse_mean_calc[i] <- mean(c(nse_bac[i], nse_q[i], nse_flux[i]))
  }
  sort(nse_mean_calc, decreasing = T) %>% enframe()
  
  run <- c(1:previous_nsims)
  nse_mean <-cbind(run, nse_mean_calc)
  colnames(nse_mean)<-c("run", "nse_mean")
  
  #should make this a 3D dataframe across the generations
  # but we will start with just saving the median scores
  
  #use the median score from the last generation to sort the keepers
  #kluge until we write a tracker csv
  all_keepers <- which(nse_mean[,2] > previous_median_score)
  n_all_keepers <- length(all_keepers)
  valid_keepers <- head(all_keepers, n = n_to_keep)
  sim_pars_keepers <- sim_pars[valid_keepers,] # this dimension is problematic for saving with everything else
  nse_mean_keepers <- nse_mean[valid_keepers,2]
  
  # find the updated unweighted kernel densities based on these new 10k simulations
  kde_next_gen <- sim_pars_keepers %>% 
    gather(key = "par", value = "parameter_range")
  
  # print the distributions
  save_kde_pdf()

  # find the median score that will be used for the next generation
  new_median_score <- median(nse_mean_keepers)
  proportion_kept <- n_all_keepers/previous_nsims

  # print results
  log_results()
  
  #fit to the normal distribution, truncate at the original range limits for each parameter.
  # fit the new distributions assuming normality from the keepers
  fit_normal_parameters()
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
  sample_truncated_normals()  

  # create tibble for next generation of monte carlo
  # pars <- create_tibble_subsequent()

  pars <- tibble(#hydrology parameters (11)
    "CN2.mgt|change = relchg"= runif(new_nsims, -0.25, 0.1),
    "GWQMN.gw|change = relchg" = runif(new_nsims, -0.5, 2),
    "ALPHA_BNK.rte|change = absval" =runif(new_nsims, 0, 1),
    "CH_K2.rte|change = absval" = runif(new_nsims, 0, 500),
    "CH_N2.rte|change = absval" = runif(new_nsims, 0, 0.3),
    "TRNSRCH.bsn|change = absval" = runif(new_nsims, 0, 1),
    "CH_N1.sub|change = absval" = runif(new_nsims, 0.01, 30),
    "CH_K1.sub|change = absval" = runif(new_nsims, 0, 300),
    "RCHRG_DP.gw|change = absval" = runif(new_nsims, 0, 1),
    "SFTMP.bsn|change = absval"= runif(new_nsims, -5, 5),
    "SMTMP.bsn|change = absval"= runif(new_nsims, -5, 5),
    #tile drainage and sediments (3)
    "DEP_IMP.hru|change = absval"= runif(new_nsims, 0, 6000),
    "DDRAIN.mgt|change = absval"= runif(new_nsims, 0, 2000),
    "GDRAIN.mgt|change = absval"= runif(new_nsims, 0, 100),
    #bacteria submodel (4)
    "BACTKDQ.bsn|change = absval" = runif(new_nsims, 0, 500),
    "BACT_SWF.bsn|change = absval" = runif(new_nsims, 0, 1),
    "THBACT.bsn|change = absval"= runif(new_nsims, 0, 10),
    "WDPRCH.bsn|change = absval"= runif(new_nsims, 0, 1)
  )
  # run swat
  bac_cal1 <- run_swat_red_cedar(base_dir, pars)

  # post-process swat with logging
  post_process_swat()
}

