options("width"=132)
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
library(xts)


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
  base_dir <- file.path("/work", "OVERFLOW", "RCR", "sim55")
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
    base_dir <- file.path("/work", "OVERFLOW", "stp")
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
#dim(bac_obs)
#View(bac_obs)
#dim(q_obs)

# create weekly average output for concentration observations
obs_data_xts <- as.xts(bac_obs$bacteria,order.by=as.Date(bac_obs$date))
bac_obs_weekly_temp <- as.data.frame(apply.weekly(obs_data_xts, mean))
bac_obs_weekly <- as.data.frame(cbind(rownames(bac_obs_weekly_temp), bac_obs_weekly_temp$V1))
colnames(bac_obs_weekly) <- c("date","bacteria")
save(bac_obs_weekly, file = file.path(base_dir, "bac_obs_w.RData"))
bac_sample_dates <- bac_obs$date

#create weekly average output for flow on days with observed concentrations
flow_obs <- right_join(q_obs, bac_obs, by="date") #needed to reduce the number fo flow observations
obs_flow_xts <- as.xts(flow_obs$discharge,order.by=as.Date(flow_obs$date))
flow_obs_weekly_temp <- as.data.frame(apply.weekly(obs_flow_xts, mean))
flow_obs_weekly <- as.data.frame(cbind(rownames(flow_obs_weekly_temp), flow_obs_weekly_temp$V1))
colnames(flow_obs_weekly) <- c("date","discharge")
save(flow_obs_weekly, file = file.path(base_dir, "q_obs_w.RData"))
q_sample_dates <- flow_obs$date

#create weekly flux output
obs_flux_xts <- as.xts(flux_obs$flux,order.by=as.Date(flux_obs$date))
flux_obs_weekly_temp <- as.data.frame(apply.weekly(obs_flux_xts, mean))
flux_obs_weekly <- as.data.frame(cbind(rownames(flux_obs_weekly_temp), flux_obs_weekly_temp$V1))
colnames(flux_obs_weekly) <- c("date","flux")
save(flux_obs_weekly, file = file.path(base_dir, "flux_obs_w.RData"))
flux_sample_dates <- flux_obs$date

# preset the generations to be simulated
# stargen = 0 means starting from scratch
startgen <- 0
ngens <- 40

# decide what time frame/interval to optimize on 
# opt_time_interval <- "daily"
opt_time_interval <- "weekly"

# should concentrations be logged (modified nash sutcliffe) 
# opt_conc_transform <- "none"
opt_conc_transform <- "logged"

# decide what metric to optimize on
opt_nse <- "mean"
#opt_nse <- "conc"
# opt_nse <- "flux"
# opt_nse <- "flow"

# every generation will have n_to_keep accepted particles
# the median score of these n_to_keep particles will be used as the cutoff for the next generation
# the first generation will be the top n_to_keep of nsims_todo simulations
# start the loop here
for(iter in startgen:ngens){
  print(paste("*********** begin generation", iter, "**********************"))
  print(paste("optimizing based on", opt_nse))
  #
  #number to keep each generation
  nsims_todo <- 1000
  n_to_keep <- 200
  if(iter==0){
    pars_tibble <- create_tibble_initial(nsims_todo)
    # create dataframe to persistently store stats
    generation_stats <- create_generation_stats(startgen, ngens, opt_nse, n_to_keep)
    #will be NA for generation 0 because not needed, this routine will clear out the vector though
    cutoff_score <- get_cutoff_score(iter, generation_stats)      
  }else{
    #load current generation stats
    load_generation_stats(iter, data_in_dir)
    #get nsims for this iteration
    nsims_todo <- generation_stats$nsims[iter+1]
    #load parameter_input_sims
    load_parameter_input_sims(iter, data_in_dir)
    #assign parameters to input tibble
    pars_tibble <- create_next_sim_tibble(nsims_todo, parameter_input_sims)
  }
  # run the swat simulations for this iteration and save daily output
  bac_cal_output <- simulate_generation_next(iter, nsims_todo, swat_path, base_dir, pars_tibble)  
  # extract only the observed days from the simulated daily output
  #which(bac_cal_output$simulation$bac_out$date==sample_dates)
  #dim(bac_cal_output$simulation$bac_out)
  #View(bac_cal_output$simulation$bac_out)
  #bac_obs$date
  bac_sims_sample_dates <- right_join(bac_cal_output$simulation$bac_out,bac_obs, by="date")
  #dim(bac_sims_sample_dates)
  #View(bac_obs)
  #dim(bac_obs)
  # then create weekly average output for sims
  bac_sims_data <- as.xts(bac_sims_sample_dates$bacteria,order.by=as.Date(bac_sims_sample_dates$date))
  #dim(bac_sims_data)
  bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_data,mean))
  colnames(bac_sims_weekly) <- c("date", "bacteria")
  #class(bac_sims_weekly)

  # save output to disk
  save_bac_cal_output(iter, bac_cal_output)
  #get parameter names and values
  sim_pars <- bac_cal_output$parameter$values
  # calculate various nses for daily data
  nse_bac_daily <- calculate_nse_bac(iter, bac_cal_output, bac_obs)
  nse_q_daily <- calculate_nse_q(iter, bac_cal_output, q_obs)
  nse_flux_daily <- calculate_nse_flux(iter, bac_cal_output, flux_obs)
  # calculate various nses for daily data with logged concentrations
  mnse_bac_daily <- calculate_modified_nse_bac(iter, bac_cal_output, bac_obs)
  mnse_q_daily <- calculate_modified_nse_q(iter, bac_cal_output, q_obs)
  mnse_flux_daily <- calculate_modified_nse_flux(iter, bac_cal_output, flux_obs)
  # calculate various nses for weekly data
  nse_bac_weekly <- calculate_nse_bac(iter, bac_sims_weekly, bac_obs_weekly)
  nse_q_weekly <- calculate_nse_q(iter, bac_sim_weekly, flow_obs_weekly)
  nse_flux_weekly <- calculate_nse_flux(iter, bac_sim_weekly, flux_obs_weekly)
  # calculate various nses for weekly data with logged concentrations
  mnse_bac_weekly <- calculate_modified_nse_bac(iter, bac_sims_weekly, bac_obs_weekly)
  mnse_q_weekly <- calculate_modified_nse_q(iter, bac_sim_weekly, flow_obs_weekly)
  mnse_flux_weekly <- calculate_modified_nse_flux(iter, bac_sim_weekly, flux_obs_weekly)
  # calculate nse means
  nse_mean_daily <- calculate_nse_mean(iter, nse_bac, nse_q, nse_flux)
  nse_mean_weekly <- calculate_nse_mean(iter, nse_bac_weekly, nse_q_weekly, nse_flux_weekly)
  # calculate modified nse means
  mnse_mean_daily <- calculate_modified_nse_mean(iter, mnse_bac, mnse_q, mnse_flux) 
  mnse_mean_weekly <- calculate_modified_nse_mean(iter, mnse_bac_weekly, mnse_q_weekly, mnse_flux_weekly) 
  # get cutoff score
  if(iter==0){
    # find the 80th percentile of target nse, top (2000 of 10000)
    # daily, none
    if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="none"){
      this_cutoff_score <- quantile(nse_mean, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_bac, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_q, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flux, probs=0.8)
    # daily, logged
    } else if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="logged"){
      this_cutoff_score <- quantile(nse_mean, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_bac, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_q, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_flux, probs=0.8)  
    # weekly, none
    } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="none"){
      this_cutoff_score <- quantile(nse_mean, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_bac, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_q, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flux, probs=0.8)  
    # weekly, logged
    } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="logged"){
      this_cutoff_score <- quantile(nse_mean, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_bac, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_q, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
      this_cutoff_score <- quantile(nse_flux, probs=0.8)  
    }
  }else{  
    #use the cutoff score from the last generation to sort the keepers
    this_cutoff_score <- get_cutoff_score(iter, generation_stats)
  }
  #determine the first Xk to keep, combine keeper nses w parameters into a df    
  # daily, none
  if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="none"){
    all_keepers <- which(nse_mean_daily > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="none") {
    all_keepers <- which(nse_bac_daily > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="none") {
    all_keepers <- which(nse_q_daily > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flux_daily > this_cutoff_score)
  # daily, logged
  } else if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="logged"){
    all_keepers <- which(mnse_bac_daily > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_bac_daily > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_q_daily > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_flux_daily > this_cutoff_score) 
  # weekly, none
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="none"){
    all_keepers <- which(nse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_q_weekly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flux_weekly > this_cutoff_score)  
  # weekly, logged
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="logged"){
    all_keepers <- which(mnse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_q_weekly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="logged") {
    all_keepers <- which(mnse_flux_weekly > this_cutoff_score)  
  }
  n_all_keepers <- length(all_keepers)
  proportion_kept <- n_all_keepers/nsims_todo
  print(paste("we had", n_all_keepers, "of", nsims_todo, "simulations that had a better cutoff score of", this_cutoff_score))
  valid_keepers <- head(all_keepers, n = n_to_keep)
  keeper <- array(data="reject", dim=nsims_todo)
  keeper[all_keepers] <- "not_kept"
  keeper[valid_keepers] <- "kept"
  keeper <- as.factor(keeper)
  #create data.frame for the keepers with parameters
  # daily, none
  if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="none"){
    # concatenate output into one big dataframe
    nses_w_parameters_all <- cbind(keeper, nse_bac_daily, nse_flux_daily, nse_q_daily, 
                               nse_mean_daily, sim_pars)
    nse_conc_keepers <- nse_bac[valid_keepers]
    nse_flow_keepers <- nse_q[valid_keepers]
    nse_flux_keepers <- nse_flux[valid_keepers]
    nse_mean_keepers <- nse_mean[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
    # plot nses versus each other
    plot_bac_v_flow_pdf(iter, nses_w_parameters_all)
    #save nses_parameters
    save_nses_parameters(iter, data_in_dir, nses_w_parameters)
  # daily, logged
  } else if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="logged"){
    
  # weekly, none
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="none"){
    
  # weekly, logged
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="logged"){
    

  # save concentration time series output to an .RData file 
  # for later sensitivity analyses TODO
  # delete the local bac_cal file TODO
  # calculate cutoff score
  if(opt_nse=="mean"){
    next_cutoff_score <- median(nse_mean_keepers)
  } else if(opt_nse=="conc") {
    next_cutoff_score <- median(nse_conc_keepers)
  } else if(opt_nse=="flow") {
    next_cutoff_score <- median(nse_flow_keepers)
  } else if(opt_nse=="flux") {
    next_cutoff_score <- median(nse_flux_keepers)
  }
  # save next mean_nse score for future use
  update_cutoff_score(iter, generation_stats, next_cutoff_score)
  # log results
  log_results(iter, this_cutoff_score, n_all_keepers, nsims_todo, nse_mean_keepers,
              nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, next_cutoff_score)
  # update and save parameter inputs
  fitted_parameter_list <- fit_normal_parameters(sim_pars[valid_keepers,])
  save_fitted_parameter_list(iter, data_in, fitted_parameter_list)
  ## find the updated unweighted kernel densities based on these new 5k simulations
  kde_next_gen <- sim_pars[valid_keepers,] %>% 
    gather(key = "par", value = "parameter_range")
  ## print the distributions
  plot_kde_pdf(iter, kde_next_gen)  
  # calculate nsims for next generation
  proportion_kept <- n_all_keepers/nsims_todo
  next_nsims <- calculate_next_nsims(n_to_keep, proportion_kept)
  # sample from input distributions for next generation
  parameter_input_sims <- sample_truncated_normals(iter, next_nsims, fitted_parameter_list)
  # save parameter inputs for next round of simulations
  save_parameter_input_sims(iter, data_in_dir, parameter_input_sims)
  # save stats for next generation
  generation_stats <- update_generation_stats(iter, generation_stats, next_nsims, 
                                              max(nse_mean_keepers), max(nse_conc_keepers), max(nse_flow_keepers), max(nse_flux_keepers), 
                                              proportion_kept, next_cutoff_score)
  save_generation_stats(iter, data_in_dir, generation_stats)
  print(paste("*********** end generation", iter, "**********************"))
}

