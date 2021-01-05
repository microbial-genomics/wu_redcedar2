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

# every generation will have n_to_keep accepted particles
# the median score of these n_to_keep particles will be used as the cutoff for the next generation
# the first generation will be the top n_to_keep of nsims_todo simulations
# start the loop here
for(iter in startgen:ngens){
  print(paste("*********** begin generation", iter, "**********************"))
  #
  #number to keep each generation
  if(iter==0){
    nsims_todo <- 10000
    n_to_keep <- 2000
    pars_tibble <- create_tibble_initial(nsims_todo)
    # create dataframe to persistently store stats
    generation_stats <- create_generation_stats(startgen, ngens, n_to_keep)
    #will be NA for generation 0 because not needed, this routine will clear out the vector though
    cutoff_mean_nse_score <- get_cutoff_mean_nse_score(iter, generation_stats)      
  }else{
    n_to_keep <- 2000 
    #load current generation stats
    load_generation_stats(iter, data_in_dir)
    #get nsims for this iteration
    nsims_todo <- generation_stats$nsims[iter+1]
    #load parameter_input_sims
    load_parameter_input_sims(iter, data_in_dir)
    #assign parameters to input tibble
    pars_tibble <- create_next_sim_tibble(nsims_todo, parameter_input_sims)
  }
  # run the swat simulations for this iteration
  bac_cal_output <- simulate_generation_next(iter, nsims_todo, swat_path, base_dir, pars_tibble)  
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
  # find the top 5k
  if(iter==0){
    # find the 75th percentile of nse_mean
    this_cutoff_mean_nse_score <- quantile(nse_mean, probs=0.75)
  }else{  
    #use the mean_nse score from the last generation to sort the keepers
    this_cutoff_mean_nse_score <- get_cutoff_mean_nse_score(iter, generation_stats)
  }
  #determine the first 5k to keep, combine keeper nses w parameters into a df
  all_keepers <- which(nse_mean > this_cutoff_mean_nse_score)
  n_all_keepers <- length(all_keepers)
  proportion_kept <- n_all_keepers/nsims_todo
  print(paste("we had", n_all_keepers, "of", nsims_todo, "simulations that had a better mean_nse score of", this_cutoff_mean_nse_score))
  valid_keepers <- head(all_keepers, n = n_to_keep)
  nses_w_parameters <- cbind(nse_bac[valid_keepers], nse_flux[valid_keepers], nse_q[valid_keepers], 
                             nse_mean[valid_keepers], sim_pars[valid_keepers,])
  ## find the updated unweighted kernel densities based on these new 5k simulations
  #kde_next_gen <- sim_pars_keepers %>% 
  #  gather(key = "par", value = "parameter_range")
  ## print the distributions
  #save_kde_pdf()  
  #save nses_parameters
  save_nses_parameters(iter, data_in_dir, nses_w_parameters)
  # save concentration time series output to an .RData file 
  # for later sensitivity analyses TODO
  # delete the local bac_cal file TODO
  # calculate mean_nse score
  next_cutoff_mean_nse_score <- median(nses_w_parameters$nse_mean)
  # save next mean_nse score for future use
  update_cutoff_mean_nse_score(iter, generation_stats, next_cutoff_mean_nse_score)
  # log results
  log_results(iter, this_cutoff_mean_nse_score, n_all_keepers, nsims_todo, nse_mean,
              nse_bac, nse_q, nse_flux, next_cutoff_mean_nse_score)
  # update and save parameter inputs
  fitted_parameter_list <- fit_normal_parameters(sim_pars)
  save_fitted_parameter_list(iter, data_in, fitted_parameter_list)
  # calculate nsims for next generation
  proportion_kept <- n_all_keepers/nsims_todo
  next_nsims <- calculate_next_nsims(n_to_keep, proportion_kept)
  # sample from input distributions for next generation
  parameter_input_sims <- sample_truncated_normals(iter, next_nsims, fitted_parameter_list)
  # save parameter inputs for next round of simulations
  save_parameter_input_sims(iter, data_in_dir, parameter_input_sims)
  # save stats for next generation
  generation_stats <- update_generation_stats(iter, generation_stats, next_nsims, max(nse_mean), 
                                              proportion_kept, next_cutoff_mean_nse_score)
  save_generation_stats(iter, data_in_dir, generation_stats)
  print(paste("*********** end generation", iter, "**********************"))
}

