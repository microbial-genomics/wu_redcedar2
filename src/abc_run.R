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
# source support functions
print("load support functions")
src_dir <- "/work/OVERFLOW/RCR/stp/MSU"
source(file.path(src_dir, "abc_functions.R"))

#set paths for local machine or hpc
# we are dumping everything in root directory on hpc
set_working_paths()

#load outside data
load_observations()

# preset the generations to be simulated
# stargen = 0 means starting from scratch
startgen <- 0
ngens <- 22

# if we are starting at zero we will wipe out the previous_median_score_vector
median_filename <- file.path(data_in_dir, "previous_median_scores.csv")
if(startgen==0){
  previous_median_score <- rep(NA, ngens+1)
  write.csv(previous_median_score, file = median_filename)
} else {
  # if starting at gen greater than zero we will read it from the existing file
  previous_median_score <- read.csv(median_filename)
}

swat_path <- "/work/OVERFLOW/RCR/stp/MSU"

## start the loop here
for(iter in startgen:ngens){

  # first (zeroeth) generation  
  if(iter==0){
    # every generation will have 5000 accepted particles
    # the median score of these 5000 particles will be used as the cutoff for the next generation
    nsims <- 5000
    pars_initial <- create_tibble_initial(nsims)
    #simulate_generation_zero(nsims, swat_path, base_dir, pars_initial)
    # run the initial set of swat simulations
    print(paste("About to run generation 0 with", nsims, "simulations"))
    swat_output0 <- run_swat_red_cedar(swat_path, pars_initial)
    
    #save the simulations
    save_file <- file.path(data_in_dir, "rcr_swat_output0.RData")
    save(swat_output0, file = save_file)
    previous_median_score[iter] <- -1000000
  }
  
  ### subsequent runs
  #load in last set of simulations
  # list with parameter and simulation elements
  load_previous_swat_simulations()
  
  #set variables
  previous_nsims <- ncol(bac_cal1$simulation$bac_out) - 1
  print(paste("the last generation ", iter, " had ", previous_nsims, " sims"))
  n_to_keep <- 5000 #number to keep each generation
  
  #load the simulated concentrations, flows, inputs for last simulations
  sim_bac<- bac_cal1$simulation$bac_out
  sim_q <- bac_cal1$simulation$q_out
  sim_pars <- bac_cal1$parameter$values
  
  # merge simulated data with above observations, calculate nses
  nse_simulations_v_observations()
    
  #calculate average nse of conc, flow and flux for all sims
  calculate_mean_nses()
  
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
  
  #reset nsims based on acceptance frequency of last generation
  new_nsims <- max(10000, round(n_to_keep/proportion_kept)*2)
  print(paste("next round we will do", new_nsims, "simulations"))
  
  #hard coding truncated parameters based on values from Sensitivity.R
  sample_truncated_normals()  

  # create tibble for next generation of monte carlo
  pars <- create_tibble_subsequent()
  
  # run swat
  bac_cal1 <- run_swat_red_cedar(base_dir, pars)

  # post-process swat with logging
  post_process_swat()
}
