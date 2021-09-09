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


i#SWATPlusR needs to be installed via devtools
# can require manual installation of tidy, etc packages with reboots
# devtools::install_github("chrisschuerz/SWATplusR")
library(SWATplusR)

# the swat executable swat2012_rev670 needs to be copied to the run directory
# and you must have exe privileges on it

set.seed(42)

#setup directory structure
#set paths for local machine or hpc
# we are dumping everything in root directory on hpc
huiyun <- TRUE #Huiyun set to true when you are running this code
print("load support functions")
if(huiyun){
  base_dir <- file.path("/work", "OVERFLOW", "RCR", "sim56")
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
  }###Huiyun's local machine
  else{
    base_dir <- file.path("/PRIV", "SSO", "wu_redcedar2")
    data_in_dir <- base_dir
    graphics_dir <- base_dir
    src_dir <- file.path(base_dir, "src")
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

#remove highest point on 2006/9/18 for bac_obs and flux_obs 
#bac_obs<-bac_obs[-85,]
#flux_obs <-flux_obs[-85,]
# create weekly average output for concentration observations
bac_obs_daily <- bac_obs$bacteria #335#removed the highest number on 2006/9/18
obs_data_xts <- as.xts(bac_obs$bacteria,order.by=as.Date(bac_obs$date))
bac_obs_weekly <- as.data.frame(apply.weekly(obs_data_xts, mean)) #204
# create monthly average output for concentration observations
bac_obs_monthly <- as.data.frame(apply.monthly(obs_data_xts, mean))#51


#create weekly average output for flow on days with observed concentrations
flow_obs <- right_join(q_obs, bac_obs, by="date") #needed to reduce the number of flow observations
flow_obs_daily <- flow_obs$discharge #335
obs_flow_xts <- as.xts(flow_obs$discharge,order.by=as.Date(flow_obs$date))
flow_obs_weekly <- as.data.frame(apply.weekly(obs_flow_xts, mean)) #204
flow_obs_monthly <- as.data.frame(apply.monthly(obs_flow_xts, mean)) #51

#create weekly flux output
flux_obs_daily <- bac_obs_daily * flow_obs_daily * 10^4 #335
obs_flux_xts <- as.xts(flux_obs_daily, order.by=as.Date(flux_obs$date))
flux_obs_weekly <- as.data.frame(apply.weekly(obs_flux_xts, mean)) #204

#create monthly flux output
flux_obs_monthly <- as.data.frame(apply.monthly(obs_flux_xts, mean)) #51

#check dates
# daily
bac_daily_dates <- bac_obs$date; flow_daily_dates <- flow_obs$date; flux_daily_dates <- flux_obs$date
# bac_daily_dates == flow_daily_dates; flow_daily_dates == flux_daily_dates
#weekly
bac_weekly_dates <- rownames(bac_obs_weekly); flow_weekly_dates <- rownames(flow_obs_weekly); flux_weekly_dates <- rownames(flux_obs_weekly)
# bac_weekly_dates == flow_weekly_dates; flow_weekly_dates == flux_weekly_dates
#monthly
bac_monthly_dates <- rownames(bac_obs_monthly); flow_monthly_dates <- rownames(flow_obs_monthly); flux_mothly_dates <- rownames(flux_obs_monthly)

#save the weekly observation outputs
save(bac_obs_weekly, file = file.path(base_dir, "bac_obs_w.RData"))
save(flow_obs_weekly, file = file.path(base_dir, "q_obs_w.RData"))
save(flux_obs_weekly, file = file.path(base_dir, "flux_obs_w.RData"))
save(bac_obs_monthly, file = file.path(base_dir, "bac_obs_m.RData"))
save(flow_obs_monthly, file = file.path(base_dir, "q_obs_m.RData"))
save(flux_obs_monthly, file = file.path(base_dir, "flux_obs_m.RData"))

# preset the generations to be simulated
# stargen = 0 means starting from scratch
startgen <- 0
ngens <- 40

# decide what time frame/interval to optimize on 
# opt_time_interval <- "daily"
#opt_time_interval <- "weekly"
opt_time_interval <- "monthly"

# should concentrations be modified (modified nash sutcliffe) 
# opt_conc_transform <- "none"
opt_conc_transform <- "modified"

# decide what metric to optimize on
opt_nse <- "mean"
#opt_nse <- "conc"
#opt_nse <- "flux"
#opt_nse <- "flow"

# every generation will have n_to_keep accepted particles
# the median score of these n_to_keep particles will be used as the cutoff for the next generation
# the first generation will be the top n_to_keep of nsims_todo simulations
# start the loop here
for(iter in startgen:ngens){
  print(paste("*********** begin generation", iter, "**********************"))
  print(paste("optimizing based on", opt_time_interval, opt_conc_transform, "(transformation)", opt_nse))
  #
  #number to keep each generation
  nsims_todo <- 100
  n_to_keep <- 20
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
  # save output to disk
  save_bac_cal_output(iter, bac_cal_output)
  
  ###### get parameter names and values
  sim_pars <- bac_cal_output$parameter$values
  
  ###### subset simulated bacteria data to observed days and average by week
  # extract only the observed days from the simulated daily output
  bac_sims_all_days <- bac_cal_output$simulation$bac_out # [3865,2215]
  #dim(bac_sims_all_days)
  # adds date and bacteria fields also
  bac_sims_daily_temp <- right_join(bac_sims_all_days, bac_obs, by="date") #[336,2216]
  #dim(bac_sims_daily_temp)
  bac_sims_daily <- bac_sims_daily_temp[,-which(colnames(bac_sims_daily_temp)=="bacteria")] #[336,2215]
  #dim(bac_sims_daily)
  #colnames(bac_sims_daily)[1]
  # then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
  #head(colnames(bac_sims_daily)) #  date 
  nsim_cols <- ncol(bac_sims_daily) #2215 date + sims field
  bac_sims_daily_data <- as.xts(bac_sims_daily[2:nsim_cols],order.by=as.Date(bac_sims_daily$date)) #[336,2214]
  bac_sims_daily <- bac_sims_daily[,-1] #[336,2215]
  #dim(bac_sims_daily) #[336,2214]
  bac_sims_weekly <- as.data.frame(apply.weekly(bac_sims_daily_data,mean)) #[204,2214]
  bac_sims_monthly <- as.data.frame(apply.monthly(bac_sims_daily_data,mean)) #51
  

  ###### subset simulated flow to observed days and average by week
  # extract only the observed days from the simulated daily output
  bac_flows_all_days <- bac_cal_output$simulation$q_out # [3865,2215]
  #dim(bac_flows_all_days)
  # adds date and flow fields also
  bac_flows_daily_temp <- right_join(bac_flows_all_days, flow_obs, by="date") #[336,2217]
  #dim(bac_flows_daily_temp)
  bac_flows_daily <- bac_flows_daily_temp[,-which((colnames(bac_flows_daily_temp)=="bacteria" | 
                                                     colnames(bac_flows_daily_temp)=="discharge"))] #[336,2215]
  #dim(bac_flows_daily)
  # then reduce daily simulated observations to weekly averages for each of the sims #[336,2215]
  #head(colnames(bac_flows_daily)) #  date 
  nsim_cols <- ncol(bac_flows_daily) #2215 date + sims field
  bac_flows_daily_data <- as.xts(bac_flows_daily[2:nsim_cols],order.by=as.Date(bac_flows_daily$date)) #[336,2214]
  bac_flows_daily <- bac_flows_daily[,-1] #[336,2215]
  #dim(bac_flows_daily) #[336,2214]
  bac_flows_weekly <- as.data.frame(apply.weekly(bac_flows_daily_data,mean)) #[203,2214]
  #dim(bac_flows_weekly)
  bac_flows_monthly <- as.data.frame(apply.monthly(bac_flows_daily_data,mean)) #51
  ###### calculate simulated flux data for observed days and average by week
  #dim(bac_sims_weekly)
  #dim(bac_flows_weekly)
  bac_fluxes_weekly <- bac_sims_weekly * bac_flows_weekly * 10^4
  #dim(bac_fluxes_weekly) #[204,2214]
  bac_fluxes_monthly <- bac_sims_monthly * bac_flows_monthly * 10^4
  # View various data states
  #View(bac_sims_weekly)
  #View(bac_obs_weekly)
  #View(bac_obs_monthly)
  
  #View(bac_flows_weekly)
  #View(bac_flows_monthly)
  #View(flux_obs_weekly) 
  
  #View(bac_fluxes_weekly)
  #View(bac_fluxes_monthly)
  #View(flux_obs_weekly)  
  
  ###### calculate various nses for daily data
  nse_bac_daily <- calculate_nse_bac_daily(iter, bac_cal_output, bac_obs)
  nse_flow_daily <- calculate_nse_flow_daily(iter, bac_cal_output, q_obs)
  nse_flux_daily <- calculate_nse_flux_daily(iter, bac_cal_output, flux_obs)
  # calculate various modified nses for daily data
  mnse_bac_daily <- calculate_mnse_bac_daily(iter, bac_cal_output, bac_obs)
  mnse_flow_daily <- calculate_mnse_flow_daily(iter, bac_cal_output, q_obs)
  mnse_flux_daily <- calculate_mnse_flux_daily(iter, bac_cal_output, flux_obs)
  # calculate various nses for weekly data
  nse_bac_weekly <- calculate_nse_bac_weekly(iter, bac_sims_weekly, bac_obs_weekly)
  nse_flow_weekly <- calculate_nse_flow_weekly(iter, bac_flows_weekly, flow_obs_weekly)
  nse_flux_weekly <- calculate_nse_flux_weekly(iter, bac_fluxes_weekly, flux_obs_weekly)
  # calculate various modified nses for weekly data
  mnse_bac_weekly <- calculate_mnse_bac_weekly(iter, bac_sims_weekly, bac_obs_weekly)
  mnse_flow_weekly <- calculate_mnse_flow_weekly(iter, bac_flows_weekly, flow_obs_weekly)
  mnse_flux_weekly <- calculate_mnse_flux_weekly(iter, bac_fluxes_weekly, flux_obs_weekly) #?? not working
  # calculate various nses for monthly data
  nse_bac_monthly <- calculate_nse_bac_monthly(iter, bac_sims_monthly, bac_obs_monthly)
  nse_flow_monthly <- calculate_nse_flow_monthly(iter, bac_flows_monthly, flow_obs_monthly)
  nse_flux_monthly <- calculate_nse_flux_monthly(iter, bac_fluxes_monthly, flux_obs_monthly)
  # calculate various modified nses for monthly data
  mnse_bac_monthly <- calculate_mnse_bac_monthly(iter, bac_sims_monthly, bac_obs_monthly)
  mnse_flow_monthly <- calculate_mnse_flow_monthly(iter, bac_flows_monthly, flow_obs_monthly)
  mnse_flux_monthly <- calculate_mnse_flux_monthly(iter, bac_fluxes_monthly, flux_obs_monthly) #?? not working
  ######### calculate means of nses
  # calculate nse means
  print("NSE mean, daily")
  nse_mean_daily <- calculate_nse_mean(iter, nse_bac_daily, nse_flow_daily, nse_flux_daily)
  print("NSE mean, weekly")
  nse_mean_weekly <- calculate_nse_mean(iter, nse_bac_weekly, nse_flow_weekly, nse_flux_weekly)
  print("NSE mean, monthly")
  nse_mean_monthly <- calculate_nse_mean(iter, nse_bac_monthly, nse_flow_monthly, nse_flux_monthly)
  # calculate modified nse means
  print("mNSE, daily")
  mnse_mean_daily <- calculate_mnse_mean(iter, mnse_bac_daily, mnse_flow_daily, mnse_flux_daily) 
  print("mNSE, weekly")
  mnse_mean_weekly <- calculate_mnse_mean(iter, mnse_bac_weekly, mnse_flow_weekly, mnse_flux_weekly)
  print("mNSE, monthly")
  mnse_mean_monthly <- calculate_mnse_mean(iter, mnse_bac_monthly, mnse_flow_monthly, mnse_flux_monthly)
  
  ###### get cutoff score
  #####need to double check, this function seems not be used (iter==0)...
  if(iter==0){
    # find the 80th percentile of target nse, top (2000 of 10000)
    print(paste0("generation 0, using the 80th percentile of the scores as the cutoff for the next generation"))
    # daily, none
    if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="none"){
      this_cutoff_score <- quantile(nse_mean_daily, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_bac_daily, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flow_daily, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flux_daily, probs=0.8)
    # daily, modified
    } else if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="modified"){
      this_cutoff_score <- quantile(mnse_mean_daily, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_bac_daily, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_flow_daily, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_flux_daily, probs=0.8)  
    # weekly, none
    } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="none"){
      this_cutoff_score <- quantile(nse_mean_weekly, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_bac_weekly, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flow_weekly, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
      this_cutoff_score <- quantile(nse_flux_weekly, probs=0.8)  
    # weekly, modified
    } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="modified"){
      this_cutoff_score <- quantile(mnse_mean_weekly, probs=0.8)
    } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_bac_weekly, probs=0.8)
    } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_flow_weekly, probs=0.8)
    } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
      this_cutoff_score <- quantile(mnse_flux_weekly, probs=0.8)  
    }
    # monthly, none
  } else if(opt_nse=="mean" && opt_time_interval=="monthly" && opt_conc_transform=="none"){
    this_cutoff_score <- quantile(nse_mean_monthly, probs=0.8)
  } else if(opt_nse=="conc" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    this_cutoff_score <- quantile(nse_bac_monthly, probs=0.8)
  } else if(opt_nse=="flow" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    this_cutoff_score <- quantile(nse_flow_monthly, probs=0.8)
  } else if(opt_nse=="flux" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    this_cutoff_score <- quantile(nse_flux_monthly, probs=0.8)  
    # monthly, modified
  } else if(opt_nse=="mean" && opt_time_interval=="monthly" && opt_conc_transform=="modified"){
    this_cutoff_score <- quantile(mnse_mean_monthly, probs=0.8)
  } else if(opt_nse=="conc" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    this_cutoff_score <- quantile(mnse_bac_monthly, probs=0.8)
  } else if(opt_nse=="flow" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    this_cutoff_score <- quantile(mnse_flow_monthly, probs=0.8)
  } else if(opt_nse=="flux" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    this_cutoff_score <- quantile(mnse_flux_monthly, probs=0.8)  
  }
    print(paste("generation 0 is based on the", opt_nse, "output;", opt_time_interval, "interval;",
          opt_conc_transform, "change to nse."))
    print(paste("Cutoff score to be applied to generation 0 =", this_cutoff_score))
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
    all_keepers <- which(nse_flow_daily > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flux_daily > this_cutoff_score)
  # daily, modified
  } else if(opt_nse=="mean" && opt_time_interval=="daily" && opt_conc_transform=="modified"){
    all_keepers <- which(mnse_mean_daily > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_bac_daily > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flow_daily > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="daily" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flux_daily > this_cutoff_score) 
  # weekly, none
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="none"){
    all_keepers <- which(nse_mean_weekly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flow_weekly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flux_weekly > this_cutoff_score)  
  # weekly, modified
  } else if(opt_nse=="mean" && opt_time_interval=="weekly" && opt_conc_transform=="modified"){
    all_keepers <- which(mnse_mean_weekly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_bac_weekly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flow_weekly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="weekly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flux_weekly > this_cutoff_score)  
#monthly, none 
  }else if(opt_nse=="mean" && opt_time_interval=="monthly" && opt_conc_transform=="none"){
    all_keepers <- which(nse_mean_monthly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_bac_monthly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flow_monthly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="monthly" && opt_conc_transform=="none") {
    all_keepers <- which(nse_flux_monthly > this_cutoff_score)  
    # monthly, modified
  } else if(opt_nse=="mean" && opt_time_interval=="monthly" && opt_conc_transform=="modified"){
    all_keepers <- which(mnse_mean_monthly > this_cutoff_score)
  } else if(opt_nse=="conc" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_bac_monthly > this_cutoff_score)
  } else if(opt_nse=="flow" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flow_monthly > this_cutoff_score)
  } else if(opt_nse=="flux" && opt_time_interval=="monthly" && opt_conc_transform=="modified") {
    all_keepers <- which(mnse_flux_monthly > this_cutoff_score)  
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
  # concatenate output into one big dataframe
  # daily, none
  if(opt_time_interval=="daily" && opt_conc_transform=="none"){
    nses_w_parameters_all <- cbind(keeper, nse_bac_daily, nse_flux_daily, nse_flow_daily, 
                               nse_mean_daily, sim_pars)
    nse_conc_keepers <- nse_bac_daily[valid_keepers]
    nse_flow_keepers <- nse_flow_daily[valid_keepers]
    nse_flux_keepers <- nse_flux_daily[valid_keepers]
    nse_mean_keepers <- nse_mean_daily[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
  # daily, modified
  } else if(opt_time_interval=="daily" && opt_conc_transform=="modified"){
    nses_w_parameters_all <- cbind(keeper, mnse_bac_daily, mnse_flux_daily, mnse_flow_daily, 
                                   mnse_mean_daily, sim_pars)
    nse_conc_keepers <- mnse_bac_daily[valid_keepers]
    nse_flow_keepers <- mnse_flow_daily[valid_keepers]
    nse_flux_keepers <- mnse_flux_daily[valid_keepers]
    nse_mean_keepers <- mnse_mean_daily[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
  # weekly, none
  } else if(opt_time_interval=="weekly" && opt_conc_transform=="none"){
    nses_w_parameters_all <- cbind(keeper, nse_bac_weekly, nse_flux_weekly, nse_flow_weekly, 
                                   nse_mean_weekly, sim_pars)
    nse_conc_keepers <- nse_bac_weekly[valid_keepers]
    nse_flow_keepers <- nse_flow_weekly[valid_keepers]
    nse_flux_keepers <- nse_flux_weekly[valid_keepers]
    nse_mean_keepers <- nse_mean_weekly[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
  # weekly, modified
  } else if(opt_time_interval=="weekly" && opt_conc_transform=="modified"){
    nses_w_parameters_all <- cbind(keeper, mnse_bac_weekly, mnse_flux_weekly, mnse_flow_weekly, 
                                   mnse_mean_weekly, sim_pars)
    nse_conc_keepers <- mnse_bac_weekly[valid_keepers]
    nse_flow_keepers <- mnse_flow_weekly[valid_keepers]
    nse_flux_keepers <- mnse_flux_weekly[valid_keepers]
    nse_mean_keepers <- mnse_mean_weekly[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
    # monthly, none
  } else if(opt_time_interval=="monthly" && opt_conc_transform=="none"){
    nses_w_parameters_all <- cbind(keeper, nse_bac_monthly, nse_flux_monthly, nse_flow_monthly, 
                                   nse_mean_monthly, sim_pars)
    nse_conc_keepers <- nse_bac_monthly[valid_keepers]
    nse_flow_keepers <- nse_flow_monthly[valid_keepers]
    nse_flux_keepers <- nse_flux_monthly[valid_keepers]
    nse_mean_keepers <- nse_mean_monthly[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
    # monthly, modified
  } else if(opt_time_interval=="monthly" && opt_conc_transform=="modified"){
    nses_w_parameters_all <- cbind(keeper, mnse_bac_monthly, mnse_flux_monthly, mnse_flow_monthly, 
                                   mnse_mean_monthly, sim_pars)
    nse_conc_keepers <- mnse_bac_monthly[valid_keepers]
    nse_flow_keepers <- mnse_flow_monthly[valid_keepers]
    nse_flux_keepers <- mnse_flux_monthly[valid_keepers]
    nse_mean_keepers <- mnse_mean_monthly[valid_keepers]
    nses_w_parameters <- cbind(nse_conc_keepers, nse_flow_keepers, nse_flux_keepers, 
                               nse_mean_keepers, sim_pars[valid_keepers,])
  }

  # plot nses versus each other
  # plot_bac_v_flow_pdf(iter, nses_w_parameters_all) #TEMPORARILY DISABLING
  
  #save nses_parameters
  save_nses_parameters(iter, data_in_dir, nses_w_parameters)
  # save concentration time series output to an .RData file 
  # for later sensitivity analyses TODO
  # delete the local bac_cal file TODO
  # calculate cutoff score
  # USING MEDIAN!!!

  # find relevant cutoff and save for next generation
  # nse_*_keepers saves appropriate calcualtion
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

