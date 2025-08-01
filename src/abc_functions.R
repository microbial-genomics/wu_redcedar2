options("width"=132)

calculate_next_nsims <- function(n_to_keep, proportion_kept){
  new_nsims <- (round(n_to_keep/proportion_kept)*2)+500
  print(paste("next round we will do", new_nsims, "simulations"))
  return(new_nsims)
}

### bacteria nash-sutcliffes (6)
calculate_mnse_bac_daily <- function(iter, bac_cal_output, bac_obs){
  #load the simulated concentrations for last simulations
  sim_bac <- bac_cal_output$simulation$bac_out
  # merge simulated and observed bacteria concentrations, calculate nses for all sims
  nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
    dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
    map_dbl(., ~mNSE(.x, bac_obs$bacteria))##the same with log???
  print(paste("range of all daily bacteria modified nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

calculate_mnse_bac_weekly <- function(iter, bac_sims_weekly, bac_obs_weekly){
  nse_bac <- mapply(mNSE, bac_sims_weekly, bac_obs_weekly)
  print(paste("range of all weekly bacteria modified nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

calculate_mnse_bac_monthly <- function(iter, bac_sims_monthly, bac_obs_monthly){
  nse_bac <- mapply(mNSE, bac_sims_monthly, bac_obs_monthly)
  print(paste("range of all monthly bacteria modified nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

calculate_nse_bac_daily <- function(iter, bac_cal_output, bac_obs){
  #load the simulated concentrations for last simulations
  sim_bac <- bac_cal_output$simulation$bac_out
  # merge simulated and observed bacteria concentrations, calculate nses for all sims
  nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
    dplyr::select(-date) %>% dplyr::select(-bacteria) %>%
    map_dbl(., ~NSE(.x, bac_obs$bacteria))
  print(paste("range of daily bacteria nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

calculate_nse_bac_weekly <- function(iter, bac_sims_weekly, bac_obs_weekly){
  nse_bac <- mapply(NSE, bac_sims_weekly, bac_obs_weekly)
  print(paste("range of all weekly bacteria nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

calculate_nse_bac_monthly <- function(iter, bac_sims_monthly, bac_obs_monthly){
  nse_bac <- mapply(NSE, bac_sims_monthly, bac_obs_monthly)
  print(paste("range of all monthly bacteria nse is (", round(min(nse_bac),4), ",", round(max(nse_bac),4), ") for generation", iter))
  return(nse_bac)
}

### flow nash-sutcliffes (6)
calculate_mnse_flow_daily <- function(iter, bac_cal_output, q_obs){
  #load the simulated flows, inputs for last simulations
  sim_q <- bac_cal_output$simulation$q_out
  # merge simulated and observed flows, calculate nses for all sims
  nse_q <- right_join(sim_q,q_obs,by="date") %>%
    dplyr::select(-date) %>% dplyr::select(-discharge) %>%
    map_dbl(., ~mNSE(.x, q_obs$discharge))
  print(paste("range of all overall flow modified nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

calculate_mnse_flow_weekly <- function(iter, bac_flows_weekly, flow_obs_weekly){
  nse_q <- mapply(mNSE, bac_flows_weekly, flow_obs_weekly)
  print(paste("range of all weekly flow modified nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

calculate_mnse_flow_monthly <- function(iter, bac_flows_monthly, flow_obs_monthly){
  nse_q <- mapply(mNSE, bac_flows_monthly, flow_obs_monthly)
  print(paste("range of all monthly flow modified nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

calculate_nse_flow_daily <- function(iter, bac_cal_output, q_obs){
  #load the simulated flows, inputs for last simulations
  sim_q <- bac_cal_output$simulation$q_out
  # merge simulated and observed flows, calculate nses for all sims
  nse_q <- right_join(sim_q,q_obs,by="date") %>%
    dplyr::select(-date) %>% dplyr::select(-discharge) %>%
    map_dbl(., ~NSE(.x, q_obs$discharge))
  print(paste("range of all overall flow nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

calculate_nse_flow_weekly <- function(iter, bac_flows_weekly, flow_obs_weekly){
  nse_q <- mapply(NSE, bac_flows_weekly, flow_obs_weekly)
  print(paste("range of all weekly flow nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

calculate_nse_flow_monthly <- function(iter, bac_flows_monthly, flow_obs_monthly){
  nse_q <- mapply(NSE, bac_flows_monthly, flow_obs_monthly)
  print(paste("range of all monthly flow modified nse is (", round(min(nse_q),4), ",", round(max(nse_q),4), ") for generation", iter))
  return(nse_q)
}

### flux nash-sutcliffes (6)
calculate_mnse_flux_daily <- function(iter, bac_cal_output, flux_obs){  
  #load the simulated concentrations for last simulations
  sim_bac <- bac_cal_output$simulation$bac_out
  sim_q <- bac_cal_output$simulation$q_out
  # calculate the simulated fluxes from concs and flows for all sims
  # flux_sim <- sim_bac[,c(-1)]*
  #   sim_q[, c(-1)]*10^4
  date <-bac_cal_output$simulation$bac_out$date
  flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
  sim_flux<- cbind(date, flux_sim)
  #merge simulated and observed fluxes, calculate nses for all sims
  nse_flux <-  right_join(sim_flux, flux_obs, by = "date") %>%
    dplyr::select(-date) %>%dplyr::select(-flux) %>%
    map_dbl(., ~mNSE(.x, flux_obs$flux))
  print(paste("range of all overall flux modified nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)
}

calculate_mnse_flux_weekly <- function(iter, bac_fluxes_weekly, flux_obs_weekly){
  nse_flux <- mapply(mNSE, bac_fluxes_weekly, flux_obs_weekly)
  print(paste("range of all weekly bacteria flux modified nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)##changed nse_bac to nse_flux
}

calculate_mnse_flux_monthly <- function(iter, bac_fluxes_monthly, flux_obs_monthly){
  nse_flux <- mapply(mNSE, bac_fluxes_monthly, flux_obs_monthly)
  print(paste("range of all monthly bacteria flux modified nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)##changed nse_bac to nse_flux
}

calculate_nse_flux_daily <- function(iter, bac_cal_output, flux_obs){  
  #load the simulated concentrations for last simulations
  sim_bac <- bac_cal_output$simulation$bac_out
  sim_q <- bac_cal_output$simulation$q_out
  # calculate the simulated fluxes from concs and flows for all sims
  # flux_sim <- sim_bac[,c(-1)]*
  #   sim_q[, c(-1)]*10^4
  date <-bac_cal_output$simulation$bac_out$date
  flux_sim <-sim_bac[,c(-1)]* sim_q[, c(-1)]*10^4
  sim_flux<- cbind(date, flux_sim)
  #merge simulated and observed fluxes, calculate nses for all sims
  nse_flux <-  right_join(sim_flux, flux_obs, by = "date") %>%
    dplyr::select(-date) %>%dplyr::select(-flux) %>%
    map_dbl(., ~NSE(.x, flux_obs$flux))
  print(paste("range of all overall flux nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)
}

calculate_nse_flux_weekly <- function(iter, bac_fluxes_weekly, flux_obs_weekly){
  nse_flux <- mapply(NSE, bac_fluxes_weekly, flux_obs_weekly)##this one might be confusing nse_bac?
  print(paste("range of all weekly bacteria flux nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)
}

calculate_nse_flux_monthly <- function(iter, bac_fluxes_monthly, flux_obs_monthly){
  nse_flux <- mapply(NSE, bac_fluxes_monthly, flux_obs_monthly)
  print(paste("range of all monthly bacteria flux modified nse is (", round(min(nse_flux),4), ",", round(max(nse_flux),4), ") for generation", iter))
  return(nse_flux)##changed nse_bac to nse_flux
}

### mean of (conc, flow, flux) nash-sutcliffes (2)
#function does not care if daily or weekly or monthly
calculate_mnse_mean <- function(iter, mnse_bac, mnse_q, mnse_flux){
  # calculate nse means
  mnse_mean <- rowMeans(cbind(mnse_bac, mnse_q, mnse_flux))
  print(paste("range of all overall modified mean nse is (", round(min(mnse_mean),4), ",", round(max(mnse_mean),4), ") for generation", iter))
  return(mnse_mean)
}

calculate_nse_mean <- function(iter, nse_bac, nse_q, nse_flux){
  # calculate nse means
  nse_mean <- rowMeans(cbind(nse_bac, nse_q, nse_flux))
  print(paste("range of all overall mean nse is (", round(min(nse_mean),4), ",", round(max(nse_mean),4), ") for generation", iter))
  return(nse_mean)
}

##
create_next_sim_tibble <- function(nsims_todo, parameter_input_sims){
  
  if(nsims_todo > nrow(parameter_input_sims)){
    print(paste("major problem, more sims requested than parameter inputs"))
  }else{
    next_tibble <- tibble(
      "CN2.mgt|change = relchg"= parameter_input_sims[1:nsims_todo,1], #$CN2,initial SCS runoff curve number for mosture condition II
      "OV_N.hru|change = relchg" = parameter_input_sims[1:nsims_todo,2],
      "GW_DELAY.gw|change = relchg" = parameter_input_sims[1:nsims_todo, 3],
      "GWQMN.gw|change = relchg" = parameter_input_sims[1:nsims_todo,4], #$GWQMN, Threshold depth of water in the shallow aquifer required for return flow to occruemm (H2O)
      "SLSUBBSN.hru|change = relchg" = parameter_input_sims[1:nsims_todo,5],
      "ALPHA_BNK.rte|change = absval" = parameter_input_sims[1:nsims_todo,6], #$ALPHA_BNK, baseflow apha factor for bank storage (days)
      "CH_K2.rte|change = absval" = parameter_input_sims[1:nsims_todo,7], #$CH_K2,Effective hydraulic conductivity in main channel alluvium (mm/hr)
      "CH_N2.rte|change = absval" = parameter_input_sims[1:nsims_todo,8], #$CH_N2,Manning's n value for the main channels
      "TRNSRCH.bsn|change = absval" = parameter_input_sims[1:nsims_todo,9], #$TRNSRCH,Fraction of trasmission losses from main channel that enter deep aquifer
      "CH_K1.sub|change = absval" = parameter_input_sims[1:nsims_todo, 10], #$CH_K1, Effective hydraulic conductivity in tributary channel alluvium (mm/hr)
      "GW_REVAP.gw|change = absval" = parameter_input_sims[1:nsims_todo, 11],
      "RCHRG_DP.gw|change = absval" = parameter_input_sims[1:nsims_todo,12], #RCHRG_DP, Deep aquifer percolation fraction 
      "SFTMP.bsn|change = absval"= parameter_input_sims[1:nsims_todo,13], #SFTMP, snowfall temperature C
      "SMTMP.bsn|change = absval"= parameter_input_sims[1:nsims_todo,14], #SMTMP, snowmelt temperature C
      "TIMP.bsn|change = absval"= parameter_input_sims[1:nsims_todo,15],
      "DEP_IMP.hru|change = absval"= parameter_input_sims[1:nsims_todo,16], #DEP_IMP, Depth to impervious layer in soil profile (mm)
      "DDRAIN.mgt|change = absval"= parameter_input_sims[1:nsims_todo,17], #DDRAIN, Depth to subsurface drian (mm)
      "GDRAIN.mgt|change = absval"= parameter_input_sims[1:nsims_todo,18], #GDRAIN, Drain tile lab time (hours)
      "THBACT.bsn|change = absval"= parameter_input_sims[1:nsims_todo,19], #THBACT, Temperature adjsument factor for bacteria die-off/growth
      "WDPRCH.bsn|change = absval"= parameter_input_sims[1:nsims_todo,20]) #WDPRCH, Die-off factor for persistent bacteria in streas  at 20C (1/day) 
    return(next_tibble)
  }
}

create_generation_stats <- function(startgen, ngens, opt_nse, nsims){
  n <- ngens + 1
  df <- data.frame(generation = integer(n),
                   nsims = integer(n),
                   proportion_kept = double(n),
                   opt_variable = character(n),
                   best_mean_nse = double(n),
                   best_bac_nse = double(n),
                   best_q_nse = double(n),
                   best_flux_nse = double(n),
                   cutoff_score = double(n),
                   stringsAsFactors=FALSE)
  df$generation <- seq(startgen, ngens)
  df$nsims[1] <- nsims
  df$proportion_kept[1] <- 1.0
  df$opt_variable[1] <- opt_nse
  df$cutoff_score[1] <- -1e+12 #accept everything first time through
  print(paste("The cutoff score for generation -1 will be -1e+12 so that we accept everything in generation 0."))
  return(df)
}

create_tibble_initial <- function(nsims){
  pars_tibble <- tibble(#hydrology parameters (15)
    "CN2.mgt|change = relchg"= runif(nsims, -0.25, 0.1), #
    "OV_N.hru|change = relchg" = runif(nsims, 0.01, 0.4),
    "GW_DELAY.gw|change = relchg" = runif(nsims, -0.75,4),
    "GWQMN.gw|change = relchg" = runif(nsims, -0.5, 0.5), #
    "SLSUBBSN.hru|change = relchg" = runif(nsims, -0.5, 1),
    "ALPHA_BNK.rte|change = absval" =runif(nsims, 0.5, 1), #large for flat recessions, and small for steep recessions
    "CH_K2.rte|change = absval" = runif(nsims, 0, 50), # changed from(0,250) 
    "CH_N2.rte|change = absval" = runif(nsims, 0.05, 0.15), # changed from (0,0.1) 
    "TRNSRCH.bsn|change = absval" = runif(nsims, 0, 0.3), # default is 0.00
    "CH_K1.sub|change = absval" = runif(nsims, 0, 120), #For prennial streams with continuous groundwater contribution, the effective conductivity will be zero.
    "GW_REVAP.gw|change = absval" = runif(nsims, 0.02, 0.2),
    "RCHRG_DP.gw|change = absval" = runif(nsims, 0, 0.5), # default range is (0,1)
    "SFTMP.bsn|change = absval"= runif(nsims, -5, 5), # changed from (-2,2),default is 1.0
    "SMTMP.bsn|change = absval"= runif(nsims, -5, 5), # change from (-2,2), default is 0.5
    "TIMP.bsn|change = absval"= runif(nsims, 0.01, 1),
    
    #tile drainage and sediments (3)
    "DEP_IMP.hru|change = absval"= runif(nsims, 2000, 6000), # "generation 21 DEP_IMP 4010.712 31.167"
    "DDRAIN.mgt|change = absval"= runif(nsims, 500, 1500), # "generation 21 DDRAIN 1105.285 107.537"
    "GDRAIN.mgt|change = absval"= runif(nsims, 0, 50), # "generation 21 GDRAIN 24.686 11.354"
    #bacteria submodel (2)
    "THBACT.bsn|change = absval"= runif(nsims, 0, 2), # default value 1.07
    "WDPRCH.bsn|change = absval"= runif(nsims, 0, 1) # 
  )
}


fit_normal_parameters <- function(sim_pars_keepers){
  fitted_CN2 <- fitdist(sim_pars_keepers$CN2, "norm")
  fitted_OV_N <- fitdist(sim_pars_keepers$OV_N, "norm")
  fitted_GW_DELAY <- fitdist(sim_pars_keepers$GW_DELAY, "norm")
  fitted_GWQMN <- fitdist(sim_pars_keepers$GWQMN, "norm")
  fitted_SLSUBBSN <- fitdist(sim_pars_keepers$SLSUBBSN, "norm")
  fitted_ALPHA_BNK <- fitdist(sim_pars_keepers$ALPHA_BNK, "norm")
  fitted_CH_K2 <- fitdist(sim_pars_keepers$CH_K2, "norm")
  fitted_CH_N2 <- fitdist(sim_pars_keepers$CH_N2, "norm")
  fitted_TRNSRCH <- fitdist(sim_pars_keepers$TRNSRCH, "norm")
  fitted_CH_K1 <- fitdist(sim_pars_keepers$CH_K1, "norm")
  fitted_GW_REVAP <- fitdist(sim_pars_keepers$GW_REVAP, "norm")
  fitted_RCHRG_DP <- fitdist(sim_pars_keepers$RCHRG_DP, "norm")
  fitted_SFTMP <- fitdist(sim_pars_keepers$SFTMP, "norm")
  fitted_SMTMP <- fitdist(sim_pars_keepers$SMTMP, "norm")
  fitted_TIMP <- fitdist(sim_pars_keepers$TIMP, "norm")
  fitted_DEP_IMP <- fitdist(sim_pars_keepers$DEP_IMP, "norm")
  fitted_DDRAIN <- fitdist(sim_pars_keepers$DDRAIN, "norm")
  fitted_GDRAIN <- fitdist(sim_pars_keepers$GDRAIN, "norm")
  fitted_THBACT <- fitdist(sim_pars_keepers$THBACT, "norm")
  fitted_WDPRCH <- fitdist(sim_pars_keepers$WDPRCH, "norm")
  return(list(fitted_CN2, fitted_OV_N, fitted_GW_DELAY, fitted_GWQMN, fitted_SLSUBBSN, fitted_ALPHA_BNK, fitted_CH_K2, fitted_CH_N2,
              fitted_TRNSRCH, fitted_CH_K1, fitted_GW_REVAP, fitted_RCHRG_DP, fitted_SFTMP, fitted_SMTMP, fitted_TIMP, fitted_DEP_IMP, fitted_DDRAIN, fitted_GDRAIN, fitted_THBACT, fitted_WDPRCH))
}

get_cutoff_score <- function(iter, generation_stats){
  # load and retrieve previously saved cutoff
  return_cutoff_score <- generation_stats$cutoff_score[iter+1] #get the cutoff that was saved at the last generation
  print(paste("cutoff score from the last generation to be used for current generation ", iter, " = ", return_cutoff_score))
  return(return_cutoff_score)
}

load_parameter_input_sims <- function(iter, data_in_dir){
  parameter_input_sims_filename <- file.path(data_in_dir, paste("parameter_input_sims", iter-1, ".RData", sep=""))
  load(file = parameter_input_sims_filename, .GlobalEnv)
  print(paste("parameter inputs from generation", iter-1, "loaded from file:", parameter_input_sims_filename))  
}

load_previous_swat_simulations <- function(iter, data_in_dir){
  bac_cal_filename <- paste('bac_cal', iter-1, '.RData', sep="")
  rdata_file_in <- file.path(data_in_dir, bac_cal_filename)
  print(paste("loading data file: ", rdata_file_in))
  load(file = rdata_file_in, .GlobalEnv)
}

load_generation_stats <- function(iter, data_in_dir){
  generation_stats_filename <- file.path(data_in_dir,paste("generation_stats",iter-1,".RData",sep=""))
  load(file= generation_stats_filename, .GlobalEnv)
  print(paste("generation stats for generation", iter-1, "loaded from", generation_stats_filename))
}

load_observations <- function(){
  load(file= file.path(data_in_dir,'bac_obs.RData'), .GlobalEnv)
  load(file = file.path(data_in_dir,'flux_obs.RData'), .GlobalEnv)
  #load(file = file.path(data_in_dir,'pcp_obs.RData'), .GlobalEnv)
  #load(file = file.path(data_in_dir,'pcp_obs2.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'q_obs.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'q_obs_w.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'bac_obs_w.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'flux_obs_w.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'q_obs_m.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'bac_obs_m.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'flux_obs_m.RData'), .GlobalEnv)
}

log_results <- function(iter, this_cutoff_score, n_all_keepers, previous_nsims, nse_mean_keepers,
                        nse_bac_keepers, nse_q_keepers, nse_flux_keepers, next_cutoff_score){
  print(paste("Generation ", iter))
  print(paste("cutoff score from the last generation was:", this_cutoff_score))
  proportion_kept <- n_all_keepers/previous_nsims
  print(paste("generation x:",n_all_keepers, "of", format(previous_nsims,scientific=F), "simulations kept; proportion kept =", round(proportion_kept,4)))
  print(paste("range of kept bacteria scores is (", round(min(nse_bac_keepers),4), ",", round(max(nse_bac_keepers),4), ") for generation", iter))
  print(paste("median of kept bacteria scores is (", round(median(nse_bac_keepers),4), ") for generation", iter))
  print(paste("range of kept flow scores is (", round(min(nse_q_keepers),4), ",", round(max(nse_q_keepers),4), ") for generation", iter))
  print(paste("median of kept flow scores is (", round(median(nse_q_keepers),4), ") for generation", iter))
  print(paste("range of kept flux scores is (", round(min(nse_flux_keepers),4), ",", round(max(nse_flux_keepers),4), ") for generation", iter))
  print(paste("median of kept flux scores is (", round(median(nse_flux_keepers),4), ") for generation", iter))
  print(paste("range of kept mean scores is (", round(min(nse_mean_keepers),4), ",", round(max(nse_mean_keepers),4), ") for generation", iter))
  print(paste("median of kept mean scores is (", round(median(nse_mean_keepers),4), ") for generation", iter))
  print(paste("cutoff score to be used for the next generation is:", round(next_cutoff_score,4)))
}

plot_bac_v_flow_pdf <- function(iter, nses_w_parameters_all, nse_mean_daily){
  ggplot(nses_w_parameters_all %>% arrange(nse_mean_daily)) +
    geom_point(aes(x=nse_q, y=nse_mean_daily, color=keeper, alpha=0.5))
  nse_plot_filename <- paste("nse_bac_v_flow_gen", iter, ".pdf", sep="")
  ggsave(file.path(graphics_dir, nse_plot_filename))  
}

plot_kde_pdf <- function(iter, kde_next_gen){
  ggplot(data = kde_next_gen) +
    geom_density(aes(x = parameter_range)) +
    facet_wrap(.~par, nrow=5, scales = "free") +
    theme_bw()
  density_plot_filename <- paste("kde_mcabc_gen", iter, ".pdf", sep="")
  ggsave(file.path(graphics_dir, density_plot_filename))  
}

run_swat_red_cedar <- function(swat_path, swat_parameters){
  run_swat2012(project_path = swat_path,
               output = list(q_out = define_output(file = "rch",
                                                   variable = "FLOW_OUT",
                                                   unit = 4),
                             bac_out = define_output(file = "rch",
                                                     variable = "BACTP_OUT",
                                                     unit = 4)),
               parameter = swat_parameters,
               start_date = "2000-01-01",
               end_date = "2014-07-31",
               years_skip = 4,
               n_thread = 32)
}

sample_truncated_normals <- function(iter, new_nsims, fitted_parameter_list){
  fitted_CN2 <- fitted_parameter_list[[1]]
  CN2_mean <- fitted_CN2$estimate[1]
  CN2_sd <- fitted_CN2$estimate[2]
  print(paste("generation", iter, "CN2", round(CN2_mean,3), round(CN2_sd,3)))
  CN2 <- rtruncnorm(new_nsims, -0.25, 0.1, mean = CN2_mean, sd =  CN2_sd)
  
  fitted_OV_N <- fitted_parameter_list[[2]]
  OV_N_mean <- fitted_OV_N$estimate[1]
  OV_N_sd <- fitted_OV_N$estimate[2]
  print(paste("generation", iter, "OV_N", round(OV_N_mean,3), round(OV_N_sd,3)))
  OV_N <- rtruncnorm(new_nsims, 0.01, 0.4, mean =OV_N_mean, sd = OV_N_sd)
  
  fitted_GW_DELAY <- fitted_parameter_list[[3]]
  GW_DELAY_mean <- fitted_GW_DELAY$estimate[1]
  GW_DELAY_sd <- fitted_GW_DELAY$estimate[2]
  print(paste("generation", iter, "GW_DELAY", round(GW_DELAY_mean,3), round(GW_DELAY_sd,3)))
  GW_DELAY <- rtruncnorm(new_nsims, -0.75, 4, mean =GW_DELAY_mean, sd = GW_DELAY_sd)
  
  
  fitted_GWQMN <- fitted_parameter_list[[4]]
  GWQMN_mean <- fitted_GWQMN$estimate[1]
  GWQMN_sd <- fitted_GWQMN$estimate[2]
  print(paste("generation", iter, "GWQMN", round(GWQMN_mean,3), round(GWQMN_sd,3)))
  GWQMN <- rtruncnorm(new_nsims, -0.5, 0.5, mean =GWQMN_mean, sd = GWQMN_sd)
  
  
  fitted_SLSUBBSN <- fitted_parameter_list[[5]]
  SLSUBBSN_mean <- fitted_SLSUBBSN$estimate[1]
  SLSUBBSN_sd <- fitted_SLSUBBSN$estimate[2]
  print(paste("generation", iter, "SLSUBBSN", round(SLSUBBSN_mean,3), round(SLSUBBSN_sd,3)))
  SLSUBBSN <- rtruncnorm(new_nsims, -0.5, 1, mean =SLSUBBSN_mean, sd = SLSUBBSN_sd)
  
  # fitted_ALPHA_BNK
  # Baseflow alpha factor for bank storage
  # estimate  Std. Error
  # mean 0.4911617 0.003369632
  # sd   0.2382690 0.002382501
  fitted_ALPHA_BNK <- fitted_parameter_list[[6]]
  ALPHA_BNK_mean <- fitted_ALPHA_BNK$estimate[1]
  ALPHA_BNK_sd <- fitted_ALPHA_BNK$estimate[2]
  print(paste("generation", iter, "ALPHA_BNK", round(ALPHA_BNK_mean,3), round(ALPHA_BNK_sd,3)))
  ALPHA_BNK <- rtruncnorm(new_nsims, 0.5, 1, mean = ALPHA_BNK_mean, sd = ALPHA_BNK_sd)
  
  # fitted_CH_K2
  # Effective hydraulic conductivity in main channel alluvium
  # estimate Std. Error
  # mean 27.90991  0.1635768
  # sd   11.56662  0.1156662
  fitted_CH_K2 <- fitted_parameter_list[[7]]
  CH_K2_mean <- fitted_CH_K2$estimate[1]
  CH_K2_sd <- fitted_CH_K2$estimate[2]
  print(paste("generation", iter, "CH_K2", round(CH_K2_mean,3), round(CH_K2_sd,3)))
  CH_K2 <- rtruncnorm(new_nsims, 0, 50, mean = CH_K2_mean, sd = CH_K2_sd)
  
  # fitted_CH_N2
  # estimate   Std. Error
  # mean 0.10415042 0.0003312391
  # sd   0.02342214 0.0002323030
  fitted_CH_N2 <- fitted_parameter_list[[8]]
  CH_N2_mean <- fitted_CH_N2$estimate[1]
  CH_N2_sd <- fitted_CH_N2$estimate[2]
  print(paste("generation", iter, "CH_N2", round(CH_N2_mean,3), round(CH_N2_sd,3)))
  CH_N2 <- rtruncnorm(new_nsims, 0.05,0.15, mean = CH_N2_mean, sd = CH_N2_sd)
  
  # fitted_TRNSRCH
  # estimate   Std. Error
  # mean 0.17602745 0.0009794462
  # sd   0.06925731 0.0006919234
  fitted_TRNSRCH <- fitted_parameter_list[[9]]
  TRNSRCH_mean <- fitted_TRNSRCH$estimate[1]
  TRNSRCH_sd <- fitted_TRNSRCH$estimate[2]
  print(paste("generation", iter, "TRNSRCH", round(TRNSRCH_mean,3), round(TRNSRCH_sd,3)))
  TRNSRCH <- rtruncnorm(new_nsims, 0, 0.3, mean = TRNSRCH_mean, sd = TRNSRCH_sd)
  
  
  # fitted_CH_K1
  # estimate Std. Error
  # mean 158.93894  0.9599238
  # sd    67.87685  0.6787685
  fitted_CH_K1 <- fitted_parameter_list[[10]]
  CH_K1_mean <- fitted_CH_K1$estimate[1]
  CH_K1_sd <- fitted_CH_K1$estimate[2]
  print(paste("generation", iter, "CH_K1", round(CH_K1_mean,3), round(CH_K1_sd,3)))
  CH_K1 <- rtruncnorm(new_nsims, 0, 120, mean = CH_K1_mean, sd = CH_K1_sd)
  
  
  fitted_GW_REVAP <- fitted_parameter_list[[11]]
  GW_REVAP_mean <- fitted_GW_REVAP$estimate[1]
  GW_REVAP_sd <- fitted_GW_REVAP$estimate[2]
  print(paste("generation", iter, "GW_REVAP", round(GW_REVAP_mean,3), round(GW_REVAP_sd,3)))
  GW_REVAP <- rtruncnorm(new_nsims, 0.02,0.2, mean =GW_REVAP_mean, sd = GW_REVAP_sd)
  
  # fitted_RCHRG_DP
  # estimate  Std. Error
  # mean 0.4808134 0.003439768
  # sd   0.2432284 0.002432099
  fitted_RCHRG_DP <- fitted_parameter_list[[12]]
  RCHRG_DP_mean <- fitted_RCHRG_DP$estimate[1]
  RCHRG_DP_sd <- fitted_RCHRG_DP$estimate[2]
  print(paste("generation", iter, "RCHRG_DP", round(RCHRG_DP_mean,3), round(RCHRG_DP_sd,3)))
  RCHRG_DP <- rtruncnorm(new_nsims, 0, 0.5, mean = RCHRG_DP_mean, sd = RCHRG_DP_sd)
  
  # fitted_SFTMP
  # estimate Std. Error
  # mean 0.1230878 0.03301559
  # sd   2.3345547 0.02334553
  fitted_SFTMP <- fitted_parameter_list[[13]]
  SFTMP_mean <- fitted_SFTMP$estimate[1]
  SFTMP_sd <- fitted_SFTMP$estimate[2]
  print(paste("generation", iter, "SFTMP", round(SFTMP_mean,3), round(SFTMP_sd,3)))
  SFTMP <- rtruncnorm(new_nsims, -5, 5, mean = SFTMP_mean, sd = SFTMP_sd)
  
  # fitted_SMTMP
  # estimate Std. Error
  # mean 0.1653056 0.03339191
  # sd   2.3611645 0.02361163
  fitted_SMTMP <- fitted_parameter_list[[14]]
  SMTMP_mean <- fitted_SMTMP$estimate[1]
  SMTMP_sd <- fitted_SMTMP$estimate[2]
  print(paste("generation", iter, "SMTMP", round(SMTMP_mean,3), round(SMTMP_sd,3)))
  SMTMP <- rtruncnorm(new_nsims, -5, 5, mean = SMTMP_mean, sd = SMTMP_sd)
  
  
  fitted_TIMP <- fitted_parameter_list[[15]]
  TIMP_mean <- fitted_TIMP$estimate[1]
  TIMP_sd <- fitted_TIMP$estimate[2]
  print(paste("generation", iter, "TIMP", round(TIMP_mean,3), round(TIMP_sd,3)))
  TIMP <- rtruncnorm(new_nsims, 0.01, 1, mean =TIMP_mean, sd = TIMP_sd)
  
  # fitted_DEP_IMP
  # estimate Std. Error
  # mean 3667.422   18.83298
  # sd   1331.808   13.31908
  fitted_DEP_IMP <- fitted_parameter_list[[16]]
  DEP_IMP_mean <- fitted_DEP_IMP$estimate[1]
  DEP_IMP_sd <- fitted_DEP_IMP$estimate[2]
  print(paste("generation", iter, "DEP_IMP", round(DEP_IMP_mean,3), round(DEP_IMP_sd,3)))
  DEP_IMP <- rtruncnorm(new_nsims, 2000, 6000, mean = DEP_IMP_mean, sd = DEP_IMP_sd)
  
  # fitted_DDRAIN
  # estimate Std. Error
  # mean 902.7242   6.710620
  # sd   474.4889   4.744931
  fitted_DDRAIN <- fitted_parameter_list[[17]]
  DDRAIN_mean <- fitted_DDRAIN$estimate[1]
  DDRAIN_sd <- fitted_DDRAIN$estimate[2]
  print(paste("generation", iter, "DDRAIN", round(DDRAIN_mean,3), round(DDRAIN_sd,3)))
  DDRAIN <- rtruncnorm(new_nsims, 500, 1500, mean = DDRAIN_mean, sd = DDRAIN_sd)
  
  # fitted_GDRAIN
  # estimate Std. Error
  # mean 53.36850  0.3311553
  # sd   23.41622  0.2341622
  fitted_GDRAIN <- fitted_parameter_list[[18]]
  GDRAIN_mean <- fitted_GDRAIN$estimate[1]
  GDRAIN_sd <- fitted_GDRAIN$estimate[2]
  print(paste("generation", iter, "GDRAIN", round(GDRAIN_mean,3), round(GDRAIN_sd,3)))
  GDRAIN <- rtruncnorm(new_nsims, 0, 50, mean = GDRAIN_mean, sd = GDRAIN_sd)
  
  
  # fitted_THBACT
  # estimate Std. Error
  # mean 1.384695 0.02367388
  # sd   1.673996 0.01673994
  fitted_THBACT <- fitted_parameter_list[[19]]
  THBACT_mean <- fitted_THBACT$estimate[1]
  THBACT_sd <- fitted_THBACT$estimate[2]
  print(paste("generation", iter, "THBACT", round(THBACT_mean,3), round(THBACT_sd,3)))
  THBACT <- rtruncnorm(new_nsims, 0, 2, mean = THBACT_mean, sd = THBACT_sd)
  
  # fitted_WDPRCH
  # estimate  Std. Error
  # mean 0.5643253 0.003093901
  # sd   0.2187718 0.002187513
  fitted_WDPRCH <- fitted_parameter_list[[20]]
  WDPRCH_mean <- fitted_WDPRCH$estimate[1]
  WDPRCH_sd <- fitted_WDPRCH$estimate[2]
  print(paste("generation", iter, "WDPRCH", round(WDPRCH_mean,3), round(WDPRCH_sd,3)))
  WDPRCH <- rtruncnorm(new_nsims, 0, 1, mean = WDPRCH_mean, sd = WDPRCH_sd)
  
  inputs_df <- as.matrix(cbind(CN2, OV_N,GW_DELAY, GWQMN,SLSUBBSN, ALPHA_BNK, CH_K2, CH_N2, TRNSRCH, CH_K1, GW_REVAP, RCHRG_DP, SFTMP, SMTMP, TIMP, DEP_IMP, DDRAIN, GDRAIN, THBACT, WDPRCH))
  colnames(inputs_df)
  return(inputs_df)
}

save_generation_stats <- function(iter, data_in_dir, generation_stats){
  generation_stats_filename <- file.path(data_in_dir,paste("generation_stats",iter,".RData",sep=""))
  save(generation_stats, file= generation_stats_filename)
  print(paste("generation stats for generation ", iter, " saved to ", generation_stats_filename))
}

save_nses_parameters <- function(iter, data_in_dir, nses_parameters){
  nses_parameters_filename <- file.path(data_in_dir,paste("nses_parameters",iter,".RData",sep=""))
  save(nses_parameters, file= nses_parameters_filename)
  print(paste("nses and parameter inputs for generation ", iter, " saved to ", nses_parameters_filename))
}

simulate_generation_zero <- function(nsims, swat_path, base_dir, pars_initial){
  # run the initial set of swat simulations
  print(paste("About to run generation 0 with", nsims, "simulations"))
  # instead of this we could make a while loop and send batch jobs 
  # with 1000 inputs until we get required nuumber of winners
  bac_cal_output <- run_swat_red_cedar(swat_path, pars_initial)
  print(paste("swat runs finished for generation ", iter))
  return(bac_cal_output)
}

simulate_generation_next <- function(iter, nsims, swat_path, base_dir, pars_subsequent){
  # run the initial set of swat simulations
  print(paste("About to run generation", iter, "with", nsims, "simulations"))
  bac_cal_output <- run_swat_red_cedar(swat_path, pars_subsequent)
  print(paste("swat runs finished for generation ", iter))
  return(bac_cal_output)
}

save_bac_cal_output <- function(iter, bac_cal_output){
  rdata_out_filename <- paste('bac_cal', iter, '.RData', sep = "")
  rdata_file_out <- file.path(data_in_dir, rdata_out_filename)
  save(bac_cal_output, file=rdata_file_out)
  print(paste("rdata file for generation ", iter, " saved to ", rdata_file_out))
}

save_fitted_parameter_list <- function(iter, data_in, fitted_parameter_list){
  fitted_parameter_filename <- file.path(data_in_dir,paste("fitted_parameters", iter, ".RData", sep=""))
  save(fitted_parameter_list, file = fitted_parameter_filename)
  print(paste("fitted posterior parameters from generation", iter, "saved to file:", fitted_parameter_filename))
}

save_parameter_input_sims <- function(iter, data_in_dir, parameter_input_sims){
  parameter_input_sims_filename <- file.path(data_in_dir,paste("parameter_input_sims", iter, ".RData", sep=""))
  save(parameter_input_sims, file = parameter_input_sims_filename)
  print(paste("parameter input sims from generation", iter, "saved to file:", parameter_input_sims_filename))  
}

update_generation_stats <- function(iter, generation_stats, next_nsims, 
                                    max_mean_nse, max_bac_nse, max_q_nse, max_flux_nse,
                                    proportion_kept, new_cutoff){
  generation_stats$nsims[iter+2] <- next_nsims
  generation_stats$best_mean_nse[iter+1] <- max_mean_nse
  generation_stats$best_bac_nse[iter+1] <- max_bac_nse
  generation_stats$best_q_nse[iter+1] <- max_q_nse
  generation_stats$best_flux_nse[iter+1] <- max_flux_nse
  generation_stats$proportion_kept[iter+1] <- proportion_kept
  generation_stats$cutoff_score[iter+2] <- new_cutoff
  return(generation_stats)
}

update_cutoff_score <- function(iter, generation_stats, cutoff_score){
  generation_stats$cutoff_score[iter + 2] <- cutoff_score # +1 because zero-based and another +1 because stored in next generation slot
  print(paste("cutoff score created in this generation ", iter, "=", round(cutoff_score,6), "and will be used in the next generation =", (iter+1)))
  print(paste("cutoffs so far:", round(generation_stats$cutoff_score,6)))
}

