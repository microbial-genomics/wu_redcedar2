
load_observations <- function(){
  load(file= file.path(data_in_dir,'bac_obs.RData'), .GlobalEnv)
  load(file = file.path(data_in_dir,'flux_obs.RData'), .GlobalEnv)
  load(file = file.path(data_in_dir,'pcp_obs.RData'), .GlobalEnv)
  load(file = file.path(data_in_dir,'pcp_obs2.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'q_obs.RData'), .GlobalEnv)
  load(file= file.path(data_in_dir, 'q_obs2.RData'), .GlobalEnv)
}

get_cutoff_median_score <- function(iter){
  # initialize if starting from generation 0
  # if we are starting at zero we need to wipe out the previous_median_score_vector
  # otherwise we will save it after every loop
  cutoff_median_filename <- file.path(data_in_dir,"cutoff_median_scores.RData")
  if(iter==0){
    cutoff_median_score <- rep(NA, 40) # if we are starting over at gen 0 we need to delete past scores
    save(cutoff_median_score, file = cutoff_median_filename)
    return_cutoff_median <- NA # we do not need to save this because median score is not used for gen 0
  }else{
    # if starting at gen greater than zero we will 
    # calculate it based on the output of the previous generation,
    # but we load it to make sure the vector exists when we save it later
    load(file= cutoff_median_filename, .GlobalEnv)
    return_cutoff_median <- cutoff_median_score[iter] #get the cutoff that was saved at the last generation
  }
  return(return_cutoff_median)
}

create_tibble_initial <- function(nsims){
  pars_tibble <- tibble(#hydrology parameters (11)
                        "CN2.mgt|change = relchg"= runif(nsims, -0.25, 0.1),
                        "GWQMN.gw|change = relchg" = runif(nsims, -0.5, 2),
                        "ALPHA_BNK.rte|change = absval" =runif(nsims, 0, 1),
                        "CH_K2.rte|change = absval" = runif(nsims, 0, 500),
                        "CH_N2.rte|change = absval" = runif(nsims, 0, 0.3),
                        "TRNSRCH.bsn|change = absval" = runif(nsims, 0, 1),
                        "CH_N1.sub|change = absval" = runif(nsims, 0.01, 30),
                        "CH_K1.sub|change = absval" = runif(nsims, 0, 300),
                        "RCHRG_DP.gw|change = absval" = runif(nsims, 0, 1),
                        "SFTMP.bsn|change = absval"= runif(nsims, -5, 5),
                        "SMTMP.bsn|change = absval"= runif(nsims, -5, 5),
                        #tile drainage and sediments (3)
                        "DEP_IMP.hru|change = absval"= runif(nsims, 0, 6000),
                        "DDRAIN.mgt|change = absval"= runif(nsims, 0, 2000),
                        "GDRAIN.mgt|change = absval"= runif(nsims, 0, 100),
                        #bacteria submodel (4)
                        "BACTKDQ.bsn|change = absval" = runif(nsims, 0, 500),
                        "BACT_SWF.bsn|change = absval" = runif(nsims, 0, 1),
                        "THBACT.bsn|change = absval"= runif(nsims, 0, 10),
                        "WDPRCH.bsn|change = absval"= runif(nsims, 0, 1)
  )
}

#create_tibble_subsequent(){
#   tibble(
#    "CN2.mgt|change = relchg"= CN2,
#    "GWQMN.gw|change = relchg" = GWQMN,
#    "ALPHA_BNK.rte|change = absval" = ALPHA_BNK,
#    "CH_K2.rte|change = absval" = CH_K2,
#    "CH_N2.rte|change = absval" = CH_N2,
#    "TRNSRCH.bsn|change = absval" = TRNSRCH,
#    "CH_N1.sub|change = absval" = CH_N1,
#    "CH_K1.sub|change = absval" = CH_K1,
#    "RCHRG_DP.gw|change = absval" = RCHRG_DP,
#    "SFTMP.bsn|change = absval"= SFTMP,
#    "SMTMP.bsn|change = absval"= SMTMP,
#    "DEP_IMP.hru|change = absval"= DEP_IMP,
#    "DDRAIN.mgt|change = absval"= DDRAIN,
#    "GDRAIN.mgt|change = absval"= GDRAIN,
#    "BACTKDQ.bsn|change = absval" = BACTKDQ,
#    "BACT_SWF.bsn|change = absval" = BACT_SWF,
#    "THBACT.bsn|change = absval"= THBACT,
#    "WDPRCH.bsn|change = absval"= WDPRCH)
#}

simulate_generation_zero <- function(nsims, swat_path, base_dir, pars_initial){
  # run the initial set of swat simulations
  print(paste("About to run generation 0 with", nsims, "simulations"))
  bac_cal_output <- run_swat_red_cedar(swat_path, pars_initial)
  print(paste("swat runs finished for generation ", iter))
  return(bac_cal_output)
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
               start_date = "2011-01-01",
               end_date = "2013-12-31",
               years_skip = 2,
               n_thread = 32)
}

save_bac_cal_output <- function(iter, bac_cal_output){
  rdata_out_filename <- paste('bac_cal', iter, '.RData', sep = "")
  rdata_file_out <- file.path(data_in_dir, rdata_out_filename)
  save(bac_cal1, file=rdata_file_out)
  print(paste("rdata file for generation ", iter, " saved to ", rdata_file_out))
}

post_process_swat <- function(iter){
  previous_median_score[iter] = new_median_score
  write.csv(previous_median_score, file = median_filename)
  print(paste("the median score cutoff used to create particles for generation ", iter, " was ", previous_median_score))
}

calculate_mean_nses <- function(){
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
}

save_kde_pdf <- function(){
  ggplot(data = kde_next_gen) +
    geom_density(aes(x = parameter_range)) +
    facet_wrap(.~par, nrow=5, scales = "free") +
    theme_bw()
  density_plot_filename <- paste("kde_mcabc_gen", iter, ".pdf", sep="")
  ggsave(file.path(graphics_dir, density_plot_filename))  
}



log_results <- function(){
  print(paste("Generation ", iter-1))
  print(paste("median score for the last generation was:", previous_median_score))
  print(paste("generation x:",n_all_keepers, "of", format(previous_nsims,scientific=F), " simulations kept; proportion kept =", round(proportion_kept,4)))
  print(paste("best kept mean nse for this generation is:", max(round(nse_mean_keepers,4))))
  print(paste("best bacteria nse for this generation is:", max(round(nse_bac,4))))
  print(paste("best flow nse for this generation is:", max(round(nse_q,4))))
  print(paste("best flux nse for this generation is:", max(round(nse_flux,4))))
  print(paste("median mean nse score for this generation is:", round(new_median_score,4)))
}

load_previous_swat_simulations <- function(){
  bac_cal_filename <- paste('bac_cal', iter-1, '.RData', sep="")
  rdata_file_in <- file.path(data_in_dir, bac_cal_filename)
  print(paste("loading data file: ", rdata_file_in))
  load(file = rdata_file_in)
  
  #determine number of simulations last time
  dim(bac_cal1$simulation$bac_out)
}

nse_simulations_v_observations <- function(bac_cal_output, bac_obs, q_obs, flux_obs){
  #load the simulated concentrations, flows, inputs for last simulations
  sim_bac <- bac_cal_output$simulation$bac_out
  sim_q <- bac_cal_output$simulation$q_out
  sim_pars <- bac_cal_output$parameter$values
  
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
}

fit_normal_parameters <- function(sim_pars_keepers){
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
}

sample_truncated_normals <- function(){
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
  
}
  
