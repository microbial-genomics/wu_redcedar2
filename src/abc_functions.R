

create_tibble_initial <- function(nsims){
  pars_tibble <- tibble(#hydrology parameters (11)
                        "CN2.mgt|change = relchg"= runif(nsims,-0.25,0.1),
                        "GWQMN.gw|change = relchg" = runif(nsims,-0.5, 2),
                        "ALPHA_BNK.rte|change = absval" =c(0, 1),
                        "CH_K2.rte|change = absval" = runif(nsims, 0, 500),
                        "CH_N2.rte|change = absval" = runif(nsims, 0, 0.3),
                        "TRNSRCH.bsn|change = absval" = runif(nsims, 0, 1),
                        "CH_N1.sub|change = absval" = runif(nsims, 0.01, 30),
                        "CH_K1.sub|change = absval" = runif(nsims, 0, 300),
                        "RCHRG_DP.gw|change = absval" = runif(nsims, 0, 1),
                        "SFTMP.bsn|change = absval"= c(-5, 5),
                        "SMTMP.bsn|change = absval"= c(-5,5),
                        #tile drainage and sediments (3)
                        "DEP_IMP.hru|change = absval"= c(0,6000),
                        "DDRAIN.mgt|change = absval"= c(0, 2000),
                        "GDRAIN.mgt|change = absval"= c(0, 100),
                        #bacteria submodel (4)
                        "BACTKDQ.bsn|change = absval" = c(0, 500),
                        "BACT_SWF.bsn|change = absval" = c(0, 1),
                        "THBACT.bsn|change = absval"= c(0, 10),
                        "WDPRCH.bsn|change = absval"= c(0, 1)
  )
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

simulate_generation_zero <- function(nsims, pars_initial){
  # run the initial set of swat simulations
  print(paste("About to run generation 0 with", nsims, "simulations"))
  swat_output0 <- run_swat_red_cedar(swat_path, pars_initial)
  
  #save the simulations
  save_file <- paste(base_dir,"rcr_swat_output0.RData",sep="")
  save(swat_output0, file=save_file)
  return(swat_output0)
}