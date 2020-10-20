library(SWATplusR)
library(tibble)


if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  base_dir <- "c:/git/wu_redcedar2/"
}

par_set22_50000 <- tibble("CN2.mgt|change = relchg"= runif(50000,-0.25,0.1),
		"SOL_K(1).sol|change = relchg" = runif(50000,-0.8,0.8),
		"SOL_AWC(1).sol|change = relchg" = runif(50000,-0.8,2),
		"OV_N.hru|change = relchg" = runif(50000,-0.8,2),
		"ALPHA_BF.gw|change = relchg" = runif(50000,-0.3,0.3),
		"GW_DELAY.gw|change = relchg" = runif(50000,-15,10),
		"GWQMN.gw|change = relchg" = runif(50000,-0.5, 2),
		"HRU_SLP.hru|change = relchg" = runif(50000,-15,10),
		"SLSUBBSN.hru|change = relchg" = runif(50000,-0.5, 3),	
                "ALPHA_BNK.rte|change = absval" = runif(50000, 0, 1),
		"CH_K2.rte|change = absval" = runif(50000, 0, 500),
		"CH_N2.rte|change = absval" = runif(50000, 0, 0.3),
		"ESCO.bsn |change = absval" = runif(50000, 0, 1),
		"EPCO.bsn|change = absval" = runif(50000, 0, 1),
		"TRNSRCH.bsn|change = absval" = runif(50000, 0, 1),
		"SURLAG.bsn|change = absval" = runif(50000, 1, 24),
		"CH_N1.sub|change = absval" = runif(50000, 0.01, 30),
		"CH_K1.sub|change = absval" = runif(50000, 0, 300),
		"REVAPMN.gw |change = absval" = runif(50000, 0, 1000),
		"GW_REVAP.gw|change = absval" = runif(50000, 0.02, 0.2),
		"RCHRG_DP.gw|change = absval" = runif(50000, 0, 1),
		"GW_SPYLD.gw|change = absval" = runif(50000, 0, 0.4))

q_sim_50000 <- run_swat2012(project_path = "/work/OVERFLOW/stp/MSU",
                      output = define_output(file = "rch",
                                             variable = "FLOW_OUT",
                                             unit = c(4,22,27)),
                                             parameter = par_set22_50000,
                                             save_path = "result50000",
                                             return_output = FALSE,
		                                         n_thread =8)
