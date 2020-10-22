library(SWATplusR)
library(tibble)


if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  base_dir <- "c:/git/wu_redcedar2/"
}

nsims=100

par_set22_50000 <- tibble("CN2.mgt|change = relchg"= runif(nsims,-0.25,0.1),
		"SOL_K(1).sol|change = relchg" = runif(nsims,-0.8,0.8),
		"SOL_AWC(1).sol|change = relchg" = runif(nsims,-0.8,2),
		"OV_N.hru|change = relchg" = runif(nsims,-0.8,2),
		"ALPHA_BF.gw|change = relchg" = runif(nsims,-0.3,0.3),
		"GW_DELAY.gw|change = relchg" = runif(nsims,-15,10),
		"GWQMN.gw|change = relchg" = runif(nsims,-0.5, 2),
		"HRU_SLP.hru|change = relchg" = runif(nsims,-15,10),
		"SLSUBBSN.hru|change = relchg" = runif(nsims,-0.5, 3),	
                "ALPHA_BNK.rte|change = absval" = runif(nsims, 0, 1),
		"CH_K2.rte|change = absval" = runif(nsims, 0, 500),
		"CH_N2.rte|change = absval" = runif(nsims, 0, 0.3),
		"ESCO.bsn |change = absval" = runif(nsims, 0, 1),
		"EPCO.bsn|change = absval" = runif(nsims, 0, 1),
		"TRNSRCH.bsn|change = absval" = runif(nsims, 0, 1),
		"SURLAG.bsn|change = absval" = runif(nsims, 1, 24),
		"CH_N1.sub|change = absval" = runif(nsims, 0.01, 30),
		"CH_K1.sub|change = absval" = runif(nsims, 0, 300),
		"REVAPMN.gw |change = absval" = runif(nsims, 0, 1000),
		"GW_REVAP.gw|change = absval" = runif(nsims, 0.02, 0.2),
		"RCHRG_DP.gw|change = absval" = runif(nsims, 0, 1),
		"GW_SPYLD.gw|change = absval" = runif(nsims, 0, 0.4))

path <- "/work/OVERFLOW/stp/MSU"

swat_output0 <- run_swat2012(project_path = path,
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

save_file <- paste(path,"/swat_output0.RData",sep="")
save(swat_output0, file=save_file)


