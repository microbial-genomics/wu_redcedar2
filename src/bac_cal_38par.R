library(SWATplusR)
library(dplyr)
library(dygraphs)
library(ggplot2)
library(lubridate)
library(mapview)
library(plotly)
library(sf)
library(tibble)
library(tidyr)
library(purrr)
library(lhs)
library(hydroGOF)
library(fast)
library(forcats)
library(lubridate)
library(sensitivity)





load(file='/work/OVERFLOW/RCR/MSU/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')


path <- "/work/OVERFLOW/RCR/calibration/MSU"


par_bound <- tibble("CN2.mgt|change = relchg"= c(-0.3,0.3),
               "SOL_K(1).sol|change = relchg" = c(-0.8,0.8),
               "SOL_AWC(1).sol|change = relchg" = c(-0.8,2),
               "OV_N.hru|change = relchg" = c(-0.8,2),
               "ALPHA_BF.gw|change = relchg" = c(-0.3,0.3),
               "GW_DELAY.gw|change = relchg" = c(-0.75,4),
              "GWQMN.gw|change = relchg" = c(-0.5,2),
               "HRU_SLP.hru|change = absval" = c(0,1),
              "SLSUBBSN.hru|change = relchg" = c(-0.5, 1),
               "ALPHA_BNK.rte|change = absval" =c(0, 1),
                "CH_K2.rte|change = absval" = c(0,50),
                "CH_N2.rte|change = absval" = c(0.05, 0.15),
                "ESCO.bsn |change = absval" = c(0, 1),
                "EPCO.bsn|change = absval" = c(0, 1),
                "TRNSRCH.bsn|change = absval" = c(0,0.3),
                "SURLAG.bsn|change = absval" = c(1, 24),
                "CH_N1.sub|change = absval" = c(0.05, 0.15),
                "CH_K1.sub|change = absval" = c(0, 300),
               "REVAPMN.gw |change = absval" = c(0, 1000),
               "GW_REVAP.gw|change = absval" = c(0.02, 0.2),
               "RCHRG_DP.gw|change = absval" = c(0, 1),
              "GW_SPYLD.gw|change = absval" = c(0, 0.4),
             "SFTMP.bsn|change = absval"= c(-5, 5),
              "SMTMP.bsn|change = absval"= c(-5,5),
             "SMFMX.bsn|change = absval"= c(0, 20),
             "SMFMN.bsn|change = absval"= c(0, 20),
              "TIMP.bsn|change = absval"= c(0.01, 1),

                "BACTKDQ.bsn|change = absval" = c(0, 500),
                "BACTMX.bsn|change = absval" = c(7, 20),
               # "BACTKDDBfert.dat|change = absval"= c(0, 1),
               #"BACTKDDB(44).fert|change = absval"= c(0, 1),
            #   "BACTKDDB(44)fert.dat|change = absval"= c(0, 1),
               # "BACTKDDBfert.dat|change = absval"= c(0, 1),
                "BACT_SWF.bsn|change = absval" = c(0, 1),
#                "BCNST.bsn|change = absval" = c(0, 100000),
                "CFRT_KG.mgt|change = relchg" = c(0, 500),
               "FRT_SURFACE.mgt|change = absval"= c(0, 1),
               "THBACT.bsn|change = absval"= c(0, 10),
                "WDPRCH.bsn|change = absval"= c(0, 1),
                "WDPQ.bsn|change = absval"= c(0, 1),
                #"WGPQ.bsn|change = absval"= c(0, 1),
               "WDPS.bsn|change = absval"= c(0, 1),
               #"WGPS.bsn|change = absval"= c(0, 1),
               "WOF_P.bsn|change = absval"= c(0, 1),
               "WDPRES.bsn|change = absval"= c(0, 1))

n_sample <- 50000
par_runif <- map_df(par_bound, ~ runif(n_sample, .x[1], .x[2]))
par_runif

bac_cal1 <- run_swat2012(project_path = path,
                       output = list(q_out = define_output(file = "rch",
                                     variable = "FLOW_OUT",
                                     unit = 4),
                       		     bac_out = define_output(file = "rch",
                                     variable = "BACTP_OUT",
                                     unit = 4)),
                       parameter = par_runif,
                       start_date = "2011-01-01",
                        end_date = "2013-12-31",
                       years_skip = 2,
                       n_thread = 32)


  save(bac_cal1, file='/work/OVERFLOW/RCR/calibration/MSU/bac_cal6.RData')

