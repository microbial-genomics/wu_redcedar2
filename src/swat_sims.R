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


par_bound <- tibble(
              "CN2.mgt|change = relchg"= c(-0.3,0.3),
               "SOL_K(1).sol|change = relchg" = c(-0.8,0.8),
               "SOL_AWC(1).sol|change = relchg" = c(-0.8,2),
               "OV_N.hru|change = relchg" = c(-0.8,2),
               "ALPHA_BF.gw|change = relchg" = c(-0.3,0.3),
               "GW_DELAY.gw|change = relchg" = c(-0.75,4),
              "GWQMN.gw|change = relchg" = c(-0.5, 2),
              #"HRU_SLP.hru|change = relchg" = c(-0.25,2), # delete it further
              "SLSUBBSN.hru|change = relchg" = c(-0.5, 1),
               "ALPHA_BNK.rte|change = absval" =c(0, 1),
                "CH_K2.rte|change = absval" = c(0, 500),
                "CH_N2.rte|change = absval" = c(0.05, 0.15),
                "ESCO.bsn |change = absval" = c(0, 1),
                "EPCO.bsn|change = absval" = c(0, 1),
               "TRNSRCH.bsn|change = absval" = c(0, 0.3),
                "SURLAG.bsn|change = absval" = c(1, 24),
                "CH_N1.sub|change = absval" = c(0.05, 0.15),
                "CH_K1.sub|change = absval" = c(0, 300),
              "REVAPMN.gw |change = absval" = c(0, 1000), 
             "GW_REVAP.gw|change = absval" = c(0.02, 0.2),
               "RCHRG_DP.gw|change = absval" = c(0, 1),
              "GW_SPYLD.gw|change = absval" = c(0, 0.4),
             # "SNOEB.sub|change = absval"= c(0, 999999),
             "SFTMP.bsn|change = absval"= c(-20, 10),
              "SMTMP.bsn|change = absval"= c(-5,20),
             "SMFMX.bsn|change = absval"= c(0, 20),
             "SMFMN.bsn|change = absval"= c(0, 20),
              "TIMP.bsn|change = absval"= c(0, 1),

                "BACTKDQ.bsn|change = absval" = c(0, 500),
                "BACTMX.bsn|change = absval" = c(7, 20),
                "BACT_SWF.bsn|change = absval" = c(0, 1),

                "WDPRCH.bsn|change = absval"= c(0, 1),
                "WDPQ.bsn|change = absval"= c(0, 1),
                #"WGPQ.bsn|change = absval"= c(0, 1),
               "WDPS.bsn|change = absval"= c(0, 1),
               #"WGPS.bsn|change = absval"= c(0, 1),
               #"WOF_P.bsn|change = absval"= c(0, 1),
               "WDPRES.bsn|change = absval"= c(0, 1)
               )

n_sample <- 10000
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

nse_cal1 <- right_join(bac_cal1$simulation$bac_out,bac_obs,by="date") %>%
  select(-date) %>% select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs$bacteria))

  sort(nse_cal1, decreasing = T) %>% enframe()

col_idx <- sapply(bac_cal1$simulation$bac_out, function(x)all(!is.infinite(x) & !is.na(x)))
#coi_idx<-bac_cal1$simulation$bac_out[Reduce(`&`, lapply(bac_cal1$simulation$bac_out, is.finite)),]
valid_bac_simulation <- bac_cal1$simulation$bac_out[ ,col_idx]%>%colnames()
#nse_valid_simulation <- nse_cal1[valid_bac_simulation]
#sort_nse_bac<- sort(nse_valid_simulation, decreasing = T) %>% enframe()

idx=which(colnames(bac_cal1$simulation$bac_out)%in%valid_bac_simulation)


bac_out_valid=bac_cal1$simulation$bac_out[,idx]

q_out_valid=bac_cal1$simulation$q_out[,idx]
parameter_valid=bac_cal1$parameter$values[idx2,]


  save(bac_cal1, file='/work/OVERFLOW/RCR/MSU/bac_cal4.RData')

