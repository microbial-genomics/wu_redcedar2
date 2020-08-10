########### Preparation

#module load intel/19.0.5
#module load R/3.6.2
#module load geos/3.8.0
#module load gdal-2.4.3/intel-19.0
#module load proj-5.2.0/intel-19.0
#module load udunits-2.2.26/intel-19.0
#Rscript visualizaation.R

################ Load R library
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

load(file='/work/OVERFLOW/RCR/test2/q_obs.RData')

######## Parameter sampling



par_bound <- tibble("CN2.mgt|change = relchg"= c(-0.25,0.1),
                "SOL_K(1).sol|change = relchg" = c(-0.8,0.8),
                "SOL_AWC(1).sol|change = relchg" = c(-0.8,2),
                "OV_N.hru|change = relchg" = c(-0.8,2),
                "ALPHA_BF.gw|change = relchg" = c(-0.3,0.3),
                "GW_DELAY.gw|change = relchg" = c(-15,10),
                "GWQMN.gw|change = relchg" = c(-0.5, 2),
                "HRU_SLP.hru|change = relchg" = c(-15,10),
                "SLSUBBSN.hru|change = relchg" = c(-0.5, 3),
                "ALPHA_BNK.rte|change = absval" =c(0, 1),
                "CH_K2.rte|change = absval" = c(0, 500),
                "CH_N2.rte|change = absval" = c(0, 0.3),
                "ESCO.bsn |change = absval" = c(0, 1),
                "EPCO.bsn|change = absval" = c(0, 1),
                "TRNSRCH.bsn|change = absval" = c(0, 1),
                "SURLAG.bsn|change = absval" = c(1, 24),
                "CH_N1.sub|change = absval" = c(0.01, 30),
                "CH_K1.sub|change = absval" = c(0, 300),
                "REVAPMN.gw |change = absval" = c(0, 1000),
                "GW_REVAP.gw|change = absval" = c(0.02, 0.2),
                "RCHRG_DP.gw|change = absval" = c(0, 1),
                "GW_SPYLD.gw|change = absval" = c(0, 0.4))

par_bound

########## Random sampling with runif()
n_sample <- 50000
par_runif <- map_df(par_bound, ~ runif(n_sample, .x[1], .x[2]))

par_runif


######### Random sampling with lhs
n_sample <- 50000
n_par <- ncol(par_bound)

par_iter1 <- randomLHS(n = n_sample, k = n_par) %>%
  as_tibble(.) %>%
  map2_df(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1])) %>%
  set_names(names(par_bound))

par_iter1

################Model calibration
path_2012 <- "/work/OVERFLOW/RCR/test2"


q_iter3 <- run_swat2012(project_path = path_2012,
                        output = list(q_out = define_output(file = "rch",
                                      variable = "FLOW_OUT",
                                      unit = 4)),
                        parameter = par_iter1,
                        start_date = "2011-01-01",
                         end_date = "2013-12-31",
                        years_skip = 2,
                        n_thread = 8)

save(q_iter3, file='/work/OVERFLOW/RCR/test2/q_iter3.RData')

########## Model evaluation
nse_iter2 <- q_iter3$simulation$q_out %>% select(-date) %>% map_dbl(., ~NSE(.x/8.64, q_obs$discharge))
sort(nse_iter2, decreasing = T) %>% enframe()

save(nse_iter2, file='/work/OVERFLOW/RCR/test2/nse_iter2.RData')








