#Parameter sensitivity analysis

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

load(file='/work/OVERFLOW/RCR/sensitivity/msu_10000/q_obs.RData')
load(file='/work/OVERFLOW/RCR/sensitivity/msu_10000/q_iter2.RData')

path_2012 <- "/work/OVERFLOW/RCR/sensitivity/msu_10000"


q_obs 

par_names <-   c("CN2.mgt|change = relchg"
                "SOL_AWC(1).sol|change = relchg",
                "OV_N.hru|change = absval",
                "ALPHA_BF.gw|change = absval",
                "GW_DELAY.gw|change = absval",
                "GWQMN.gw|change = absval",
                "HRU_SLP.hru|change = relchg",
                "SLSUBBSN.hru|change = relchg",
                "ALPHA_BNK.rte|change = absval",
                "CH_K2.rte|change = absval",
                "CH_N2.rte|change = absval",
                "ESCO.bsn |change = absval",
                "EPCO.bsn|change = absval",
                "TRNSRCH.bsn|change = absval",
                "SURLAG.bsn|change = absval",
                "CH_N1.sub|change = absval",
                "CH_K1.sub|change = absval",
                "REVAPMN.gw |change = absval",
                "GW_REVAP.gw|change = absval",
                "RCHRG_DP.gw|change = absval",
                "GW_SPYLD.gw|change = absval")
                       

par_fast <- fast_parameters(
  minimum = c(-0.3, -0.8, -0.8, 0.004,  0,  0,    0,   -15, -0.5, 0,   0,  0,   0, 0, 0, 1,  0.01,  0,   0,    0.02, 0, 0),
  maximum = c(0.3,  0.8,    2,   0.8,   1,  500,  5000, 10,    3, 1, 500,  0.3, 1, 1, 1, 24, 30,    300, 1000, 0.2,  1, 0.4),
    names = par_names) %>%
  as_tibble()



par_fast

q_fast <- run_swat2012(project_path = path_2012,
            output = list(q_sim = define_output(file = "rch",
                                      variable = "FLOW_OUT",
                                      unit = 4)),
                        parameter = par_fast,
                        start_date = "2011-01-01",
                        end_date = "2013-12-31",
                        years_skip = 2,
                        n_thread = 8)

 save(q_fast, file='/work/OVERFLOW/RCR/sensitivity/result_msu_10000/q_fast.RData')


nse_fast <- q_fast$simulation$q_sim %>%
  select(-date) %>%
  map_dbl(., ~NSE(.x, q_obs$discharge))

nse_fast <- q_out_valid %>%
  select(-date) %>%
  map_dbl(., ~NSE(.x, q_obs$discharge))

nse_fast <- right_join(bac_out_valid,bac_obs,by="date") %>%
  select(-date) %>%
  map_dbl(., ~NSE(.x, bac_obs$bacteria))

sens_fast <- sensitivity(nse_fast, 40)

result_fast <- tibble(parameter = q_fast$parameter$definition$par_name,
                      fast      = sens_fast) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., fast))


result_fast <- tibble(parameter = bac_cal1$parameter$definition$par_name,
                      fast      = sens_fast) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., fast))


  ggplot(data = result_fast) +
  geom_col(aes(x = parameter, y = fast)) +
  xlab("Parameter") +
  ylab("Sensitivity") +
  coord_flip() +
  theme_bw()

  swat_sobol <- function(par, obs) {
  q_sim <- run_swatplus(project_path = proj_path,
                        output =list(q_sim = define_output(file = "channel",
                                     variable = "flo_out",
                                     unit = 1)),
                        parameter = par,
                        start_date = "2000-01-01",
                        end_date = "2007-12-31",
                        years_skip = 3, n_thread = 4,
                        add_date = FALSE)
  nse_q <- map_dbl(q_sim$simulation$q_sim, ~ NSE(.x, obs))
  return(nse_q)
}


par_bound <- tibble("cn2.hru | change = abschg" = c(-15, 10),
                    "lat_ttime.hru | change = absval" = c(0.5, 50),
                    "lat_len.hru | change = absval" = c(10, 100),
                    "k.sol | change = pctchg" = c(-50, 50),
                    "z.sol | change = pctchg" = c(-50, 50),
                    "esco.hru | change = absval" = c(0, 1),
                    "epco.hru | change = absval" = c(0, 1))
n_par  <- 7
n_samp <- 500
x1 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))
x2 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))


  sens_sobol <- sobol(model = swat_nse, X1 = x1, X2 = x2, 
                    obs = q_obs$discharge, nboot = 100)

  plot_sobol <- sens_sobol$S %>%
  mutate(parameter = rownames(.)) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., original))
ggplot(data = plot_sobol) +
  geom_pointrange(aes(x = parameter, y = original ,
                      ymin = `min. c.i.`, ymax = `max. c.i.`)) +
  coord_flip() +
  xlab("Parameter") +
  ylab("Sensitivity") +
  theme_bw()

rm(list = ls())

############################################################################

library(dplyr)
library(fast)
library(forcats)
library(ggplot2)
library(hydroGOF)
library(lubridate)
library(purrr)
library(SWATplusR)
library(sensitivity)
library(tidyr)



load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
path_2012 <- "/work/OVERFLOW/RCR/sensitivity/msu_sobol"


swat_sobol <- function(par, obs) {
  q_sim <- run_swat2012(project_path = path_2012,
                        output =list(q_sim = define_output(file = "rch",
                                     variable = "FLOW_OUT",
                                     unit = 4)),
                        parameter = par,
                        start_date = "2011-01-01",
                        end_date = "2013-12-31",
                        years_skip = 2,
                        n_thread = 8,
                        add_date = FALSE)
  nse_q <- map_dbl(q_sim$simulation$q_sim, ~ NSE(.x, obs))
  return(nse_q)
}


par_bound <- tibble("CN2.mgt|change = relchg"= c(-0.3,0.3),
                "SOL_K(1).sol|change = relchg" = c(-0.8,0.8),
                "SOL_AWC(1).sol|change = relchg" = c(-0.8,2),
                "OV_N.hru|change = absval" = c(0.004,0.8),
                "ALPHA_BF.gw|change = absval" = c(0,1),
                "GW_DELAY.gw|change = absval" = c(0,500),
                "GWQMN.gw|change = absval" = c(0, 5000),
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

n_par  <- 22

n_samp <- 500
x1 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))
x2 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))


sens_sobol <- sobol(model = swat_sobol, X1 = x1, X2 = x2, 
                    obs = q_obs$discharge, nboot = 100)

save(sens_sobol, file='/work/OVERFLOW/RCR/sensitivity/result_msu_10000/sens_sobol.RData')


