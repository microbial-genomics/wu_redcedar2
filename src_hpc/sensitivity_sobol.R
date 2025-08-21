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
library(lhs)
library(tibble)




load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')


path <- "/work/OVERFLOW/RCR/MSU"

swat_sobol <- function(par, obs) {
  q_bac_sim <- run_swat2012(project_path = path,
                        output =list(q_out = define_output(file = "rch",
                                     variable = "FLOW_OUT",
                                     unit = 4),
                                     bac_out = define_output(file = "rch",
                                     variable = "BACTP_OUT",
                                     unit = 4)),
                        parameter = par,
                        start_date = "2011-01-01",
                        end_date = "2013-12-31",
                        years_skip = 2,
                        add_date = TRUE,
                        n_thread = 32)
  save(q_bac_sim, file='/work/OVERFLOW/RCR/MSU/q_bac_sim.RData')
  #all_days <- seq(as.Date("2013-1-1"), as.Date("2013-12-31"),by="day")
  # temp <- bac_out%>%add_column(date=all_days,.before="run_00001")
  ###temp2<-right_join(temp,bac_obs,by="date")
  ###check two table date column name 
  #colnames(temp)

  temp3<-right_join(q_bac_sim$simulation$bac_out,bac_obs,by="date")%>%select(-date) %>%select(-bacteria)
  #colnames(temp)
  save(temp3, file='/work/OVERFLOW/RCR/MSU/temp3.RData')
  nse_bac <-map_dbl(temp3, ~NSE(.x,obs)) 
  #nse_bac <-map_dbl(temp, ~NSE(.x,bac_obs$bacteria)) 
  #sort(nse_bac, decreasing = T) %>% enframe()

  #

  return(nse_bac)
}

#################################################################
par_bound <- tibble("CN2.mgt|change = relchg"= c(-0.3,0.3),
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
              "GW_SPYLD.gw|change = absval" = c(0, 0.4),
             # "SNOEB.sub|change = absval"= c(0, 999999),
             "SFTMP.bsn|change = absval"= c(-20, 10),
              "SMTMP.bsn|change = absval"= c(-5,20),
             "SMFMX.bsn|change = absval"= c(0, 20),
             "SMFMN.bsn|change = absval"= c(0, 20),
              "TIMP.bsn|change = absval"= c(0, 1),

                "BACTKDQ.bsn|change = absval" = c(0, 500),
                "BACTMX.bsn|change = absval" = c(7, 20),
               # "BACTKDDBfert.dat|change = absval"= c(0, 1),
               #"BACTKDDB(44).fert|change = absval"= c(0, 1),
            #   "BACTKDDB(44)fert.dat|change = absval"= c(0, 1),
               # "BACTKDDBfert.dat|change = absval"= c(0, 1),
                "BACT_SWF.bsn|change = absval" = c(0, 1),
               # "BCNST.bsn|change = absval" = c(0, 100000),
                "CFRT_KG.mgt|change = absval" = c(0, 500),
                "FRT_SURFACE.mgt|change = absval"= c(0, 1),
                "THBACT.bsn|change = absval"= c(0, 10),
                "WDPRCH.bsn|change = absval"= c(0, 1),
                "WDPQ.bsn|change = absval"= c(0, 1),
                "WGPQ.bsn|change = absval"= c(0, 1),
               "WDPS.bsn|change = absval"= c(0, 1),
               "WGPS.bsn|change = absval"= c(0, 1),
               "WOF_P.bsn|change = absval"= c(0, 1),
               "WDPRES.bsn|change = absval"= c(0, 1))

n_par  <- 40

n_samp <- 10

x1 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))
x2 <- data.frame(matrix(runif(n_par * n_samp), nrow = n_samp)) %>%
  set_names(., names(par_bound)) %>%
  map2_dfc(., par_bound, ~ (.x * (.y[2] - .y[1]) + .y[1]))


sens_sobol_bac <- sobol(model = swat_sobol, X1 = x1, X2 = x2, 
                    obs = bac_obs$bacteria, nboot = 100)

save(sens_sobol_bac, file='/work/OVERFLOW/RCR/MSU/sens_sobol_bac.RData')

plot_sobol <- sens_sobol_bac$S %>%
  mutate(parameter = rownames(.)) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., original))


ggplot(data = plot_sobol) +
  geom_pointrange(aes(x = parameter, y = original ,
                      ymin = `min. c.i.`, ymax = `max. c.i.`)) +
  coord_flip() +
  xlab("Parameter") +
  ylab("Sensitivity") +
  theme_bw()

ggsave('bac_sobol2.pdf')

