
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
library(forcats)
library(lubridate)
library(sensitivity)





load(file='/work/OVERFLOW/RCR/sim56-sensitivity/q_obs.RData')
load(file='/work/OVERFLOW/RCR/sim56-sensitivity/bac_obs.RData')


path <- "/work/OVERFLOW/RCR/sim56-sensitivity"


par_bound <- tibble(
  #Hydrology
              "CN2.mgt|change = relchg"= c(-0.3,0.3),
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

#tile drainage and sediments
              "DEP_IMP.hru|change = absval"= c(0,6000),
              "DDRAIN.mgt|change = absval"= c(0, 2000),
              "TDRAIN.mgt|change = absval"= c(0, 72),
              "GDRAIN.mgt|change = absval"= c(0, 100),
              "SPCON.bsn|change = absval"= c(0.0001, 0.01),
              "SPEXP.bsn|change = absval"= c(1, 2),
              "PRF_BSN.bsn|change = absval"= c(0.5, 2),
              "ADJ_PKR.bsn|change = absval"= c(0.5, 2),
              


                "BACTKDQ.bsn|change = absval" = c(0, 500),
                "BACTMX.bsn|change = absval" = c(7, 20),
                "BACT_SWF.bsn|change = absval" = c(0, 1),
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

#45 parameters


n_sample <- 5000
# par_runif <- map_df(par_bound, ~ runif(n_sample, .x[1], .x[2]))
# par_runif
# 
# bac_cal1 <- run_swat2012(project_path = path,
#                        output = list(q_out = define_output(file = "rch",
#                                      variable = "FLOW_OUT",
#                                      unit = 4),
#                                      bac_out = define_output(file = "rch",
#                                      variable = "BACTP_OUT",
#                                      unit = 4)),
#                        parameter = par_runif,
#                        start_date = "2000-01-01",
#                         end_date = "2014-07-31",
#                        years_skip = 4,
#                        n_thread = 32)
# 
# save(bac_cal1, file = "/work/OVERFLOW/RCR/sim56-sensitivity/bac_sensitivity.RData")



library(dplyr)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(miceadds)
library(sensitivity)
library(vioplot)
library(xts)


# huiyun laptop
  rcdir <- path.expand("/work/OVERFLOW/RCR/sim56-sensitivity/")


# 
# # Variable names for data
# rcdir_data_in <- paste(rcdir,'data_in/',sep='')
# rcdir_data_out <- paste(rcdir,'data_out/',sep='')
# rcdir_graphics <- paste(rcdir,'graphics/',sep='')
# rcdir_src <- paste(rcdir,'src/',sep='')

# load data_out files

load(file = '/work/OVERFLOW/RCR/sim56-sensitivity/bac_sensitivity.RData')

sim_parameters <- bac_cal1$parameter$value


sim_bac_concs <- bac_cal1$simulation$bac_out

sim_flows <- bac_cal1$simulation$q_out 



# file info
#View(sim_parameters)
#View(sim_bac_concs)
#View(sim_flows)
dim(sim_parameters)
colnames(sim_parameters)
dim(sim_bac_concs)
colnames(sim_bac_concs)
dim(sim_flows)
colnames(sim_flows)

#fix problematic colnames
colnames(sim_parameters)[2] <- "SOL_K"
colnames(sim_parameters)[3] <- "SOL_AWC"

#simulation is for the 365 days of 2013
#drop date column for sim_bac_concs
sim_bac_concs2 <- sim_bac_concs[,-c(1)]
sim_dates <- sim_bac_concs[,1]
#View(sim_bac_concs2)
dim(sim_bac_concs2)
#transpose bacteria concentration data frame
sim_bac_concs3 <- t(sim_bac_concs2)
dim(sim_bac_concs3)
# pcc on row medians
sim_bac_medians <- rowMedians(sim_bac_concs3)
pcc(sim_parameters, sim_bac_medians)
# pcc on row maxs
sim_bac_maxs <- rowMaxs(sim_bac_concs3)
pcc(sim_parameters, sim_bac_maxs)


# drop dates for sim_flows
sim_flows2 <- sim_flows[,-c(1)]
dim(sim_flows2)
#transpose flows data frame
sim_flows3 <- t(sim_flows2)
dim(sim_flows3)
# pcc on row medians
sim_flows_medians <- rowMedians(sim_flows3)
pcc(sim_parameters, sim_flows_medians)
# pcc on row maxs
sim_flows_maxs <- rowMaxs(sim_flows3)
pcc(sim_parameters, sim_flows_maxs)

# flux
sim_flux3 <- sim_bac_concs3*sim_flows3*10^4
dim(sim_flux3)
#pcc on row medians
sim_flux_medians <- rowMedians(sim_flux3)
pcc(sim_parameters, sim_flux_medians)
#pcc on row maxs
sim_flux_maxs <- rowMaxs(sim_flux3)
pcc(sim_parameters, sim_flux_maxs)


# daily pcc for bacteria and flow
bac_pcc <- matrix(data=NA, nrow=3865, ncol=45)
flows_pcc <- matrix(data=NA, nrow=3865, ncol=45)
flux_pcc <- matrix(data=NA, nrow=3865, ncol=45)
for(i in 1:3865){
  print(i)
  daily_bac_pcc <- pcc(sim_parameters, sim_bac_concs3[,i])
  daily_flows_pcc <- pcc(sim_parameters, sim_flows3[,i])
  daily_flux_pcc <- pcc(sim_parameters, sim_flux3[,i])
  length(bac_pcc[i,])
  length(t(daily_bac_pcc$PCC))
  bac_pcc[i,] <- t(daily_bac_pcc$PCC)
  flows_pcc[i,] <- t(daily_flows_pcc$PCC)
  flux_pcc[i,] <- t(daily_flux_pcc$PCC)
}


save(bac_pcc, file = "/work/OVERFLOW/RCR/sim56-sensitivity/bac_pcc.RData")

save(flows_pcc, file = "/work/OVERFLOW/RCR/sim56-sensitivity/flows_pcc.RData")

save(flux_pcc, file = "/work/OVERFLOW/RCR/sim56-sensitivity/flux_pcc.RData")

#print violin plot
mklab <- function(y_var){
  if(y_var){
    names(mf)[response]
  } else {
    paste(names(mf)[-response], collapse = " : ")
  }
}
dim(bac_pcc)
colnames(bac_pcc) <- colnames(sim_parameters)
colnames(flows_pcc) <- colnames(sim_parameters)
colnames(flux_pcc) <- colnames(sim_parameters)
pdf(paste("pcc_violin_50000_45.pdf",sep=""),width=55,height=30,onefile=TRUE)
  vioplot(bac_pcc)
  vioplot(flows_pcc)
  vioplot(flux_pcc)
dev.off()


#create time series
bac_dates <- ts(sim_dates, start=c(2004, 1), end=c(2014, 3865), frequency=3865)


#simple ggplot
dim(bac_pcc)
pdf(paste(rcdir_graphics,"pcc_ts_50000_45.pdf",sep=""),width=11,height=8, onefile=TRUE)
  for(i in 1:45){
    data <- data.frame(
      day = as.Date("2014-01-01") - 0:3864,
      bac_value = bac_pcc[,i],
      flows_value = flows_pcc[,i],
      flux_value = flux_pcc[,i]
    )
    p <- ggplot(data, aes(x=day)) +
      geom_line(aes(y=bac_value), color = "darkred") + 
      geom_line(aes(y=flows_value), color = "steelblue", linetype="twodash") +
      geom_line(aes(y=flux_value), color = "orange") +
      xlab("") +
      ylab("pcc sensitivity") +
      ylim(-1,1)
    print(p + ggtitle(colnames(bac_pcc)[i]))
  }
dev.off()


