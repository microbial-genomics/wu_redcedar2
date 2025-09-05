library(sensitivity)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(miceadds)
#library(dtwclust)
library(vioplot)
library(xts)
#library(SWATplusR)
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
library(ggtext)
library(RColorBrewer)
library(cowplot)
#library(devtools)
#library(dplyr)
#library(dygraphs)
##library(fast) #error
#library(forcats)
#library(ggplot2)
#library(hydroGOF)
#library(lhs)
#library(lubridate)
#library(mapview)
#library(plotly)
#library(purrr)
#library(sensitivity)
#library(sf)
#library(tibble)
#library(tidyr)
#library(fitdistrplus)
#library(truncnorm)
#library(xts)
#require(data.table) # load it

## setup root directory path
if(Sys.info()[4]=="LZ26TPURUCKE-2"){
  # tom epa windows
  of_root <- file.path("c:", "Users", "tpurucke", "git", "wu_redcedar2")
}

data_in_dir <- file.path(of_root, "data_in")
graphics_dir <- file.path(of_root, "graphics")
src_dir <- file.path(of_root, "src")
hpc_data <- file.path(of_root, "hpc_data", "sim56-flux-weekly")
hpc_data_sensitivity <- file.path(of_root, "hpc_data", "sim56-sensitivity")

# simulation sets 0, 3, 6, 9, 11
# load parameter fits but in list format
# 18 parameter fits for each variable
# set 0
load(file.path(hpc_data, "fitted_parameters0.RData"))
fitted_parameter_list0 <- fitted_parameter_list
str(fitted_parameter_list)
for(i in 1:18){
  print(fitted_parameter_list0[[i]]$estimate)
}

# set 3
load(file.path(hpc_data, "fitted_parameters3.RData"))
fitted_parameter_list3 <- fitted_parameter_list
for(i in 1:18){
  print(fitted_parameter_list3[[i]]$estimate)
}

# set 6
load(file.path(hpc_data, "fitted_parameters6.RData"))
fitted_parameter_list6 <- fitted_parameter_list
for(i in 1:18){
  print(fitted_parameter_list6[[i]]$estimate)
}

# set 9
load(file.path(hpc_data, "fitted_parameters9.RData"))
fitted_parameter_list9 <- fitted_parameter_list
for(i in 1:18){
  print(fitted_parameter_list9[[i]]$estimate)
}

# set 11
load(file.path(hpc_data, "fitted_parameters11.RData"))
fitted_parameter_list11 <- fitted_parameter_list
for(i in 1:18){
  print(fitted_parameter_list11[[i]]$estimate)
}

# this loads everything for the generation stats
load(file.path(hpc_data, "generation_stats11.RData"))
generation_stats <- generation_stats[1:11,]
#View(generation_stats)

# nse scores for each generation
load(file.path(hpc_data, "nses_parameters0.RData"))
nses_parameters0 <- nses_parameters
dim(nses_parameters0)
parameter_names <- colnames(nses_parameters0)[5:22]
parameter_names

load(file.path(hpc_data, "nses_parameters3.RData"))
nses_parameters3 <- nses_parameters
dim(nses_parameters3)

load(file.path(hpc_data, "nses_parameters6.RData"))
nses_parameters6 <- nses_parameters
dim(nses_parameters6)

load(file.path(hpc_data, "nses_parameters9.RData"))
nses_parameters9 <- nses_parameters
dim(nses_parameters9)

load(file.path(hpc_data, "nses_parameters11.RData"))
nses_parameters11 <- nses_parameters
dim(nses_parameters11)
# View(nses_parameters11)

# parameter inputs for each generation
load(file.path(hpc_data, "parameter_input_sims0.RData"))
parameter_inputs0 <- parameter_input_sims
dim(parameter_inputs0)
simulated_parameter_list <- colnames(parameter_inputs0)

load(file.path(hpc_data, "parameter_input_sims3.RData"))
parameter_inputs3 <- parameter_input_sims
dim(parameter_inputs3)

load(file.path(hpc_data, "parameter_input_sims6.RData"))
parameter_inputs6 <- parameter_input_sims
dim(parameter_inputs6)

load(file.path(hpc_data, "parameter_input_sims9.RData"))
parameter_inputs9 <- parameter_input_sims
dim(parameter_inputs9)

load(file.path(hpc_data, "parameter_input_sims11.RData"))
parameter_inputs11 <- parameter_input_sims
dim(parameter_inputs11)

# list of output
load(file.path(hpc_data, "bac_cal0.RData"))
bac_cal_output0 <- bac_cal_output
dim(bac_cal_output0)

load(file.path(hpc_data, "bac_cal3.RData"))
bac_cal_output3 <- bac_cal_output
dim(bac_cal_output3)

load(file.path(hpc_data, "bac_cal6.RData"))
bac_cal_output6 <- bac_cal_output
dim(bac_cal_output6)

load(file.path(hpc_data, "bac_cal9.RData"))
bac_cal_output9 <- bac_cal_output
str(bac_cal_output9)

load(file.path(hpc_data, "bac_cal11.RData"))
bac_cal_output11 <- bac_cal_output
str(bac_cal_output11)

# these are the initial parameters (prior distributions) for the sensitivity analysis
par_bound <- tibble(
  #Hydrology
  "CN2.mgt|change = relchg"= c(-0.3,0.3),
  "SOL_K(1).sol|change = relchg" = c(-0.8,0.8),
  "SOL_AWC(1).sol|change = relchg" = c(-0.8,2),
  "OV_N.hru|change = relchg" = c(-0.8,2),
  "ALPHA_BF.gw|change = relchg" = c(-0.3,0.3),
  "GW_DELAY.gw|change = relchg" = c(-0.75,4),
  "GWQMN.gw|change = relchg" = c(-0.5,2),
  # "HRU_SLP.hru|change = absval" = c(0,1),
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

par_bound
dim(par_bound)
colnames(par_bound)

# empirical observations
# empirical flow observations
load(file.path(hpc_data_sensitivity, "q_obs.RData"))
dim(q_obs)
colnames(q_obs)

# empirical bacteria concentrations
load(file.path(hpc_data_sensitivity, "bac_obs.RData"))
dim(bac_obs)
colnames(bac_obs)

# rainfall
load(file.path(hpc_data_sensitivity, "pcp_obs.RData"))
dim(pcp_obs)
colnames(pcp_obs)


