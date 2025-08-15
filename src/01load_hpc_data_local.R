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

# parameter inputs for each generation
load(file.path(hpc_data, "parameter_input_sims0.RData"))
parameter_inputs0 <- parameter_input_sims
dim(parameter_inputs0)

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


