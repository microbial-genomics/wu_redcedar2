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

# 0, 1, 5, 10, 15, 19
# load parameters but in list format
# 18 parameter fits for each variable
load(file.path(hpc_data, "fitted_parameters0.RData"))
fitted_parameter_list0 <- fitted_parameter_list
str(fitted_parameter_list)
fitted_parameter_list0[[2]]$estimate
fitted_parameter_list0[[5]]$estimate
fitted_parameter_list0[[8]]$estimate
fitted_parameter_list0[[18]]$estimate

load(file.path(hpc_data, "fitted_parameters1.RData"))
fitted_parameter_list1 <- fitted_parameter_list
fitted_parameter_list1[[2]]$estimate

load(file.path(hpc_data, "fitted_parameters4.RData"))
fitted_parameter_list4 <- fitted_parameter_list
fitted_parameter_list4[[2]]$estimate

load(file.path(hpc_data, "fitted_parameters7.RData"))
fitted_parameter_list7 <- fitted_parameter_list
fitted_parameter_list7[[2]]$estimate

load(file.path(hpc_data, "fitted_parameters10.RData"))
fitted_parameter_list10 <- fitted_parameter_list
fitted_parameter_list10[[2]]$estimate

load(file.path(hpc_data, "fitted_parameters13.RData"))
fitted_parameter_list13 <- fitted_parameter_list
fitted_parameter_list13[[2]]$estimate

# this loads everything for the generation stats
load(file.path(hpc_data, "generation_stats13.RData"))
generation_stats <- generation_stats[1:14,]

# nse scores for each generation
load(file.path(hpc_data, "nses_parameters0.RData"))
nses_parameters0 <- nses_parameters
dim(nses_parameters0)
parameter_names <- colnames(nses_parameters0)[5:22]
parameter_names

load(file.path(hpc_data, "nses_parameters1.RData"))
nses_parameters1 <- nses_parameters
dim(nses_parameters1)

load(file.path(hpc_data, "nses_parameters4.RData"))
nses_parameters4 <- nses_parameters
dim(nses_parameters4)

load(file.path(hpc_data, "nses_parameters7.RData"))
nses_parameters7 <- nses_parameters
dim(nses_parameters7)

load(file.path(hpc_data, "nses_parameters10.RData"))
nses_parameters10 <- nses_parameters
dim(nses_parameters10)

load(file.path(hpc_data, "nses_parameters13.RData"))
nses_parameters13 <- nses_parameters
dim(nses_parameters13)

# parameter inputs for each generation
load(file.path(hpc_data, "parameter_input_sims0.RData"))
parameter_inputs0 <- parameter_input_sims
dim(parameter_inputs0)

load(file.path(hpc_data, "parameter_input_sims1.RData"))
parameter_inputs1 <- parameter_input_sims
dim(parameter_inputs1)

load(file.path(hpc_data, "parameter_input_sims4.RData"))
parameter_inputs4 <- parameter_input_sims
dim(parameter_inputs4)

load(file.path(hpc_data, "parameter_input_sims7.RData"))
parameter_inputs7 <- parameter_input_sims
dim(parameter_inputs7)

load(file.path(hpc_data, "parameter_input_sims10.RData"))
parameter_inputs10 <- parameter_input_sims
dim(parameter_inputs10)

load(file.path(hpc_data, "parameter_input_sims13.RData"))
parameter_inputs13 <- parameter_input_sims
dim(parameter_inputs13)

# list of output
load(file.path(hpc_data, "bac_cal0.RData"))
bac_cal_output0 <- bac_cal_output
dim(bac_cal_output0)

load(file.path(hpc_data, "bac_cal1.RData"))
bac_cal_output1 <- bac_cal_output
dim(bac_cal_output1)

load(file.path(hpc_data, "bac_cal4.RData"))
bac_cal_output4 <- bac_cal_output
dim(bac_cal_output4)

load(file.path(hpc_data, "bac_cal7.RData"))
bac_cal_output7 <- bac_cal_output
str(bac_cal_output7)

load(file.path(hpc_data, "bac_cal10.RData"))
bac_cal_output10 <- bac_cal_output
dim(bac_cal_output10)

load(file.path(hpc_data, "bac_cal13.RData"))
bac_cal_output13 <- bac_cal_output
str(bac_cal_output13)

