library(SWATplusR)
library(tibble)

# specify base directory
if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  base_dir <- file.path("c:", "git", "wu_redcedar2") #laptop
  src_dir <- file.path(base_dir, "src")
  data_in_dir <- file.path(base_dir,"data_in")
}else{
  base_dir <- "/work/OVERFLOW/stp/MSU/" #hpc
  src_dir <- base_dir
  data_in_dir <- base_dir
  swat_path <- base_dir
}

# source support functions
source(file.path(src_dir,"abc_functions.R"))

# every generation will have 5000 accepted particles
# the median score of these candidates will be used as the cutoff for the next generation
nsims=5000
pars_initial <- create_tibble_initial(nsims)
simulate_generation_zero(nsims, swat_path, base_dir, pars_initial)
