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

# run the initial set of swat simulations
swat_output0 <- run_swat_red_cedar(swat_path, pars_initial)

save_file <- paste(base_dir,"rcr_swat_output0.RData",sep="")
save(swat_output0, file=save_file)


