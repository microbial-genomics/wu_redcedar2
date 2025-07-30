
## setup root directory path
if(Sys.info()[4]=="LZ26TPURUCKE-2"){
  # tom epa windows
  of_root <- file.path("c:", "Users", "tpurucke", "git", "wu_redcedar2")
}

data_in_dir <- file.path(of_root, "data_in")
graphics_dir <- file.path(of_root, "graphics")
src_dir <- file.path(of_root, "src")
hpc_data <- file.path(of_root, "hpc_data")

# load parameters but in list format
load(file.path(hpc_data, "fitted_parameters0.RData"))
str(fitted_parameter_list)

# this loads everything for the generation stats
load(file.path(hpc_data, "generation_stats18.RData"))
generation_stats <- generation_stats[1:19,]

# nse scores for each generation
load(file.path(hpc_data, "nses_parameters0.RData"))
nses_parameters0 <- nses_parameters
dim(nses_parameters0)

load(file.path(hpc_data, "nses_parameters6.RData"))
nses_parameters6 <- nses_parameters
dim(nses_parameters6)

load(file.path(hpc_data, "nses_parameters12.RData"))
nses_parameters12 <- nses_parameters
dim(nses_parameters12)

load(file.path(hpc_data, "nses_parameters18.RData"))
nses_parameters18 <- nses_parameters
dim(nses_parameters18)

# parameter inputs for each generation
load(file.path(hpc_data, "parameter_input_sims0.RData"))
parameter_inputs0 <- parameter_input_sims
dim(parameter_inputs0)

load(file.path(hpc_data, "parameter_input_sims6.RData"))
parameter_inputs6 <- parameter_input_sims
dim(parameter_inputs6)

load(file.path(hpc_data, "parameter_input_sims12.RData"))
parameter_inputs12 <- parameter_input_sims
dim(parameter_inputs12)

load(file.path(hpc_data, "parameter_input_sims18.RData"))
parameter_inputs18 <- parameter_input_sims
dim(parameter_inputs18)

# list of output
load(file.path(hpc_data, "bac_cal0.RData"))
bac_cal_output0 <- bac_cal_output
dim(bac_cal_output0)



plot_kde_pdf <- function(iter, kde_next_gen){
  ggplot(data = kde_next_gen) +
    geom_density(aes(x = parameter_range)) +
    facet_wrap(.~par, nrow=5, scales = "free") +
    theme_bw()
  density_plot_filename <- paste("kde_mcabc_gen", iter, ".pdf", sep="")
  ggsave(file.path(graphics_dir, density_plot_filename))  
}



## find the updated unweighted kernel densities based on these new 5k simulations
kde_next_gen <- sim_pars[valid_keepers,] %>% 
  gather(key = "par", value = "parameter_range")

plot_kde_pdf(iter, kde_next_gen)