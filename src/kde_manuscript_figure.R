
## setup root directory path
if(Sys.info()[4]=="LZ26TPURUCKE-2"){
  # tom epa windows
  of_root <- file.path("c:", "git", "wu_redcedar2")
}

data_in_dir <- file.path(of_root, "data_in")
graphics_dir <- file.path(of_root, "graphics")
src_dir <- file.path(of_root, "src")
hpc_data <- file.path(of_root, "hpc_data")

# load initials
load("data/myfile.RData")


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