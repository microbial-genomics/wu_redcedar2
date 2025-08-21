# create data set of distribution moments for kde figure
simulated_parameter_list

# change to NAs
moments_matrix1 <- matrix(NA, nrow = 6, ncol=18)   # shape: 6 x 18
moments_matrix2 <- matrix(NA, nrow = 6, ncol=18)   # shape: 6 x 18

colnames(moments_matrix1) <- simulated_parameter_list
colnames(moments_matrix2) <- simulated_parameter_list

#parameter_names
# manually input original priors for sensitivity analysis
# these are mins and maxes of a uniform for the initial sensitivity analysis
# CN2.mgt|change = relchg= c(-0.3,0.3)
moments_matrix1[1,1] <- -0.3
moments_matrix2[1,1] <- 0.3
# GWQMN.gw|change = relchg = c(-0.5,2)
moments_matrix1[1,2] <- -0.5
moments_matrix2[1,2] <- 2
# ALPHA_BNK.rte|change = absval =c(0, 1) #large for flat recessions, and small for steep recessions
moments_matrix1[1,3] <- 0
moments_matrix2[1,3] <- 1
# CH_K2.rte|change = absval = c(0,50), # changed from(0,250) 
moments_matrix1[1,4] <- 0
moments_matrix2[1,4] <- 50
# CH_N2.rte|change = absval = c(0.05, 0.15), # changed from (0,0.1) 
moments_matrix1[1,5] <- 0.05
moments_matrix2[1,5] <- 0.15
# TRNSRCH.bsn|change = absval = c(0,0.3) # default is 0.00
moments_matrix1[1,6] <- 0
moments_matrix2[1,6] <- 0.3
# CH_N1.sub|change = absval = c(0.05, 0.15)
moments_matrix1[1,7] <- 0.05
moments_matrix2[1,7] <- 0.15
# CH_K1.sub|change = absval = c(0, 300) #For perennial streams with continuous groundwater contribution, the effective conductivity will be zero.
moments_matrix1[1,8] <- 0
moments_matrix2[1,8] <- 300
# RCHRG_DP.gw|change = absval = c(0, 1) # default range is (0,1)
moments_matrix1[1,9] <- 0
moments_matrix2[1,9] <- 1
# SFTMP.bsn|change = absval= c(-5, 5) # changed from (-2,2),default is 1.0
moments_matrix1[1,10] <- -5
moments_matrix2[1,10] <- 5
# SMTMP.bsn|change = absval= c(-5,5) # changed from (-2,2), default is 1.0
moments_matrix1[1,11] <- -5
moments_matrix2[1,11] <- 5
# DEP_IMP.hru|change = absval= c(0,6000) # "generation 21 DEP_IMP 4010.712 31.167"
moments_matrix1[1,12] <- 0
moments_matrix2[1,12] <- 6000
# DDRAIN.mgt|change = absval= c(0, 2000) # "generation 21 DDRAIN 1105.285 107.537"
moments_matrix1[1,13] <- 0
moments_matrix2[1,13] <- 2000
# GDRAIN.mgt|change = absval= c(0, 100) # "generation 21 GDRAIN 24.686 11.354"
moments_matrix1[1,14] <- 0
moments_matrix2[1,14] <-100
# BACTKDQ.bsn|change = absval = c(0, 500)
moments_matrix1[1,15] <- 0
moments_matrix2[1,15] <- 500
# BACT_SWF.bsn|change = absval = c(0, 1)
moments_matrix1[1,16] <- 0
moments_matrix2[1,16] <- 1
# THBACT.bsn|change = absval= c(0, 10) # default value 1.07
moments_matrix1[1,17] <- 0
moments_matrix2[1,17] <- 10
# WDPRCH.bsn|change = absval= c(0, 1) # 
moments_matrix1[1,18] <- 0 
moments_matrix2[1,18] <- 1

moments_matrix1
moments_matrix2
  
# assign values of the posterior fitted distribution
# after each generation to display
fitted_parameter_list0
for(i in 1:18){
  moments_matrix1[2,i] <- fitted_parameter_list0[[i]]$estimate[[1]]
  moments_matrix2[2,i]<- fitted_parameter_list0[[i]]$estimate[[2]]
  moments_matrix1[3,i] <- fitted_parameter_list3[[i]]$estimate[[1]]
  moments_matrix2[3,i]<- fitted_parameter_list3[[i]]$estimate[[2]]
  moments_matrix1[4,i] <- fitted_parameter_list6[[i]]$estimate[[1]]
  moments_matrix2[4,i]<- fitted_parameter_list6[[i]]$estimate[[2]]
  moments_matrix1[5,i] <- fitted_parameter_list9[[i]]$estimate[[1]]
  moments_matrix2[5,i]<- fitted_parameter_list9[[i]]$estimate[[2]]
  moments_matrix1[6,i] <- fitted_parameter_list11[[i]]$estimate[[1]]
  moments_matrix2[6,i]<- fitted_parameter_list11[[i]]$estimate[[2]]
}
moments_matrix1
moments_matrix2

# Assign level and panel names
levels <- c("gen0", "gen1", "gen4", "gen7", "gen10", "gen12")
#levels <- factor(levels, levels = c("gen0", "gen3", "gen6", "gen9", "gen11"))
#levels <- fct_relevel(levels, "gen0", "gen3", "gen6", "gen9", "gen11")

rownames(moments_matrix1) <- levels
rownames(moments_matrix2) <- levels

# stopping here 8/21/2025

# Convert to tidy long format
means_long <- as.data.frame(means_matrix) %>%
  mutate(level = rownames(means_matrix)) %>%
  pivot_longer(-level, names_to = "panel", values_to = "mean")

sds_long <- as.data.frame(sds_matrix) %>%
  mutate(level = rownames(sds_matrix)) %>%
  pivot_longer(-level, names_to = "panel", values_to = "sd")

mins_long <- as.data.frame(mins_matrix) %>%
  mutate(level = rownames(mins_matrix)) %>%
  pivot_longer(-level, names_to = "panel", values_to = "mins")

maxs_long <- as.data.frame(maxs_matrix) %>%
  mutate(level = rownames(maxs_matrix)) %>%
  pivot_longer(-level, names_to = "panel", values_to = "maxs")

parameters2 <- left_join(means_long, sds_long, by = c("level", "panel"))
parameters1 <- left_join(parameters2, mins_long, by = c("level", "panel"))
parameters <- as.data.frame(left_join(parameters1, maxs_long, by = c("level", "panel")))

parameters
colnames(parameters)



# Suppose your long data is in a dataframe called 'df'
df_curves <- parameters %>%
  group_by(panel, level) %>%
  group_modify(
    ~{
      x = seq(.x$mins, .x$maxs, length.out = 200)
      data.frame(
        x = x,
        density = dnorm(x, mean = .x$mean, sd = .x$sd),
        mean = .x$mean,
        sd = .x$sd,
        min = .x$mins,
        max = .x$maxs
      )
    }
  ) %>%
  ungroup()

colnames(df_curves)
View(df_curves)

# since these are overlapping kde plots can't use scales = free
# need to manually specify the limits
scales_xranges <- list(
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)), #
  scale_x_continuous(limits = c(X, X)) #
)


ggplot(df_curves, aes(x = x, y = density, color = level)) +
  geom_line(size = 1) +
  facet_wrap(~ panel, ncol = 3, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Posterior Densities by Parameter (Truncated to Initial Uniform Prior)",
    x = "Value",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1")


library(patchwork)
plots <- lapply(
  split(df_curves, df_curves$panel),
  function(d) {
    ggplot(d, aes(x = x, y = density, color = level)) +
      geom_line(size = 1) +
      xlim(min, max) +
      theme_minimal() +
      labs(x = "Value", y = "Density")
  }
)
wrap_plots(plots, ncol = 6)





# Create a grid of x-values covering the full range expected in all densities
x_min <- min(parameters$mean - 3 * parameters$sd)
x_max <- max(parameters$mean + 3 * parameters$sd)
x_vals <- seq(x_min, x_max, length.out = 200)

# For each panel and level, calculate the density
plot_data <- parameters %>%
  group_by(panel, level) %>%
  do({
    data.frame(
      x = x_vals,
      density = dnorm(x_vals, mean = .$mean, sd = .$sd),
      level = .$level,
      panel = .$panel
    )
  }) %>%
  ungroup()

#View(plot_data)



ggplot(plot_data, aes(x = x, y = density, color = level)) +
  geom_line(size = 1) +
  facet_wrap(~ panel, ncol = 6, scales='free') +
  theme_minimal() +
  labs(
    title = "Updated Posterior Densities of Sensitive Parameters",
    x = "Value", y = "Density"
  ) +
  scale_color_brewer(palette = "Set1")




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
