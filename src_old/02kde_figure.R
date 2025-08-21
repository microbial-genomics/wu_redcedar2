# Example: Replace these with your real data

# change to NAs
means_matrix <- matrix(runif(6 * 18, 0, 10), nrow = 6)   # shape: 6 x 18
dim(means_matrix)
sds_matrix   <- matrix(runif(6 * 18, 0.5, 2), nrow = 6)  # shape: 6 x 18
mins_matrix   <- matrix(runif(6 * 18, 0.5, 2), nrow = 6)  # shape: 6 x 18
maxs_matrix   <- matrix(runif(6 * 18, 0.5, 2), nrow = 6)  # shape: 6 x 18

parameter_names

#"CN2.mgt|change = relchg"= runif(nsims, -0.3, 0.3), #
mins_matrix[1:4,1] <- -0.25
maxs_matrix[1:4,1] <- 0.1
#"GWQMN.gw|change = relchg" = runif(nsims, -0.5, 2), #
mins_matrix[1:4,2] <- -0.5
maxs_matrix[1:4,2] <- 0.5
#"ALPHA_BNK.rte|change = absval" =runif(nsims, 0, 1), #large for flat recessions, and small for steep recessions
mins_matrix[1:4,3] <- 0.5
maxs_matrix[1:4,3] <- 1
#"CH_K2.rte|change = absval" = runif(nsims, 0, 50), # changed from(0,250) 
mins_matrix[1:4,4] <- 0
maxs_matrix[1:4,4] <- 50
#"CH_N2.rte|change = absval" = runif(nsims, 0.05, 0.15), # changed from (0,0.1) 
mins_matrix[1:4,5] <- 0
maxs_matrix[1:4,5] <- 0.1
# TRNSRCH.bsn|change = absval" = runif(nsims, 0, 0.3), # default is 0.00
mins_matrix[1:4,6] <- 0
maxs_matrix[1:4,6] <- 0.3
#CH_N1
mins_matrix[1:4,7] <- 0
maxs_matrix[1:4,7] <- 0.15
#"CH_K1.sub|change = absval" = runif(nsims, 0, 120), #For prennial streams with continuous groundwater contribution, the effective conductivity will be zero.
mins_matrix[1:4,8] <- 0
maxs_matrix[1:4,8] <- 300
#"RCHRG_DP.gw|change = absval" = runif(nsims, 0, 0.5), # default range is (0,1)
mins_matrix[1:4,9] <- 0
maxs_matrix[1:4,9] <- 1
#"SFTMP.bsn|change = absval"= runif(nsims, -5, 5), # changed from (-2,2),default is 1.0
mins_matrix[1:4,10] <- -5
maxs_matrix[1:4,10] <- 5
#"SMTMP.bsn|change = absval"= runif(nsims, -5, 5), # changed from (-2,2),default is 1.0
mins_matrix[1:4,11] <- -5
maxs_matrix[1:4,11] <-5
#"DEP_IMP.hru|change = absval"= runif(nsims, 2000, 6000), # "generation 21 DEP_IMP 4010.712 31.167"
mins_matrix[1:4,12] <- 0
maxs_matrix[1:4,12] <- 6000
#"DDRAIN.mgt|change = absval"= runif(nsims, 500, 1500), # "generation 21 DDRAIN 1105.285 107.537"
mins_matrix[1:4,13] <- 0
maxs_matrix[1:4,13] <- 2000
#"GDRAIN.mgt|change = absval"= runif(nsims, 0, 50), # "generation 21 GDRAIN 24.686 11.354"
mins_matrix[1:4,14] <- 0
maxs_matrix[1:4,14] <- 72
#BACKTKDQ
mins_matrix[1:4,15] <- 0
maxs_matrix[1:4,15] <- 500
#BACT_SWF
mins_matrix[1:4,16] <- 0
maxs_matrix[1:4,16] <- 1
#"THBACT.bsn|change = absval"= runif(nsims, 0, 2), # default value 1.07
mins_matrix[1:4,17] <- 0
maxs_matrix[1:4,17] <- 10
#"WDPRCH.bsn|change = absval"= runif(nsims, 0, 1) # 
mins_matrix[1:4,18] <- 0 
maxs_matrix[1:4,18] <- 1
  
# not carried forward?
#"OV_N.hru|change = relchg" = runif(nsims, 0.01, 0.4),
#"GW_DELAY.gw|change = relchg" = runif(nsims, -0.75,4),
#"SLSUBBSN.hru|change = relchg" = runif(nsims, -0.5, 1),
#"GW_REVAP.gw|change = absval" = runif(nsims, 0.02, 0.2),
#"TIMP.bsn|change = absval"= runif(nsims, 0.01, 1),

# assign values
for(i in 1:18){
  means_matrix[1,i] <- fitted_parameter_list0[[i]]$estimate[[1]]
  sds_matrix[1,i]<- fitted_parameter_list0[[i]]$estimate[[2]]
  means_matrix[2,i] <- fitted_parameter_list1[[i]]$estimate[[1]]
  sds_matrix[2,i]<- fitted_parameter_list1[[i]]$estimate[[2]]
  means_matrix[3,i] <- fitted_parameter_list4[[i]]$estimate[[1]]
  sds_matrix[3,i]<- fitted_parameter_list4[[i]]$estimate[[2]]
  means_matrix[4,i] <- fitted_parameter_list7[[i]]$estimate[[1]]
  sds_matrix[4,i]<- fitted_parameter_list7[[i]]$estimate[[2]]
  means_matrix[5,i] <- fitted_parameter_list10[[i]]$estimate[[1]]
  sds_matrix[5,i]<- fitted_parameter_list10[[i]]$estimate[[2]]
  means_matrix[6,i] <- fitted_parameter_list13[[i]]$estimate[[1]]
  sds_matrix[6,i]<- fitted_parameter_list13[[i]]$estimate[[2]]
}

# Assign level and panel names
levels <- c("gen0", "gen1", "gen4", "gen7", "gen10", "gen13")
levels <- factor(levels, levels = c("gen0", "gen1", "gen4", "gen7", "gen10", "gen13"))
levels <- fct_relevel(levels, "gen0", "gen1", "gen4", "gen7", "gen10", "gen13")


panels <- parameter_names
rownames(means_matrix) <- levels
rownames(sds_matrix) <- levels
rownames(mins_matrix) <- levels
rownames(maxs_matrix) <- levels
colnames(means_matrix) <- panels
colnames(sds_matrix) <- panels
colnames(mins_matrix) <- panels
colnames(maxs_matrix) <- panels

means_matrix
sds_matrix
mins_matrix
maxs_matrix



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

ggplot(df_curves, aes(x = x, y = density, color = level)) +
  geom_line(size = 1) +
  facet_wrap(~ panel, ncol = 6, scales = "free") +
  theme_minimal() +
  labs(
    title = "Posterior Densities by Parameter (Truncated to Initial Uniform Prior)",
    x = "Value",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1")








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
