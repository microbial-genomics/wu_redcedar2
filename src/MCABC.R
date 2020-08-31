
module load intel/19.0.5
module load R/3.6.2
module load geos/3.8.0
module load gdal-2.4.3/intel-19.0
module load proj-5.2.0/intel-19.0
module load udunits-2.2.26/intel-19.0
module load imagemagick


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
library(fast)
library(forcats)
library(lubridate)
library(sensitivity)
library(fitdistrplus)


load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file='/work/OVERFLOW/RCR/sensitivity/msu_sobol/q_obs.RData')
load(file='/work/OVERFLOW/RCR/MSU/bac_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/q_obs2.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/flux_obs.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs.RData')
load(file = '/work/OVERFLOW/RCR/calibration/MSU/pcp_obs2.RData')

load(file = '/work/OVERFLOW/RCR/calibration/MSU/bac_cal12.RData')




# steps for the Approximate Bayesian Computation with sequential Monte Carlo simulation (ABC-MC):##############

#Step 1 You have already done 50k simulations from the uniform priors. We are going to take the first 10k (not the best 10k, do not sort) of these as the first generation of simulations. We are simply throwing away the last 40k and not using them for the official simulations for the manuscript. This is necessary to honor the assumptions behind ABC-MC.

sim_bac<- bac_cal1$simulation$bac_out[,1:10001]
sim_q <- bac_cal1$simulation$q_out[,1:10001]
sim_pars <- bac_cal1$parameter$values[1:10000,]

#Step 2 Calculate the average of the 3 nses for each of the 10k simulations: nse_average = mean(nse_conc, nse_flow, nse_ flux)

nse_bac <- right_join(sim_bac,bac_obs,by="date")%>%
  select(-date) %>% select(-bacteria) %>%
  map_dbl(., ~NSE(.x, bac_obs$bacteria))

 sort(nse_bac, decreasing = T) %>% enframe()

nse_q <- right_join(sim_q,q_obs,by="date") %>%
  select(-date) %>% select(-discharge) %>%
  map_dbl(., ~NSE(.x, q_obs$discharge))
 
 sort(nse_q, decreasing = T) %>% enframe()


flux_sim <- sim_bac[c(97:167),c(-1)]*
            sim_q  [c(97:167), c(-1)]*10^4
nse_flux <- flux_sim %>%
     map_dbl(., ~NSE(.x, flux_obs[,1]))
 
 sort(nse_flux, decreasing = T) %>% enframe()


nse_average <- matrix(data=NA, nrow=10000, ncol=1)
	for(i in 1:10000){
		print(i)
  nse_average[i] <- mean(c(nse_bac[i], nse_q[i], nse_flux[i]))
					}
sort(nse_average, decreasing = T) %>% enframe()


#Step3.	Take the top 5k simulations based on the nse_average, discard the other 5K simulations.
#run<-colnames(sim_bac)

#run <- run[2:10001]
#nse_average <- nse_average
#nse_mean <-cbind(run,nse_average)
#top5k_sim <- nse_mean[order(-nse_average,run),][1:5000,]
#valid_sims <-c(top5k[,1])
#sim_bac_2 <-sim_bac[valid_sims]
#sim_q_2 <-sim_q[valid_sims]


run<-c(1:10000)

nse_average <- nse_average
nse_mean <-cbind(run,nse_average)
top5k_par <- nse_mean[order(-nse_average,run),][1:2500,]
median(top5k_par[,2])


valid_pars <-c(top5k_par[,1])
sim_pars_2 <- sim_pars[valid_pars,]

##Step 4 Use the 5k set of inputs associated with these top 5k simulations. Calculate the unweighted kernel densities using the kde package and fit to a normal distribution, truncate at the range limits for each parameter.


kde_mcabc <- sim_pars_2 %>% 
  #filter(nse > -10) %>%
  gather(key = "par", value = "parameter_range")

     ggplot(data = kde_mcabc) +
   geom_density(aes(x = parameter_range)) +
   facet_wrap(.~par, nrow=5, scales = "free") +
    theme_bw()
ggsave("/home/hwu/wu_redcedar2/graphics/kde_mcabc.pdf")
###


 ####.	

	
#5.	Find the median_average_nse of these 5k nse_averages, this will be the average of the 2500 and 2501st highest average_nse.
#6.	Use the truncated normal distributions from 4) to set up the next round of simulations. Simulate with these inputs. For each simulation, calculate the average_nse. Only keep the simulation if the average_nse is higher that the median_average_nse from the last set of simulations.
#7.	Keep simulating until we get 5k new simulations with an average_nse greater than the median_average_nse calculated for the previous set of 5k simulation results.
#8.	Repeat step 4 to calculate the updated unweighted kernel densities and fit to the normal distribution.
#9.	Now repeat step 5 for the new set of 5k simulations. Calculate the updated median_average_nse.
#10.	Repeat steps 6)-9) over and over until the median_average_nse fails to improve by X% versus the median_average_nse from the last generation. We have not explicitly defined what this percentage is just yet. It will probably be something like 1% or less.



