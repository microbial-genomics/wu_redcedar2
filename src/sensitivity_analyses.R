library(dplyr)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(miceadds)
library(sensitivity)
library(vioplot)
library(xts)

# tom epa laptop
if(Sys.info()[4] == "LZ2626UTPURUCKE"){
  rcdir <- path.expand("C:/git/wu_redcedar/")
}

# huiyun laptop
if(Sys.info()[4] == "computer_name"){
  rcdir <- path.expand("C:/dir1/dir2/")
}

# Variable names for data
rcdir_data_in <- paste(rcdir,'data_in/',sep='')
rcdir_data_out <- paste(rcdir,'data_out/',sep='')
rcdir_graphics <- paste(rcdir,'graphics/',sep='')
rcdir_src <- paste(rcdir,'src/',sep='')

# load data_out files
parameter_valid_file <- paste(rcdir_data_out,"parameter_valid.RData", sep="")
file.exists(parameter_valid_file)
load.Rdata(parameter_valid_file, "sim_parameters")

bac_out_file <- paste(rcdir_data_out,"bac_out_valid.RData", sep="")
file.exists(bac_out_file)
load.Rdata(bac_out_file, "sim_bac_concs")

flow_out_file <- paste(rcdir_data_out,"q_out_valid2.RData", sep="")
file.exists(flow_out_file)
load.Rdata(flow_out_file, "sim_flows")

# file info
#View(sim_parameters)
#View(sim_bac_concs)
#View(sim_flows)
dim(sim_parameters)
colnames(sim_parameters)
dim(sim_bac_concs)
colnames(sim_bac_concs)
dim(sim_flows)
colnames(sim_flows)

#fix problematic colnames
colnames(sim_parameters)[2] <- "SOL_K"
colnames(sim_parameters)[3] <- "SOL_AWC"

#simulation is for the 365 days of 2013
#drop date column for sim_bac_concs
sim_bac_concs2 <- sim_bac_concs[,-c(1)]
sim_dates <- sim_bac_concs[,1]
#View(sim_bac_concs2)
dim(sim_bac_concs2)
#transpose bacteria concentration data frame
sim_bac_concs3 <- t(sim_bac_concs2)
dim(sim_bac_concs3)
# pcc on row medians
sim_bac_medians <- rowMedians(sim_bac_concs3)
pcc(sim_parameters, sim_bac_medians)
# pcc on row maxs
sim_bac_maxs <- rowMaxs(sim_bac_concs3)
pcc(sim_parameters, sim_bac_maxs)

# drop dates for sim_flows
sim_flows2 <- sim_flows[,-c(1)]
dim(sim_flows2)
#transpose flows data frame
sim_flows3 <- t(sim_flows2)
dim(sim_flows3)
# pcc on row medians
sim_flows_medians <- rowMedians(sim_flows3)
pcc(sim_parameters, sim_flows_medians)
# pcc on row maxs
sim_flows_maxs <- rowMaxs(sim_flows3)
pcc(sim_parameters, sim_flows_maxs)


# daily pcc for bacteria and flow
bac_pcc <- matrix(data=NA, nrow=365, ncol=40)
flows_pcc <- matrix(data=NA, nrow=365, ncol=40)
for(i in 1:365){
  print(i)
  daily_bac_pcc <- pcc(sim_parameters, sim_bac_concs3[,i])
  daily_flows_pcc <- pcc(sim_parameters, sim_flows3[,i])
  length(bac_pcc[i,])
  length(t(daily_bac_pcc$PCC))
  bac_pcc[i,] <- t(daily_bac_pcc$PCC)
  flows_pcc[i,] <- t(daily_flows_pcc$PCC)
}

#print violin plot
mklab <- function(y_var){
  if(y_var){
    names(mf)[response]
  } else {
    paste(names(mf)[-response], collapse = " : ")
  }
}
dim(bac_pcc)
colnames(bac_pcc) <- colnames(sim_parameters)
colnames(flows_pcc) <- colnames(sim_parameters)
pdf(paste(rcdir_graphics,"pcc_violin.pdf",sep=""),width=55,height=30,onefile=TRUE)
  vioplot(bac_pcc)
  vioplot(flows_pcc)
dev.off()


#create time series
bac_dates <- ts(sim_dates, start=c(2013, 1), end=c(2013, 365), frequency=365)


#simple ggplot
dim(bac_pcc)
pdf(paste(rcdir_graphics,"pcc_ts.pdf",sep=""),width=11,height=8, onefile=TRUE)
  for(i in 1:40){
    data <- data.frame(
      day = as.Date("2014-01-01") - 0:364,
      bac_value = bac_pcc[,i],
      flows_value = flows_pcc[,i]
    )
    p <- ggplot(data, aes(x=day)) +
      geom_line(aes(y=bac_value), color = "darkred") + 
      geom_line(aes(y=flows_value), color = "steelblue", linetype="twodash") +
      xlab("") +
      ylab("pcc sensitivity") +
      ylim(-0.2,0.1)
    print(p + ggtitle(colnames(bac_pcc)[i]))
  }
dev.off()
