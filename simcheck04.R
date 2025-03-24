##############################################
### simcheck04.R: Run a few iterations.    ###
##############################################

# Load relevant libraries and source simcheck02.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(foreach) # for using foreach loop
source("simcheck02.R")

# Set seed
set.seed(12345678)

# run 3 repetitions
results <- foreach( i = 1:3, .combine="rbind") %do% {
  dataframe <- gendata( nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)) )
  GLMmod(dataframe,greatmod="time+X+Xt")
}

results_Obs <- foreach( i = 1:3, .combine="rbind") %do% {
  dataframe_Obs <- gendata_Obs(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2)
  Obsmod(dataframe=dataframe_Obs)
}

results_Par <- foreach( i = 1:3, .combine="rbind") %do% {
  dataframe_Par <- gendata_Par(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2,sigma_mu2=0.5)
  Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")
}

# view results
results
results_Obs 
results_Par