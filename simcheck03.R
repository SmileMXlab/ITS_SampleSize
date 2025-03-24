################################################################
### simcheck03.R: Simulate and analyse a single data set.    ###
################################################################

# Load relevant libraries and source simcheck02.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
source("simcheck02.R")

# Set seed
set.seed(12345678)


# generate a single large data set
dataframe <- gendata(nsmp=100000,pre_nsmp=50000,beta=c(log(2),log(2),log(1.5),log(4)))
dataframe_Obs <- gendata_Obs(nsmp=100000,pre_nsmp=50000,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2)
dataframe_Par <- gendata_Par(nsmp=100000,pre_nsmp=50000,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2,sigma_mu2=0.5)


# summarise the data
summary(dataframe)
summary(dataframe_Obs)
summary(dataframe_Par)


# analyse the data
GLMmod(dataframe=dataframe,greatmod="time+X+Xt")

Obsmod(dataframe=dataframe_Obs)

Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")


