#######################################################################
### simcheck06.R: Make it easy to re-create any simulated data set. ###
#######################################################################

# Load relevant libraries and source simcheck05.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(foreach) # for using foreach loop
source("simcheck05.R")

# Set seed
set.seed(12345678)

# run 3 repetitions
results <- foreach( i = 1:3) %do% {  
  dataframe <- gendata( nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)) )
  resul<-GLMmod(dataframe,greatmod="time+X+Xt")
  attr(resul, "seed")<-.Random.seed 
  resul
}

results_Obs <- foreach( i = 1:3) %do% {  
  dataframe_Obs <- gendata_Obs(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2)
  resul_Obs<-Obsmod(dataframe=dataframe_Obs)
  attr(resul_Obs, "seed")<-.Random.seed 
  resul_Obs
}

results_Par <- foreach( i = 1:3) %do% {  
  dataframe_Par <- gendata_Par(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2,sigma_mu2=0.5)
  resul_Par<-Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")
  attr(resul_Par, "seed")<-.Random.seed 
  resul_Par
}


# view stored results for 3rd repetition
results[[3]]
results_Obs[[3]]
results_Par[[3]]


# reconstruct data set for 3rd repetition and check it gives same results
.Random.seed <- attr(results[[2]], "seed") # Need the seed status after running 2 repetitions
dataframe <- gendata( nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)) )
resul<-GLMmod(dataframe,greatmod="time+X+Xt")
resul            # can verify that these results are the same

.Random.seed <- attr(results_Obs[[2]], "seed") # Need the seed status after running 2 repetitions
dataframe <- gendata_Obs(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2)
resul_Obs<-Obsmod(dataframe=dataframe_Obs)
resul_Obs            # can verify that these results are the same

.Random.seed <- attr(results_Par[[2]], "seed") # Need the seed status after running 2 repetitions
dataframe <- gendata_Par(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),log(1.5),log(4)),rho=0.2,sigma_mu2=0.5)
resul_Par<-Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")
resul_Par           # can verify that these results are the same



