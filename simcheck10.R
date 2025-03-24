################################################
### simcheck11.R: Check Monte Carlo errors.  ###
################################################

# Load relevant libraries and source simcheck05.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(foreach)
library(rsimsum) 
source("simcheck05.R")

# Set seed
set.seed(12345678)

# run simulation study
results <- foreach( i = 1:1000) %do% {    
  dataframe <- gendata( nsmp=100,pre_nsmp=50,beta=c(log(2),log(2),0,0) )
  resul<-GLMmod(dataframe,greatmod="time+X+Xt")
  attr(resul, "seed")<-.Random.seed  
  resul
}
rresults<-bind_rows(results, .id = "rep")
summary(rresults)


results_Obs <- foreach( i = 1:1000) %do% {    
  dataframe_Obs <- gendata_Obs(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),0,0),rho=0.2)
  resul_Obs<- Obsmod(dataframe=dataframe_Obs)
  attr(resul_Obs, "seed")<-.Random.seed  
  resul_Obs
}
rresults_Obs<-bind_rows(results_Obs, .id = "rep")
summary(rresults_Obs)


results_Par <- foreach( i = 1:1000) %do% {
  dataframe_Par <- gendata_Par(nsmp=500,pre_nsmp=250,beta=c(log(2),log(2),0,0),rho=0.2,sigma_mu2=0.5)
  resul_Par <-Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")
  attr(resul_Par, "seed")<-.Random.seed  
  resul_Par
}
rresults_Par<-bind_rows(results_Par, .id = "rep")
summary(rresults_Par)
Sys.time()


sX <- simsum(data = rresults, estvarname = "biasX", true = 0, se = "seX")
sXt <- simsum(data = rresults, estvarname = "biasXt", true = 0, se = "seXt")
summary(sX)
summary(sXt)

sX_Obs <- simsum(data = rresults_Obs, estvarname = "biasX", true = 0, se = "seX")
sXt_Obs <- simsum(data = rresults_Obs, estvarname = "biasXt", true = 0, se = "seXt")
summary(sX_Obs)
summary(sXt_Obs)

sX_Par <- simsum(data = rresults_Par, estvarname = "biasX", true = 0, se = "seX")
sXt_Par <- simsum(data = rresults_Par, estvarname = "biasXt", true = 0, se = "seXt")
summary(sX_Par)
summary(sXt_Par)

