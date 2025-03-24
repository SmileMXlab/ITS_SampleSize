##############################################
### simcheck07.R: Assess method failures.  ###
##############################################
  
# Load relevant libraries and source simcheck05.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(foreach) # for using foreach loop
source("simcheck05.R")

# Set seed
set.seed(12345678)

# run simulation study
# we increase reps to 100 for illustrative purposes
results <- foreach( i = 1:100) %do% {    
  dataframe <- gendata( nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0) )
  resul<-GLMmod(dataframe,greatmod="time+X+Xt")
  attr(resul, "seed")<-.Random.seed  
  resul
}
rresults<-bind_rows(results, .id = "rep")
summary(rresults)


results_Obs <- foreach( i = 1:100) %do% {    
  dataframe_Obs <- gendata_Obs(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.8)
  resul_Obs<- Obsmod(dataframe=dataframe_Obs)
  attr(resul_Obs, "seed")<-.Random.seed  
  resul_Obs
}
rresults_Obs<-bind_rows(results_Obs, .id = "rep")
summary(rresults_Obs)


results_Par <- foreach( i = 1:100) %do% {
  dataframe_Par <- gendata_Par(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.2,sigma_mu2=0.5)
  resul_Par <-Parmod(dataframe=dataframe_Par,greatmod="time+X+Xt")
  attr(resul_Par, "seed")<-.Random.seed  
  resul_Par
}
rresults_Par<-bind_rows(results_Par, .id = "rep")
summary(rresults_Par)



# Inspect results for method failures
rresults[is.na(rresults$biasX),]
rresults[rresults$rep==20,]

rresults_Obs[is.na(rresults_Obs$biasX),]
rresults_Obs[rresults_Obs$rep==72,]

rresults_Par[is.na(rresults_Par$biasX),]
rresults_Par[rresults_Par$rep==95,]


# Reconstruct one data set with method failures
.Random.seed <- attr(results[[19]], "seed") # Need the seed status after running 109 repetitions
dataframe <- gendata( nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0) )
dataframe
mod_GLM <- glm(y~time+X+Xt,data=dataframe,family=poisson)
summary(mod_GLM)


.Random.seed <- attr(results_Obs[[71]], "seed") # Need the seed status after running 109 repetitions
dataframe_Obs <- gendata_Obs(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.8)
dataframe_Obs
y <- dataframe_Obs[,1]
Xtdata <- dataframe_Obs[,-1]
mod_INGARCH <- tsglm(y,model=list(past_obs=1),xreg=Xtdata,link="log", distr="poisson")
summary(mod_INGARCH)


.Random.seed <- attr(results_Par[[94]], "seed") # Need the seed status after running 109 repetitions
dataframe_Par<- gendata_Par(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.2,sigma_mu2=0.5)
dataframe_Par
mod_GLM <- glm(y~time+X+Xt,data=dataframe_Par,family=poisson)
summary(mod_GLM)



# Save results
save(results,rresults, file="simcheck07_results.RData")
save(results_Obs,rresults_Obs, file="simcheck07_results_Obs.RData")
save(results_Par,rresults_Par, file="simcheck07_results_Par.RData")




