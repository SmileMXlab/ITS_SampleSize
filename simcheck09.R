############################################
### simcheck09.R: Investigate outliers.  ###
############################################

# Load relevant libraries and source simcheck05.R
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(foreach) # for using foreach loop
library(dplyr) # To transform list fo data frame in single data frame
source("simcheck05.R")

load("simcheck07_results.RData")
load("simcheck07_results_Obs.RData")
load("simcheck07_results_Par.RData")

# Find where very large se occurs
large.ses<-unique(rresults[rresults$seX>100,"rep"])
large.ses_Obs<-unique(rresults_Obs[rresults_Obs$seX>100,"rep"])
large.ses_Par<-unique(rresults_Par[rresults_Par$seX>100,"rep"])

# let's pick out the second rep
rresults[rresults$rep==2,]
rresults_Obs[rresults_Obs$rep==4,]
rresults_Par[rresults_Par$rep==3,]

# Reconstruct this data set with method failures
.Random.seed <- attr(results[[1]], "seed") # Need the seed status after running 1 repetitions
dataframe <- gendata( nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0) )

.Random.seed <- attr(results_Obs[[3]], "seed") # Need the seed status after running 1 repetitions
dataframe_Obs <- gendata_Obs(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.8)

.Random.seed <- attr(results_Par[[2]], "seed") # Need the seed status after running 1 repetitions
dataframe_Par <- gendata_Par(nsmp=100,pre_nsmp=50,beta=c(log(0.1),log(0.1),0,0),rho=0.2,sigma_mu2=0.5)

# and explore it
dataframe
mod_GLM <- glm(y~time+X+Xt,data=dataframe,family=poisson)
summary(mod_GLM)

dataframe_Obs
y <- dataframe_Obs[,1]
Xtdata <- dataframe_Obs[,-1]
mod_INGARCH <- tsglm(y,model=list(past_obs=1),xreg=Xtdata,link="log", distr="poisson")
summary(mod_INGARCH)

dataframe_Par
mod_GLM <- glm(y~time+X+Xt,data=dataframe_Par,family=poisson)
Xmatrix <- as.matrix(dataframe_Par[,-1])
dummy<-rhofun(dataframe$y,Xmatrix,coef(mod_GLM))
dummyUB<-rhofunUB(dataframe$y,Xmatrix,coef(mod_GLM),dummy$GEEsigma,dummy$GEErho)
vovmax<-correctcov(Xmatrix,coef(mod_GLM),dummyUB$GEEsigmaUB,dummyUB$GEErhomax)
vovmax$sebeta