############################################
### simcheck08.R: Check for  outliers.   ###
############################################

# Load relevant libraries, source simcheck05.R and load simcheck08_results.RData
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   
library(ggplot2)
library(tidyr)

source("simcheck05.R")
load("simcheck07_results.RData")
load("simcheck07_results_Obs.RData")
load("simcheck07_results_Par.RData")


# Find where very large se occurs
plotdata<-data.frame(rep=rep(1:100,2),Method="GLM",Type=c(rep("X",100),rep("Xt",100)),
                     Bias=c(rresults$biasX,rresults$biasXt),
                     Se=c(rresults$seX,rresults$seXt),
                     P=c(rresults$PX,rresults$PXt))

plotdata_Obs<-data.frame(rep=rep(1:100,2),Method="INGARCH",Type=c(rep("X",100),rep("Xt",100)),
                         Bias=c(rresults_Obs$biasX,rresults_Obs$biasXt),
                         Se=c(rresults_Obs$seX,rresults_Obs$seXt),
                         P=c(rresults_Obs$PX,rresults_Obs$PXt))

plotdata_Par<-data.frame(rep=rep(1:100,2),Method="MSRC",Type=c(rep("X",100),rep("Xt",100)),
                         Bias=c(rresults_Par$biasX,rresults_Par$biasXt),
                         Se=c(rresults_Par$seX,rresults_Par$seXt),
                         P=c(rresults_Par$PX,rresults_Par$PXt))


large.ses<-unique(plotdata[plotdata$Se>100,"rep"])
large.ses_Obs<-unique(plotdata_Obs[plotdata_Obs$Se>100,"rep"])
large.ses_Par<-unique(plotdata_Par[plotdata_Par$Se>100,"rep"])


# Scatterplot of SE against estimate
ggplot(plotdata,aes(x=Bias,y=Se)) + 
  geom_point(size = 2) + 
  facet_wrap(~Type, ncol = 2) + 
  labs(x="Point estimate", y="Standard error estimate")

ggplot(plotdata_Obs,aes(x=Bias,y=Se)) + 
  geom_point(size = 2) + 
  facet_wrap(~Type, ncol = 2) + 
  labs(x="Point estimate", y="Standard error estimate")

ggplot(plotdata_Par,aes(x=Bias,y=Se)) + 
  geom_point(size = 2) + 
  facet_wrap(~Type, ncol = 2) + 
  labs(x="Point estimate", y="Standard error estimate")
