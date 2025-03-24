####################################################################################################
### example.R: Sample size estimation of the impact of the EPI on hepatitis A incidence          ###
####################################################################################################
# Load relevant libraries
library(openxlsx)
library(tsModel)
library(tscount)
library(lubridate)


#-----------------------------------------------------------------------------------------------
#                                    Determination of parameters
#-----------------------------------------------------------------------------------------------
## Load the data
data <- read.xlsx("HAVdata.xlsx")
data$Date<-as.Date(as.numeric(data$Date),origin="1899-12-30")
data$month<-as.numeric(format(data$Date, "%m"))
data$time<-rep(1:180)/180
data$Xt<-data$Intervention*(data$time*180-43)/180 
dim(data)
head(data)
str(data)

## GLM
mod1 <- glm(Frequency ~ time + Intervention + Xt + harmonic(month,1,12), family=poisson, data)
summary(mod1)
pacf(residuals(mod1,type="deviance"),main="Residuals over time")

## INGARCH model
Xtdata <- data.frame(model.matrix(mod1))[,-1]
colnames(Xtdata)[4:5]<-c("sin","cos")
mod2 <- tsglm(data$Frequency, model=list(past_obs=1), xreg=Xtdata, link="log", distr="poisson" )
summary(mod2)
pacf(residuals(mod2,type="pearson"),main="Residuals over time")


#-----------------------------------------------------------------------------------------------
#                                     Sample size estimation
#-----------------------------------------------------------------------------------------------
Obspower<-function(nsmp,pre_nsmp,beta,rho,sim){
  nstart<-50
  X      <- c(rep(0,pre_nsmp), rep(1,nsmp-pre_nsmp))
  time   <- c(c(1:nsmp)/nsmp)
  Xt	   <- X*(time*nsmp-pre_nsmp)/nsmp
  season <- harmonic(as.numeric(format(seq.Date(from=as.Date("2008-07-01")-months(pre_nsmp), 
                           to=as.Date("2008-07-01")+months(nsmp-pre_nsmp-1),by="month"),"%m")),1,12)
  sin    <- season[,1]
  cos    <- season[,2]
  Xtdata <- model.matrix( ~ time + X + Xt + sin + cos)
  mu0	   <- c(rep(beta[1],nstart),apply(Xtdata, 1, function(s){sum(s*beta)}))
  mu1	   <- exp(mu0)
  Xtdata <- data.frame(time,X,Xt,sin=sin,cos=cos)
  
  matrixy_OAR<-NULL;
  bias_INGARCHX<-NULL;bias_INGARCHXt<-NULL
  rho_INGARCH<-NULL
  powerX_INGARCH<-powerXt_INGARCH<-NULL
  modwarn<-NULL
  for(i in 1:sim){
    
    ntotal<-nsmp+nstart
    #=============== generate data ===============#
    ### AR(1)
    y_OAR <- mu_OAR <- rep(NA, ntotal + 1)
    mu_OAR[1]<- exp(0)                                                 
    y_OAR[1] = rpois(1,exp(beta[1]/(1-rho)))
    for (i in 2:(ntotal+1)){
      ylag = rho*log(y_OAR[i-1]+1)
      mu_OAR[i] = mu0[i-1] + ylag
      y_OAR[i] = rpois(1,exp(mu_OAR[i])) 						
    }
    matrixy_OAR<-rbind(matrixy_OAR,y_OAR)
    
    #=============== fit and test ===============#
    mod_INGARCHw<- tryCatch({ model_INGARCH<-tsglm(y_OAR[-(1:51)], model=list(past_obs=1), xreg=Xtdata, 
                                               link="log", distr="poisson" ) },warning = function(w){1})
    if (is.list(mod_INGARCHw)){
      coef_INGARCH<-coef(model_INGARCH)
      vov_INGARCH<-vcov(model_INGARCH)
      se_INGARCH<-sqrt(diag(vov_INGARCH))
      bias_INGARCHX<-cbind(bias_INGARCHX,coef_INGARCH[4])
      bias_INGARCHXt<-cbind(bias_INGARCHXt,coef_INGARCH[5])
      rho_INGARCH<-cbind(rho_INGARCH,coef_INGARCH[2])
      
      # Level change/Trend change
      powerX_INGARCH<-cbind(powerX_INGARCH,ifelse(( 2*(1-pnorm(abs(coef_INGARCH[4]/se_INGARCH[4] )))<0.05 ),1,0)) 
      powerXt_INGARCH<-cbind(powerXt_INGARCH,ifelse(( 2*(1-pnorm(abs(coef_INGARCH[5]/se_INGARCH[5] )))<0.05 ),1,0)) 
      
      modwarn<-cbind(modwarn,0) } else {modwarn<-cbind(modwarn,1)}
  }
  
  #------------------------------Monte Carlo simulation
  #Level change
  powerX_MC<-ifelse(( 2*(1-pnorm(abs(bias_INGARCHX/sd(bias_INGARCHX) )))<0.05 ),1,0)
  
  #Trend change
  powerXt_MC<-ifelse(( 2*(1-pnorm(abs(bias_INGARCHXt/sd(bias_INGARCHXt) )))<0.05 ),1,0)
    
  
  ##===========Presentation of results
  NSum5<-colSums(matrixy_OAR[which(modwarn!=1),])/(sim-sum(modwarn))
  Tvalue1<-beta[3]
  Tvalue2<-beta[4]
  pointest5<-data.frame(N=nsmp,beta_int=beta[1],beta_t=beta[2],beta_X=beta[3],beta_Xt=beta[4],rho=rho,
                        biasX=mean(bias_INGARCHX)-Tvalue1,seX=sd(bias_INGARCHX),MSEX=mean((bias_INGARCHX-Tvalue1)^2),
                        biasXt=mean(bias_INGARCHXt)-Tvalue2,seXt=sd(bias_INGARCHXt),MSEXt=mean((bias_INGARCHXt-Tvalue2)^2),
                        AlphaMCX=sum(powerX_MC)/(sim-sum(modwarn)),
                        AlphaMCXt=sum(powerXt_MC)/(sim-sum(modwarn)),
                        AlphaINGARCHX  =sum(powerX_INGARCH)/(sim-sum(modwarn)),
                        AlphaINGARCHXt =sum(powerXt_INGARCH)/(sim-sum(modwarn)),
                        before_N=mean(NSum5[1:(nsmp/2)]),after_N=mean(NSum5[(nsmp/2+1):(nsmp)]) , modwarn=sum(modwarn))
  
  return(pointest5)
}


#========================================================================================
ObsN<-function(pold,beta,rho,sim,itmax,method,type){
set.seed(202305)

	ITS1<-Obspower(nsmp=pold,pre_nsmp=round(pold*0.25),beta=beta,rho=rho,sim=sim)
	if (method == "INGARCH") {
		if (type == "X") {
			power0 <- ITS1$AlphaINGARCHX
		} else {
			power0 <- ITS1$AlphaINGARCHXt
		}
	} else {
		if (type == "X") {
			power0 <- ITS1$AlphaMCX
		} else {
			power0 <- ITS1$AlphaMCXt
		}
	}

	NiterObs<-function(powerold){
	   iterdObs = data.frame(Method=method,type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],rho=rho,
                               N=pold,Power=powerold,k=1,index=1)
	   k<-2
	   index1<-1
	   while(k<itmax){
	   if(powerold==1) {powerold<-powerold-0.000001}
		pnew<-ceiling( (qnorm(0.025)+qnorm(1-0.8))^2 / (qnorm(0.025)+qnorm(1-powerold))^2 *pold)
		pnew<-ifelse(pnew%%2==0,pnew,pnew+1)
		if (method == "INGARCH") {
			if (type == "X") {
				powernew<-Obspower(nsmp=pnew,pre_nsmp=round(pnew*0.25),beta,rho=rho,sim=sim)$AlphaINGARCHX
			} else {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=round(pnew*0.25),beta,rho=rho,sim=sim)$AlphaINGARCHXt
			}
		} else {
			if (type == "X") {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=round(pnew*0.25),beta,rho=rho,sim=sim)$AlphaMCX
			} else {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=round(pnew*0.25),beta,rho=rho,sim=sim)$AlphaMCXt
			}
		}
		iterdObs = rbind(iterdObs, data.frame(Method=method,type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],rho=rho,
                                                  N=pnew,Power=powernew,k=k,index=index1) )


	   print(iterdObs)
	   print(Sys.time())



		if(is.na(powernew)){index1<-0;break};
		if(abs(pnew-pold)<1 & powernew>0.8){index1<-0;break};
		k<-k+1
		pold<-pnew
		powerold<-powernew
		}
	return(iterdObs)
	}
	NObs<-NiterObs(power0)

	return(NObs)
}


Sys.time()
ObsN(pold=400,beta=c(2.51,0.47,-0.15,-1.40,-0.08,-0.24),rho=0.37,sim=1000,itmax=1000,method="INGARCH",type="X")
Sys.time()
ObsN(pold=400,beta=c(2.51,0.47,-0.15,-1.40,-0.08,-0.24),rho=0.37,sim=1000,itmax=1000,method="INGARCH",type="Xt")
Sys.time()
ObsN(pold=400,beta=c(2.51,0.47,-0.15,-1.40,-0.08,-0.24),rho=0.37,sim=1000,itmax=1000,method="MC",type="X")
Sys.time()
ObsN(pold=400,beta=c(2.51,0.47,-0.15,-1.40,-0.08,-0.24),rho=0.37,sim=1000,itmax=1000,method="MC",type="Xt")
Sys.time()




