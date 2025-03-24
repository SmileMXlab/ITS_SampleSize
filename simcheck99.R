################################################################################################
### simcheck99.R: A complete successful simulation study.                                    ###
###               Codes are presented for three cases: ①no autocorrelation;                 ###
###                                                    ②observation-driven autocorrelation; ###
###                                                    ③parameter-driven autocorrelation    ###
################################################################################################

# Load relevant libraries
library(MASS)
library(Epi)
library(dplyr)
library(tscount)   


#-----------------------------------------------------------------------------------------------
#                                Non-autocorrelation scenario
#-----------------------------------------------------------------------------------------------
# First, we write a function to generate a single data set:

gendata<-function(nsmp,pre_nsmp,beta){

	time   <- c(1:nsmp)/nsmp
	X      <- c(rep(0,pre_nsmp), rep(1,nsmp-pre_nsmp))
	Xt     <- X*(time*nsmp-pre_nsmp)/nsmp	
	Xtdata <- model.matrix( ~ time + X + Xt )
	mu0    <- apply(Xtdata, 1, function(s){sum(s*beta)})
	mu1    <- exp(mu0)

	#=============== Generate count outcome ===============#
	y      <- rpois(nsmp,mu1) 
      data   <- data.frame(y=y,time=time,X=X,Xt=Xt)

	return(data)
}

# Second, we provide a function to analyse observed data.set:

GLMmod<-function(dataframe,greatmod){

	mod_GLMw<- tryCatch({ mod_GLM<-eval(parse(text=paste( "glm(y~",greatmod,",data=dataframe,family=poisson)",sep="")))  },warning = function(w){1})
	if (is.list(mod_GLMw)) {
		dummy<-ci.lin(mod_GLM,subset=c("X"),Exp=T)
      	res <- data.frame(method="GLM",biasX=dummy[1,1],biasXt=dummy[2,1],seX=dummy[1,2],seXt=dummy[2,2],PX=dummy[1,4],PXt=dummy[2,4])
	} else {
		res<-data.frame(method="GLM",biasX=NA,biasXt=NA,seX=NA,seXt=NA,PX=NA,PXt=NA)
	}

	return(res)
}

# Third, we provide a function to calculate power:

GLMpower<-function(nsmp,pre_nsmp,beta,greatmod,sim){

	matrixy<-NULL
	result<-NULL
	for(i in 1:sim){
		dataframe<-gendata(nsmp=nsmp,pre_nsmp=pre_nsmp,beta=beta)
		matrixy<-rbind(matrixy,dataframe$y)
		result<-rbind(result,GLMmod(dataframe=dataframe,greatmod=greatmod))
		}
	warncount<-sum(is.na(result[, 2]))
	NSum<-colSums(matrixy,na.rm = TRUE)/(sim-warncount) 
	pointest<-data.frame(Method="GLM",N=nsmp,beta0=beta[1],beta_t=beta[2],beta_X=beta[3],beta_Xt=beta[4],
		biasX=mean(result$biasX,na.rm = TRUE)-beta[3],
		biasXt=mean(result$biasXt,na.rm = TRUE)-beta[4],
		AlphaGLMX  =sum(result$PX<0.05,na.rm = TRUE)/(sim-warncount),
		AlphaGLMXt =sum(result$PXt<0.05,na.rm = TRUE)/(sim-warncount),
		before_N=mean(NSum[1:(nsmp/2)]),
		after_N=mean(NSum[(nsmp/2+1):(nsmp)]),modwarn=warncount)

	return(pointest)
}

# Final, we write a function to iterate the sample size:

GLMN<-function(pold,beta,sim,itmax,type){

	ITS1<-GLMpower(nsmp=pold,pre_nsmp=pold/2,beta=beta,greatmod="time+X+Xt",sim=sim)
      if (type=="X") {power0<- ITS1$AlphaGLMX} else {power0<- ITS1$AlphaGLMXt}

	NiterGLM<-function(powerold){
	   iterdGLM = data.frame(Method="GLM",type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],
                               N=pold,Power=powerold,k=1,index=1)
	   k<-2
	   index1<-1
	   while(k<itmax){
	   if(powerold==1) {powerold<-powerold-0.000001}
		pnew<-ceiling( (qnorm(0.025)+qnorm(1-0.8))^2 / (qnorm(0.025)+qnorm(1-powerold))^2 *pold)
		pnew<-ifelse(pnew%%2==0,pnew,pnew+1)
		if (type=="X"){
		   powernew<-GLMpower(nsmp=pnew,pre_nsmp=pnew/2,beta,greatmod="time+X+Xt",sim=sim)$AlphaGLMX
		} else {
		   powernew<-GLMpower(nsmp=pnew,pre_nsmp=pnew/2,beta,greatmod="time+X+Xt",sim=sim)$AlphaGLMXt
		}
		iterdGLM = rbind(iterdGLM, data.frame(Method="GLM",type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],
                                                  N=pnew,Power=powernew,k=k,index=index1) )
		if(is.na(powernew)){index1<-0;break};
		if(abs(pnew-pold)<1 & powernew>0.8){index1<-0;break};
		k<-k+1
		pold<-pnew
		powerold<-powernew
		}
	return(iterdGLM)
	}
	NGLM<-NiterGLM(power0)

	return(NGLM)
}






#-----------------------------------------------------------------------------------------------
#                              Observation-driven autocorrelation case
#-----------------------------------------------------------------------------------------------
# First, we write a function to generate a single data set:

gendata_Obs<-function(nsmp,pre_nsmp,beta,rho){

	nstart <- 50
	time   <- c(1:nsmp)/nsmp
	X      <- c(rep(0,pre_nsmp), rep(1,nsmp-pre_nsmp))
	Xt     <- X*(time*nsmp-pre_nsmp)/nsmp	
	Xtdata <- model.matrix( ~ time + X + Xt )
	mu0    <- c(rep(beta[1],nstart),apply(Xtdata, 1, function(s){sum(s*beta)}))
	mu1    <- exp(mu0)

	#=============== Generate count outcome ===============#
	ntotal    <- nsmp+nstart
	y_OAR     <- mu_OAR <- rep(NA, ntotal + 1)
	mu_OAR[1] <- exp(0)                                              
	y_OAR[1]  <- rpois(1,exp(beta[1]/(1-rho)))
	for (i in 2:(ntotal+1)){
		ylag = rho*log(y_OAR[i-1]+1)
		mu_OAR[i] = mu0[i-1] + ylag
		y_OAR[i] = rpois(1,exp(mu_OAR[i])) 					
	}
      data   <- data.frame(y=y_OAR[-(1:51)],time=time,X=X,Xt=Xt)

	return(data)
}

# Second, we provide a function to analyse observed data.set:

Obsmod<-function(dataframe){

	y <- dataframe[,1]
	Xtdata <- dataframe[,-1]
	mod_INGARCHw<- tryCatch({ mod_INGARCH<-tsglm(y, model=list(past_obs=1), xreg=Xtdata, 
                                                 link="log", distr="poisson" ) },warning = function(w){1})
	if (is.list(mod_INGARCHw)){
	  vov_INGARCHe<-tryCatch({ dummy<- ci.lin(mod_INGARCH,subset=c("X"),Exp=T) },error = function(e){1}) 
	  if(is.matrix(vov_INGARCHe)){
      		res <- data.frame(method="INGARCH",biasX=dummy[1,1],biasXt=dummy[2,1],seX=dummy[1,2],seXt=dummy[2,2],PX=dummy[1,4],PXt=dummy[2,4])
		} else {
			res<-data.frame(method="INGARCH",biasX=NA,biasXt=NA,seX=NA,seXt=NA,PX=NA,PXt=NA)
		}
	} else {
		res<-data.frame(method="INGARCH",biasX=NA,biasXt=NA,seX=NA,seXt=NA,PX=NA,PXt=NA)
	}

	return(res)
}

# Third, we provide a function to calculate power:

Obspower<-function(nsmp,pre_nsmp,beta,rho,sim){

	matrixy<-NULL
	result<-NULL
	for(i in 1:sim){
		dataframe<-gendata_Obs(nsmp=nsmp,pre_nsmp=pre_nsmp,beta=beta,rho=rho)
		matrixy<-rbind(matrixy,dataframe$y)
		result<-rbind(result,Obsmod(dataframe=dataframe))
		}
	warncount<-sum(is.na(result[, 2]))
	NSum<-colSums(matrixy,na.rm = TRUE)/(sim-warncount)  
	pointest<-data.frame(Method="INGARCH",N=nsmp,beta0=beta[1],beta_t=beta[2],beta_X=beta[3],beta_Xt=beta[4],
		biasX=mean(result$biasX,na.rm = TRUE)-beta[3],
		biasXt=mean(result$biasXt,na.rm = TRUE)-beta[4],
		AlphaINGARCHX  =sum(result$PX<0.05,na.rm = TRUE)/(sim-warncount) ,
		AlphaINGARCHXt =sum(result$PXt<0.05,na.rm = TRUE)/(sim-warncount) ,
		AlphaMCX   =sum( ifelse(( 2*(1-pnorm(abs(result$biasX/sd(result$biasX,na.rm = TRUE) )))<0.05 ),1,0),na.rm = TRUE )/(sim-warncount) ,
		AlphaMCXt  =sum( ifelse(( 2*(1-pnorm(abs(result$biasXt/sd(result$biasXt,na.rm = TRUE) )))<0.05 ),1,0),na.rm = TRUE )/(sim-warncount) ,
		before_N=mean(NSum[1:(nsmp/2)]),
		after_N=mean(NSum[(nsmp/2+1):(nsmp)]),modwarn=warncount)

	return(pointest)
}

# Final, we write a function to iterate the sample size:

ObsN<-function(pold,beta,rho,sim,itmax,method,type){

	ITS1<-Obspower(nsmp=pold,pre_nsmp=pold/2,beta=beta,rho=rho,sim=sim)
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
				powernew<-Obspower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sim=sim)$AlphaINGARCHX
			} else {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sim=sim)$AlphaINGARCHXt
			}
		} else {
			if (type == "X") {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sim=sim)$AlphaMCX
			} else {
				powernew <- Obspower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sim=sim)$AlphaMCXt
			}
		}
		iterdObs = rbind(iterdObs, data.frame(Method=method,type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],rho=rho,
                                                  N=pnew,Power=powernew,k=k,index=index1) )
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






#-----------------------------------------------------------------------------------------------
#                                   Parameter-driven autocorrelation scenario
#-----------------------------------------------------------------------------------------------
# First, we write a function to generate a single data set:

gendata_Par<-function(nsmp,pre_nsmp,beta,rho,sigma_mu2){

	time   <- c(1:nsmp)/nsmp
	X      <- c(rep(0,pre_nsmp), rep(1,nsmp-pre_nsmp))
	Xt     <- X*(time*nsmp-pre_nsmp)/nsmp	
	Xtdata <- model.matrix( ~ time + X + Xt )
	mu0    <- apply(Xtdata, 1, function(s){sum(s*beta)})
	mu1    <- exp(mu0)
	sigma_ep2<-sigma_mu2*(1-rho^2)

	#=============== Generate count outcome ===============#
	u5<-rep(NA,nsmp+1)
	y<-rep(NA,nsmp)
	u5[1]<-rnorm(1,-sigma_mu2/2,sqrt(sigma_mu2))
	for (j in 1:nsmp){ 
		u5[j+1]<-rho*(u5[j]+sigma_mu2/2)+rnorm(1,0,sqrt(sigma_ep2))-sigma_mu2/2
		y[j]<-rpois(1,exp(mu0[j]+u5[j+1])) }
      data   <- data.frame(y=y,intercept=rep(1,nsmp),time=time,X=X,Xt=Xt)

	return(data)
}

# Second, we provide functions to analyse observed data.set:

correctcov<-function(data,beta,sigma,rho){
   n<-dim(data)[1]
   m<-dim(data)[2]
   SigmaI<-matrix(0,m,m)
   SigmaII<-matrix(0,m,m)

   for (p in 1:n){ 
      dummy_SigmaI <- (data[p,]%*%t(data[p,]))*c(exp(t(data[p,])%*%beta))
      SigmaI<-SigmaI+dummy_SigmaI
      }

   for (p in 1:n){ 
      for (q in 1:n){
         dummy_SigmaII <- (data[p,]%*%t(data[q,]))*c(exp((t(data[p,])+t(data[q,]))%*%beta))*(sigma*(rho^abs(p-q)))
         SigmaII<-SigmaII+dummy_SigmaII
         }}

   varbeta<-solve(SigmaI)+solve(SigmaI)%*%SigmaII%*%solve(SigmaI)
   sebeta<-sqrt(diag(varbeta))
   list(varbeta=varbeta,sebeta=sebeta)
}

rhofun<-function(yy,Xdata,beta){
	nsmp<-length(yy)
	GEEmu<-exp(Xdata%*%beta)
	GEEsigma2<-sum((yy-GEEmu)^2-GEEmu)/sum(GEEmu^2)
	if(GEEsigma2<0) {GEEsigma2<-1e-4}

	num<-den<-rep(NA,nsmp-1)
	for(m in 1:(nsmp-1)){
		num[m]<-(yy[m]-GEEmu[m])*(yy[m+1]-GEEmu[m+1])
		den[m]<-(GEEmu[m])*(GEEmu[m+1])
					}
	GEErho1<-sum(num,na.rm = T)/(GEEsigma2*sum(den))
	if(GEErho1>1) {GEErho1<-0.99}
	if(GEErho1<0) {GEErho1<-0}
	list(GEEsigma=GEEsigma2,GEErho=GEErho1)
}

correctcovUB<-function(data,beta,sigma,rho,L=15){
   n<-dim(data)[1]
   m<-dim(data)[2]
   SigmaI<-matrix(0,m,m)
   SigmaII<-matrix(0,m,m)
   mu<-exp(data%*%beta)

   for (p in 1:n){ 
      dummy_SigmaI <- (data[p,]%*%t(data[p,]))*mu[p]
      SigmaI<-SigmaI+dummy_SigmaI
      }

   for (p in (-L):L){ 
      for (q in max(1-p,1):min(n-p,n)){
         dummy_SigmaII <- (data[q,]%*%t(data[q+p,]))*mu[q]*mu[q+p]*(sigma*(rho^abs(p)))
         SigmaII<-SigmaII+dummy_SigmaII
         }}

   G<-solve(SigmaI)+solve(SigmaI)%*%SigmaII%*%solve(SigmaI)
   return(G)
}

rhofunUB<-function(yy,Xdata,beta,sigmapre,rhopre){
	nsmp<-length(yy)

	## σUB
	GEEmu<-exp(Xdata%*%beta)
	G<-correctcovUB(Xdata,beta=beta,sigmapre,rhopre,L=15)
	dummy1<-exp(-2*Xdata%*%G%*%t(Xdata))
	dummy2<-exp(2*Xdata%*%G%*%t(Xdata))
	dummy3<-exp(Xdata%*%G%*%t(Xdata)/2)
	sigmaUB<-sum((yy-GEEmu)^2+GEEmu^2*diag(dummy1)*diag(dummy2-2*dummy3+1)-GEEmu)/sum(GEEmu^2*diag(dummy1))
	if(is.na(sigmaUB) | sigmaUB<0) {sigmaUB<-1e-4}

	## ρUB & ρMSRC
	rhoUBvector<-rhoUBvectorC<-VrhoUBvector<-Prhovector<-NULL;rankUB<-1
	while(rankUB>0 & rankUB<6){
	g<-gamaUB1<-gamaUB2<-VrhoUB1<-VrhoUB2<-rep(NA,nsmp-rankUB)
	for (k in 1:(nsmp-rankUB)){
		g[k]<-exp(-(Xdata[k,]+Xdata[k+rankUB,])%*%G%*%(Xdata[k,]+Xdata[k+rankUB,])/2)
		gamaUB1[k]<-GEEmu[k]*GEEmu[k+rankUB]*g[k]
		gamaUB2[k]<-(yy[k]-GEEmu[k])*(yy[k+rankUB]-GEEmu[k+rankUB])+GEEmu[k]*GEEmu[k+rankUB]*g[k]*(1-diag(dummy3)[k]-diag(dummy3)[k+rankUB]+1/g[k])
		VrhoUB1[k]<-GEEmu[k]*GEEmu[k+rankUB]
		VrhoUB2[k]<-GEEmu[k]^2*GEEmu[k+rankUB]^2*(1+1/GEEmu[k]/sigmaUB)*(1+1/GEEmu[k+rankUB]/sigmaUB)
	}
	gamaUB<-sum(gamaUB2)/sum(gamaUB1)
	rhoUB<-gamaUB/sigmaUB
	if(is.na(rhoUB)) {rhoUB<-0}
	if(rhoUB>1) {rhoUB<-0.99}
	rhoUBvector<-rbind(rhoUBvector,rhoUB)

	VrhoUB<-sum(VrhoUB2)/(sum(VrhoUB1)^2)
	VrhoUBvector<-rbind(VrhoUBvector,VrhoUB)
	Prho<-ifelse(( 2*(1-pnorm( abs(rhoUB/sqrt(VrhoUB)) ))<0.01 ),1,0)  
	Prhovector<-rbind(Prhovector,Prho)
	if(rhoUB<0) {rhoUBvectorC<-rbind(rhoUBvectorC,NA)} else {rhoUBvectorC<-rbind(rhoUBvectorC,(rhoUB)^(1/rankUB))}
	rankUB=rankUB+1
	}

	if(rhoUBvector[1]<0) {rhoUBvector[1]<-0}
	rhoUBvectorC[which(is.na(rhoUBvectorC))]<-0
	if (sum(Prhovector==1)==0) {order=1} else {order=c(which(Prhovector==1))}

	list(GEEsigmaUB=sigmaUB,GEErhoUB=rhoUBvector[1],GEErhomax=max(rhoUBvectorC[order],na.rm = T),order=mean(order) )
}

Parmod<-function(dataframe,greatmod){

	mod_GLMw<- tryCatch({ mod_GLM<-eval(parse(text=paste( "glm(y~",greatmod,",data=dataframe,family=poisson)",sep=""))) },warning = function(w){1})
	if (is.list(mod_GLMw)) {
		Xmatrix <- as.matrix(dataframe[,-1])
		dummy<-rhofun(dataframe$y,Xmatrix,coef(mod_GLM))
		dummyUB<-rhofunUB(dataframe$y,Xmatrix,coef(mod_GLM),dummy$GEEsigma,dummy$GEErho)
		vovmax<-correctcov(Xmatrix,coef(mod_GLM),dummyUB$GEEsigmaUB,dummyUB$GEErhomax)
		pvalue<-2*(1-pnorm(abs(coef(mod_GLM)/vovmax$sebeta)))
      	res <- data.frame(method="MSRC",biasX=coef(mod_GLM)[3],biasXt=coef(mod_GLM)[4],seX=vovmax$sebeta[3],seXt=vovmax$sebeta[4],PX=pvalue[3],PXt=pvalue[4])
	} else {
		res<-data.frame(method="MSRC",biasX=NA,biasXt=NA,seX=NA,seXt=NA,PX=NA,PXt=NA)
	}

	return(res)
}

# Third, we provide a function to calculate power:

Parpower<-function(nsmp,pre_nsmp,beta,rho,sigma_mu2,greatmod,sim){

	matrixy<-NULL
	result<-NULL
	for(i in 1:sim){
		dataframe<-gendata_Par(nsmp=nsmp,pre_nsmp=pre_nsmp,beta=beta,rho=rho,sigma_mu2=sigma_mu2)
		matrixy<-rbind(matrixy,dataframe$y)
		result<-rbind(result,Parmod(dataframe=dataframe,greatmod=greatmod))
		}
	warncount<-sum(is.na(result[, 2]))
	NSum<-colSums(matrixy,na.rm = TRUE)/(sim-warncount) 
	pointest<-data.frame(Method="MSRC",N=nsmp,beta0=beta[1],beta_t=beta[2],beta_X=beta[3],beta_Xt=beta[4],rho=rho,sigma_mu2=sigma_mu2,
		biasX=mean(result$biasX,na.rm = TRUE)-beta[3],
		biasXt=mean(result$biasXt,na.rm = TRUE)-beta[4],
		AlphaMSRCX  =sum(result$PX<0.05,na.rm = TRUE)/(sim-warncount),
		AlphaMSRCXt =sum(result$PXt<0.05,na.rm = TRUE)/(sim-warncount),
		before_N=mean(NSum[1:(nsmp/2)]),
		after_N=mean(NSum[(nsmp/2+1):(nsmp)]),modwarn=warncount)

	return(pointest)
}

# Final, we write a function to iterate the sample size:

ParN<-function(pold,beta,rho,sigma_mu2,sim,itmax,type){

	ITS1<-Parpower(nsmp=pold,pre_nsmp=pold/2,beta=beta,rho=rho,sigma_mu2=sigma_mu2,greatmod="time+X+Xt",sim=sim)
      if (type=="X") {power0<- ITS1$AlphaMSRCX} else {power0<- ITS1$AlphaMSRCXt}

	NiterPar<-function(powerold){
	   iterdPar = data.frame(Method="MSRC",type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],
                               rho=rho,sigma_mu2=sigma_mu2,N=pold,Power=powerold,k=1,index=1)
	   k<-2
	   index1<-1
	   while(k<itmax){
	   if(powerold==1) {powerold<-powerold-0.000001}
		pnew<-ceiling( (qnorm(0.025)+qnorm(1-0.8))^2 / (qnorm(0.025)+qnorm(1-powerold))^2 *pold)
		pnew<-ifelse(pnew%%2==0,pnew,pnew+1)
		if (type=="X"){
		   powernew<-Parpower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sigma_mu2=sigma_mu2,greatmod="time+X+Xt",sim=sim)$AlphaMSRCX
		} else {
		   powernew<-Parpower(nsmp=pnew,pre_nsmp=pnew/2,beta,rho=rho,sigma_mu2=sigma_mu2,greatmod="time+X+Xt",sim=sim)$AlphaMSRCXt
		}
		iterdPar = rbind(iterdPar, data.frame(Method="MSRC",type=type,bet_inter=beta[1],bet_t=beta[2],bet_X=beta[3],bet_Xt=beta[4],
                                                  rho=rho,sigma_mu2=sigma_mu2,N=pnew,Power=powernew,k=k,index=index1) )
		if(is.na(powernew)){index1<-0;break};
		if(abs(pnew-pold)<1 & powernew>0.8){index1<-0;break};
		k<-k+1
		pold<-pnew
		powerold<-powernew
		}
	return(iterdPar)
	}
	NPar<-NiterPar(power0)

	return(NPar)
}






Sys.time()
set.seed(202305)
GLMN(pold=400,beta=c(log(2),log(2),log(1.5),0),sim=1000,itmax=1000,type="X") 
Sys.time()
set.seed(202305)
GLMN(pold=400,beta=c(log(2),log(2),0,log(4)),sim=1000,itmax=1000,type="Xt")  
Sys.time()
set.seed(202305)
ObsN(pold=400,beta=c(log(2),log(2),log(1.5),0),rho=0.6,sim=1000,itmax=1000,method="INGARCH",type="X")  
Sys.time()
set.seed(202305)
ObsN(pold=400,beta=c(log(2),log(2),0,log(4)),rho=0.6,sim=1000,itmax=1000,method="INGARCH",type="Xt")    
Sys.time()
set.seed(202305)
ParN(pold=400,beta=c(log(2),log(2),log(1.5),0),rho=0.6,sigma_mu2=0.5,sim=1000,itmax=1000,type="X")     
Sys.time()
set.seed(202305)
ParN(pold=400,beta=c(log(2),log(2),0,log(4)),rho=0.2,sigma_mu2=0.5,sim=1000,itmax=1000,type="Xt")     
Sys.time()



