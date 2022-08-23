###############################################################################
###############################################################################
###       Zero-inflated negative binomial mixed models for predicting number 
###		  of wildfires
###
###       Created on Wednesday May 18 2022
###       Directory  /Users/mariabugalloporto/Desktop/Code_Folder (Local Directory)                                               
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(glmmTMB)
options(warn=-1)

set.seed(1605)

########################################
######## PART 1 
########################################
# Load data
datos=read.csv("MODEL_DATA/prepared_data.csv")
datos$anho_numeric=datos$anho
datos$anho=factor(datos$anho); datos$mes=factor(datos$mes)
datos$idprovincia=factor(datos$idprovincia); datos$nombre=factor(datos$nombre)
group_anho3=cut(datos$anho_numeric, breaks=c(2001, 2006, 2012, 2015), labels=c('1', '2', '3'))  
attach(datos)

# Fitting data: 2002-2014; Test data: 2015
datos$group_anho3=group_anho3
datos_break=datos[anho_numeric<2015,]
fit_binomialZIP1prov_mes<-glmmTMB(numero_de_incendios~1+e+hr+n_llu+n_nie+np_300+nt_00
	+nw_55+nw_91+q_mar+q_max+q_min+ta_max+ta_min+tm_mes+tm_min+group_anho3+(1|mes+idprovincia), 
	ziformula=~1+hr+np_300+ta_max+group_anho3+(1|mes+idprovincia), 
	data=datos_break, family=nbinom2)

aux2015=data.frame(anho, mes, idprovincia, e, hr, n_llu, n_nie, np_300, nt_00, nw_55, nw_91, q_mar, q_max, 
	q_min, ta_max, ta_min, tm_mes, tm_min, group_anho3)[7801:8400,]	
rownames(aux2015)=seq(1,600,by=1)
	
pred2015=predict(fit_binomialZIP1prov_mes,  newdata=aux2015, type="response")	
	
sigma_11=sqrt(VarCorr(fit_binomialZIP1prov_mes)[[2]][[1]][[1]]) 
sigma_12=sqrt(VarCorr(fit_binomialZIP1prov_mes)[[2]][[2]][[1]]) 
sigma_21=sqrt(VarCorr(fit_binomialZIP1prov_mes)[[1]][[1]][[1]]) 
sigma_22=sqrt(VarCorr(fit_binomialZIP1prov_mes)[[1]][[2]][[1]]) 

# p_ijk* == predict(fit_binomialZIP1prov_mes, type="zprob")
prob=function(u_11,u_12){	exp( fixef(fit_binomialZIP1prov_mes)$zi[1]+fixef(fit_binomialZIP1prov_mes)$zi[2]*hr+
	fixef(fit_binomialZIP1prov_mes)$zi[3]*np_300+fixef(fit_binomialZIP1prov_mes)$zi[4]*ta_max+
	fixef(fit_binomialZIP1prov_mes)$zi[5]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$zi[6]*(group_anho3==3)
	+sigma_11*rep(u_11, 14,each=50)+rep(u_12,14*12)*sigma_12 ) * 
	(1+ exp( fixef(fit_binomialZIP1prov_mes)$zi[1]+fixef(fit_binomialZIP1prov_mes)$zi[2]*hr+
	fixef(fit_binomialZIP1prov_mes)$zi[3]*np_300+fixef(fit_binomialZIP1prov_mes)$zi[4]*ta_max+
	fixef(fit_binomialZIP1prov_mes)$zi[5]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$zi[6]*(group_anho3==3)
	+sigma_11*rep(u_11, 14,each=50)+rep(u_12,14*12)*sigma_12 ) )^(-1)	}	

# lambda_ijk* == predict(fit_binomialZIP1prov_mes, type="conditional")
lambda=function(u_21,u_22){ exp(fixef(fit_binomialZIP1prov_mes)$cond[1]+fixef(fit_binomialZIP1prov_mes)$cond[2]*e+
	fixef(fit_binomialZIP1prov_mes)$cond[3]*hr+fixef(fit_binomialZIP1prov_mes)$cond[4]*n_llu+
	fixef(fit_binomialZIP1prov_mes)$cond[5]*n_nie+fixef(fit_binomialZIP1prov_mes)$cond[6]*np_300+
	fixef(fit_binomialZIP1prov_mes)$cond[7]*nt_00+fixef(fit_binomialZIP1prov_mes)$cond[8]*nw_55+
	fixef(fit_binomialZIP1prov_mes)$cond[9]*nw_91+fixef(fit_binomialZIP1prov_mes)$cond[10]*q_mar+
	fixef(fit_binomialZIP1prov_mes)$cond[11]*q_max+fixef(fit_binomialZIP1prov_mes)$cond[12]*q_min+
	fixef(fit_binomialZIP1prov_mes)$cond[13]*ta_max+fixef(fit_binomialZIP1prov_mes)$cond[14]*ta_min+
	fixef(fit_binomialZIP1prov_mes)$cond[15]*tm_mes+fixef(fit_binomialZIP1prov_mes)$cond[16]*tm_min+
	fixef(fit_binomialZIP1prov_mes)$cond[17]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$cond[18]*(group_anho3==3)
	+rep(u_21, 14,each=50)*sigma_21+rep(u_22,14*12)*sigma_22)	}
	
########################################
######## PART 2: Bootstrap resampling
########################################	

R=500 # number of simulations
g1<-12 # months
g2<-50 # provinces
nobs=7800

y_boot=rep(NA,nobs+600)
nparam=18+6+1+4
est_boot_parametros=matrix(ncol=nparam, nrow=R)
se_boot_parametros=matrix(ncol=nparam, nrow=R)

mu_PIr=matrix(ncol=600, nrow=R)
mu_ijk=matrix(ncol=nobs+600, nrow=R)

y_boot_ijk=matrix(ncol=nobs+600, nrow=R)
var_boot_ijk=rep(NA,600)
IC_percentilbasico=matrix(ncol=600, nrow=R)

t <- proc.time() # time

for(r in 1:R){
	print(r)
	u_11<-rnorm(g1,mean=0,sd=1)  
	u_12<-rnorm(g2,mean=0,sd=1)  
	
	u_21<-rnorm(g1,mean=0,sd=1) 
	u_22<-rnorm(g2,mean=0,sd=1)  
	
	proba=prob(u_11, u_12)
	
	lambdaa=lambda(u_21, u_22)
	
	z<-rbinom(nobs+600,size=1, prob=proba)   # z_ijk* Be(p)
	
	y_boot[z==1]=0	# si z_ijk*=1 -> y_ijk*=0
	y_boot[z==0]=rpois(length(lambdaa[z==0]),lambda=lambdaa[z==0]) # si z_ijk*=0 -> y_ijk*\in BN(lambda_ijk*)
	y_boot_ijk[r,]=y_boot
	
	mu_ijk[r,]=(1-proba)*lambdaa
	
	datos_break$y_boot=y_boot[1:7800]
	
	try(fit_binomialZIP1prov_mes_boot<-glmmTMB(y_boot~1+e+hr+n_llu+n_nie+np_300+nt_00
	+nw_55+nw_91+q_mar+q_max+q_min+ta_max+ta_min+tm_mes+tm_min+group_anho3+(1|mes+idprovincia), 
	ziformula=~1+hr+np_300+ta_max+group_anho3+(1|mes+idprovincia), 
	data=datos_break, family=nbinom2), silent=FALSE)

	est_boot_parametros[r,]=fit_binomialZIP1prov_mes_boot$sdr$par.fixed
	mu_PIr[r,]=predict(fit_binomialZIP1prov_mes_boot, newdata=aux2015, type="response")		
  
	# Basic percentile method
	IC_percentilbasico[r,]=(y_boot_ijk[r,7801:8400]-pred2015)
}		

proc.time()-t # time	

print('###### NEGATIVE BINOMIAL PARAMETERS ###### ')
quantile(est_boot_parametros[,1],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,2],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,3],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,4],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,5],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,6],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,7],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,8],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,9],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,10],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,11],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,12],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,13],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,14],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,15],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,16],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,17],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,18],probs=c(0.025, 0.975))

print('###### BERNOULLI PARAMETERS ###### ')
quantile(est_boot_parametros[,19],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,20],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,21],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,22],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,23],probs=c(0.025, 0.975))
quantile(est_boot_parametros[,24],probs=c(0.025, 0.975))

print('###### NEGATIVE BINOMIAL GAMMA PARAMETER ###### ')
quantile(est_boot_parametros[,25],probs=c(0.025, 0.975))

print('###### NEGATIVE BINOMIAL PARAMETER month ###### ')
quantile(est_boot_parametros[,26],probs=c(0.025, 0.975))

print('###### NEGATIVE BINOMIAL PARAMETER province ###### ')
quantile(est_boot_parametros[,27],probs=c(0.025, 0.975))

print('###### BERNOULLI PARAMETER month ###### ')
quantile(est_boot_parametros[,28],probs=c(0.025, 0.975))

print('###### BERNOULLI PARAMETER province ###### ')
quantile(est_boot_parametros[,29],probs=c(0.025, 0.975))

# RMSE and RRMSE predictions
mu_ijk2015=mu_ijk[,7801:8400]

rmse_muPIr=sqrt(colMeans( (mu_PIr-mu_ijk2015)^2 ,na.rm=T))
rrmse_muPIr=100*rmse_muPIr/pred2015

# Variance Boot
y_boot_ijk2015=y_boot_ijk[, 7801:8400]
y_mean_ijk=t(replicate(expr=colMeans(y_boot_ijk2015, na.rm=T), n=R))

var_boot_ijk=1/(R-1)*R*colMeans( (y_boot_ijk2015-y_mean_ijk)^2 ,na.rm=T)

df=data.frame(rmse_muPIr, rrmse_muPIr, var_boot_ijk)
colnames(df)=c('rmse_muPIr', 'rrmse_muPIr', 'var_boot_ijk')
write.csv(df,"RRMSE_var_Fires500.csv", row.names = FALSE)

# Basic percentile method
IC_05LB=rep(NA, 600)
IC_05UB=rep(NA, 600)
IC_01LB=rep(NA, 600)
IC_01UB=rep(NA, 600)
for (i in 1:600){
  IC_05LB[i]=quantile(IC_percentilbasico[,i], probs=c(0.05/2))
  IC_05UB[i]=quantile(IC_percentilbasico[,i], probs=c(1-0.05/2))
  IC_01LB[i]=quantile(IC_percentilbasico[,i], probs=c(0.01/2))
  IC_01UB[i]=quantile(IC_percentilbasico[,i], probs=c(1-0.01/2))
}

dff=data.frame(IC_05LB ,IC_05UB, IC_01LB, IC_01UB)
colnames(dff)=c('IC_05LB', 'IC_05UB', 'IC_01LB', 'IC_01uB')
write.csv(dff,"ICperbasico_Fires500.csv", row.names = FALSE)

