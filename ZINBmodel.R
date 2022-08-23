###############################################################################
###############################################################################
###       Zero-inflated negative binomial mixed models for predicting number 
###		  of wildfires
###
###       Created on Wednesday April 20 2022
###       Directory  /Users/mariabugalloporto/Desktop/Code_Folder (Local Directory)                                                 
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(readxl)
library(mefa)
library(glmmTMB)

########################################
######## PART 1 
########################################
# Load data
datos=read.csv("MODEL_DATA/prepared_data.csv")
datos$anho_numeric=datos$anho
datos$anho=factor(datos$anho)
datos$mes=factor(datos$mes)
datos$idprovincia=factor(datos$idprovincia)
datos$nombre=factor(datos$nombre)
attach(datos)

# Fire count by year, month and province
datff=data.frame(aggregate(numero_de_incendios, by=list(anho_numeric, mes, idprovincia), sum))
colnames(datff)=c('anho', 'mes', 'idprovincia', 'recuento')

zeros=datff[datff$recuento==0,]
zeros$cont=1

data.frame(aggregate((numero_de_incendios==0)*1, by=list(idprovincia), sum)) # zeros by province
data.frame(aggregate(zeros$cont, by=list(zeros$mes), sum)) # zeros by month
data.frame(aggregate(zeros$cont, by=list(zeros$anho), sum)) # zeros by year

data.frame(aggregate(numero_de_incendios, by=list(idprovincia), sum)) # wildfires by province
data.frame(aggregate(numero_de_incendios, by=list(mes),sum)) # wildfires by month
data.frame(aggregate(numero_de_incendios, by=list(anho), sum)) # wildfires by year

group_anho3=cut(datos$anho_numeric, breaks=c(2001, 2006, 2012, 2015), labels=c('1', '2', '3'))  

########################################
######## PART 2 : Basic Poisson and Negative Binomial models
########################################

# Let us fit a Poisson model (without random effects) with all available auxiliary variables and 
# recursively eliminate those that do not help to improve the AIC

fit_poisson<-glm(numero_de_incendios~1+e+hr+n_fog+n_gra+n_llu+n_nie+np_001+np_010+np_100+np_300+
	nt_00+nt_30+n_tor+nw_55+nw_91+p_max+p_mes+q_mar+q_max+q_med+q_min+ta_max+ta_min+ti_max+tm_max+
	tm_mes+tm_min+ts_min+w_med+tasa_de_paro+group_anho3,
	data=datos, family=poisson)

library(RcmdrMisc)	
stepwise(fit_poisson, direction='backward', criterion='AIC') # best AIC	435094
# Does not remove anything

library(MASS)		
fit_binomial<-glm.nb(numero_de_incendios~1+e+hr+n_fog+n_gra+n_llu+n_nie+np_001+np_010+np_100+np_300+
	nt_00+nt_30+n_tor+nw_55+nw_91+p_max+p_mes+q_mar+q_max+q_med+q_min+ta_max+ta_min+ti_max+tm_max+
	tm_mes+tm_min+ts_min+w_med+tasa_de_paro+group_anho3, data=datos)

stepwise(fit_binomial, direction='backward', criterion='AIC') # best 65409	
 
fit_binomial1<-glm.nb(numero_de_incendios~1+e+hr+n_fog+n_gra+n_llu+n_nie+np_001+np_010+np_300+nt_00+nt_30+n_tor+
	nw_55+nw_91+q_mar+q_max+q_min+ta_max+ta_min+tm_mes+tm_min+ts_min+w_med+tasa_de_paro+group_anho3,  data=datos)
summary(fit_binomial1)

# Continue with fit_binomial1

########################################
######## PART 3
########################################

# Negative binomial distribution: quadratic parameterization. Variance=mu*(1+mu/phi) = mu+mu^2/phi

# Model fit
# The first step was with hr, np_001, np_010, np_300, nt_00, ta_max, ta_min, tm_mes, tm_min, group_anho in the Bernoulli model 
datos_break=datos[anho_numeric<2015,]
group_anho3_break=group_anho3[anho_numeric<2015]
fit_binomialZIP1prov_mes<-glmmTMB(numero_de_incendios~1+e+hr+n_llu+n_nie+np_300+nt_00
	+nw_55+nw_91+q_mar+q_max+q_min+ta_max+ta_min+tm_mes+tm_min+group_anho3_break+(1|mes+idprovincia), 
	ziformula=~1+hr+np_300+ta_max+group_anho3_break+(1|mes+idprovincia), 
	data=datos_break, family=nbinom2)	# AIC  55796.8

summary(fit_binomialZIP1prov_mes)

round(confint(fit_binomialZIP1prov_mes),4)

# CONDITIONAL MODEL delete: 
	# n_gra, n_fog, np_100, nt_30,	n_tor, p_max, p_mes, q_med, ti_max, tm_max, ts_min, w_med, tasa_de_paro
# Zero inflated MODEL delete:
	# np_001, np_010, nt_00, ta_min, tm_min, tm_mes
	

########################################
######## PART 4 : Check the behaviour of the fitted function
########################################

sigma_11 = sqrt(VarCorr(fit_binomialZIP1prov_mes)[[2]][[1]][[1]]) 
sigma_12 = sqrt(VarCorr(fit_binomialZIP1prov_mes)[[2]][[2]][[1]]) 
sigma_21 = sqrt(VarCorr(fit_binomialZIP1prov_mes)[[1]][[1]][[1]]) 
sigma_22 = sqrt(VarCorr(fit_binomialZIP1prov_mes)[[1]][[2]][[1]]) 

# p_ijk* == predict(fit_binomialZIP1prov_mes, type="zprob")
prob = function(u_11,u_12){	
	exp( fixef(fit_binomialZIP1prov_mes)$zi[1]+fixef(fit_binomialZIP1prov_mes)$zi[2]*hr+
	fixef(fit_binomialZIP1prov_mes)$zi[3]*np_300+fixef(fit_binomialZIP1prov_mes)$zi[4]*ta_max+
	fixef(fit_binomialZIP1prov_mes)$zi[5]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$zi[6]*(group_anho3==3)
	+sigma_11*rep(u_11, 14,each=50)+rep(u_12,14*12)*sigma_12 ) * 
	(1+ exp( fixef(fit_binomialZIP1prov_mes)$zi[1]+fixef(fit_binomialZIP1prov_mes)$zi[2]*hr+
	fixef(fit_binomialZIP1prov_mes)$zi[3]*np_300+fixef(fit_binomialZIP1prov_mes)$zi[4]*ta_max+
	fixef(fit_binomialZIP1prov_mes)$zi[5]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$zi[6]*(group_anho3==3)
	+sigma_11*rep(u_11, 14,each=50)+rep(u_12,14*12)*sigma_12 ) )^(-1)	
	}	

# lambda_ijk* == predict(fit_binomialZIP1prov_mes, type="conditional")
lambda = function(u_21,u_22){ 
	exp(fixef(fit_binomialZIP1prov_mes)$cond[1]+fixef(fit_binomialZIP1prov_mes)$cond[2]*e+
	fixef(fit_binomialZIP1prov_mes)$cond[3]*hr+fixef(fit_binomialZIP1prov_mes)$cond[4]*n_llu+
	fixef(fit_binomialZIP1prov_mes)$cond[5]*n_nie+fixef(fit_binomialZIP1prov_mes)$cond[6]*np_300+
	fixef(fit_binomialZIP1prov_mes)$cond[7]*nt_00+fixef(fit_binomialZIP1prov_mes)$cond[8]*nw_55+
	fixef(fit_binomialZIP1prov_mes)$cond[9]*nw_91+fixef(fit_binomialZIP1prov_mes)$cond[10]*q_mar+
	fixef(fit_binomialZIP1prov_mes)$cond[11]*q_max+fixef(fit_binomialZIP1prov_mes)$cond[12]*q_min+
	fixef(fit_binomialZIP1prov_mes)$cond[13]*ta_max+fixef(fit_binomialZIP1prov_mes)$cond[14]*ta_min+
	fixef(fit_binomialZIP1prov_mes)$cond[15]*tm_mes+fixef(fit_binomialZIP1prov_mes)$cond[16]*tm_min+
	fixef(fit_binomialZIP1prov_mes)$cond[17]*(group_anho3==2)+fixef(fit_binomialZIP1prov_mes)$cond[18]*(group_anho3==3)
	+rep(u_21, 14,each=50)*sigma_21+rep(u_22,14*12)*sigma_22)	
	}
	
u_11 = ranef(fit_binomialZIP1prov_mes)$zi$mes/sigma_11 
u_12 = ranef(fit_binomialZIP1prov_mes)$zi$idprovincia/sigma_12 
u_21= ranef(fit_binomialZIP1prov_mes)$cond$mes/sigma_21 
u_22 = ranef(fit_binomialZIP1prov_mes)$cond$idprovincia/sigma_22

predic_manual = (1-prob(u_11,u_12))*lambda(u_21,u_22)
sum((predic_manual[1:7800,1]-fitted(fit_binomialZIP1prov_mes))^2) 


########################################
######## PART 5 : Residuals
########################################

# Raw residuals
res.bt=(numero_de_incendios[anho_numeric<2015]-fitted(fit_binomialZIP1prov_mes)[anho_numeric<2015])

# Standardized residuals
res.st=(res.bt-mean(res.bt))/sd(res.bt)

# Standardized residuals vs domain indexes
plot(res.st, pch=20, xlab='', ylab='', main='Domains', cex.main=1.95, axes = FALSE, ylim=c(-9,20))
axis(1, cex.axis=1.8); axis(2, cex.axis=1.8)
abline(h=0, lty=2, col="red", lwd=3)

# Standardized residuals vs predicted values
predic=fitted(fit_binomialZIP1prov_mes)
plot(predic,res.st, main='Predicted values', pch=20, las=0, ylim=c(-9,20), 
	cex.main=1.95, axes = FALSE, xlab='', ylab='')
axis(1, cex.axis=1.8); axis(2, cex.axis=1.8)
abline(h=0, lty=2, col="red", lwd=3)

# Standardized residuals vs log predicted values
predic_l=predict(fit_binomialZIP1prov_mes)
plot(predic_l,res.st, xlab='', pch=20, las=0,
	ylab='', main='Log predicted values', cex.main=1.95, axes = FALSE, ylim=c(-9,20))
axis(1, cex.axis=1.8); axis(2, cex.axis=1.8)
abline(h=0, lty=2, col="red")

# Boxplots
par(mfrow=c(1,3))
boxplot(res.st~anho[anho_numeric<2015], xlab='', ylab='', main='Year', cex.main=1.5, axes = FALSE, ylim=c(-9,20))
axis(1, cex.axis=1.5, at=seq(1,13,length=5), labels=c('2002', '2006', '2008', '2012', '2014')); axis(2, cex.axis=1.5)

boxplot(res.st~mes[anho_numeric<2015], xlab='', ylab='', main='Month', cex.main=1.5, axes = FALSE, xlim=c(0,12), ylim=c(-9,20))
axis(1, cex.axis=1.5); axis(2, cex.axis=1.5)

boxplot(res.st~idprovincia[anho_numeric<2015], xlab='', ylab='', main='Province', cex.main=1.5, axes = FALSE, ylim=c(-9,20))
axis(1, cex.axis=1.5); axis(2, cex.axis=1.5)

	
########################################
######## PART 7 : Predictions
########################################

# 1. For a specific province 
j=32 # Ourense	
fitted_prov=fitted(fit_binomialZIP1prov_mes)[seq(j, 50*(12*14), by=50)]
plot(fitted_prov, type='l', ylab='N. de incendios', xlab='Anho', ylim=c(0,1100), axes=F, cex.lab=1.3)
points(seq(1, 12*14, by=12),fitted_prov[seq(1, 12*14, by=12)], col='green', pch=20) 
points(seq(8, 12*14, by=12),fitted_prov[seq(8, 12*14, by=12)], col='red', pch=20) 
points(numero_de_incendios[seq(j, 50*(12*14), by=50)], pch=4, col='blue', type='l')
axis(1, cex.axis=1.3, at=seq(1,12*14+1,length=8), labels=c('2002', '2004', '2006', '2008', '2010', '2012', '2014', '2016'))
axis(2, cex.axis=1.3)

# 2. All predictions for 2015
fitted_futuro=predic_manual[7801:8400,1]
plot(fitted_futuro, type='l', ylab='Wildfires', xlab='Month', ylim=c(0,600), axes=F, cex.lab=1.3, col='blue')
points(numero_de_incendios[7801:8400], type='l')
axis(1, cex.axis=1.2, at=c(seq(1,50*11,length=6), 600), labels=c('January', 'March', 'May', 'July', 'Sept.', 'Nov.', ' '))
axis(2, cex.axis=1.2)
legend(1,550, c('Observed', 'Plug-In'), col=c('black', 'blue'), lty=c(1,1), pt.cex = 1.35, cex=1.3, bty = "n", lwd=c(2,2))


########################################
######## PART 6: 2015 Map information
########################################

fitted_futuro=predic_manual[7801:8400,1] # year 2015
province=1:50

junio2015=fitted_futuro[datos$mes[7801:8400]==6]
julio2015=fitted_futuro[datos$mes[7801:8400]==7]
agosto2015=fitted_futuro[datos$mes[7801:8400]==8]
septiembre2015=fitted_futuro[datos$mes[7801:8400]==9]

used=data.frame(datos$nombre[1:50], province, junio2015, julio2015, agosto2015, septiembre2015)
colnames(used)=c('nombre', 'province', 'junio2015', 'julio2015', 'agosto2015', 'septiembre2015')
write.table(used, "pred_summer2015.txt", sep=";", col.names = TRUE, row.names = FALSE)

write.table(data.frame(fitted_futuro), 'pred2015.txt', sep=";", 
	col.names = TRUE, row.names = FALSE)

# Absolute and relative squared prediction error (RSPE) by month
error_abs_provincia=rep(NA,50)
error_relativo_provincia=rep(NA,50)
for (i in 1:50){
	aux=sum(numero_de_incendios[7801:8400][idprovincia[7801:8400]==i])
	error_abs_provincia[i]=
	sum((numero_de_incendios[7801:8400][idprovincia[7801:8400]==i]-fitted_futuro[idprovincia[7801:8400]==i])^2)
	error_relativo_provincia[i]=sqrt(error_abs_provincia[i])/aux
}

write.table(error_relativo_provincia,"RSPE2015.txt", sep=";",col.names = TRUE, row.names = FALSE)

# Absolute and relative squared prediction error (RSPE) error by province
error_abs_mes=rep(NA,12)
error_relativo_mes=rep(NA,12)
for (j in 1:12){
	aux=sum(numero_de_incendios[7801:8400][mes[7801:8400]==j])
	error_abs_mes[j]=
	sum((numero_de_incendios[7801:8400][mes[7801:8400]==j]-fitted_futuro[mes[7801:8400]==j])^2)
	error_relativo_mes[j]=sqrt(error_abs_mes[j])/aux
}

