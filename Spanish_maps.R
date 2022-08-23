###############################################################################
###############################################################################
###       Zero-inflated negative binomial mixed models for predicting number 
###		  of wildfires
###
###       Created on Tuesday April 30 2022
###       Directory  /Users/mariabugalloporto/Desktop/Code_Folder (Local Directory)                                                 
###                                                     
###       @author: mariabugalloporto
           
rm(list=ls())           
library(maptools)
library(RColorBrewer)
library(mefa)

########################################
######## PART 1 : Main functions
########################################

GroupClassification <- function(data,datacompare,intervals)
{
   n = length(data)
   group = matrix(0,nrow=n) 
   ninterv = length(intervals)
   for (i in 1:n)
   {
      for (j in 1:ninterv)
         if (datacompare[i]<intervals[j])
         {
            group[i]= intervals[j]
            break
         }
   }
   result = split(data,group)
   return (result)
}
 
PrintSpainMap <- function(pathmap,datos,colors,titlemap,textlegend,eliminarprov)
{

   m <- matrix(c(1,1,1,2),2,2)
   layout(m, widths=c(1.5, 1), heights=c(1.5, 1), respect=F)

   xName <- readShapePoly(pathmap, IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))     

   xName$datos <- NA
   for (i in 1:length(colors))
      xName$datos[datos[[i]]] <- colors[i]
 
   xSC <- xName[xName$ESP_PROV_I < 35 | xName$ESP_PROV_I >38 | xName$ESP_PROV_I==36 | xName$ESP_PROV_I ==37,]
   plot(xSC,  xlab="",  col=xSC$datos, axes=F)
  title(titlemap, line=-1.2, cex.main=1.7)
  legend( "topright", textlegend, col = 1, pt.bg = colors, pch=21, bty = "n", cex=1.2, pt.cex = 2.5)  #  cex=1.3

   xC <- xName[xName$ESP_PROV_I == 35 | xName$ESP_PROV_I == 38,]
   plot(xC,  xlab="",  col=xC$datos)
   box()
}

########################################
######## PART 2 : Paper maps
########################################
pathmap    <- "SPAIN_MAP/esp_prov.shp"

##############################################################
# Observed values
datos=read.csv("MODEL_DATA/prepared_data.csv")
numero_de_incendios=datos$numero_de_incendios
dom=datos$idprovincia[1:50]
estML=numero_de_incendios[datos$anho==2015 & datos$mes==7]

# Intervals  
intervals_prop <- c(0,12,50,100,180, 230, Inf)  

# Colors
colorsprop <- c("white", brewer.pal(8,"YlOrRd")[c(2,4,5,6,8)] )
 
# Legend
legend_prop <- expression("<= 12", ">12 & <= 50", ">50 & <=100", ">100 & <=180", ">180 & <=230", '>230') 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"Wildfires July 2015",legend_prop,eliminarprov)


##############################################################
# Preditions
file_est   <- "pred_summer2015.txt"  
c_dom     <- 2     # domain column number
c_estML   <- 6     # column number
datoss <- read.table(file=file_est, header=TRUE, sep=';', dec='.')
nprov <- nrow(datoss)
dom     <- datoss[,c_dom]
estML   <- round(datoss[,c_estML],0) # integers

intervals_prop <- c(0,12,50,100, Inf)  

colorsprop <- c("white", brewer.pal(8,"YlOrRd")[c(2,4,5)] )
 
legend_prop <- expression("<= 12", ">12 & <= 50", ">50 & <=100", ">100 & <=180") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"Forecast September 2015",legend_prop,eliminarprov)


##############################################################
# RRMSE
file_est   <- "RRMSE_var_Fires500.csv" 
datosss <- read.csv(file=file_est, header=TRUE, sep=',', dec='.')$rrmse_muPIr
messs=rep(seq(1,12), each=50)[1:600]
provinceee=rep(seq(1,50),12)[1:600]
datossss=data.frame(messs, provinceee, datosss)

julio=datossss[messs==7,]
agosto=datossss[messs==8,]
septiembre=datossss[messs==9,]
selec=septiembre

###
dom=selec[,2]
estML=selec[,3]/100

intervals_prop = c(0,0.08, 0.16, 0.25, 0.40, 0.50, Inf)

colorsprop <- c("white", brewer.pal(8,"YlOrRd")[c(2,3,4,5,8)] ) 

legend_prop = expression("under 8%", "8 - 16 %", "16 - 25 %", "25 - 40 %", " 40 - 50 %", 'over 50 %')

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"RRMSE September 2015",legend_prop,eliminarprov)


##############################################################
# RSPE
dom=read.csv("MODEL_DATA/prepared_data.csv")$idprovincia[1:50]
file_est   <- "RSPE2015.txt"  
estML   <- read.table(file=file_est, header=TRUE, sep=';', dec='.')[,1]

intervals_prop <- c(0,0.10, 0.15, 0.225, 0.30, 0.4, Inf) 

colorsprop <- c("white", brewer.pal(8,"YlOrRd")[c(2,4,5,6,8)] ) 
 
legend_prop <- expression("under 10%", "10 - 15 %", "15 - 22.5 %", "22.5 - 30 %", " 30 - 40 %", '40 - 45 %')

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"RSPE 2015",legend_prop,eliminarprov)


##############################################################
# Confidence intervals
incendios=read.table("pred_summer2015.txt", sep=';',header=T)
messs=rep(seq(1,12), each=50)[1:600]

real_incendios=read.csv("MODEL_DATA/prepared_data.csv")
real_incendios$anho=factor(real_incendios$anho)
real_incendios2015=real_incendios$numero_de_incendios[real_incendios$anho==2015]
fitted_futuro=read.table("pred2015.txt", sep=';',header=T)

sqrtVar_yijk=sqrt(read.csv('RRMSE_var_Fires500.csv', header=TRUE, sep=',', dec='.')$var_boot_ijk)
estML=matrix(nrow=12, ncol=50)
for (i in 1:12){
	ics=data.frame(1:50, fitted_futuro[messs==i,]+qnorm(0.05/2)*sqrtVar_yijk[messs==i],
	fitted_futuro[messs==i,]-qnorm(0.05/2)*sqrtVar_yijk[messs==i], real_incendios2015[messs==i])	
	colnames(ics)=c('Province', 'Lower', 'Upper', 'Real')
	estML[i,ics$Lower<=ics$Real & ics$Upper>=ics$Real]=1 # Inside
	estML[i,ics$Lower>ics$Real]=0 # Wrong IC Lower Bound
	estML[i,ics$Upper<ics$Real]=0 # Wrong IC Upper Bound
}
estML=colMeans(estML)*100
 
intervals_prop <- c(8.3, 16.6, 58.3,  90, 100, Inf) 

colorsprop <- c(brewer.pal(8,"YlGn")[c(2,3,4,5,7)] ) 
 
legend_prop <- expression("8.3 - 16.6 %", "16.6 - 58.3 %", "58.3 - 83.3 %", '83.3 - 91.6%', '91.6 - 100%')

result = GroupClassification(1:50,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"Coverage 2015",legend_prop,eliminarprov)


##### Average coverage probabilities per month
### RRMSE_var_Fires500 aggregated
estML=matrix(nrow=50, ncol=12)
provinceee=rep(seq(1,50),12)

for (j in 1:50){
	ics=data.frame(1:12, fitted_futuro[provinceee==j,]+qnorm(0.05/2)*sqrtVar_yijk[provinceee==j],
	fitted_futuro[provinceee==j,]-qnorm(0.05/2)*sqrtVar_yijk[provinceee==j], real_incendios2015[provinceee==j])	
	colnames(ics)=c('Month', 'Lower', 'Upper', 'Real')
	estML[j,ics$Lower<=ics$Real & ics$Upper>=ics$Real]=1 # Inside
	estML[j,ics$Lower>ics$Real]=0 # Wrong IC Lower Bound
	estML[j,ics$Upper<ics$Real]=0 # Wrong IC Upper Bound
}
estML=colMeans(estML)*100



