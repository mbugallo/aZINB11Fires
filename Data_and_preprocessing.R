###############################################################################
###############################################################################
###       Zero-inflated negative binomial mixed models for predicting number 
###		  of wildfires
###
###       Created on Tuesday April 19 2022
###       Directory  /Users/mariabugalloporto/Desktop/Code_Folder (Local Directory)                                                 
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(readxl)
library(mefa)
library(dplyr)

########################################
######## PART 1 : Construct the dataset
########################################

# Target variable
n_incendiosi=as.data.frame(read_excel('fires_year_month_province.xlsx'))
	## Monthly and provincial data (50, excluding Ceuta and Melilla), from 01/1999 to 12/2015
	
# Auxiliary variables
cov_AEMETi=read.csv("AEMET_year_month_province.csv",header=T, sep=';') 
	## Monthly and provincial data (52, including Ceuta and Melilla), from 01/2000 to 12/2015
	
	# Unemployment Rate (no. of unemployed/active population)*100
economia=read.csv("unrate_year_month_province.csv")
	## Monthly and provincial data (52, including Ceuta and Melilla), from 01/2002 to 12/2015

## Initial pre-processing
# 1. By removing Ceuta and Melilla from 'cov_AEMETi', the province ids of both data.frames matches
cov_AEMET=cov_AEMETi[cov_AEMETi$idprovincia!=51,]
cov_AEMET=cov_AEMET[cov_AEMET$idprovincia!=52,] 

# 2. Delete the province column in 'cov_AEMET' (we already have its id)
cov_AEMET=cov_AEMET[, -4] 

# 3. Removing 1999, the years of 'n_incendiosi' and 'cov_AEMETi' match
n_incendios=n_incendiosi[n_incendiosi$anho!=1999, ]

# 4. Create a joint data.frame
datosii=merge(n_incendios, cov_AEMET, by = c("anho", "mes", "idprovincia"))
datosi=merge(datosii, economia, by = c("anho", "mes", "idprovincia"))
datos=datosi[order(datosi$anho, datosi$mes, datosi$idprovincia),]

write.csv(datos,"MODEL_DATA/raw_data.csv", row.names = FALSE)


########################################
######## PART 2 : Data Analysis
########################################

rm(list=ls())
datos=read.csv("MODEL_DATA/raw_data.csv")

# There are many variables for which most of the data are -9999 or NA
# These cases were analysed manually and the variables in question will not be incorporated

# i. When there are' many' -9999 or NA (more than 250), these variables are deleted
datosv2 = datos %>% select(-c(evap, glo, inso, n_cub, n_des, n_nub, nv_0050, nv_0100, nv_1000,
		p_sol, ts_10, ts_20, ts_50, w_rec))

# ii. When there are 'few' -9999 or NA, they are substituted by the median value of the remaining values
datosv2$e[datosv2$e==-9999]=median(datosv2$e[datosv2$e!=-9999])
datosv2$hr[datosv2$hr==-9999]=median(datosv2$hr[datosv2$hr!=-9999])

# Days of...
datosv2$n_fog[datosv2$n_fog==-9999]=median(datosv2$n_fog[datosv2$n_fog!=-9999])
datosv2$n_gra[datosv2$n_gra==-9999]=median(datosv2$n_gra[datosv2$n_gra!=-9999])
datosv2$n_llu[datosv2$n_llu==-9999]=median(datosv2$n_llu[datosv2$n_llu!=-9999])
datosv2$n_nie[datosv2$n_nie==-9999]=median(datosv2$n_nie[datosv2$n_nie!=-9999])
datosv2$np_001[datosv2$np_001==-9999]=median(datosv2$np_001[datosv2$np_001!=-9999])
datosv2$np_010[datosv2$np_010==-9999]=median(datosv2$np_010[datosv2$np_010!=-9999])
datosv2$np_100[datosv2$np_100==-9999]=median(datosv2$np_100[datosv2$np_100!=-9999])
datosv2$np_300[datosv2$np_300==-9999]=median(datosv2$np_300[datosv2$np_300!=-9999])
datosv2$nt_00[datosv2$nt_00==-9999]=median(datosv2$nt_00[datosv2$nt_00!=-9999])
datosv2$nt_30[datosv2$nt_30==-9999]=median(datosv2$nt_30[datosv2$nt_30!=-9999])
datosv2$n_tor[datosv2$n_tor==-9999]=median(datosv2$n_tor[datosv2$n_tor!=-9999])
datosv2$nw_55[datosv2$nw_55==-9999]=median(datosv2$nw_55[datosv2$nw_55!=-9999])
datosv2$nw_91[datosv2$nw_91==-9999]=median(datosv2$nw_91[datosv2$nw_91!=-9999])

# Converted into percentage (rate) of days of...
datosv2$n_fog=100/30*datosv2$n_fog
datosv2$n_gra=100/30*datosv2$n_gra
datosv2$n_llu=100/30*datosv2$n_llu
datosv2$n_nie=100/30*datosv2$n_nie
datosv2$np_001=100/30*datosv2$np_001
datosv2$np_010=100/30*datosv2$np_010
datosv2$np_100=100/30*datosv2$np_100
datosv2$np_300=100/30*datosv2$np_300
datosv2$nt_00=100/30*datosv2$nt_00
datosv2$nt_30=100/30*datosv2$nt_30
datosv2$n_tor=100/30*datosv2$n_tor
datosv2$nw_55=100/30*datosv2$nw_55
datosv2$nw_91=100/30*datosv2$nw_91

datosv2$p_max[datosv2$p_max==-9999]=median(datosv2$p_max[datosv2$p_max!=-9999])
datosv2$p_mes[datosv2$p_mes==-9999]=median(datosv2$p_mes[datosv2$p_mes!=-9999])

# Very large scale (hPa pressure == 100Pa)
datosv2$q_mar[datosv2$q_mar==-9999]=median(datosv2$q_mar[datosv2$q_mar!=-9999])
datosv2$q_max[datosv2$q_max==-9999]=median(datosv2$q_max[datosv2$q_max!=-9999])
datosv2$q_med[datosv2$q_med==-9999]=median(datosv2$q_med[datosv2$q_med!=-9999])
datosv2$q_min[datosv2$q_min==-9999]=median(datosv2$q_min[datosv2$q_min!=-9999])

# We divide by 10 and work in larger units (KPa == 1000Pa)
datosv2$q_mar=datosv2$q_mar/10
datosv2$q_max=datosv2$q_max/10
datosv2$q_med=datosv2$q_med/10
datosv2$q_min=datosv2$q_min/10

datosv2$ta_max[datosv2$ta_max==-9999]=median(datosv2$ta_max[datosv2$ta_max!=-9999])
datosv2$ta_min[datosv2$ta_min==-9999]=median(datosv2$ta_min[datosv2$ta_min!=-9999])
datosv2$ti_max[datosv2$ti_max==-9999]=median(datosv2$ti_max[datosv2$ti_max!=-9999])
datosv2$tm_max[datosv2$tm_max==-9999]=median(datosv2$tm_max[datosv2$tm_max!=-9999])
datosv2$tm_mes[datosv2$tm_mes==-9999]=median(datosv2$tm_mes[datosv2$tm_mes!=-9999])
datosv2$tm_min[datosv2$tm_min==-9999]=median(datosv2$tm_min[datosv2$tm_min!=-9999])
datosv2$ts_min[datosv2$ts_min==-9999]=median(datosv2$ts_min[datosv2$ts_min!=-9999])
datosv2$w_med[datosv2$w_med==-9999]=median(datosv2$w_med[datosv2$w_med!=-9999])

# NOTES:
# Area id: year:month:idprovince (anho:mes:idprovincia)

# As we only have the unemployment rate from 2002 onwards, we will only consider 
# data from that date onwards
datosv2 = datosv2[datosv2$anho>2001,]

write.csv(datosv2,"MODEL_DATA/prepared_data.csv", row.names = FALSE)

