################################################################################################################
################################################################################################################
################################################################################################################

#This script does initial data processing for early warning signals analysis of changing ecological resilience
#at Hubbard Brook Experimental Forest, New Hampshire, USA.

#Data inputs are from the Original Data folder. Data outputs are preprocessed data that have been reformatted so that 
#all data sources are aggregated at annual timesteps. They are located in the Preprocessed Data folder.

#code developed by A Contosta
#most recent version 5/16/2022

#changes include code editing and cleaning for manuscript currently in press in Environmental Research Letters 

################################################################################################################
################################################################################################################
################################################################################################################

#Initial set-up

#call libraries
library(data.table)
library(earlywarnings)
library(forecast)
library(smoother)
library(stats)
library(stringr)
library(tseries)
library(zoo)

#set working directory and import files
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Original Data")

#########
#climate#
#########

#atemp = daily air temperature, prec = daily total precipitation;

atemp = read.table("HBEF_air_temp_daily_1957-2019.csv",
                   head = T, sep = ",")
prec = read.table("dailyWatershedPrecip1956-2020.csv",
                  head = T, sep = ",")

#########################
#precipitation chemistry#
#########################

#ws6_prcpc = watershed 6 annual volume weighted precipitation chemistry
#ws6_prcpf = watershed 6 annual precipitation solute flux

ws6_prcpc = read.table("ws6_precip_WY_volwt_conc.csv", head = T, sep = ",", 
                      na.strings = c("-888.88", "-888.888"))
ws6_prcpf = read.table("ws6_precip_WY_flux_gHa.csv", head = T, sep = ",", 
                       na.strings = c("-888.88", "-888.888"))

#############################################
#bird abundance, lepidopteran biomass,      #
#microbial biomass,                         #
#and basal area increment                   #
#############################################

#bird = bird abundance, mic = microbial biomass and activity, lep = lepidopteran biomass
#tcor = tree woody increment

bird = read.table("hb_bird.txt", head = T, sep = ",", na.strings = c("t", "tr"))
mic = read.table("micbio.csv", head = T, sep = ",", na.strings = c("-9999.99", "-9999.99.00"))
lep = read.table("lep_per_1000.csv", head = T, sep = ",", na.strings = "NA")
tcor = tcor = read.table("bai_nupert.csv", head = T, sep = ",")

##################
#stream chemistry#
##################

#call the concentration files in the specified directory, add column identifying watershed
#and merge into one file
ws1_volwt = read.table("ws1_stream_WY_volwt_conc.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws1_volwt$ws = 1

ws2_volwt =  read.table("ws2_stream_WY_volwt_conc.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws2_volwt$ws = 2

ws4_volwt =  read.table("ws4_stream_WY_volwt_conc.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws4_volwt$ws = 4

ws5_volwt =  read.table("ws5_stream_WY_volwt_conc.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws5_volwt$ws = 5

ws6_volwt =  read.table("ws6_stream_WY_volwt_conc.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws6_volwt$ws = 6

volwt = data.frame(rbind(ws1_volwt, ws2_volwt, ws4_volwt, ws5_volwt, ws6_volwt))

#call the flux files in the specified directory, add column identifying watershed
#and merge into one file. Remove all values < -888 (na.strings not grabbing all of them)
ws1_flux =  read.table("ws1_stream_WY_flux_gHa.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws1_flux$ws = 1
ws1_flux[1:ncol(ws1_flux)] = lapply(ws1_flux[1:ncol(ws1_flux)], function(x) ifelse(x < -888, NA, x))

ws2_flux =  read.table("ws2_stream_WY_flux_gHa.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws2_flux$ws = 2
ws2_flux[1:ncol(ws2_flux)] = lapply(ws2_flux[1:ncol(ws2_flux)], function(x) ifelse(x < -888, NA, x))

ws4_flux =  read.table("ws4_stream_WY_flux_gHa.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws4_flux$ws = 4
ws4_flux[1:ncol(ws4_flux)] = lapply(ws4_flux[1:ncol(ws4_flux)], function(x) ifelse(x < -888, NA, x))

ws5_flux =  read.table("ws5_stream_WY_flux_gHa.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws5_flux$ws = 5
ws5_flux[1:ncol(ws5_flux)] = lapply(ws5_flux[1:ncol(ws5_flux)], function(x) ifelse(x < -888, NA, x))

ws6_flux =  read.table("ws6_stream_WY_flux_gHa.csv", head = T, sep = ",", na.strings = c("-888.888"))
ws6_flux$ws = 6
ws6_flux[1:ncol(ws6_flux)] = lapply(ws6_flux[1:ncol(ws6_flux)], function(x) ifelse(x < -888, NA, x))

flux = data.frame(rbind(ws1_flux, ws2_flux, ws4_flux, ws5_flux, ws6_flux))

################################################################################################################
################################################################################################################
################################################################################################################

#data pre-processing

#########
#climate#
#########

#####air temperature

#subset atemp to only include weather station #1
atemp.1 = atemp[atemp$STA == "STA1", ]

#convert atemp.1$date to a posix-compliant date object
atemp.1$Date = as.Date(strptime(atemp.1$date, "%Y-%m-%d"))

#create columns for year and month
atemp.1$year = as.numeric(strftime(atemp.1$Date, "%Y"))
atemp.1$month = as.numeric(strftime(atemp.1$Date, "%m"))

#make new column for WaterYear (to align with stream chemistry data)
atemp.1$WaterYear = ifelse(atemp.1$month < 6, atemp.1$year - 1, atemp.1$year)

#calculate annual average temperature 
atemp.2 = data.table(atemp.1)
antemp = atemp.2[ , list(TMEAN = mean(AVE)), by = WaterYear]

#####precipitation

#subset prec to only include watershed 3
prec.1 = prec[prec$watershed == "W3", ]

#convert prec.1$date to a posix-compliant date object
prec.1$Date = as.Date(strptime(prec.1$DATE, "%Y-%m-%d"))

#create columns for year and month
prec.1$year = as.numeric(strftime(prec.1$Date, "%Y"))
prec.1$month = as.numeric(strftime(prec.1$Date, "%m"))

#make new column for WaterYear (to align with stream chemistry data)
prec.1$WaterYear = ifelse(prec.1$month < 6, prec.1$year - 1, prec.1$year)

#calculate annual total precipitation
prec.2 = data.table(prec.1)
anprec = prec.2[ , list(PRCP = sum(Precip)), by = WaterYear]

#####merge atemp and prec datasets

#merge antemp with anprec
anclim = merge(antemp, anprec, by.x = "WaterYear", by.y = "WaterYear")

#write table for downstream analysis (for comparing trends in MAT and total Precip to other variables)
write.table(anclim, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\anclim.csv", sep = ",", row.names = F, col.names = T)

#########################
#precipitation chemistry#
#########################

#calculate acid neutralizing capacity by subtracting the sum of anions (SO4 + NO3 + Cl)
#from the sum of cations (Ca + Mg + Na + K + NH4)

#note that units of concentration (microg/L) have to be converted
#to meq / L as follows: (mass x valence) / molecular weight
#the 1000 is a multiplier for converting mg to micrograms

ws6_prcpc$ANC = (((ws6_prcpc$volwt_Ca * 2) / 40.078) * 1000 + 
                   ((ws6_prcpc$volwt_Mg * 1) / 24.305) * 1000  +   
                   ((ws6_prcpc$volwt_Na * 1) / 22.990) * 1000 +
                   ((ws6_prcpc$volwt_K * 1) / 39.0983)  * 1000  +  
                   ((ws6_prcpc$volwt_NH4 * 1) / 18.038)  * 1000) -
                (((ws6_prcpc$volwt_SO4 * 2) / 96.07 * 1000) +  
                   ((ws6_prcpc$volwt_NO3 * 1) / 62.0049 * 1000) + 
                   ((ws6_prcpc$volwt_Cl * 1) / 35.45 * 1000))

#subset to only include WaterYear, pH, and ANC
ws6_prcpc.1 = ws6_prcpc[ , c("WaterYear", "volwt_pH", "ANC")]

#merge with precipitation solute fluxes (ws6_prcpf)
ws6_prcp = merge(ws6_prcpc.1, ws6_prcpf, by.x = "WaterYear", by.y = "WaterYear")

#Calculate total inorganic N deposition
ws6_prcp$totN = ws6_prcp$NH4_flux + ws6_prcp$NO3_flux

#export ws6_prcp for downstream analysis of how changes in prcp chemistry might be related to other vars
write.table(ws6_prcp, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws6annual.csv", sep = ",", col.names = T, row.names = F)

##################
#stream chemistry#
##################

#omit records that come before experimental watersheds return to baseline. This is based on Fakrei et al. (2020)
#and is cited in Methods. 

#omit records from W1 after 2018 so that all watershed records end at the same time

ws1_volwt = ws1_volwt[ws1_volwt$WaterYear < 2019, ]
ws2_volwt = ws2_volwt[ws2_volwt$WaterYear >= 1972, ]
ws4_volwt = ws4_volwt[ws4_volwt$WaterYear >= 1976, ]
ws5_volwt = ws5_volwt[ws5_volwt$WaterYear >= 1988, ]

ws1_flux = ws1_flux[ws1_flux$WaterYear < 2019, ]
ws2_flux = ws2_flux[ws2_flux$WaterYear >= 1972, ]
ws4_flux = ws4_flux[ws4_flux$WaterYear >= 1976, ]
ws5_flux = ws5_flux[ws5_flux$WaterYear >= 1988, ]

#export for downstream analysis
write.table(ws1_volwt, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws1_volwt.csv", sep = ",", col.names = T, row.names = F)
write.table(ws2_volwt, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws2_volwt.csv", sep = ",", col.names = T, row.names = F)
write.table(ws4_volwt, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws4_volwt.csv", sep = ",", col.names = T, row.names = F)
write.table(ws5_volwt, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws5_volwt.csv", sep = ",", col.names = T, row.names = F)
write.table(ws6_volwt, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws6_volwt.csv", sep = ",", col.names = T, row.names = F)

write.table(ws1_flux, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws1_flux.csv", sep = ",", col.names = T, row.names = F)
write.table(ws2_flux, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws2_flux.csv", sep = ",", col.names = T, row.names = F)
write.table(ws4_flux, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws4_flux.csv", sep = ",", col.names = T, row.names = F)
write.table(ws5_flux, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws5_flux.csv", sep = ",", col.names = T, row.names = F)
write.table(ws6_flux, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\ws6_flux.csv", sep = ",", col.names = T, row.names = F)

################
#bird abundance#
################

#reformat bird file
newbird = melt(bird)

#extract American Redstart
amre = newbird[newbird$Bird.Species == "American Redstart", ]

#extract black-throated blue warbler
btbw = newbird[newbird$Bird.Species == "Bk_thr_Blue Warbler", ]

#extract black-throated green warbler
btgw = newbird[newbird$Bird.Species == "Bk_thr_Green Warbler", ]

#extract least flycatcher
lefl = newbird[newbird$Bird.Species == "Least Flycatcher", ]

#extract ovenbird
oven = newbird[newbird$Bird.Species == "Ovenbird", ]

#extract red-eyed vireo
revi = newbird[newbird$Bird.Species == "Red_eyed Vireo", ]

#calculate total abundance per year 

######
#amre#
######

annamre = aggregate(amre$value, by = list(amre$variable), sum, na.rm = T)
names(annamre) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annamre$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annamre data.frame
annamre$year = sps$year

######
#btbw#
######

annbtbw = aggregate(btbw$value, by = list(btbw$variable), sum, na.rm = T)
names(annbtbw) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annbtbw$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annbtbw data.frame
annbtbw$year = sps$year

######
#btgw#
######

annbtgw = aggregate(btgw$value, by = list(btgw$variable), sum, na.rm = T)
names(annbtgw) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annbtgw$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annbtgw data.frame
annbtgw$year = sps$year

######
#lefl#
######

annlefl = aggregate(lefl$value, by = list(lefl$variable), sum, na.rm = T)
names(annlefl) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annlefl$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annlefl data.frame
annlefl$year = sps$year

######
#oven#
######

annoven = aggregate(oven$value, by = list(oven$variable), sum, na.rm = T)
names(annoven) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annoven$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annoven data.frame
annoven$year = sps$year

######
#revi#
######

annrevi = aggregate(revi$value, by = list(revi$variable), sum, na.rm = T)
names(annrevi) = c("Xyear", "totbird")

#make year column by separating "X" character from year
sps <- data.frame(do.call(rbind, str_split(annrevi$Xyear, "X")))
names(sps) <- c("X", "year")

#make sps$year into integer since the string split call forces it into being a factor class
sps$year = as.integer(sps$year)

#add sps$year column to annrevi data.frame
annrevi$year = sps$year

#export for downstream analysis
write.table(annamre, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annamre.csv", sep = ",", col.names = T, row.names = F)
write.table(annbtbw, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annbtbw.csv", sep = ",", col.names = T, row.names = F)
write.table(annbtgw, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annbtgw.csv", sep = ",", col.names = T, row.names = F)
write.table(annlefl, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annlefl.csv", sep = ",", col.names = T, row.names = F)
write.table(annoven, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annoven.csv", sep = ",", col.names = T, row.names = F)
write.table(annrevi, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\annrevi.csv", sep = ",", col.names = T, row.names = F)

###################
#microbial biomass#
###################

#subset microbial biomass data to include only reference data at low-elevations in organic horizons
micsub = mic[mic$Treatment == "BearBrook" & mic$Se == "SU" & mic$El == "L" & mic$Hor != "Min", ]

#add an ID column for pasting Year and SAM
micsub$ID = paste(micsub$Year, micsub$Plot, sep = "")

#make into data.table for additional data reduction
micsub.1 = data.table(micsub)

#calculate total microbial biomass carbon (BIOC) and nitrogen
#(BION) within each year (across the Oi/Oe and Oa/A)
micsub.2 = micsub.1[ , list(Year = unique(Year),
                            BIOC = sum(BIOC, na.rm = T),
                            BION = sum (BION, na.rm = T)),
                     by = ID]

#calculate annual average microbial biomass carbon (BIOC) and nitrogen (BION)
micyear = micsub.2[ , list(BIOC = mean(BIOC, na.rm = T),
                           BION = mean(BION, na.rm = T)), by = Year]

#export for downstream analysis
write.table(micyear, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\micyear.csv", sep = ",", col.names = T, row.names = F)


######################
#lepidopteran biomass#
######################

#median lepidopteran biomass by year
lepyear = aggregate(lep$BiomassPer1000, by = list(lep$Year), median, na.rm = T)
names(lepyear) = c("Year", "biomass")

#export for downstream analysis
write.table(lepyear, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data\\lepyear.csv", sep = ",", col.names = T, row.names = F)

################################################################################################################
################################################################################################################
################################################################################################################

#use the augmented Dickey-Fuller test to evaluate stationarity of time series
#p-values greater than 0.05 indicates accept null hypothesis of non-stationarity
#and therefore data will have to be detrended prior to calculating EWS

#create empty dataframe for holding results the ADF test for each of the variables
#in the dataset

#first create list of variables in dataset
vars = c("TMEAN", "PRCP", "ws6_ANC", "ws6_totN", "ws6_NO3", "ws6_NH4", "ws6_SO4", "ws6_pH",
         "amre", "btbw", "btgw", "lefl", "oven", "revi", "lep",
         "BIOC", "BION",
         "tcor.ACSA", "tcor.FAGR",
         "ws1.volwt_Ca", "ws2.volwt_Ca", "ws4.volwt_Ca", "ws5.volwt_Ca", "ws6.volwt_Ca",
         "ws1.volwt_NO3", "ws2.volwt_NO3", "ws4.volwt_NO3", "ws5.volwt_NO3", "ws6.volwt_NO3",
         "ws1.flux_Ca", "ws2.flux_Ca", "ws4.flux_Ca", "ws5.flux_Ca", "ws6.flux_Ca",
         "ws1.flux_NO3", "ws2.flux_NO3", "ws4.flux_NO3", "ws5.flux_NO3", "ws6.flux_NO3")


adfvals = data.frame(matrix(vector(), nrow = length(vars), ncol = 2))
names(adfvals) = c("vars", "pval")

adfvals$vars = vars

#####################################
#climate (temperature and precip)
#####################################

adf.TMEAN = adf.test(anclim$TMEAN)
adf.PRCP = adf.test(anclim$PRCP)

adfvals[adfvals$vars == "TMEAN", ]$pval = as.numeric(adf.TMEAN[4])
adfvals[adfvals$vars == "PRCP", ]$pval = as.numeric(adf.PRCP[4])

#####################################################
#precip chemistry (ANC, totN, NO3, NH4, SO4, and pH)#
#####################################################

adf.ws6_ANC = adf.test(ws6_prcp$ANC)
adf.ws6_totN = adf.test(ws6_prcp$totN)
adf.ws6_NO3 = adf.test(ws6_prcp$NO3_flux)
adf.ws6_NH4 = adf.test(ws6_prcp$NH4_flux)
adf.ws6_SO4 = adf.test(ws6_prcp$SO4_flux)
adf.ws6_pH = adf.test(ws6_prcp$volwt_pH)

adfvals[adfvals$vars == "ws6_ANC", ]$pval = as.numeric(adf.ws6_ANC[4])
adfvals[adfvals$vars == "ws6_totN", ]$pval = as.numeric(adf.ws6_totN[4])
adfvals[adfvals$vars == "ws6_NO3", ]$pval = as.numeric(adf.ws6_NO3[4])
adfvals[adfvals$vars == "ws6_NH4", ]$pval = as.numeric(adf.ws6_NH4[4])
adfvals[adfvals$vars == "ws6_SO4", ]$pval = as.numeric(adf.ws6_SO4[4])
adfvals[adfvals$vars == "ws6_pH", ]$pval = as.numeric(adf.ws6_pH[4])

#####################################
#birds (American redstart,
#black-throated blue warbler,
#black-throated green warbler, 
#least flycatcher
#ovenbird, and
#red-eyed vireo#
#####################################

adf.amre = adf.test(annamre$totbird)
adf.btbw = adf.test(annbtbw$totbird)
adf.btgw = adf.test(annbtgw$totbird)
adf.lefl = adf.test(annlefl$totbird)
adf.oven = adf.test(annoven$totbird)
adf.revi = adf.test(annrevi$totbird)

adfvals[adfvals$vars == "amre", ]$pval = as.numeric(adf.amre[4])
adfvals[adfvals$vars == "btbw", ]$pval = as.numeric(adf.btbw[4])
adfvals[adfvals$vars == "btgw", ]$pval = as.numeric(adf.btgw[4])
adfvals[adfvals$vars == "lefl", ]$pval = as.numeric(adf.lefl[4])
adfvals[adfvals$vars == "oven", ]$pval = as.numeric(adf.oven[4])
adfvals[adfvals$vars == "revi", ]$pval = as.numeric(adf.revi[4])

###################
#microbial biomass#
###################

adf.BIOC = adf.test(micyear$BIOC)
adf.BION = adf.test(micyear$BION)

adfvals[adfvals$vars == "BIOC", ]$pval = as.numeric(adf.BIOC[4])
adfvals[adfvals$vars == "BION", ]$pval = as.numeric(adf.BION[4])

######################
#lepidopteran biomass#
######################

adf.lep = adf.test(lepyear$biomass)

adfvals[adfvals$vars == "lep", ]$pval = as.numeric(adf.lep[4])

######################
#basal area increment#
######################

adf.ACSA = adf.test(tcor$ACSA_mean)
adf.FAGR = adf.test(tcor$FAGR_mean)

adfvals[adfvals$vars == "tcor.ACSA", ]$pval = as.numeric(adf.ACSA[4])
adfvals[adfvals$vars == "tcor.FAGR", ]$pval = as.numeric(adf.FAGR[4])

######################
#stream concentration#
######################

#calcium

adf.ws1.volwt_Ca = adf.test(ws1_volwt$volwt_Ca)
adf.ws2.volwt_Ca = adf.test(ws2_volwt$volwt_Ca)
adf.ws4.volwt_Ca = adf.test(ws4_volwt$volwt_Ca)
adf.ws5.volwt_Ca = adf.test(ws5_volwt$volwt_Ca)
adf.ws6.volwt_Ca = adf.test(ws6_volwt$volwt_Ca)

adfvals[adfvals$vars == "ws1.volwt_Ca", ]$pval = as.numeric(adf.ws1.volwt_Ca[4])
adfvals[adfvals$vars == "ws2.volwt_Ca", ]$pval = as.numeric(adf.ws2.volwt_Ca[4])
adfvals[adfvals$vars == "ws4.volwt_Ca", ]$pval = as.numeric(adf.ws4.volwt_Ca[4])
adfvals[adfvals$vars == "ws5.volwt_Ca", ]$pval = as.numeric(adf.ws5.volwt_Ca[4])
adfvals[adfvals$vars == "ws6.volwt_Ca", ]$pval = as.numeric(adf.ws6.volwt_Ca[4])

#nitrogen
#omit rows where NO3 is NA (adf.test will not run with NA values, which occur at the beginning of
#the record)
ws_1_cc = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]
ws_2_cc = ws2_volwt[complete.cases(ws2_volwt$volwt_NO3), ]
ws_4_cc = ws4_volwt[complete.cases(ws4_volwt$volwt_NO3), ]
ws_5_cc = ws5_volwt[complete.cases(ws5_volwt$volwt_NO3), ]
ws_6_cc = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]

adf.ws1.volwt_NO3 = adf.test(ws_1_cc$volwt_NO3)
adf.ws2.volwt_NO3 = adf.test(ws_2_cc$volwt_NO3)
adf.ws4.volwt_NO3 = adf.test(ws_4_cc$volwt_NO3)
adf.ws5.volwt_NO3 = adf.test(ws_5_cc$volwt_NO3)
adf.ws6.volwt_NO3 = adf.test(ws_6_cc$volwt_NO3)

adfvals[adfvals$vars == "ws1.volwt_NO3", ]$pval = as.numeric(adf.ws1.volwt_NO3[4])
adfvals[adfvals$vars == "ws2.volwt_NO3", ]$pval = as.numeric(adf.ws2.volwt_NO3[4])
adfvals[adfvals$vars == "ws4.volwt_NO3", ]$pval = as.numeric(adf.ws4.volwt_NO3[4])
adfvals[adfvals$vars == "ws5.volwt_NO3", ]$pval = as.numeric(adf.ws5.volwt_NO3[4])
adfvals[adfvals$vars == "ws6.volwt_NO3", ]$pval = as.numeric(adf.ws6.volwt_NO3[4])

#############
#stream flux#
#############

#calcium

adf.ws1.flux_Ca = adf.test(ws1_flux$Ca_flux)
adf.ws2.flux_Ca = adf.test(ws2_flux$Ca_flux)
adf.ws4.flux_Ca = adf.test(ws4_flux$Ca_flux)
adf.ws5.flux_Ca = adf.test(ws5_flux$Ca_flux)
adf.ws6.flux_Ca = adf.test(ws6_flux$Ca_flux)

adfvals[adfvals$vars == "ws1.flux_Ca", ]$pval = as.numeric(adf.ws1.flux_Ca[4])
adfvals[adfvals$vars == "ws2.flux_Ca", ]$pval = as.numeric(adf.ws2.flux_Ca[4])
adfvals[adfvals$vars == "ws4.flux_Ca", ]$pval = as.numeric(adf.ws4.flux_Ca[4])
adfvals[adfvals$vars == "ws5.flux_Ca", ]$pval = as.numeric(adf.ws5.flux_Ca[4])
adfvals[adfvals$vars == "ws6.flux_Ca", ]$pval = as.numeric(adf.ws6.flux_Ca[4])

#nitrogen
#omit rows where NO3 is NA (adf.test will not run with NA values, which occur at the beginning of
#the record)
ws_1_cc = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]
ws_2_cc = ws2_flux[complete.cases(ws2_flux$NO3_flux), ]
ws_4_cc = ws4_flux[complete.cases(ws4_flux$NO3_flux), ]
ws_5_cc = ws5_flux[complete.cases(ws5_flux$NO3_flux), ]
ws_6_cc = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]

adf.ws1.flux_NO3 = adf.test(ws_1_cc$NO3_flux)
adf.ws2.flux_NO3 = adf.test(ws_2_cc$NO3_flux)
adf.ws4.flux_NO3 = adf.test(ws_4_cc$NO3_flux)
adf.ws5.flux_NO3 = adf.test(ws_5_cc$NO3_flux)
adf.ws6.flux_NO3 = adf.test(ws_6_cc$NO3_flux)

adfvals[adfvals$vars == "ws1.flux_NO3", ]$pval = as.numeric(adf.ws1.flux_NO3[4])
adfvals[adfvals$vars == "ws2.flux_NO3", ]$pval = as.numeric(adf.ws2.flux_NO3[4])
adfvals[adfvals$vars == "ws4.flux_NO3", ]$pval = as.numeric(adf.ws4.flux_NO3[4])
adfvals[adfvals$vars == "ws5.flux_NO3", ]$pval = as.numeric(adf.ws5.flux_NO3[4])
adfvals[adfvals$vars == "ws6.flux_NO3", ]$pval = as.numeric(adf.ws6.flux_NO3[4])

#note that only six of 39 variables have p-values <0.05. This indicates non-stationarity
#in the time series, and data will have to be detrended

#write adfvals table
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data")
write.table(adfvals, "adfstats.csv", row.names = F, col.names = T, sep = ",")

################################################################################################################
################################################################################################################
################################################################################################################


