################################################################################################################
################################################################################################################
################################################################################################################

#This script calculates generic early warning signals using the earlywarnings package
#on data from the Hubbard Brook Experimental Forest to examine changes in ecological resilience.

#code developed by Dakos et al. (2015) and adapted by A Contosta
#most recent version 5/16/2022

#changes include code editing and cleaning in preparation for submitting manuscript revisions for peer review

################################################################################################################
################################################################################################################
################################################################################################################


#Initial set-up

#call libraries
library(dplyr)
library(ggplot2)
library(earlywarnings)
library(reshape2)

#set working directory and import files
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Preprocessed Data")

#############################################
#bird abundance, lepidopteran biomass,      #
#microbial biomass,                         #
#and basal area increment                   #
#############################################

#bird abundance. annamre = annual american redstart; annbtbw = annual black-throated blue warbler; 
#annbtgw = annual black-throated green warbler; annlefl = annual least flycatcher;
#annoven = annual ovenbird; annrevi = annual red-eyed vireo
annamre = read.table("annamre.csv", head = T, sep = ",")
annbtbw = read.table("annbtbw.csv", head = T, sep = ",")
annbtgw = read.table("annbtgw.csv", head = T, sep = ",")
annlefl = read.table("annlefl.csv", head = T, sep = ",")
annoven = read.table("annoven.csv", head = T, sep = ",")
annrevi = read.table("annrevi.csv", head = T, sep = ",")

#lepidopteran biomass
lepyear = read.table("lepyear.csv", head = T, sep = ",")

#microbial biomass
micyear = read.table("micyear.csv", head = T, sep = ",")

#basal area increment
baiyear = read.table("tcor.csv", head = T, sep = ",")

##################
#stream chemistry#
##################

#call the concentration files
ws1_volwt =  read.table("ws1_volwt.csv", head = T, sep = ",", na.strings = c("NA", "-888.888"))
ws2_volwt =  read.table("ws2_volwt.csv", head = T, sep = ",", na.strings = c("NA", "-888.888"))
ws4_volwt =  read.table("ws4_volwt.csv", head = T, sep = ",", na.strings = c("NA", "-888.888"))
ws5_volwt =  read.table("ws5_volwt.csv", head = T, sep = ",", na.strings = c("NA", "-888.888"))
ws6_volwt =  read.table("ws6_volwt.csv", head = T, sep = ",", na.strings = c("NA", "-888.888"))

#call the flux files
ws1_flux =  read.table("ws1_flux.csv", head = T, sep = ",")
ws2_flux =  read.table("ws2_flux.csv", head = T, sep = ",")
ws4_flux =  read.table("ws4_flux.csv", head = T, sep = ",")
ws5_flux =  read.table("ws5_flux.csv", head = T, sep = ",")
ws6_flux =  read.table("ws6_flux.csv", head = T, sep = ",")

################################################################################################################
################################################################################################################
################################################################################################################

#####################################
#####################################
#bird abundance                     #
#####################################
#####################################

######
#AMRE#
######

#add column with timeindex for merging with generic_ews output
annamre$timeindex = c(1:nrow(annamre))

#subset to include just year and time index
annamre.sub = annamre[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annamre.pdf")

annamre.ews = generic_ews(annamre$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annamre.ews output with actual years in timeseries
annamre.ewst = merge(annamre.sub, annamre.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annamre.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annamre_EWS.csv",
            col.names = T, row.names = F, sep = ",")

######
#BTBW#
######

#add column with timeindex for merging with generic_ews output
annbtbw$timeindex = c(1:nrow(annbtbw))

#subset to include just year and time index
annbtbw.sub = annbtbw[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annbtbw.pdf")

annbtbw.ews = generic_ews(annbtbw$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annbtbw.ews output with actual years in timeseries
annbtbw.ewst = merge(annbtbw.sub, annbtbw.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annbtbw.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annbtbw_EWS.csv",
            col.names = T, row.names = F, sep = ",")

######
#BTGW#
######

#add column with timeindex for merging with generic_ews output
annbtgw$timeindex = c(1:nrow(annbtgw))

#subset to include just year and time index
annbtgw.sub = annbtgw[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annbtgw.pdf")

annbtgw.ews = generic_ews(annbtgw$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annbtgw.ews output with actual years in timeseries
annbtgw.ewst = merge(annbtgw.sub, annbtgw.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annbtgw.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annbtgw_EWS.csv",
            col.names = T, row.names = F, sep = ",")

######
#LEFL#
######

#add column with timeindex for merging with generic_ews output
annlefl$timeindex = c(1:nrow(annlefl))

#subset to include just year and time index
annlefl.sub = annlefl[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annlefl.pdf")

annlefl.ews = generic_ews(annlefl$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annlefl.ews output with actual years in timeseries
annlefl.ewst = merge(annlefl.sub, annlefl.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annlefl.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annlefl_EWS.csv",
            col.names = T, row.names = F, sep = ",")

######
#OVEN#
######

#add column with timeindex for merging with generic_ews output
annoven$timeindex = c(1:nrow(annoven))

#subset to include just year and time index
annoven.sub = annoven[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annoven.pdf")

annoven.ews = generic_ews(annoven$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annoven.ews output with actual years in timeseries
annoven.ewst = merge(annoven.sub, annoven.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annoven.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annoven_EWS.csv",
            col.names = T, row.names = F, sep = ",")


######
#REVI#
######

#add column with timeindex for merging with generic_ews output
annrevi$timeindex = c(1:nrow(annrevi))

#subset to include just year and time index
annrevi.sub = annrevi[, c("year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_annrevi.pdf")

annrevi.ews = generic_ews(annrevi$totbird, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge annrevi.ews output with actual years in timeseries
annrevi.ewst = merge(annrevi.sub, annrevi.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(annrevi.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\annrevi_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#####################################
#####################################
#lepidopteran biomass               #
#####################################
#####################################

#add column with timeindex for merging with generic_ews output
lepyear$timeindex = c(1:nrow(lepyear))

#subset to include just year and time index
lepyear.sub = lepyear[, c("Year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_lepyear.pdf")

lepyear.ews = generic_ews(lepyear$biomass, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge lepyear.ews output with actual years in timeseries
lepyear.ewst = merge(lepyear.sub, lepyear.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(lepyear.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\lepyear_EWS.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#Microbial Biomass                  #
#####################################
#####################################

#####################
#Microbial Biomass C#
#####################

#add column with timeindex for merging with generic_ews output
micyear$timeindex = c(1:nrow(micyear))

#subset to include just year and time index
BIOC.sub = micyear[, c("Year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_BIOC.pdf")

BIOC.ews = generic_ews(micyear$BIOC, winsize = 50, detrending = "gaussian",
                          bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge micyear.ews output with actual years in timeseries
BIOC.ewst = merge(BIOC.sub, BIOC.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(BIOC.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\BIOC_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#####################
#Microbial Biomass N#
#####################

#add column with timeindex for merging with generic_ews output
micyear$timeindex = c(1:nrow(micyear))

#subset to include just year and time index
BION.sub = micyear[, c("Year", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_BION.pdf")

BION.ews = generic_ews(micyear$BION, winsize = 50, detrending = "gaussian",
                       bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge micyear.ews output with actual years in timeseries
BION.ewst = merge(BION.sub, BION.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(BION.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\BION_EWS.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#basal area increment               #
#####################################
#####################################

######
#ACSA#
######

#add column with timeindex for merging with generic_ews output
baiyear$timeindex = c(1:nrow(baiyear))

#subset to include just year and time index
ACSA.sub = baiyear[, c("YEAR", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ACSA.pdf")

ACSA.ews = generic_ews(baiyear$ACSA_mean, winsize = 50, detrending = "gaussian",
                       bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge baiyear.ews output with actual years in timeseries
ACSA.ewst = merge(ACSA.sub, ACSA.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ACSA.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ACSA_EWS.csv",
            col.names = T, row.names = F, sep = ",")


######
#FAGR#
######

#add column with timeindex for merging with generic_ews output
baiyear$timeindex = c(1:nrow(baiyear))

#subset to include just year and time index
FAGR.sub = baiyear[, c("YEAR", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_FAGR.pdf")

FAGR.ews = generic_ews(baiyear$FAGR_mean, winsize = 50, detrending = "gaussian",
                       bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge baiyear.ews output with actual years in timeseries
FAGR.ewst = merge(FAGR.sub, FAGR.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(FAGR.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\FAGR_EWS.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#stream chemistry: concentration    #
#####################################
#####################################

##################
#Ca concentration#
##################

#########
#Ca ws1#
#########

#add column with timeindex for merging with generic_ews output
ws1_volwt$timeindex = c(1:nrow(ws1_volwt))

#subset to include just year and time index
ws1_volwt_Ca.sub = ws1_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws1_volwt_Ca.pdf")

ws1_volwt_Ca.ews = generic_ews(ws1_volwt$volwt_Ca, winsize = 50, detrending = "gaussian",
                       bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws1_volwt_Ca.ewst = merge(ws1_volwt_Ca.sub, ws1_volwt_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws1_volwt_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws1_volwt_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws2#
#########

#add column with timeindex for merging with generic_ews output
ws2_volwt$timeindex = c(1:nrow(ws2_volwt))

#subset to include just year and time index
ws2_volwt_Ca.sub = ws2_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws2_volwt_Ca.pdf")

ws2_volwt_Ca.ews = generic_ews(ws2_volwt$volwt_Ca, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws2_volwt_Ca.ewst = merge(ws2_volwt_Ca.sub, ws2_volwt_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws2_volwt_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws2_volwt_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws4#
#########

#add column with timeindex for merging with generic_ews output
ws4_volwt$timeindex = c(1:nrow(ws4_volwt))

#subset to include just year and time index
ws4_volwt_Ca.sub = ws4_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws4_volwt_Ca.pdf")

ws4_volwt_Ca.ews = generic_ews(ws4_volwt$volwt_Ca, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws4_volwt_Ca.ewst = merge(ws4_volwt_Ca.sub, ws4_volwt_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws4_volwt_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws4_volwt_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws5#
#########

#add column with timeindex for merging with generic_ews output
ws5_volwt$timeindex = c(1:nrow(ws5_volwt))

#subset to include just year and time index
ws5_volwt_Ca.sub = ws5_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws5_volwt_Ca.pdf")

ws5_volwt_Ca.ews = generic_ews(ws5_volwt$volwt_Ca, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws5_volwt_Ca.ewst = merge(ws5_volwt_Ca.sub, ws5_volwt_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws5_volwt_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws5_volwt_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws6#
#########

#add column with timeindex for merging with generic_ews output
ws6_volwt$timeindex = c(1:nrow(ws6_volwt))

#subset to include just year and time index
ws6_volwt_Ca.sub = ws6_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws6_volwt_Ca.pdf")

ws6_volwt_Ca.ews = generic_ews(ws6_volwt$volwt_Ca, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws6_volwt_Ca.ewst = merge(ws6_volwt_Ca.sub, ws6_volwt_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws6_volwt_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws6_volwt_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")



##################
#NO3 concentration#
##################

#########
#NO3 ws1#
#########

#add column with timeindex for merging with generic_ews output
ws1_volwt$timeindex = c(1:nrow(ws1_volwt))

#subset to include just year and time index
#only include rows where NO3 values are not NA as the generic_ews function will not run with missing values
ws1_volwt_NO3.sub = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3) , c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws1_volwt_NO3.pdf")


ws1_volwt_cc = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]

ws1_volwt_NO3.ews = generic_ews(ws1_volwt_cc$volwt_NO3, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)



dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws1_volwt_NO3.ewst = merge(ws1_volwt_NO3.sub, ws1_volwt_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

ws1_volwt_NO3.ewst$WaterYear = c(1993:2018)

#write table to folder
write.table(ws1_volwt_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws1_volwt_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws2#
#########

#add column with timeindex for merging with generic_ews output
ws2_volwt$timeindex = c(1:nrow(ws2_volwt))

#subset to include just year and time index
ws2_volwt_NO3.sub = ws2_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws2_volwt_NO3.pdf")

ws2_volwt_NO3.ews = generic_ews(ws2_volwt$volwt_NO3, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_2_volwt.ews output with actual years in timeseries
ws2_volwt_NO3.ewst = merge(ws2_volwt_NO3.sub, ws2_volwt_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws2_volwt_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws2_volwt_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws4#
#########

#add column with timeindex for merging with generic_ews output
ws4_volwt$timeindex = c(1:nrow(ws4_volwt))

#subset to include just year and time index
ws4_volwt_NO3.sub = ws4_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws4_volwt_NO3.pdf")

ws4_volwt_NO3.ews = generic_ews(ws4_volwt$volwt_NO3, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws4_volwt_NO3.ewst = merge(ws4_volwt_NO3.sub, ws4_volwt_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws4_volwt_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws4_volwt_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws5#
#########

#add column with timeindex for merging with generic_ews output
ws5_volwt$timeindex = c(1:nrow(ws5_volwt))

#subset to include just year and time index
ws5_volwt_NO3.sub = ws5_volwt[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws5_volwt_NO3.pdf")

ws5_volwt_NO3.ews = generic_ews(ws5_volwt$volwt_NO3, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws5_volwt_NO3.ewst = merge(ws5_volwt_NO3.sub, ws5_volwt_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws5_volwt_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws5_volwt_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws6#
#########

#add column with timeindex for merging with generic_ews output
ws6_volwt$timeindex = c(1:nrow(ws6_volwt))

#subset to include just year and time index
#only include rows where NO3 values are not NA as the generic_ews function will not run with missing values
ws6_volwt_NO3.sub = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3) , c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws6_volwt_NO3.pdf")

ws6_volwt_cc = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]

ws6_volwt_NO3.ews = generic_ews(ws6_volwt_cc$volwt_NO3, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_volwt.ews output with actual years in timeseries
ws6_volwt_NO3.ewst = merge(ws6_volwt_NO3.sub, ws6_volwt_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws6_volwt_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws6_volwt_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")



#####################################
#####################################
#stream chemistry: flux    #
#####################################
#####################################

##################
#Ca flux#
##################

#########
#Ca ws1#
#########

#add column with timeindex for merging with generic_ews output
ws1_flux$timeindex = c(1:nrow(ws1_flux))

#subset to include just year and time index
ws1_flux_Ca.sub = ws1_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws1_flux_Ca.pdf")

ws1_flux_Ca.ews = generic_ews(ws1_flux$Ca_flux, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws1_flux_Ca.ewst = merge(ws1_flux_Ca.sub, ws1_flux_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws1_flux_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws1_flux_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws2#
#########

#add column with timeindex for merging with generic_ews output
ws2_flux$timeindex = c(1:nrow(ws2_flux))

#subset to include just year and time index
ws2_flux_Ca.sub = ws2_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws2_flux_Ca.pdf")

ws2_flux_Ca.ews = generic_ews(ws2_flux$Ca_flux, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws2_flux_Ca.ewst = merge(ws2_flux_Ca.sub, ws2_flux_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws2_flux_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws2_flux_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws4#
#########

#add column with timeindex for merging with generic_ews output
ws4_flux$timeindex = c(1:nrow(ws4_flux))

#subset to include just year and time index
ws4_flux_Ca.sub = ws4_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws4_flux_Ca.pdf")

ws4_flux_Ca.ews = generic_ews(ws4_flux$Ca_flux, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws4_flux_Ca.ewst = merge(ws4_flux_Ca.sub, ws4_flux_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws4_flux_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws4_flux_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws5#
#########

#add column with timeindex for merging with generic_ews output
ws5_flux$timeindex = c(1:nrow(ws5_flux))

#subset to include just year and time index
ws5_flux_Ca.sub = ws5_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws5_flux_Ca.pdf")

ws5_flux_Ca.ews = generic_ews(ws5_flux$Ca_flux, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws5_flux_Ca.ewst = merge(ws5_flux_Ca.sub, ws5_flux_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws5_flux_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws5_flux_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws6#
#########

#add column with timeindex for merging with generic_ews output
ws6_flux$timeindex = c(1:nrow(ws6_flux))

#subset to include just year and time index
ws6_flux_Ca.sub = ws6_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws6_flux_Ca.pdf")

ws6_flux_Ca.ews = generic_ews(ws6_flux$Ca_flux, winsize = 50, detrending = "gaussian",
                               bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws6_flux_Ca.ewst = merge(ws6_flux_Ca.sub, ws6_flux_Ca.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws6_flux_Ca.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws6_flux_Ca_EWS.csv",
            col.names = T, row.names = F, sep = ",")



##################
#NO3 flux#
##################

#########
#NO3 ws1#
#########

#add column with timeindex for merging with generic_ews output
ws1_flux$timeindex = c(1:nrow(ws1_flux))

#subset to include just year and time index
#only include rows where NO3 values are not NA as the generic_ews function will not run with missing values
ws1_flux_NO3.sub = ws1_flux[complete.cases(ws1_flux$NO3_flux) , c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws1_flux_NO3.pdf")


ws1_flux_cc = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]

ws1_flux_NO3.ews = generic_ews(ws1_flux_cc$NO3_flux, winsize = 50, detrending = "gaussian",
                                bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws1_flux_NO3.ewst = merge(ws1_flux_NO3.sub, ws1_flux_NO3.ews, by.x = "timeindex", by.y =  "timeindex")
ws1_flux_Ca.ewst$WaterYear = c(1993:2018)

#write table to folder
write.table(ws1_flux_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws1_flux_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws2#
#########

#add column with timeindex for merging with generic_ews output
ws2_flux$timeindex = c(1:nrow(ws2_flux))

#subset to include just year and time index
ws2_flux_NO3.sub = ws2_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws2_flux_NO3.pdf")

ws2_flux_NO3.ews = generic_ews(ws2_flux$NO3_flux, winsize = 50, detrending = "gaussian",
                                bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_2_flux.ews output with actual years in timeseries
ws2_flux_NO3.ewst = merge(ws2_flux_NO3.sub, ws2_flux_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws2_flux_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws2_flux_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws4#
#########

#add column with timeindex for merging with generic_ews output
ws4_flux$timeindex = c(1:nrow(ws4_flux))

#subset to include just year and time index
ws4_flux_NO3.sub = ws4_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws4_flux_NO3.pdf")

ws4_flux_NO3.ews = generic_ews(ws4_flux$NO3_flux, winsize = 50, detrending = "gaussian",
                                bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws4_flux_NO3.ewst = merge(ws4_flux_NO3.sub, ws4_flux_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws4_flux_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws4_flux_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws5#
#########

#add column with timeindex for merging with generic_ews output
ws5_flux$timeindex = c(1:nrow(ws5_flux))

#subset to include just year and time index
ws5_flux_NO3.sub = ws5_flux[, c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws5_flux_NO3.pdf")

ws5_flux_NO3.ews = generic_ews(ws5_flux$NO3_flux, winsize = 50, detrending = "gaussian",
                                bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws5_flux_NO3.ewst = merge(ws5_flux_NO3.sub, ws5_flux_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws5_flux_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws5_flux_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws6#
#########

#add column with timeindex for merging with generic_ews output
ws6_flux$timeindex = c(1:nrow(ws6_flux))

#subset to include just year and time index
#only include rows where NO3 values are not NA as the generic_ews function will not run with missing values
ws6_flux_NO3.sub = ws6_flux[complete.cases(ws6_flux$NO3_flux) , c("WaterYear", "timeindex")]

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
#commands for dev.copy and dev.off allow for the copying of the 
#plot made in device 4 to the folder

pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Figures\\EWS_ws6_flux_NO3.pdf")

ws6_flux_cc = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]

ws6_flux_NO3.ews = generic_ews(ws6_flux_cc$NO3_flux, winsize = 50, detrending = "gaussian",
                                bandwidth = 50, logtransform = F, interpolate = F)


dev.copy(which = dev.list()["pdf"])
dev.off(which = dev.list()["pdf"])

#merge ws_1_flux.ews output with actual years in timeseries
ws6_flux_NO3.ewst = merge(ws6_flux_NO3.sub, ws6_flux_NO3.ews, by.x = "timeindex", by.y =  "timeindex")

#write table to folder
write.table(ws6_flux_NO3.ewst, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables\\ws6_flux_NO3_EWS.csv",
            col.names = T, row.names = F, sep = ",")

