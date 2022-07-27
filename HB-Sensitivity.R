################################################################################################################
################################################################################################################
################################################################################################################

#This script performs sensitivity analysis on EWS, examining Kendall's tau and p-value statistics for change over time
#for different window lengths and smoothing bandwiths using the earlywarnings package.

#Data inputs are from the Preprocessed Data folder. 

#Outputs include diagnostic plots and related tables of Kendall's Tau and p-values that result from 
#using different smoothing bandwiths for detrending data and different window lengths for calculating 
#early warning signals.

#Outputs also include Figure_S4, which shows the distribution of p-values for change over time for each of four 
#early warning signals (standard deviation, autocorrelation, skewness, and kurtosis) over five possible bandwiths 
#for smoothing and detrending the data and over five possible sliding windows for calculating early warning signals 

#Data outputs are in the Sensitivity Analysis Figures and Sensitivity Analysis Tables folders for the results 
#of the sensivitiy analyses and in the EWS Summary Figures for Figure S4.  

#code developed by Dakos et al. (2015) and adapted by A Contosta
#most recent version 2/16/2022

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

#run sensitivity analysis for system state variables. EWS are sd, acf1, sk, kt 
#window size varies from 25% to 75% of the time series (default), 
#with an increment of 5 (incrwinsize)
#bandwith size for Gaussian detrending varies from 5 to 100% (default), 
#with an increment of 5 (incrbandwith)

#output for each state variable and EWS includes contour plots for
#Kendall tau estimates and their associated p-values for the range of both the rolling 
#window sizes and bandwidths used. A triangle indicates which combination of the 
#window size and bandwith results in the highest Kendall tau and lowest p-value

#output also includes histograms with the range of tau and p-values for the ranges of
#window sizes and bandwiths used to calculate EWS

#Kendall's tau statistics for each state variable, EWS, window size, and bandwith 
#are then compiled into a single dataframe for downstream diagnostics and analysis

#####################################
#####################################
#bird abundance                     #
#####################################
#####################################

######
#AMRE#
######

#run sensitivity analyses and write diagnostic plots to folder
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annamre.pdf")

#run analyses
annamre.sd.sa = sensitivity_ews(timeseries = annamre$totbird, indicator = "sd", 
                            detrending = 'gaussian', 
                            incrwinsize = 5, incrbandwidth = 5)

annamre.ac.sa = sensitivity_ews(timeseries = annamre$totbird, indicator = "acf1", 
                            detrending = 'gaussian', 
                            incrwinsize = 5, incrbandwidth = 5)

annamre.sk.sa = sensitivity_ews(timeseries = annamre$totbird, indicator = "sk", 
                            detrending = 'gaussian', 
                            incrwinsize = 5, incrbandwidth = 5)

annamre.kt.sa = sensitivity_ews(timeseries = annamre$totbird, indicator = "kurt", 
                            detrending = 'gaussian', 
                            incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annamre.sd.sa$bandyears = as.numeric(row.names(annamre.sd.sa)) 
annamre.ac.sa$bandyears = as.numeric(row.names(annamre.ac.sa)) 
annamre.sk.sa$bandyears = as.numeric(row.names(annamre.sk.sa)) 
annamre.kt.sa$bandyears = as.numeric(row.names(annamre.kt.sa)) 

#add column to each dataframe with the name of the EWS
annamre.sd.sa$li = "sd"
annamre.ac.sa$li = "ac"
annamre.sk.sa$li = "sk"
annamre.kt.sa$li = "kt"

#combine dataframes
annamre.sa = rbind(annamre.sd.sa, annamre.ac.sa, annamre.sk.sa, annamre.kt.sa)

#add column with state variable name
annamre.sa$var = "annamre"

#melt data frame so that columns with window lengths become rows
annamre.sa.lon = reshape2::melt(annamre.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annamre.sa.lon$bandwith = annamre.sa.lon$bandyears / nrow(annamre)
annamre.sa.lon$window = as.numeric(as.character(annamre.sa.lon$variable)) / nrow(annamre)

#rename columns
names(annamre.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annamre.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annamre_sa.csv",
            col.names = T, row.names = F, sep = ",")



######
#BTBW#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annbtbw.pdf")

#run analyses
annbtbw.sd.sa = sensitivity_ews(timeseries = annbtbw$totbird, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtbw.ac.sa = sensitivity_ews(timeseries = annbtbw$totbird, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtbw.sk.sa = sensitivity_ews(timeseries = annbtbw$totbird, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtbw.kt.sa = sensitivity_ews(timeseries = annbtbw$totbird, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annbtbw.sd.sa$bandyears = as.numeric(row.names(annbtbw.sd.sa)) 
annbtbw.ac.sa$bandyears = as.numeric(row.names(annbtbw.ac.sa)) 
annbtbw.sk.sa$bandyears = as.numeric(row.names(annbtbw.sk.sa)) 
annbtbw.kt.sa$bandyears = as.numeric(row.names(annbtbw.kt.sa)) 

#add column to each dataframe with the name of the EWS
annbtbw.sd.sa$li = "sd"
annbtbw.ac.sa$li = "ac"
annbtbw.sk.sa$li = "sk"
annbtbw.kt.sa$li = "kt"

#combine dataframes
annbtbw.sa = rbind(annbtbw.sd.sa, annbtbw.ac.sa, annbtbw.sk.sa, annbtbw.kt.sa)

#add column with state variable name
annbtbw.sa$var = "annbtbw"

#melt data frame so that columns with window lengths become rows
annbtbw.sa.lon = reshape2::melt(annbtbw.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annbtbw.sa.lon$bandwith = annbtbw.sa.lon$bandyears / nrow(annbtbw)
annbtbw.sa.lon$window = as.numeric(as.character(annbtbw.sa.lon$variable)) / nrow(annbtbw)

#rename columns
names(annbtbw.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annbtbw.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annbtbw_sa.csv",
            col.names = T, row.names = F, sep = ",")



######
#BTGW#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annbtgw.pdf")

#run analyses
annbtgw.sd.sa = sensitivity_ews(timeseries = annbtgw$totbird, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtgw.ac.sa = sensitivity_ews(timeseries = annbtgw$totbird, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtgw.sk.sa = sensitivity_ews(timeseries = annbtgw$totbird, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annbtgw.kt.sa = sensitivity_ews(timeseries = annbtgw$totbird, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annbtgw.sd.sa$bandyears = as.numeric(row.names(annbtgw.sd.sa)) 
annbtgw.ac.sa$bandyears = as.numeric(row.names(annbtgw.ac.sa)) 
annbtgw.sk.sa$bandyears = as.numeric(row.names(annbtgw.sk.sa)) 
annbtgw.kt.sa$bandyears = as.numeric(row.names(annbtgw.kt.sa)) 

#add column to each dataframe with the name of the EWS
annbtgw.sd.sa$li = "sd"
annbtgw.ac.sa$li = "ac"
annbtgw.sk.sa$li = "sk"
annbtgw.kt.sa$li = "kt"

#combine dataframes
annbtgw.sa = rbind(annbtgw.sd.sa, annbtgw.ac.sa, annbtgw.sk.sa, annbtgw.kt.sa)

#add column with state variable name
annbtgw.sa$var = "annbtgw"

#melt data frame so that columns with window lengths become rows
annbtgw.sa.lon = reshape2::melt(annbtgw.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annbtgw.sa.lon$bandwith = annbtgw.sa.lon$bandyears / nrow(annbtgw)
annbtgw.sa.lon$window = as.numeric(as.character(annbtgw.sa.lon$variable)) / nrow(annbtgw)

#rename columns
names(annbtgw.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annbtgw.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annbtgw_sa.csv",
            col.names = T, row.names = F, sep = ",")



######
#LEFL#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annlefl.pdf")

#run analyses
annlefl.sd.sa = sensitivity_ews(timeseries = annlefl$totbird, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annlefl.ac.sa = sensitivity_ews(timeseries = annlefl$totbird, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annlefl.sk.sa = sensitivity_ews(timeseries = annlefl$totbird, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annlefl.kt.sa = sensitivity_ews(timeseries = annlefl$totbird, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annlefl.sd.sa$bandyears = as.numeric(row.names(annlefl.sd.sa)) 
annlefl.ac.sa$bandyears = as.numeric(row.names(annlefl.ac.sa)) 
annlefl.sk.sa$bandyears = as.numeric(row.names(annlefl.sk.sa)) 
annlefl.kt.sa$bandyears = as.numeric(row.names(annlefl.kt.sa)) 

#add column to each dataframe with the name of the EWS
annlefl.sd.sa$li = "sd"
annlefl.ac.sa$li = "ac"
annlefl.sk.sa$li = "sk"
annlefl.kt.sa$li = "kt"

#combine dataframes
annlefl.sa = rbind(annlefl.sd.sa, annlefl.ac.sa, annlefl.sk.sa, annlefl.kt.sa)

#add column with state variable name
annlefl.sa$var = "annlefl"

#melt data frame so that columns with window lengths become rows
annlefl.sa.lon = reshape2::melt(annlefl.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annlefl.sa.lon$bandwith = annlefl.sa.lon$bandyears / nrow(annlefl)
annlefl.sa.lon$window = as.numeric(as.character(annlefl.sa.lon$variable)) / nrow(annlefl)

#rename columns
names(annlefl.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annlefl.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annlefl_sa.csv",
            col.names = T, row.names = F, sep = ",")


######
#OVEN#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annoven.pdf")

#run analyses
annoven.sd.sa = sensitivity_ews(timeseries = annoven$totbird, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annoven.ac.sa = sensitivity_ews(timeseries = annoven$totbird, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annoven.sk.sa = sensitivity_ews(timeseries = annoven$totbird, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annoven.kt.sa = sensitivity_ews(timeseries = annoven$totbird, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annoven.sd.sa$bandyears = as.numeric(row.names(annoven.sd.sa)) 
annoven.ac.sa$bandyears = as.numeric(row.names(annoven.ac.sa)) 
annoven.sk.sa$bandyears = as.numeric(row.names(annoven.sk.sa)) 
annoven.kt.sa$bandyears = as.numeric(row.names(annoven.kt.sa)) 

#add column to each dataframe with the name of the EWS
annoven.sd.sa$li = "sd"
annoven.ac.sa$li = "ac"
annoven.sk.sa$li = "sk"
annoven.kt.sa$li = "kt"

#combine dataframes
annoven.sa = rbind(annoven.sd.sa, annoven.ac.sa, annoven.sk.sa, annoven.kt.sa)

#add column with state variable name
annoven.sa$var = "annoven"

#melt data frame so that columns with window lengths become rows
annoven.sa.lon = reshape2::melt(annoven.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annoven.sa.lon$bandwith = annoven.sa.lon$bandyears / nrow(annoven)
annoven.sa.lon$window = as.numeric(as.character(annoven.sa.lon$variable)) / nrow(annoven)

#rename columns
names(annoven.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annoven.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annoven_sa.csv",
            col.names = T, row.names = F, sep = ",")


######
#REVI#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_annrevi.pdf")

#run analyses
annrevi.sd.sa = sensitivity_ews(timeseries = annrevi$totbird, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annrevi.ac.sa = sensitivity_ews(timeseries = annrevi$totbird, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annrevi.sk.sa = sensitivity_ews(timeseries = annrevi$totbird, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

annrevi.kt.sa = sensitivity_ews(timeseries = annrevi$totbird, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
annrevi.sd.sa$bandyears = as.numeric(row.names(annrevi.sd.sa)) 
annrevi.ac.sa$bandyears = as.numeric(row.names(annrevi.ac.sa)) 
annrevi.sk.sa$bandyears = as.numeric(row.names(annrevi.sk.sa)) 
annrevi.kt.sa$bandyears = as.numeric(row.names(annrevi.kt.sa)) 

#add column to each dataframe with the name of the EWS
annrevi.sd.sa$li = "sd"
annrevi.ac.sa$li = "ac"
annrevi.sk.sa$li = "sk"
annrevi.kt.sa$li = "kt"

#combine dataframes
annrevi.sa = rbind(annrevi.sd.sa, annrevi.ac.sa, annrevi.sk.sa, annrevi.kt.sa)

#add column with state variable name
annrevi.sa$var = "annrevi"

#melt data frame so that columns with window lengths become rows
annrevi.sa.lon = reshape2::melt(annrevi.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
annrevi.sa.lon$bandwith = annrevi.sa.lon$bandyears / nrow(annrevi)
annrevi.sa.lon$window = as.numeric(as.character(annrevi.sa.lon$variable)) / nrow(annrevi)

#rename columns
names(annrevi.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(annrevi.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\annrevi_sa.csv",
            col.names = T, row.names = F, sep = ",")


#####################################
#####################################
#lepidopteran biomass               #
#####################################
#####################################

#run sensitivity analyses and write diagnostic plots to folder

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_lepyear.pdf")

#run analyses
lepyear.sd.sa = sensitivity_ews(timeseries = lepyear$biomass, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

lepyear.ac.sa = sensitivity_ews(timeseries = lepyear$biomass, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

lepyear.sk.sa = sensitivity_ews(timeseries = lepyear$biomass, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

lepyear.kt.sa = sensitivity_ews(timeseries = lepyear$biomass, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()


#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
lepyear.sd.sa$bandyears = as.numeric(row.names(lepyear.sd.sa)) 
lepyear.ac.sa$bandyears = as.numeric(row.names(lepyear.ac.sa)) 
lepyear.sk.sa$bandyears = as.numeric(row.names(lepyear.sk.sa)) 
lepyear.kt.sa$bandyears = as.numeric(row.names(lepyear.kt.sa)) 

#add column to each dataframe with the name of the EWS
lepyear.sd.sa$li = "sd"
lepyear.ac.sa$li = "ac"
lepyear.sk.sa$li = "sk"
lepyear.kt.sa$li = "kt"

#combine dataframes
lepyear.sa = rbind(lepyear.sd.sa, lepyear.ac.sa, lepyear.sk.sa, lepyear.kt.sa)

#add column with state variable name
lepyear.sa$var = "lepyear"

#melt data frame so that columns with window lengths become rows
lepyear.sa.lon = reshape2::melt(lepyear.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
lepyear.sa.lon$bandwith = lepyear.sa.lon$bandyears / nrow(lepyear)
lepyear.sa.lon$window = as.numeric(as.character(lepyear.sa.lon$variable)) / nrow(lepyear)

#rename columns
names(lepyear.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(lepyear.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\lepyear_sa.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#Microbial Biomass                  #
#####################################
#####################################

#####################
#Microbial Biomass C#
#####################

#run sensitivity analyses and write diagnostic plots to folder

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_BIOC.pdf")

#run analyses
BIOC.sd.sa = sensitivity_ews(timeseries = micyear$BIOC, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BIOC.ac.sa = sensitivity_ews(timeseries = micyear$BIOC, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BIOC.sk.sa = sensitivity_ews(timeseries = micyear$BIOC, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BIOC.kt.sa = sensitivity_ews(timeseries = micyear$BIOC, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
BIOC.sd.sa$bandyears = as.numeric(row.names(BIOC.sd.sa)) 
BIOC.ac.sa$bandyears = as.numeric(row.names(BIOC.ac.sa)) 
BIOC.sk.sa$bandyears = as.numeric(row.names(BIOC.sk.sa)) 
BIOC.kt.sa$bandyears = as.numeric(row.names(BIOC.kt.sa)) 

#add column to each dataframe with the name of the EWS
BIOC.sd.sa$li = "sd"
BIOC.ac.sa$li = "ac"
BIOC.sk.sa$li = "sk"
BIOC.kt.sa$li = "kt"

#combine dataframes
BIOC.sa = rbind(BIOC.sd.sa, BIOC.ac.sa, BIOC.sk.sa, BIOC.kt.sa)

#add column with state variable name
BIOC.sa$var = "BIOC"

#melt data frame so that columns with window lengths become rows
BIOC.sa.lon = reshape2::melt(BIOC.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
BIOC.sa.lon$bandwith = BIOC.sa.lon$bandyears / nrow(micyear)
BIOC.sa.lon$window = as.numeric(as.character(BIOC.sa.lon$variable)) / nrow(micyear)

#rename columns
names(BIOC.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(BIOC.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\BIOC_sa.csv",
            col.names = T, row.names = F, sep = ",")


#####################
#Microbial Biomass N#
#####################

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_BION.pdf")

#run analyses
BION.sd.sa = sensitivity_ews(timeseries = micyear$BION, indicator = "sd", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BION.ac.sa = sensitivity_ews(timeseries = micyear$BION, indicator = "acf1", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BION.sk.sa = sensitivity_ews(timeseries = micyear$BION, indicator = "sk", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

BION.kt.sa = sensitivity_ews(timeseries = micyear$BION, indicator = "kurt", 
                                detrending = 'gaussian', 
                                incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
BION.sd.sa$bandyears = as.numeric(row.names(BION.sd.sa)) 
BION.ac.sa$bandyears = as.numeric(row.names(BION.ac.sa)) 
BION.sk.sa$bandyears = as.numeric(row.names(BION.sk.sa)) 
BION.kt.sa$bandyears = as.numeric(row.names(BION.kt.sa)) 

#add column to each dataframe with the name of the EWS
BION.sd.sa$li = "sd"
BION.ac.sa$li = "ac"
BION.sk.sa$li = "sk"
BION.kt.sa$li = "kt"

#combine dataframes
BION.sa = rbind(BION.sd.sa, BION.ac.sa, BION.sk.sa, BION.kt.sa)

#add column with state variable name
BION.sa$var = "BION"

#melt data frame so that columns with window lengths become rows
BION.sa.lon = reshape2::melt(BION.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
BION.sa.lon$bandwith = BION.sa.lon$bandyears / nrow(micyear)
BION.sa.lon$window = as.numeric(as.character(BION.sa.lon$variable)) / nrow(micyear)

#rename columns
names(BION.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(BION.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\BION_sa.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#basal area increment               #
#####################################
#####################################

######
#ACSA#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ACSA.pdf")

#run analyses
ACSA.sd.sa = sensitivity_ews(timeseries = baiyear$ACSA_mean, indicator = "sd", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ACSA.ac.sa = sensitivity_ews(timeseries = baiyear$ACSA_mean, indicator = "acf1", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ACSA.sk.sa = sensitivity_ews(timeseries = baiyear$ACSA_mean, indicator = "sk", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ACSA.kt.sa = sensitivity_ews(timeseries = baiyear$ACSA_mean, indicator = "kurt", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ACSA.sd.sa$bandyears = as.numeric(row.names(ACSA.sd.sa)) 
ACSA.ac.sa$bandyears = as.numeric(row.names(ACSA.ac.sa)) 
ACSA.sk.sa$bandyears = as.numeric(row.names(ACSA.sk.sa)) 
ACSA.kt.sa$bandyears = as.numeric(row.names(ACSA.kt.sa)) 

#add column to each dataframe with the name of the EWS
ACSA.sd.sa$li = "sd"
ACSA.ac.sa$li = "ac"
ACSA.sk.sa$li = "sk"
ACSA.kt.sa$li = "kt"

#combine dataframes
ACSA.sa = rbind(ACSA.sd.sa, ACSA.ac.sa, ACSA.sk.sa, ACSA.kt.sa)

#add column with state variable name
ACSA.sa$var = "ACSA"

#melt data frame so that columns with window lengths become rows
ACSA.sa.lon = reshape2::melt(ACSA.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ACSA.sa.lon$bandwith = ACSA.sa.lon$bandyears / nrow(baiyear)
ACSA.sa.lon$window = as.numeric(as.character(ACSA.sa.lon$variable)) / nrow(baiyear)

#rename columns
names(ACSA.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ACSA.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ACSA_sa.csv",
            col.names = T, row.names = F, sep = ",")


######
#FAGR#
######

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_FAGR.pdf")

#run analyses
FAGR.sd.sa = sensitivity_ews(timeseries = baiyear$FAGR_mean, indicator = "sd", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

FAGR.ac.sa = sensitivity_ews(timeseries = baiyear$FAGR_mean, indicator = "acf1", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

FAGR.sk.sa = sensitivity_ews(timeseries = baiyear$FAGR_mean, indicator = "sk", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

FAGR.kt.sa = sensitivity_ews(timeseries = baiyear$FAGR_mean, indicator = "kurt", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
FAGR.sd.sa$bandyears = as.numeric(row.names(FAGR.sd.sa)) 
FAGR.ac.sa$bandyears = as.numeric(row.names(FAGR.ac.sa)) 
FAGR.sk.sa$bandyears = as.numeric(row.names(FAGR.sk.sa)) 
FAGR.kt.sa$bandyears = as.numeric(row.names(FAGR.kt.sa)) 

#add column to each dataframe with the name of the EWS
FAGR.sd.sa$li = "sd"
FAGR.ac.sa$li = "ac"
FAGR.sk.sa$li = "sk"
FAGR.kt.sa$li = "kt"

#combine dataframes
FAGR.sa = rbind(FAGR.sd.sa, FAGR.ac.sa, FAGR.sk.sa, FAGR.kt.sa)

#add column with state variable name
FAGR.sa$var = "FAGR"

#melt data frame so that columns with window lengths become rows
FAGR.sa.lon = reshape2::melt(FAGR.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
FAGR.sa.lon$bandwith = FAGR.sa.lon$bandyears / nrow(baiyear)
FAGR.sa.lon$window = as.numeric(as.character(FAGR.sa.lon$variable)) / nrow(baiyear)

#rename columns
names(FAGR.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(FAGR.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\FAGR_sa.csv",
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

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws1_volwt_Ca.pdf")

#run analyses
ws1_volwt_Ca.sd.sa = sensitivity_ews(timeseries = ws1_volwt$volwt_Ca, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_Ca.ac.sa = sensitivity_ews(timeseries = ws1_volwt$volwt_Ca, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_Ca.sk.sa = sensitivity_ews(timeseries = ws1_volwt$volwt_Ca, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_Ca.kt.sa = sensitivity_ews(timeseries = ws1_volwt$volwt_Ca, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws1_volwt_Ca.sd.sa$bandyears = as.numeric(row.names(ws1_volwt_Ca.sd.sa)) 
ws1_volwt_Ca.ac.sa$bandyears = as.numeric(row.names(ws1_volwt_Ca.ac.sa)) 
ws1_volwt_Ca.sk.sa$bandyears = as.numeric(row.names(ws1_volwt_Ca.sk.sa)) 
ws1_volwt_Ca.kt.sa$bandyears = as.numeric(row.names(ws1_volwt_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws1_volwt_Ca.sd.sa$li = "sd"
ws1_volwt_Ca.ac.sa$li = "ac"
ws1_volwt_Ca.sk.sa$li = "sk"
ws1_volwt_Ca.kt.sa$li = "kt"

#combine dataframes
ws1_volwt_Ca.sa = rbind(ws1_volwt_Ca.sd.sa, ws1_volwt_Ca.ac.sa, ws1_volwt_Ca.sk.sa, ws1_volwt_Ca.kt.sa)

#add column with state variable name
ws1_volwt_Ca.sa$var = "ws1_volwt_Ca"

#melt data frame so that columns with window lengths become rows
ws1_volwt_Ca.sa.lon = reshape2::melt(ws1_volwt_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws1_volwt_Ca.sa.lon$bandwith = ws1_volwt_Ca.sa.lon$bandyears / nrow(ws1_volwt)
ws1_volwt_Ca.sa.lon$window = as.numeric(as.character(ws1_volwt_Ca.sa.lon$variable)) / nrow(ws1_volwt)

#rename columns
names(ws1_volwt_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws1_volwt_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws1_volwt_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws2#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws2_volwt_Ca.pdf")

#run analyses
ws2_volwt_Ca.sd.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_Ca, indicator = "sd", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_Ca.ac.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_Ca, indicator = "acf1", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_Ca.sk.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_Ca, indicator = "sk", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_Ca.kt.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_Ca, indicator = "kurt", 
                             detrending = 'gaussian', 
                             incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws2_volwt_Ca.sd.sa$bandyears = as.numeric(row.names(ws2_volwt_Ca.sd.sa)) 
ws2_volwt_Ca.ac.sa$bandyears = as.numeric(row.names(ws2_volwt_Ca.ac.sa)) 
ws2_volwt_Ca.sk.sa$bandyears = as.numeric(row.names(ws2_volwt_Ca.sk.sa)) 
ws2_volwt_Ca.kt.sa$bandyears = as.numeric(row.names(ws2_volwt_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws2_volwt_Ca.sd.sa$li = "sd"
ws2_volwt_Ca.ac.sa$li = "ac"
ws2_volwt_Ca.sk.sa$li = "sk"
ws2_volwt_Ca.kt.sa$li = "kt"

#combine dataframes
ws2_volwt_Ca.sa = rbind(ws2_volwt_Ca.sd.sa, ws2_volwt_Ca.ac.sa, ws2_volwt_Ca.sk.sa, ws2_volwt_Ca.kt.sa)

#add column with state variable name
ws2_volwt_Ca.sa$var = "ws2_volwt_Ca"

#melt data frame so that columns with window lengths become rows
ws2_volwt_Ca.sa.lon = reshape2::melt(ws2_volwt_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws2_volwt_Ca.sa.lon$bandwith = ws2_volwt_Ca.sa.lon$bandyears / nrow(ws2_volwt)
ws2_volwt_Ca.sa.lon$window = as.numeric(as.character(ws2_volwt_Ca.sa.lon$variable)) / nrow(ws2_volwt)

#rename columns
names(ws2_volwt_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws2_volwt_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws2_volwt_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws4#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws4_volwt_Ca.pdf")

#run analyses
ws4_volwt_Ca.sd.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_Ca, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_Ca.ac.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_Ca, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_Ca.sk.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_Ca, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_Ca.kt.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_Ca, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws4_volwt_Ca.sd.sa$bandyears = as.numeric(row.names(ws4_volwt_Ca.sd.sa)) 
ws4_volwt_Ca.ac.sa$bandyears = as.numeric(row.names(ws4_volwt_Ca.ac.sa)) 
ws4_volwt_Ca.sk.sa$bandyears = as.numeric(row.names(ws4_volwt_Ca.sk.sa)) 
ws4_volwt_Ca.kt.sa$bandyears = as.numeric(row.names(ws4_volwt_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws4_volwt_Ca.sd.sa$li = "sd"
ws4_volwt_Ca.ac.sa$li = "ac"
ws4_volwt_Ca.sk.sa$li = "sk"
ws4_volwt_Ca.kt.sa$li = "kt"

#combine dataframes
ws4_volwt_Ca.sa = rbind(ws4_volwt_Ca.sd.sa, ws4_volwt_Ca.ac.sa, ws4_volwt_Ca.sk.sa, ws4_volwt_Ca.kt.sa)

#add column with state variable name
ws4_volwt_Ca.sa$var = "ws4_volwt_Ca"

#melt data frame so that columns with window lengths become rows
ws4_volwt_Ca.sa.lon = reshape2::melt(ws4_volwt_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws4_volwt_Ca.sa.lon$bandwith = ws4_volwt_Ca.sa.lon$bandyears / nrow(ws4_volwt)
ws4_volwt_Ca.sa.lon$window = as.numeric(as.character(ws4_volwt_Ca.sa.lon$variable)) / nrow(ws4_volwt)

#rename columns
names(ws4_volwt_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws4_volwt_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws4_volwt_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")


#########
#Ca ws5#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws5_volwt_Ca.pdf")

#run analyses
ws5_volwt_Ca.sd.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_Ca, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_Ca.ac.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_Ca, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_Ca.sk.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_Ca, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_Ca.kt.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_Ca, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws5_volwt_Ca.sd.sa$bandyears = as.numeric(row.names(ws5_volwt_Ca.sd.sa)) 
ws5_volwt_Ca.ac.sa$bandyears = as.numeric(row.names(ws5_volwt_Ca.ac.sa)) 
ws5_volwt_Ca.sk.sa$bandyears = as.numeric(row.names(ws5_volwt_Ca.sk.sa)) 
ws5_volwt_Ca.kt.sa$bandyears = as.numeric(row.names(ws5_volwt_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws5_volwt_Ca.sd.sa$li = "sd"
ws5_volwt_Ca.ac.sa$li = "ac"
ws5_volwt_Ca.sk.sa$li = "sk"
ws5_volwt_Ca.kt.sa$li = "kt"

#combine dataframes
ws5_volwt_Ca.sa = rbind(ws5_volwt_Ca.sd.sa, ws5_volwt_Ca.ac.sa, ws5_volwt_Ca.sk.sa, ws5_volwt_Ca.kt.sa)

#add column with state variable name
ws5_volwt_Ca.sa$var = "ws5_volwt_Ca"

#melt data frame so that columns with window lengths become rows
ws5_volwt_Ca.sa.lon = reshape2::melt(ws5_volwt_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws5_volwt_Ca.sa.lon$bandwith = ws5_volwt_Ca.sa.lon$bandyears / nrow(ws5_volwt)
ws5_volwt_Ca.sa.lon$window = as.numeric(as.character(ws5_volwt_Ca.sa.lon$variable)) / nrow(ws5_volwt)

#rename columns
names(ws5_volwt_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws5_volwt_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws5_volwt_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws6#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws6_volwt_Ca.pdf")

#run analyses
ws6_volwt_Ca.sd.sa = sensitivity_ews(timeseries = ws6_volwt$volwt_Ca, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_Ca.ac.sa = sensitivity_ews(timeseries = ws6_volwt$volwt_Ca, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_Ca.sk.sa = sensitivity_ews(timeseries = ws6_volwt$volwt_Ca, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_Ca.kt.sa = sensitivity_ews(timeseries = ws6_volwt$volwt_Ca, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws6_volwt_Ca.sd.sa$bandyears = as.numeric(row.names(ws6_volwt_Ca.sd.sa)) 
ws6_volwt_Ca.ac.sa$bandyears = as.numeric(row.names(ws6_volwt_Ca.ac.sa)) 
ws6_volwt_Ca.sk.sa$bandyears = as.numeric(row.names(ws6_volwt_Ca.sk.sa)) 
ws6_volwt_Ca.kt.sa$bandyears = as.numeric(row.names(ws6_volwt_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws6_volwt_Ca.sd.sa$li = "sd"
ws6_volwt_Ca.ac.sa$li = "ac"
ws6_volwt_Ca.sk.sa$li = "sk"
ws6_volwt_Ca.kt.sa$li = "kt"

#combine dataframes
ws6_volwt_Ca.sa = rbind(ws6_volwt_Ca.sd.sa, ws6_volwt_Ca.ac.sa, ws6_volwt_Ca.sk.sa, ws6_volwt_Ca.kt.sa)

#add column with state variable name
ws6_volwt_Ca.sa$var = "ws6_volwt_Ca"

#melt data frame so that columns with window lengths become rows
ws6_volwt_Ca.sa.lon = reshape2::melt(ws6_volwt_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws6_volwt_Ca.sa.lon$bandwith = ws6_volwt_Ca.sa.lon$bandyears / nrow(ws6_volwt)
ws6_volwt_Ca.sa.lon$window = as.numeric(as.character(ws6_volwt_Ca.sa.lon$variable)) / nrow(ws6_volwt)

#rename columns
names(ws6_volwt_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws6_volwt_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws6_volwt_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")


##################
#NO3 concentration#
##################

#########
#NO3 ws1#
#########

#run sensitivity analyses and write diagnostic plots to folder

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws1_volwt_NO3.pdf")

#run analyses
ws1_volwt_NO3.sd.sa = sensitivity_ews(timeseries = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]$volwt_NO3, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_NO3.ac.sa = sensitivity_ews(timeseries = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]$volwt_NO3, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_NO3.sk.sa = sensitivity_ews(timeseries = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]$volwt_NO3, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_volwt_NO3.kt.sa = sensitivity_ews(timeseries = ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ]$volwt_NO3, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws1_volwt_NO3.sd.sa$bandyears = as.numeric(row.names(ws1_volwt_NO3.sd.sa)) 
ws1_volwt_NO3.ac.sa$bandyears = as.numeric(row.names(ws1_volwt_NO3.ac.sa)) 
ws1_volwt_NO3.sk.sa$bandyears = as.numeric(row.names(ws1_volwt_NO3.sk.sa)) 
ws1_volwt_NO3.kt.sa$bandyears = as.numeric(row.names(ws1_volwt_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws1_volwt_NO3.sd.sa$li = "sd"
ws1_volwt_NO3.ac.sa$li = "ac"
ws1_volwt_NO3.sk.sa$li = "sk"
ws1_volwt_NO3.kt.sa$li = "kt"

#combine dataframes
ws1_volwt_NO3.sa = rbind(ws1_volwt_NO3.sd.sa, ws1_volwt_NO3.ac.sa, ws1_volwt_NO3.sk.sa, ws1_volwt_NO3.kt.sa)

#add column with state variable name
ws1_volwt_NO3.sa$var = "ws1_volwt_NO3"

#melt data frame so that columns with window lengths become rows
ws1_volwt_NO3.sa.lon = reshape2::melt(ws1_volwt_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws1_volwt_NO3.sa.lon$bandwith = ws1_volwt_NO3.sa.lon$bandyears / nrow(ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ])
ws1_volwt_NO3.sa.lon$window = as.numeric(as.character(ws1_volwt_NO3.sa.lon$variable)) /  nrow(ws1_volwt[complete.cases(ws1_volwt$volwt_NO3), ])

#rename columns
names(ws1_volwt_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws1_volwt_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws1_volwt_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws2#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws2_volwt_NO3.pdf")

#run analyses
ws2_volwt_NO3.sd.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_NO3, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_NO3.ac.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_NO3, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_NO3.sk.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_NO3, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_volwt_NO3.kt.sa = sensitivity_ews(timeseries = ws2_volwt$volwt_NO3, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws2_volwt_NO3.sd.sa$bandyears = as.numeric(row.names(ws2_volwt_NO3.sd.sa)) 
ws2_volwt_NO3.ac.sa$bandyears = as.numeric(row.names(ws2_volwt_NO3.ac.sa)) 
ws2_volwt_NO3.sk.sa$bandyears = as.numeric(row.names(ws2_volwt_NO3.sk.sa)) 
ws2_volwt_NO3.kt.sa$bandyears = as.numeric(row.names(ws2_volwt_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws2_volwt_NO3.sd.sa$li = "sd"
ws2_volwt_NO3.ac.sa$li = "ac"
ws2_volwt_NO3.sk.sa$li = "sk"
ws2_volwt_NO3.kt.sa$li = "kt"

#combine dataframes
ws2_volwt_NO3.sa = rbind(ws2_volwt_NO3.sd.sa, ws2_volwt_NO3.ac.sa, ws2_volwt_NO3.sk.sa, ws2_volwt_NO3.kt.sa)

#add column with state variable name
ws2_volwt_NO3.sa$var = "ws2_volwt_NO3"

#melt data frame so that columns with window lengths become rows
ws2_volwt_NO3.sa.lon = reshape2::melt(ws2_volwt_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws2_volwt_NO3.sa.lon$bandwith = ws2_volwt_NO3.sa.lon$bandyears / nrow(ws2_volwt)
ws2_volwt_NO3.sa.lon$window = as.numeric(as.character(ws2_volwt_NO3.sa.lon$variable)) / nrow(ws2_volwt)

#rename columns
names(ws2_volwt_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws2_volwt_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws2_volwt_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws4#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws4_volwt_NO3.pdf")

#run analyses
ws4_volwt_NO3.sd.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_NO3, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_NO3.ac.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_NO3, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_NO3.sk.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_NO3, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_volwt_NO3.kt.sa = sensitivity_ews(timeseries = ws4_volwt$volwt_NO3, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws4_volwt_NO3.sd.sa$bandyears = as.numeric(row.names(ws4_volwt_NO3.sd.sa)) 
ws4_volwt_NO3.ac.sa$bandyears = as.numeric(row.names(ws4_volwt_NO3.ac.sa)) 
ws4_volwt_NO3.sk.sa$bandyears = as.numeric(row.names(ws4_volwt_NO3.sk.sa)) 
ws4_volwt_NO3.kt.sa$bandyears = as.numeric(row.names(ws4_volwt_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws4_volwt_NO3.sd.sa$li = "sd"
ws4_volwt_NO3.ac.sa$li = "ac"
ws4_volwt_NO3.sk.sa$li = "sk"
ws4_volwt_NO3.kt.sa$li = "kt"

#combine dataframes
ws4_volwt_NO3.sa = rbind(ws4_volwt_NO3.sd.sa, ws4_volwt_NO3.ac.sa, ws4_volwt_NO3.sk.sa, ws4_volwt_NO3.kt.sa)

#add column with state variable name
ws4_volwt_NO3.sa$var = "ws4_volwt_NO3"

#melt data frame so that columns with window lengths become rows
ws4_volwt_NO3.sa.lon = reshape2::melt(ws4_volwt_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws4_volwt_NO3.sa.lon$bandwith = ws4_volwt_NO3.sa.lon$bandyears / nrow(ws4_volwt)
ws4_volwt_NO3.sa.lon$window = as.numeric(as.character(ws4_volwt_NO3.sa.lon$variable)) / nrow(ws4_volwt)

#rename columns
names(ws4_volwt_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws4_volwt_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws4_volwt_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws5#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws5_volwt_NO3.pdf")

#run analyses
ws5_volwt_NO3.sd.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_NO3, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_NO3.ac.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_NO3, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_NO3.sk.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_NO3, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_volwt_NO3.kt.sa = sensitivity_ews(timeseries = ws5_volwt$volwt_NO3, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws5_volwt_NO3.sd.sa$bandyears = as.numeric(row.names(ws5_volwt_NO3.sd.sa)) 
ws5_volwt_NO3.ac.sa$bandyears = as.numeric(row.names(ws5_volwt_NO3.ac.sa)) 
ws5_volwt_NO3.sk.sa$bandyears = as.numeric(row.names(ws5_volwt_NO3.sk.sa)) 
ws5_volwt_NO3.kt.sa$bandyears = as.numeric(row.names(ws5_volwt_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws5_volwt_NO3.sd.sa$li = "sd"
ws5_volwt_NO3.ac.sa$li = "ac"
ws5_volwt_NO3.sk.sa$li = "sk"
ws5_volwt_NO3.kt.sa$li = "kt"

#combine dataframes
ws5_volwt_NO3.sa = rbind(ws5_volwt_NO3.sd.sa, ws5_volwt_NO3.ac.sa, ws5_volwt_NO3.sk.sa, ws5_volwt_NO3.kt.sa)

#add column with state variable name
ws5_volwt_NO3.sa$var = "ws5_volwt_NO3"

#melt data frame so that columns with window lengths become rows
ws5_volwt_NO3.sa.lon = reshape2::melt(ws5_volwt_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws5_volwt_NO3.sa.lon$bandwith = ws5_volwt_NO3.sa.lon$bandyears / nrow(ws5_volwt)
ws5_volwt_NO3.sa.lon$window = as.numeric(as.character(ws5_volwt_NO3.sa.lon$variable)) / nrow(ws5_volwt)

#rename columns
names(ws5_volwt_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws5_volwt_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws5_volwt_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws6#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws6_volwt_NO3.pdf")

#run analyses
ws6_volwt_NO3.sd.sa = sensitivity_ews(timeseries = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]$volwt_NO3, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_NO3.ac.sa = sensitivity_ews(timeseries = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]$volwt_NO3, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_NO3.sk.sa = sensitivity_ews(timeseries = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]$volwt_NO3, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_volwt_NO3.kt.sa = sensitivity_ews(timeseries = ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ]$volwt_NO3, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws6_volwt_NO3.sd.sa$bandyears = as.numeric(row.names(ws6_volwt_NO3.sd.sa)) 
ws6_volwt_NO3.ac.sa$bandyears = as.numeric(row.names(ws6_volwt_NO3.ac.sa)) 
ws6_volwt_NO3.sk.sa$bandyears = as.numeric(row.names(ws6_volwt_NO3.sk.sa)) 
ws6_volwt_NO3.kt.sa$bandyears = as.numeric(row.names(ws6_volwt_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws6_volwt_NO3.sd.sa$li = "sd"
ws6_volwt_NO3.ac.sa$li = "ac"
ws6_volwt_NO3.sk.sa$li = "sk"
ws6_volwt_NO3.kt.sa$li = "kt"

#combine dataframes
ws6_volwt_NO3.sa = rbind(ws6_volwt_NO3.sd.sa, ws6_volwt_NO3.ac.sa, ws6_volwt_NO3.sk.sa, ws6_volwt_NO3.kt.sa)

#add column with state variable name
ws6_volwt_NO3.sa$var = "ws6_volwt_NO3"

#melt data frame so that columns with window lengths become rows
ws6_volwt_NO3.sa.lon = reshape2::melt(ws6_volwt_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws6_volwt_NO3.sa.lon$bandwith = ws6_volwt_NO3.sa.lon$bandyears / nrow(ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ])
ws6_volwt_NO3.sa.lon$window = as.numeric(as.character(ws6_volwt_NO3.sa.lon$variable)) / nrow(ws6_volwt[complete.cases(ws6_volwt$volwt_NO3), ])

#rename columns
names(ws6_volwt_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws6_volwt_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws6_volwt_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#####################################
#####################################
#stream chemistry: flux             #
#####################################
#####################################

##################
#Ca flux         #
##################

#########
#Ca ws1#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws1_flux_Ca.pdf")

#run analyses
ws1_flux_Ca.sd.sa = sensitivity_ews(timeseries = ws1_flux$Ca_flux, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_flux_Ca.ac.sa = sensitivity_ews(timeseries = ws1_flux$Ca_flux, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_flux_Ca.sk.sa = sensitivity_ews(timeseries = ws1_flux$Ca_flux, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws1_flux_Ca.kt.sa = sensitivity_ews(timeseries = ws1_flux$Ca_flux, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws1_flux_Ca.sd.sa$bandyears = as.numeric(row.names(ws1_flux_Ca.sd.sa)) 
ws1_flux_Ca.ac.sa$bandyears = as.numeric(row.names(ws1_flux_Ca.ac.sa)) 
ws1_flux_Ca.sk.sa$bandyears = as.numeric(row.names(ws1_flux_Ca.sk.sa)) 
ws1_flux_Ca.kt.sa$bandyears = as.numeric(row.names(ws1_flux_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws1_flux_Ca.sd.sa$li = "sd"
ws1_flux_Ca.ac.sa$li = "ac"
ws1_flux_Ca.sk.sa$li = "sk"
ws1_flux_Ca.kt.sa$li = "kt"

#combine dataframes
ws1_flux_Ca.sa = rbind(ws1_flux_Ca.sd.sa, ws1_flux_Ca.ac.sa, ws1_flux_Ca.sk.sa, ws1_flux_Ca.kt.sa)

#add column with state variable name
ws1_flux_Ca.sa$var = "ws1_flux_Ca"

#melt data frame so that columns with window lengths become rows
ws1_flux_Ca.sa.lon = reshape2::melt(ws1_flux_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws1_flux_Ca.sa.lon$bandwith = ws1_flux_Ca.sa.lon$bandyears / nrow(ws1_volwt)
ws1_flux_Ca.sa.lon$window = as.numeric(as.character(ws1_flux_Ca.sa.lon$variable)) / nrow(ws1_volwt)

#rename columns
names(ws1_flux_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws1_flux_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws1_flux_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws2#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws2_flux_Ca.pdf")

#run analyses
ws2_flux_Ca.sd.sa = sensitivity_ews(timeseries = ws2_flux$Ca_flux, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_flux_Ca.ac.sa = sensitivity_ews(timeseries = ws2_flux$Ca_flux, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_flux_Ca.sk.sa = sensitivity_ews(timeseries = ws2_flux$Ca_flux, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws2_flux_Ca.kt.sa = sensitivity_ews(timeseries = ws2_flux$Ca_flux, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws2_flux_Ca.sd.sa$bandyears = as.numeric(row.names(ws2_flux_Ca.sd.sa)) 
ws2_flux_Ca.ac.sa$bandyears = as.numeric(row.names(ws2_flux_Ca.ac.sa)) 
ws2_flux_Ca.sk.sa$bandyears = as.numeric(row.names(ws2_flux_Ca.sk.sa)) 
ws2_flux_Ca.kt.sa$bandyears = as.numeric(row.names(ws2_flux_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws2_flux_Ca.sd.sa$li = "sd"
ws2_flux_Ca.ac.sa$li = "ac"
ws2_flux_Ca.sk.sa$li = "sk"
ws2_flux_Ca.kt.sa$li = "kt"

#combine dataframes
ws2_flux_Ca.sa = rbind(ws2_flux_Ca.sd.sa, ws2_flux_Ca.ac.sa, ws2_flux_Ca.sk.sa, ws2_flux_Ca.kt.sa)

#add column with state variable name
ws2_flux_Ca.sa$var = "ws2_flux_Ca"

#melt data frame so that columns with window lengths become rows
ws2_flux_Ca.sa.lon = reshape2::melt(ws2_flux_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws2_flux_Ca.sa.lon$bandwith = ws2_flux_Ca.sa.lon$bandyears / nrow(ws2_volwt)
ws2_flux_Ca.sa.lon$window = as.numeric(as.character(ws2_flux_Ca.sa.lon$variable)) / nrow(ws2_volwt)

#rename columns
names(ws2_flux_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws2_flux_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws2_flux_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws4#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws4_flux_Ca.pdf")

#run analyses
ws4_flux_Ca.sd.sa = sensitivity_ews(timeseries = ws4_flux$Ca_flux, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_flux_Ca.ac.sa = sensitivity_ews(timeseries = ws4_flux$Ca_flux, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_flux_Ca.sk.sa = sensitivity_ews(timeseries = ws4_flux$Ca_flux, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws4_flux_Ca.kt.sa = sensitivity_ews(timeseries = ws4_flux$Ca_flux, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws4_flux_Ca.sd.sa$bandyears = as.numeric(row.names(ws4_flux_Ca.sd.sa)) 
ws4_flux_Ca.ac.sa$bandyears = as.numeric(row.names(ws4_flux_Ca.ac.sa)) 
ws4_flux_Ca.sk.sa$bandyears = as.numeric(row.names(ws4_flux_Ca.sk.sa)) 
ws4_flux_Ca.kt.sa$bandyears = as.numeric(row.names(ws4_flux_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws4_flux_Ca.sd.sa$li = "sd"
ws4_flux_Ca.ac.sa$li = "ac"
ws4_flux_Ca.sk.sa$li = "sk"
ws4_flux_Ca.kt.sa$li = "kt"

#combine dataframes
ws4_flux_Ca.sa = rbind(ws4_flux_Ca.sd.sa, ws4_flux_Ca.ac.sa, ws4_flux_Ca.sk.sa, ws4_flux_Ca.kt.sa)

#add column with state variable name
ws4_flux_Ca.sa$var = "ws4_flux_Ca"

#melt data frame so that columns with window lengths become rows
ws4_flux_Ca.sa.lon = reshape2::melt(ws4_flux_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws4_flux_Ca.sa.lon$bandwith = ws4_flux_Ca.sa.lon$bandyears / nrow(ws4_volwt)
ws4_flux_Ca.sa.lon$window = as.numeric(as.character(ws4_flux_Ca.sa.lon$variable)) / nrow(ws4_volwt)

#rename columns
names(ws4_flux_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws4_flux_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws4_flux_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws5#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws5_flux_Ca.pdf")

#run analyses
ws5_flux_Ca.sd.sa = sensitivity_ews(timeseries = ws5_flux$Ca_flux, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_flux_Ca.ac.sa = sensitivity_ews(timeseries = ws5_flux$Ca_flux, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_flux_Ca.sk.sa = sensitivity_ews(timeseries = ws5_flux$Ca_flux, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws5_flux_Ca.kt.sa = sensitivity_ews(timeseries = ws5_flux$Ca_flux, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws5_flux_Ca.sd.sa$bandyears = as.numeric(row.names(ws5_flux_Ca.sd.sa)) 
ws5_flux_Ca.ac.sa$bandyears = as.numeric(row.names(ws5_flux_Ca.ac.sa)) 
ws5_flux_Ca.sk.sa$bandyears = as.numeric(row.names(ws5_flux_Ca.sk.sa)) 
ws5_flux_Ca.kt.sa$bandyears = as.numeric(row.names(ws5_flux_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws5_flux_Ca.sd.sa$li = "sd"
ws5_flux_Ca.ac.sa$li = "ac"
ws5_flux_Ca.sk.sa$li = "sk"
ws5_flux_Ca.kt.sa$li = "kt"

#combine dataframes
ws5_flux_Ca.sa = rbind(ws5_flux_Ca.sd.sa, ws5_flux_Ca.ac.sa, ws5_flux_Ca.sk.sa, ws5_flux_Ca.kt.sa)

#add column with state variable name
ws5_flux_Ca.sa$var = "ws5_flux_Ca"

#melt data frame so that columns with window lengths become rows
ws5_flux_Ca.sa.lon = reshape2::melt(ws5_flux_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws5_flux_Ca.sa.lon$bandwith = ws5_flux_Ca.sa.lon$bandyears / nrow(ws5_volwt)
ws5_flux_Ca.sa.lon$window = as.numeric(as.character(ws5_flux_Ca.sa.lon$variable)) / nrow(ws5_volwt)

#rename columns
names(ws5_flux_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws5_flux_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws5_flux_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#Ca ws6#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws6_flux_Ca.pdf")

#run analyses
ws6_flux_Ca.sd.sa = sensitivity_ews(timeseries = ws6_flux$Ca_flux, indicator = "sd", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_flux_Ca.ac.sa = sensitivity_ews(timeseries = ws6_flux$Ca_flux, indicator = "acf1", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_flux_Ca.sk.sa = sensitivity_ews(timeseries = ws6_flux$Ca_flux, indicator = "sk", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

ws6_flux_Ca.kt.sa = sensitivity_ews(timeseries = ws6_flux$Ca_flux, indicator = "kurt", 
                                     detrending = 'gaussian', 
                                     incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws6_flux_Ca.sd.sa$bandyears = as.numeric(row.names(ws6_flux_Ca.sd.sa)) 
ws6_flux_Ca.ac.sa$bandyears = as.numeric(row.names(ws6_flux_Ca.ac.sa)) 
ws6_flux_Ca.sk.sa$bandyears = as.numeric(row.names(ws6_flux_Ca.sk.sa)) 
ws6_flux_Ca.kt.sa$bandyears = as.numeric(row.names(ws6_flux_Ca.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws6_flux_Ca.sd.sa$li = "sd"
ws6_flux_Ca.ac.sa$li = "ac"
ws6_flux_Ca.sk.sa$li = "sk"
ws6_flux_Ca.kt.sa$li = "kt"

#combine dataframes
ws6_flux_Ca.sa = rbind(ws6_flux_Ca.sd.sa, ws6_flux_Ca.ac.sa, ws6_flux_Ca.sk.sa, ws6_flux_Ca.kt.sa)

#add column with state variable name
ws6_flux_Ca.sa$var = "ws6_flux_Ca"

#melt data frame so that columns with window lengths become rows
ws6_flux_Ca.sa.lon = reshape2::melt(ws6_flux_Ca.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws6_flux_Ca.sa.lon$bandwith = ws6_flux_Ca.sa.lon$bandyears / nrow(ws6_volwt)
ws6_flux_Ca.sa.lon$window = as.numeric(as.character(ws6_flux_Ca.sa.lon$variable)) / nrow(ws6_volwt)

#rename columns
names(ws6_flux_Ca.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws6_flux_Ca.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws6_flux_Ca_sa.csv",
            col.names = T, row.names = F, sep = ",")

##################
#NO3 flux        #  
##################

#########
#NO3 ws1#
#########

#run sensitivity analyses and write diagnostic plots to folder

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws1_flux_NO3.pdf")

#run analyses
ws1_flux_NO3.sd.sa = sensitivity_ews(timeseries = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]$NO3_flux, indicator = "sd", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws1_flux_NO3.ac.sa = sensitivity_ews(timeseries = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]$NO3_flux, indicator = "acf1", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws1_flux_NO3.sk.sa = sensitivity_ews(timeseries = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]$NO3_flux, indicator = "sk", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws1_flux_NO3.kt.sa = sensitivity_ews(timeseries = ws1_flux[complete.cases(ws1_flux$NO3_flux), ]$NO3_flux, indicator = "kurt", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws1_flux_NO3.sd.sa$bandyears = as.numeric(row.names(ws1_flux_NO3.sd.sa)) 
ws1_flux_NO3.ac.sa$bandyears = as.numeric(row.names(ws1_flux_NO3.ac.sa)) 
ws1_flux_NO3.sk.sa$bandyears = as.numeric(row.names(ws1_flux_NO3.sk.sa)) 
ws1_flux_NO3.kt.sa$bandyears = as.numeric(row.names(ws1_flux_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws1_flux_NO3.sd.sa$li = "sd"
ws1_flux_NO3.ac.sa$li = "ac"
ws1_flux_NO3.sk.sa$li = "sk"
ws1_flux_NO3.kt.sa$li = "kt"

#combine dataframes
ws1_flux_NO3.sa = rbind(ws1_flux_NO3.sd.sa, ws1_flux_NO3.ac.sa, ws1_flux_NO3.sk.sa, ws1_flux_NO3.kt.sa)

#add column with state variable name
ws1_flux_NO3.sa$var = "ws1_flux_NO3"

#melt data frame so that columns with window lengths become rows
ws1_flux_NO3.sa.lon = reshape2::melt(ws1_flux_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws1_flux_NO3.sa.lon$bandwith = ws1_flux_NO3.sa.lon$bandyears / nrow(ws1_volwt[complete.cases(ws1_flux$NO3_flux), ])
ws1_flux_NO3.sa.lon$window = as.numeric(as.character(ws1_flux_NO3.sa.lon$variable)) /  nrow(ws1_volwt[complete.cases(ws1_flux$NO3_flux), ])

#rename columns
names(ws1_flux_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws1_flux_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws1_flux_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws2#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws2_flux_NO3.pdf")

#run analyses
ws2_flux_NO3.sd.sa = sensitivity_ews(timeseries = ws2_flux$NO3_flux, indicator = "sd", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws2_flux_NO3.ac.sa = sensitivity_ews(timeseries = ws2_flux$NO3_flux, indicator = "acf1", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws2_flux_NO3.sk.sa = sensitivity_ews(timeseries = ws2_flux$NO3_flux, indicator = "sk", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws2_flux_NO3.kt.sa = sensitivity_ews(timeseries = ws2_flux$NO3_flux, indicator = "kurt", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws2_flux_NO3.sd.sa$bandyears = as.numeric(row.names(ws2_flux_NO3.sd.sa)) 
ws2_flux_NO3.ac.sa$bandyears = as.numeric(row.names(ws2_flux_NO3.ac.sa)) 
ws2_flux_NO3.sk.sa$bandyears = as.numeric(row.names(ws2_flux_NO3.sk.sa)) 
ws2_flux_NO3.kt.sa$bandyears = as.numeric(row.names(ws2_flux_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws2_flux_NO3.sd.sa$li = "sd"
ws2_flux_NO3.ac.sa$li = "ac"
ws2_flux_NO3.sk.sa$li = "sk"
ws2_flux_NO3.kt.sa$li = "kt"

#combine dataframes
ws2_flux_NO3.sa = rbind(ws2_flux_NO3.sd.sa, ws2_flux_NO3.ac.sa, ws2_flux_NO3.sk.sa, ws2_flux_NO3.kt.sa)

#add column with state variable name
ws2_flux_NO3.sa$var = "ws2_flux_NO3"

#melt data frame so that columns with window lengths become rows
ws2_flux_NO3.sa.lon = reshape2::melt(ws2_flux_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws2_flux_NO3.sa.lon$bandwith = ws2_flux_NO3.sa.lon$bandyears / nrow(ws2_volwt)
ws2_flux_NO3.sa.lon$window = as.numeric(as.character(ws2_flux_NO3.sa.lon$variable)) / nrow(ws2_volwt)

#rename columns
names(ws2_flux_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws2_flux_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws2_flux_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")


#########
#NO3 ws4#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws4_flux_NO3.pdf")

#run analyses
ws4_flux_NO3.sd.sa = sensitivity_ews(timeseries = ws4_flux$NO3_flux, indicator = "sd", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws4_flux_NO3.ac.sa = sensitivity_ews(timeseries = ws4_flux$NO3_flux, indicator = "acf1", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws4_flux_NO3.sk.sa = sensitivity_ews(timeseries = ws4_flux$NO3_flux, indicator = "sk", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws4_flux_NO3.kt.sa = sensitivity_ews(timeseries = ws4_flux$NO3_flux, indicator = "kurt", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws4_flux_NO3.sd.sa$bandyears = as.numeric(row.names(ws4_flux_NO3.sd.sa)) 
ws4_flux_NO3.ac.sa$bandyears = as.numeric(row.names(ws4_flux_NO3.ac.sa)) 
ws4_flux_NO3.sk.sa$bandyears = as.numeric(row.names(ws4_flux_NO3.sk.sa)) 
ws4_flux_NO3.kt.sa$bandyears = as.numeric(row.names(ws4_flux_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws4_flux_NO3.sd.sa$li = "sd"
ws4_flux_NO3.ac.sa$li = "ac"
ws4_flux_NO3.sk.sa$li = "sk"
ws4_flux_NO3.kt.sa$li = "kt"

#combine dataframes
ws4_flux_NO3.sa = rbind(ws4_flux_NO3.sd.sa, ws4_flux_NO3.ac.sa, ws4_flux_NO3.sk.sa, ws4_flux_NO3.kt.sa)

#add column with state variable name
ws4_flux_NO3.sa$var = "ws4_flux_NO3"

#melt data frame so that columns with window lengths become rows
ws4_flux_NO3.sa.lon = reshape2::melt(ws4_flux_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws4_flux_NO3.sa.lon$bandwith = ws4_flux_NO3.sa.lon$bandyears / nrow(ws4_volwt)
ws4_flux_NO3.sa.lon$window = as.numeric(as.character(ws4_flux_NO3.sa.lon$variable)) / nrow(ws4_volwt)

#rename columns
names(ws4_flux_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws4_flux_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws4_flux_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws5#
#########

#run sensitivity analyses and write diagnostic plots to folder

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures")

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws5_flux_NO3.pdf")

#run analyses
ws5_flux_NO3.sd.sa = sensitivity_ews(timeseries = ws5_flux$NO3_flux, indicator = "sd", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws5_flux_NO3.ac.sa = sensitivity_ews(timeseries = ws5_flux$NO3_flux, indicator = "acf1", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws5_flux_NO3.sk.sa = sensitivity_ews(timeseries = ws5_flux$NO3_flux, indicator = "sk", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws5_flux_NO3.kt.sa = sensitivity_ews(timeseries = ws5_flux$NO3_flux, indicator = "kurt", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws5_flux_NO3.sd.sa$bandyears = as.numeric(row.names(ws5_flux_NO3.sd.sa)) 
ws5_flux_NO3.ac.sa$bandyears = as.numeric(row.names(ws5_flux_NO3.ac.sa)) 
ws5_flux_NO3.sk.sa$bandyears = as.numeric(row.names(ws5_flux_NO3.sk.sa)) 
ws5_flux_NO3.kt.sa$bandyears = as.numeric(row.names(ws5_flux_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws5_flux_NO3.sd.sa$li = "sd"
ws5_flux_NO3.ac.sa$li = "ac"
ws5_flux_NO3.sk.sa$li = "sk"
ws5_flux_NO3.kt.sa$li = "kt"

#combine dataframes
ws5_flux_NO3.sa = rbind(ws5_flux_NO3.sd.sa, ws5_flux_NO3.ac.sa, ws5_flux_NO3.sk.sa, ws5_flux_NO3.kt.sa)

#add column with state variable name
ws5_flux_NO3.sa$var = "ws5_flux_NO3"

#melt data frame so that columns with window lengths become rows
ws5_flux_NO3.sa.lon = reshape2::melt(ws5_flux_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws5_flux_NO3.sa.lon$bandwith = ws5_flux_NO3.sa.lon$bandyears / nrow(ws5_volwt)
ws5_flux_NO3.sa.lon$window = as.numeric(as.character(ws5_flux_NO3.sa.lon$variable)) / nrow(ws5_volwt)

#rename columns
names(ws5_flux_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws5_flux_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws5_flux_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

#########
#NO3 ws6#
#########

#run sensitivity analyses and write diagnostic plots to folder

#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Figures\\SA_ws6_flux_NO3.pdf")

#run analyses
ws6_flux_NO3.sd.sa = sensitivity_ews(timeseries = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]$NO3_flux, indicator = "sd", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws6_flux_NO3.ac.sa = sensitivity_ews(timeseries = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]$NO3_flux, indicator = "acf1", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws6_flux_NO3.sk.sa = sensitivity_ews(timeseries = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]$NO3_flux, indicator = "sk", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

ws6_flux_NO3.kt.sa = sensitivity_ews(timeseries = ws6_flux[complete.cases(ws6_flux$NO3_flux), ]$NO3_flux, indicator = "kurt", 
                                      detrending = 'gaussian', 
                                      incrwinsize = 5, incrbandwidth = 5)

dev.off()

#create single table identifying with Kendall's tau values for each rolling window by 
#fitted bandwith combination

#add column with names of row indices. These are the numbers of years used in the bandwiths
ws6_flux_NO3.sd.sa$bandyears = as.numeric(row.names(ws6_flux_NO3.sd.sa)) 
ws6_flux_NO3.ac.sa$bandyears = as.numeric(row.names(ws6_flux_NO3.ac.sa)) 
ws6_flux_NO3.sk.sa$bandyears = as.numeric(row.names(ws6_flux_NO3.sk.sa)) 
ws6_flux_NO3.kt.sa$bandyears = as.numeric(row.names(ws6_flux_NO3.kt.sa)) 

#add column to each dataframe with the name of the EWS
ws6_flux_NO3.sd.sa$li = "sd"
ws6_flux_NO3.ac.sa$li = "ac"
ws6_flux_NO3.sk.sa$li = "sk"
ws6_flux_NO3.kt.sa$li = "kt"

#combine dataframes
ws6_flux_NO3.sa = rbind(ws6_flux_NO3.sd.sa, ws6_flux_NO3.ac.sa, ws6_flux_NO3.sk.sa, ws6_flux_NO3.kt.sa)

#add column with state variable name
ws6_flux_NO3.sa$var = "ws6_flux_NO3"

#melt data frame so that columns with window lengths become rows
ws6_flux_NO3.sa.lon = reshape2::melt(ws6_flux_NO3.sa, id = c("var", "li", "bandyears"))

#make new columns that computes window size (for EWS calculation) and bandwith (for Gaussian filtering) as a proportion of the time series
ws6_flux_NO3.sa.lon$bandwith = ws6_flux_NO3.sa.lon$bandyears / nrow(ws6_volwt[complete.cases(ws6_flux$NO3_flux), ])
ws6_flux_NO3.sa.lon$window = as.numeric(as.character(ws6_flux_NO3.sa.lon$variable)) / nrow(ws6_volwt[complete.cases(ws6_flux$NO3_flux), ])

#rename columns
names(ws6_flux_NO3.sa.lon) = c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#export table
write.table(ws6_flux_NO3.sa.lon, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Sensitivity Analysis Tables\\ws6_flux_NO3_sa.csv",
            col.names = T, row.names = F, sep = ",")

########################################################################################
########################################################################################
########################################################################################

#determine the most common window size by bandwith combination that results in the max
#tau statistic and min p-value

#set working directory 
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience")

#read in files, create new column identifying response variable, and combine into single dataframe

#call file names in Leading Indicators Data subdirectory
files.sens = list.files(path = "Sensitivity Analysis Tables\\", full.names = T)

#make empty data frame for pasting files
all.sens = data.frame(matrix(vector(), nrow = 1, ncol = 7))
names(all.sens)=c("var", "li", "bandyears", "winyears", "ktau", "bandwith", "window")

#loop through stations in files.sens
for(i in 1:length(files.sens)){#
  
  sens = read.table(files.sens[i], header = T, sep = ",")
  
  all.sens = rbind(all.sens, sens)
  
}

#remove first row (all NAs)
all.sens = all.sens[-1, ]

#for each state variable and EWS, determine the max Kendall's tau and associated 
#window size and bandwith.

maxkt = (all.sens
              %>% group_by(var, li)
              %>% slice( which.max( abs(ktau) ) )
              %>% as.data.frame
)


#add columns to all.sens and maxkt for plotting distributions of ktau values
all.sens$bin = "All"
maxkt$bin = "Max"

pl.sens = rbind(all.sens, maxkt)

#plot the density of window sizes and bandwiths that result in the highest Kendall's tau
#rearrange data so that window and bandwith are rows (long format)

maxwin = ggplot(maxkt, aes(window))+
        geom_density()
maxwin + facet_grid(. ~li)

maxband = ggplot(maxkt, aes(bandwith))+
  geom_density()
maxband + facet_grid(. ~li)

winkt = pl.sens[ , c("var", "li", "ktau", "window", "bin")]
bandkt = pl.sens[ , c("var", "li", "ktau", "bandwith", "bin")]

names(winkt) = c("var", "li", "ktau", "Duration", "Dataset")
names(bandkt) = c("var", "li", "ktau", "Duration", "Dataset")

winkt$sens = "Window"
bandkt$sens = "Bandwith"

allkt = rbind(winkt, bandkt)
allkt = allkt[allkt$Dataset == "Max", ]

#plot density of max ktau values for each ews

#redefine and rename variables for plotting in the order you want
allkt$li2 = ifelse(allkt$li == "sd", "A",
                   ifelse(allkt$li == "ac", "B",
                          ifelse(allkt$li == "sk", "C",
                                 ifelse(allkt$li == "kt", "D", NA
                   ))))

allkt$li2 = factor(allkt$li2, levels = c("A", "B", "C", "D"), 
                             labels = c("SD", "AC", "SK", "KT"))

#write plot to folder


#call pdf options
pdf.options(width= 6.5, height= 4.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\EWS Summary Figures\\Figure_S4.pdf")

plotkt = ggplot(allkt, aes(x = Duration))+
  geom_density(alpha = 0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Density")
plotkt + facet_grid(rows = vars(sens), cols = vars(li2))

dev.off()
