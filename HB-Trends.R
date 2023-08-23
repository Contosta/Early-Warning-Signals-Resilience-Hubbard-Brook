################################################################################################################
################################################################################################################
################################################################################################################

#This script calculates trends over time for early warning signals 
#(standard deviation, autocorrelation, skewness, and kurtosis)
#calculated from detrended data to examine changing ecological resilience.

#Input data are calculated early warning signals located in the Early Warning Signals Tables folder. Outputs are
#results from Mann-Kendall analysis of change over time and are written as two tables. alltrends.csv is in long format and is an
#input for further analysis in the HB-Synthesis script. tabs2.csv is in wide format and is the basis for
#Table S2 in the manuscript. Both output tables are written to the EWS Summary Tables folder.

#code developed by A Contosta
#most recent version 5/16/2023

#changes include code editing and cleaning for manuscript currently in press in Environmental Research Letters 

################################################################################################################
################################################################################################################
################################################################################################################

#Initial set-up

#call libraries
library(dplyr)
library(stringr)
library(trend)

####################################
#read in early warning signals data#
####################################

#set working directory 
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables")

#read in files, create new column identifying response variable, and combine into single dataframe

#call file names in Early Warning Signals Tables subdirectory
files.dat = list.files(path = ".", full.names = F)

#make empty data frame for pasting files
dat.all = data.frame(matrix(vector(), nrow = 1, ncol = 11))
names(dat.all) = c("timeindex", "Year", "ar1", "sd", "sk", "kurt", "cv", "returnrate",
                 "densratio", "acf1", "var")


#loop through stations in files.ews
for(i in 1:length(files.dat)){
  
  #split file name
  var.name = data.frame(do.call(rbind, str_split(files.dat, "_EWS.csv")))
  names(var.name) = c("var", "blank")

  #read in file
  dat = read.table(files.dat[i], header = T, sep = ",")
  
  #extract state variable from file name
  dat$var = var.name$var[i]
  
  #change name of second column to Year for consistency across files
  names(dat)[2] = "Year"
  
  #bind files together
  dat.all = rbind(dat.all, dat)
  
}

#remove first row (all NAs)
dat.all = dat.all[-1, ]

###########################################################################################
#calculate trends over time in standard deviation, autocorrelation, skewness, and kurtosis#
###########################################################################################

#omit duplicate variables from the beginning and the end of all.ews to determine when to 
#start / end model subsets
start.rows = !duplicated(dat.all$var)
end.rows = !duplicated(dat.all$var, fromLast = TRUE)

#extract original row indexes associated with the start and end of unique state variable
sr.ind = seq_along(dat.all$var)[start.rows]
er.ind = seq_along(dat.all$var)[end.rows]

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#create containers to hold results of trend analysis
sd.trends = data.frame(matrix(nrow = length(sr.ind), ncol = 4))
names(sd.trends) = c("var", "ews", "tau", "pval")

ar1.trends = data.frame(matrix(nrow = length(sr.ind), ncol = 4))
names(ar1.trends) = c("var", "ews", "tau", "pval")

sk.trends = data.frame(matrix(nrow = length(sr.ind), ncol = 4))
names(sk.trends) = c("var", "ews", "tau", "pval")

kt.trends = data.frame(matrix(nrow = length(sr.ind), ncol = 4))
names(kt.trends) = c("var", "ews", "tau", "pval")

#begin loop
for (h in seq(1,length(sr.ind))) {
    
    #subset dat.all for state variable [h]
    dat.sub = dat.all[sr.ind[h]:er.ind[h], ]
    
    #write name of state variable to dataframe
    sd.trends[h, "var"] = unique(dat.sub$var)
    ar1.trends[h, "var"] = unique(dat.sub$var)
    sk.trends[h, "var"] = unique(dat.sub$var)
    kt.trends[h, "var"] = unique(dat.sub$var)
    
    #calculate Kendall's tau and associated p-value and write to 
    #dataframes
    
    ####
    #sd#
    ####
    
    #make a data frame for year x sd combination 
    sd = data.frame(dat.sub$Year, dat.sub$sd)
    #remove NA values (otherwise mk.test won't work)
    sd. = sd[complete.cases(sd),]
    if(nrow(sd.) > 0 ){
      
    sd.trends[h, "tau"] = mk.test(sd.[,2])$estimates[3]
    sd.trends[h, "pval"] = mk.test(sd.[,2])$p.value
 
    }
    
    
    #####
    #ar1#
    #####
    
    #make a data frame for year x ar1 combination 
    ar1 = data.frame(dat.sub$Year, dat.sub$ar1)
    #remove NA values (otherwise mmkh won't work)
    ar1. = ar1[complete.cases(ar1),]
    if(nrow(ar1.) > 0 ){
      
      ar1.trends[h, "tau"] = mk.test(ar1.[,2])$estimates[3]
      ar1.trends[h, "pval"] = mk.test(ar1.[,2])$p.value
      
    }
    
    
    ####
    #sk#
    ####
    
    #make a data frame for year x sk combination 
    sk = data.frame(dat.sub$Year, dat.sub$sk)
    #remove NA values (otherwise mmkh won't work)
    sk. = sk[complete.cases(sk),]
    if(nrow(sk.) > 0 ){
      
      sk.trends[h, "tau"] = mk.test(sk.[,2])$estimates[3]
      sk.trends[h, "pval"] = mk.test(sk.[,2])$p.value
      
    }
    
    
    ####
    #kt#
    ####
    
    #make a data frame for year x kt combination 
    kt = data.frame(dat.sub$Year, dat.sub$kurt)
    #remove NA values (otherwise mmkh won't work)
    kt. = kt[complete.cases(kt),]
    if(nrow(kt.) > 0 ){
      
      kt.trends[h, "tau"] = mk.test(kt.[,2])$estimates[3]
      kt.trends[h, "pval"] = mk.test(kt.[,2])$p.value
      
    }
}

#add ews id to dataframe
sd.trends$ews = "sd"
ar1.trends$ews = "ar1"
sk.trends$ews = "sk"
kt.trends$ews = "kt"

#combine trend analysis dataframes in long format for downstream analysis
all.trends = rbind(sd.trends, ar1.trends, sk.trends, kt.trends)

#combine trend analysis datafrmaes in wide format for manuscript
tabs2 = cbind(sd.trends, ar1.trends, sk.trends, kt.trends)

names(tabs2)=c("var_sd", "ews_sd", "tau_sd", "pval_sd", 
               "var_ac", "ews_ac", "tau_ac", "pval_ac", 
               "var_sk", "ews_sk", "tau_sk", "pval_sk",
               "var_kt", "ews_kt", "tau_kt", "pval_kt")

tabs2.1 = tabs2 %>% mutate(across(contains("pval"), ~ p.adjust(.x, method = "fdr", n = length(.x))))


#add columns to tabs2 for formatting final manuscript table
varnames = data.frame(matrix(vector(), nrow = 31, ncol = 2))
names(varnames) = c("vars", "nums")

#create vector with variable names in the order you want them to appear in tabs2
varnames$vars = c(
          "ws1_volwt_Ca", "ws2_volwt_Ca", "ws4_volwt_Ca", "ws5_volwt_Ca", "ws6_volwt_Ca", 
          "ws1_flux_Ca", "ws2_flux_Ca",  "ws4_flux_Ca", "ws5_flux_Ca", "ws6_flux_Ca",
          "ws1_volwt_NO3", "ws2_volwt_NO3", "ws4_volwt_NO3", "ws5_volwt_NO3", "ws6_volwt_NO3", 
          "ws1_flux_NO3", "ws2_flux_NO3",  "ws4_flux_NO3", "ws5_flux_NO3", "ws6_flux_NO3",
          "BIOC", "BION", "FAGR", "ACSA", "lepyear",
          "annamre", "annbtbw", "annbtgw", "annlefl", "annoven", "annrevi")


#add index number for ordering response variables as they will appear in supplemental table
varnames$nums = c(1:31)

#merge varnames with tabs2.1
#rename the first "var" in tabs2 so that it has a unique column name for merging
names(tabs2.1)[1] = "var_sd"
tabs2.2 = merge(tabs2.1, varnames, by.x = "var_sd", by.y = "vars")

#flag where trends are significantly positive or negative for summarizing data
spsd = ifelse(tabs2.2$tau_sd > 0 & tabs2.2$pval_sd <= 0.05, 1, 0)
snsd = ifelse(tabs2.2$tau_sd < 0 & tabs2.2$pval_sd <= 0.05, 1, 0)

spac = ifelse(tabs2.2$tau_ac > 0 & tabs2.2$pval_ac <= 0.05, 1, 0)
snac = ifelse(tabs2.2$tau_ac < 0 & tabs2.2$pval_ac <= 0.05, 1, 0)

spsk = ifelse(tabs2.2$tau_sk > 0 & tabs2.2$pval_sk <= 0.05, 1, 0)
snsk = ifelse(tabs2.2$tau_sk < 0 & tabs2.2$pval_sk <= 0.05, 1, 0)

spkt = ifelse(tabs2.2$tau_kt > 0 & tabs2.2$pval_kt <= 0.05, 1, 0)
snkt = ifelse(tabs2.2$tau_kt < 0 & tabs2.2$pval_kt <= 0.05, 1, 0)

sigdir = cbind(spsd, snsd, spac, snac, spsk, snsk, spkt, snkt)

colSums(sigdir)

#export all.trends and tabs2.2
write.table(all.trends, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\EWS Summary Tables\\alltrends.csv", 
            col.names = T, row.names = F, sep = ",")

write.table(tabs2.2, "C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\EWS Summary Tables\\tabs2.csv", 
            col.names = T, row.names = F, sep = ",")

