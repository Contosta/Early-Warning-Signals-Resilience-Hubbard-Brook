################################################################################################################
################################################################################################################
################################################################################################################

#This script summarizes how early warning signals change over time and as a function of disturbance drivers. 
#It contains most of the code for producing many of the figures and tables presented in the manuscript

#Inputs include hypothesized drivers of climate, precipitation chemistry, and two climate oscillation indices located in the 
#Original Data folder. Inputs also include early warning signals as response variables, both calculated early
#warning signals for each response variable from the Early Warning Signals Tables folder and trends over time 
#in early warning signals (alltrends.csv) in the EWS Summary Tables folder. 

#Outputs include Figure_2, Figure_3, Figure_5, FigureS2, Figure_S3, and Figure_S5 in the manuscript, which are
#written to the EWS Summary Figures folder. Outputs also consist of Kendall's Tau and associated p-values
#that compare early warning signals to hypothesized drivers (EWS_regtab) that forms the basis of Tables S3-S6, as
#well as p- and r2 values derived by fitting hypothesized drivers onto nonmetric multidimensional scaling 
#(NMDS) ordinations of early warning signals. Both output tables are located in the EWS Summary Tables folder.


#code developed by A Contosta
#most recent version 07/25/22 

#changes include code editing and cleaning in preparation for submitting manuscript revisions for peer review

################################################################################################################
################################################################################################################
################################################################################################################

#Initial set-up

#call libraries
library(cowplot)
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtable)
library(Kendall)
library(mblm)
library(psych)
library(stringr)
library(vegan)

#set working directory 
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience")

##########################
#read in driver variables#
##########################

#climate
clim = read.table("Original Data\\anclim.csv", head = T, sep = ",")

#precipitation chemistry
ws6prcp = read.table("Original Data\\ws6annual.csv", head = T, sep = ",")


#######################################
#read in climate osciallation indices#
#NAO and SOI                         #
######################################

nao = read.table("Original Data\\NAO.csv", head = T, sep = ",")
soi = read.table("Original Data\\SOI.csv", head = T, sep = ",")

############################
#read in response variables#
############################

#read in files, create new column identifying response variable, and combine into single dataframe

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\Early Warning Signals Tables")

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
  
  names(dat) = c("timeindex", "Year", "ar1", "sd", "sk", "kurt", "cv", "returnrate",
                     "densratio", "acf1", "var")
  #bind files together
  dat.all = rbind(dat.all, dat)
  
}

#remove first row (all NAs)
dat.all = dat.all[-1, ]

#change wd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience")

#call all.trends table with results of trend analysis
all.trends = read.table("EWS Summary Tables\\alltrends.csv", head = T, sep = ",")

###############################################################
#adjust p-values and flag significance and direction of trends#
###############################################################

#adjust p-values to control for false-discovery rate
all.trends$padj = p.adjust(all.trends$pval, method = "fdr")

#flag adj pvalues that are significant
#when alpha = 0.05
all.trends$sig = ifelse(all.trends$padj <= 0.05, 1, 0)

#flag ews that are increasing
all.trends$posdir = ifelse(all.trends$tau > 0, 1, 0)

#flag ews that are decreasing
all.trends$negdir = ifelse(all.trends$tau < 0, 1, 0)

#flag where adj pvalues are significant AND trends are increasing
all.trends$sigposdir = ifelse(all.trends$padj <= 0.05 & all.trends$tau > 0, 1, 0)

#flag where adjpvalues are significant AND taus are decreasing
all.trends$signegdir = ifelse(all.trends$padj <= 0.05 & all.trends$tau < 0, 1, 0)

#summarize by variable where one or both of these condtions occur
all.trends.1 = data.table(all.trends)

trendsum = all.trends.1[ , list(sig = sum(sig, na.rm = T), 
                              posdir = sum(posdir, na.rm = T),
                              negdir = sum(negdir, na.rm = T),
                              sigposdir = sum(sigposdir, na.rm = T),
                              signegdir = sum(signegdir, na.rm = T)), by = var]


################################################################################################################
################################################################################################################
################################################################################################################

#Plot change over time in EWS

###########################################
#format dat.all to prepare it for graphing#
###########################################

#add column to delineate variables into groups
dat.all$cat = ifelse(dat.all$var == "annamre" | dat.all$var == "annbtbw" | dat.all$var =="annbtgw"  |
                       dat.all$var == "annlefl" | dat.all$var == "annoven" | dat.all$var == "annrevi", "bird",
                       ifelse(dat.all$var == "lepyear", "bug",
                            ifelse(dat.all$var == "ACSA" | dat.all$var == "FAGR", "bai",
                                   ifelse(dat.all$var == "BIOC" | dat.all$var == "BION", "mic",
                                          ifelse(dat.all$var == "ws1_volwt_Ca" | dat.all$var == "ws2_volwt_Ca" | dat.all$var == "ws4_volwt_Ca" |
                                                   dat.all$var == "ws5_volwt_Ca" | dat.all$var == "ws6_volwt_Ca", "volwt_Ca",
                                                 ifelse(dat.all$var == "ws1_volwt_NO3" | dat.all$var == "ws2_volwt_NO3" | dat.all$var == "ws4_volwt_NO3" |
                                                          dat.all$var == "ws5_volwt_NO3" | dat.all$var == "ws6_volwt_NO3", "volwt_NO3",
                                                        ifelse(dat.all$var == "ws1_flux_Ca" | dat.all$var == "ws2_flux_Ca" | dat.all$var == "ws4_flux_Ca" |
                                                                 dat.all$var == "ws5_flux_Ca" | dat.all$var == "ws6_flux_Ca", "flux_Ca",
                                                               ifelse(dat.all$var == "ws1_flux_NO3" | dat.all$var == "ws2_flux_NO3" | dat.all$var == "ws4_flux_NO3" |
                                                                        dat.all$var == "ws5_flux_NO3" | dat.all$var == "ws6_flux_NO3", "flux_NO3",
                                                   NA))))))))
  
#subset to only include columns used in plotting
dat.sub = dat.all[ , c("cat", "var", "Year", "sd", "ar1", "sk", "kurt")]
names(dat.sub) = c("cat", "var", "Year", "sd", "ac", "sk", "kt")

##############################
##############################
#Make Figure 2 for Manuscript#
##############################
##############################

#note that the figure has 16 panels
#the process consists of
#making individual panels with data
#making dummy panels for customizing the legend
#extracting legend parameters
#arranging panels into grid for plotting
#plotting grid with panels and legend and export to folder

#####################
#make Ca_volwt plots#
#####################

#subset dat.sub to include just Ca_volwt data 
dat.Ca_volwt = dat.sub[dat.sub$cat == "volwt_Ca", ]

Ca_volwt.sd.plot = ggplot(data = dat.Ca_volwt, aes(x = Year, y= sd, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = "Standard Deviation") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.0", "0.1", "0.2"), breaks = c(0, 0.1, 0.2)) +
  scale_size_manual(values = c(1.5,1.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(1,1,1,2,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.Ca_volwt[is.na(dat.Ca_volwt$var) == F, ]$Year) - 1), y = (0.98 * max(dat.Ca_volwt$sd, na.rm = T)), 
           label = "a", fontface = "bold")

Ca_volwt.ac.plot = ggplot(data = dat.Ca_volwt, aes(x = Year, y= ac, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = "Autocorrelation") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.3", "0.6", "0.9"), breaks = c(0.3, 0.6, 0.9)) +
  scale_size_manual(values = c(0.5,1.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(1,1,1,1,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.Ca_volwt[is.na(dat.Ca_volwt$var) == F, ]$Year) - 1), y = 1.01, 
           label = "b", fontface = "bold")

Ca_volwt.sk.plot = ggplot(data = dat.Ca_volwt, aes(x = Year, y= sk, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = "Skewness") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.4", "1.2", "2.0"), breaks = c(0.4, 1.2, 2.0)) +
  scale_size_manual(values = c(0.5,1.5,1.5,1.5,0.5)) +
  scale_linetype_manual(values = c(1,1,1,1,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.Ca_volwt[is.na(dat.Ca_volwt$var) == F, ]$Year) - 1), y = (0.98 * max(dat.Ca_volwt$sk, na.rm = T)), 
           label = "c", fontface = "bold")

Ca_volwt.kt.plot = ggplot(data = dat.Ca_volwt, aes(x = Year, y= kt, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = "Year", y = "Kurtosis") +
  theme(plot.margin = margin(0.1, 0.4, -11, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("3.0", "5.5", "8.0"), breaks = c(3.0, 5.5, 8.0)) +
  scale_x_continuous(
    labels = c("1995", "2005", "2015"), breaks = c(1995, 2005, 2015)) +
  scale_size_manual(values = c(0.5,1.5,0.5,1.5,0.5)) +
  scale_linetype_manual(values = c(1,1,1,1,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.Ca_volwt[is.na(dat.Ca_volwt$var) == F, ]$Year) - 1), y = (0.98 * max(dat.Ca_volwt$kt, na.rm = T)), 
           label = "d", fontface = "bold")

#####################
#make NO3_flux plots#
####################

#subset dat.sub to include just NO3_flux data 
dat.NO3_flux = dat.sub[dat.sub$cat == "flux_NO3", ]

NO3_flux.sd.plot = ggplot(data = dat.NO3_flux, aes(x = Year, y= sd, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -15))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("3.0", "7.5", "12.0"), breaks = c(3000, 7500, 12000)) +
  scale_size_manual(values = c(0.5,0.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(2,2,2,1,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.NO3_flux[is.na(dat.NO3_flux$var) == F, ]$Year) + 1), y = 14500, 
           label = "e", fontface = "bold")

NO3_flux.ac.plot = ggplot(data = dat.NO3_flux, aes(x = Year, y= ac, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.0", "0.3", "0.6"), breaks = c(0, 0.3, 0.6)) +
  scale_size_manual(values = c(0.5,1.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(1,2,1,2,2), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.NO3_flux[is.na(dat.NO3_flux$var) == F, ]$Year) - 1), y = (0.98 * max(dat.NO3_flux$ac, na.rm = T)), 
           label = "f", fontface = "bold")

NO3_flux.sk.plot = ggplot(data = dat.NO3_flux, aes(x = Year, y= sk, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("1.0", "2.0", "3.0"), breaks = c(1, 2, 3)) +
  scale_size_manual(values = c(0.5,1.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(1,1,1,1,1), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.NO3_flux[is.na(dat.NO3_flux$var) == F, ]$Year) - 1), y = (0.98 * max(dat.NO3_flux$sk, na.rm = T)), 
           label = "g", fontface = "bold")

NO3_flux.kt.plot = ggplot(data = dat.NO3_flux, aes(x = Year, y= kt, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, -11, -15))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("5.0", "10.0", "15.5"), breaks = c(5, 10, 15)) +
  scale_x_continuous(
    labels = c("1995", "2005", "2015"), breaks = c(1995, 2005, 2015)) +
  scale_size_manual(values = c(1.5,1.5,1.5,1.5,1.5)) +
  scale_linetype_manual(values = c(2,1,1,1,1), labels = c("WS1", "WS2", "WS4", "WS5", "WS6"))+
  scale_colour_manual(values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2"),
                      labels = c("WS1", "WS2", "WS4", "WS5", "WS6"), guide = guide_legend(ncol = 1))+
  annotate("text", x = (min(dat.NO3_flux[is.na(dat.NO3_flux$var) == F, ]$Year) - 1), y = (0.98 * max(dat.NO3_flux$kt, na.rm = T)), 
           label = "h", fontface = "bold")


#################
#make bird plots#
#################

#subset dat.sub to include just bird data for btbw, lefl, and oven
dat.bird = dat.sub[dat.sub$var == "annbtbw" | dat.sub$var == "annlefl" | dat.sub$var == "annoven", ]

birds.sd.plot = ggplot(data = dat.bird, aes(x = Year, y= sd, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "blank") +
  scale_colour_manual(values = c("dodgerblue4", "darkolivegreen4", "tomato4"),
                      labels = c("BTBW","LEFL", "OVEN"))+
  guides(color = guide_legend(ncol=1))+
  scale_y_continuous(
    labels = c("1.5", "5.5", "9.5"), breaks = c(1.5, 5.5, 9.5)) +
  scale_size_manual(values = c(0.5,1.5,1.5)) +
  scale_linetype_manual(values = c(2,2,1), labels = c("BTBW", "LEFL", "OVEN"))+
  annotate("text", x = (min(dat.bird[is.na(dat.bird$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bird$sd, na.rm = T)), 
           label = "i", fontface = "bold")

birds.ac.plot = ggplot(data = dat.bird, aes(x = Year, y= ac, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_colour_manual(values = c("dodgerblue4", "darkolivegreen4", "tomato4"),
                      labels = c("BTBW", "LEFL", "OVEN"))+
  guides(color = guide_legend(ncol=1))+
  scale_y_continuous(
    labels = c("0.2", "0.5", "0.8"), breaks = c(0.2, 0.5, 0.8)) +
  scale_size_manual(values = c(1.5,0.5,1.5)) +
  scale_linetype_manual(values = c(1,2,1), labels = c("BTBW", "LEFL", "OVEN"))+
  annotate("text", x = (min(dat.bird[is.na(dat.bird$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bird$ac, na.rm = T)), 
           label = "j", fontface = "bold")

birds.sk.plot = ggplot(data = dat.bird, aes(x = Year, y= sk, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_colour_manual(values = c("dodgerblue4","darkolivegreen4", "tomato4"),
                      labels = c("BTBW", "LEFL", "OVEN"))+
  guides(color = guide_legend(ncol=1))+
  scale_y_continuous(
    labels = c("0.3", "0.9", "1.4"), breaks = c(0.3, 0.9, 1.4)) +
  scale_size_manual(values = c(0.5,1.5,1.5)) +
  scale_linetype_manual(values = c(1,1,2), labels = c("BTBW", "LEFL", "OVEN"))+
  annotate("text", x = (min(dat.bird[is.na(dat.bird$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bird$sk, na.rm = T)), 
           label = "k", fontface = "bold")

birds.kt.plot = ggplot(data = dat.bird, aes(x = Year, y= kt, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, -11, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  scale_colour_manual(values = c("dodgerblue4", "darkolivegreen4", "tomato4"),
                      labels = c("BTBW", "LEFL", "OVEN"))+
  guides(color = guide_legend(ncol=1))+
  scale_y_continuous(
    labels = c("2.0", "4.0", "6.0"), breaks = c(2,4,6)) +
  scale_x_continuous(
    labels = c("1992", "2002", "2012"), breaks = c(1992, 2002, 2012),
    limits = c(min(dat.bird[is.na(dat.bird$var) == F, ]$Year -1, na.rm = T), 
               max(dat.bird[is.na(dat.bird$var) == F, ]$Year, na.rm = T))) +
  scale_size_manual(values = c(1.5,0.5,0.5)) +
  scale_linetype_manual(values = c(1,1,2), labels = c("BTBW", "LEFL", "OVEN"))+
  annotate("text", x = (min(dat.bird[is.na(dat.bird$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bird$kt, na.rm = T)), 
           label = "l", fontface = "bold")

#################
#make bai plots #
#################

#subset dat.sub to include just bai data 
dat.bai = dat.sub[dat.sub$cat == "bai", ]

bai.sd.plot = ggplot(data = dat.bai, aes(x = Year, y= sd, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "blank") +
  scale_y_continuous(
    labels = c("0.7", "1.3", "1.9"), breaks = c(0.7, 1.3, 1.9)) +
  scale_size_manual(values = c(1.5,0.5)) +
  scale_linetype_manual(values = c(1,1), labels = c("ACSA", "FAGR"))+
  scale_colour_manual(values = c("orangered3", "seagreen4"),
                      labels = c("ACSA", "FAGR"))+
  annotate("text", x = (min(dat.bai[is.na(dat.bai$var) == F, ]$Year) + 0.5), y = (0.99 * max(dat.bai$sd, na.rm = T)), 
           label = "m", fontface = "bold")

bai.ac.plot = ggplot(data = dat.bai, aes(x = Year, y= ac, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.0", "0.3", "0.6"), breaks = c(0, 0.3, 0.6)) +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_linetype_manual(values = c(1,1), labels = c("ACSA", "FAGR"))+
  scale_colour_manual(values = c("orangered3", "seagreen4"),
                      labels = c("Sugar Maple", "American Beech"))+
  annotate("text", x = (min(dat.bai[is.na(dat.bai$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bai$ac, na.rm = T)), 
           label = "n", fontface = "bold")

bai.sk.plot = ggplot(data = dat.bai, aes(x = Year, y= sk, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, 0, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("0.2", "0.6", "1.0"), breaks = c(0.2, 0.6, 1.0)) +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_linetype_manual(values = c(2,1), labels = c("ACSA", "FAGR"))+
  scale_colour_manual(values = c("orangered3", "seagreen4"),
                      labels = c("Sugar Maple", "American Beech"))+
  annotate("text", x = (min(dat.bai[is.na(dat.bai$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bai$sk, na.rm = T)), 
           label = "o", fontface = "bold")

bai.kt.plot = ggplot(data = dat.bai, aes(x = Year, y= kt, color = var)) +
  geom_line(aes(linetype = var, size = var))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  labs(x = " ", y = " ") +
  theme(plot.margin = margin(0.1, 0.4, -11, -10))+
  theme(axis.title.y=element_text(face = "bold"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(
    labels = c("2.2", "3.2", "4.2"), breaks = c(2.2, 3.2, 4.2)) +
  scale_x_continuous(
    labels = c("1985", "1995", "2005"), breaks = c(1985, 1995, 2005)) +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_linetype_manual(values = c(2,1), labels = c("ACSA", "FAGR"))+
  scale_colour_manual(values = c("orangered3", "seagreen4"),
                      labels = c("Sugar Maple", "American Beech"))+
  annotate("text", x = (min(dat.bai[is.na(dat.bai$var) == F, ]$Year) - 1), y = (0.98 * max(dat.bai$kt, na.rm = T)), 
           label = "p", fontface = "bold")


############################################################
#make dummy dataframes and panels                          #    
#for customizing the legend                                #  
############################################################

#########################################################
#make dummy dataframe with all ecosystem state variables#
#########################################################

vars = c("WS 1",  "WS 2", "WS 4", "WS 5", "WS 6", "NA", "BTBW", "LEFL", "OVEN",  "Sugar Maple", 
         "American Beech")
years = c(1960:2010)

df = data.frame(matrix(vector(), nrow = length(vars) * length(years), ncol = 3))
names(df) = c("vars", "years", "kt")

df$vars = rep(vars, each = length(years))
df$years = years
df$kt = runif(length(df)) 

#factorize and relevel the early warning signals names so they
#show up in the right order on the plot
df$vars = factor(df$vars, levels = unique(df$vars), 
                 labels =  unique(df$vars))

#############################################################
#create dummy dataframe for trend significance and direction#
#############################################################

tsig = c("Significant Increase", "Significant Decrease", 
         "Non-Signficant Increase", "Non-Signficant Decrease", "NA", "NA")
years = c(1960:2010)

dfsig = data.frame(matrix(vector(), nrow = length(tsig) * length(years), ncol = 3))
names(dfsig) = c("tsig", "years", "kt")

dfsig$tsig = rep(tsig, each = length(years))
dfsig$years = years
dfsig$kt = runif(length(tsig)) 

#factorize and relevel the early warning signals names so they
#show up in the right order on the plot
dfsig$tsig = factor(dfsig$tsig, levels = unique(dfsig$tsig), 
                 labels =  unique(dfsig$tsig))


###################
#make dummy panels#
###################

#make df plot
df.plot = ggplot(data = df, aes(x = years, y= kt, color = vars))+
  geom_line(aes(size = vars))+
  theme(legend.position = "top")+
  theme(legend.title = element_text(
    face = "bold"))+
  theme(legend.key = element_rect(fill = NA, color = NA))+
  theme(legend.key.width = unit(1, "line")) +
  scale_colour_manual(name = "Ecosystem State Variable",
                      values = c("firebrick1", "deepskyblue1", "olivedrab3", 
                                 "mediumpurple", "goldenrod2", "white",
                                 "dodgerblue4", "darkolivegreen4", "tomato4", 
                                 "orangered3", "seagreen4"
                      ),
                      labels = c("Watershed 1",  "Watershed 2", "Watershed 4", "Watershed 5", "Watershed 6", " ",
                                 "Black-Throated Blue Warbler", "Least Flycatcher", "Ovenbird", 
                                 "Sugar Maple", "American Beech"
                      ))+
  scale_size_manual(name = "Ecosystem State Variable", values = rep(1.5, 11), 
                    labels = c("Watershed 1",  "Watershed 2", "Watershed 4", "Watershed 5", "Watershed 6", " ",
                               "Black-Throated Blue Warbler", "Least Flycatcher", "Ovenbird", 
                               "Sugar Maple", "American Beech"))+
  
  guides(colour = guide_legend(title.vjust = 1, ncol=2, title.position="top"))


#make signficant and direction plot
dfsig.plot = ggplot(data = dfsig, aes(x = years, y= kt, color = tsig))+
  geom_line(aes(linetype = tsig, size = tsig))+
  theme(legend.title = element_text(
    face = "bold"))+
  theme(legend.direction = "vertical")+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(fill = NA, color = NA))+
  theme(legend.key.width = unit(5, "line")) +
  scale_colour_manual(name = "Trend Signficance and Direction",
                      values = c("gray30", "gray30", "gray30", "gray30", "white", "white"),
                      labels = c("Significant Increase", "Significant Decrease",
                                 "Non-Significant Increase", "Non-Significant Decrease", " ", " "))+
  scale_linetype_manual(name = "Trend Signficance and Direction",
                        values = c(2,1,2,1,1,1), 
                        labels = c("Significant Increase", "Significant Decrease",
                                   "Non-Significant Increase", "Non-Significant Decrease", " ", " "))+
  scale_size_manual(name = "Trend Signficance and Direction", values = c(1.5, 1.5, 0.5, 0.5, 1.5, 1.5), 
                    labels = c("Significant Increase", "Significant Decrease",
                               "Non-Significant Increase", "Non-Significant Decrease", " ", " "))+
  guides(colour = guide_legend(title.vjust = 1, ncol=1, title.position="top"))

############################################
#extract legend parameters from dummy plots#
############################################

# Create user-defined function, which extracts legends from ggplots

extract_df <- function(df.plot) {
  step1 <- ggplot_gtable(ggplot_build(df.plot))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}


extract_dfsig <- function(dfsig.plot) {
  step1 <- ggplot_gtable(ggplot_build(dfsig.plot))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}



# Apply user-defined function to extract legend
#write function
extract_legend<-function(ggplot_object){
  tmp <- ggplot_gtable(ggplot_build(ggplot_object))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

df_legend = extract_legend(df.plot)
dfsig_legend = extract_legend(dfsig.plot)

##############################################
#create text grobs for labeling panel columns#
##############################################

Ca_conc = text_grob("Ca Concentration", face = "bold", hjust = 0.4)
NO3_flux = text_grob(expression(paste(bold("NO"[3]*" Flux"))), face = "bold", hjust = 0.25)
bird = text_grob("Bird Abundance", face = "bold", hjust = 0.35)
bai = text_grob("Basal Area Increment", face = "bold", hjust = 0.45)

########################################
#arrange panels into grid for plotting #
#and export to folder                  #
########################################

#add height onto dfsig_legend so that it aligns with df_legend in plot
dfsig_legend_2 = gtable_add_rows(dfsig_legend,  unit(0.61, "cm"), pos = -1)

#call pdf options
pdf.options(width= 6.5, height= 8.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="EWS Summary Figures//Figure_2.pdf")

grid.arrange(

  arrangeGrob(
  Ca_conc, NO3_flux, bird, bai, ncol = 4, nrow = 1),
  
  arrangeGrob(
    Ca_volwt.sd.plot, NO3_flux.sd.plot, birds.sd.plot, bai.sd.plot, 
    Ca_volwt.ac.plot, NO3_flux.ac.plot, birds.ac.plot, bai.ac.plot, 
    Ca_volwt.sk.plot, NO3_flux.sk.plot, birds.sk.plot, bai.sk.plot, 
    Ca_volwt.kt.plot, NO3_flux.kt.plot, birds.kt.plot, bai.kt.plot, 
    ncol = 4, nrow = 4, bottom = text_grob("Year", face = "bold", vjust = 1),
    padding = unit(1, "cm")), 
  

  arrangeGrob(df_legend, dfsig_legend_2, 
              nrow = 1, ncol = 2, as.table = F), 
  
  left = text_grob(" ", rot = 90, hjust= -12.5),
  heights = c(0.5, 10, 3)
  
)
dev.off()

##################################################################################################
##################################################################################################
##################################################################################################

#Examine relationships between hypothesized drivers and early warning signal
#response variables

###################################################################
#examine relationship between decadal climate oscillations and EWS#
###################################################################

#combine nao and soi
clos = merge(nao, soi, by.x = "Year", by.y = "YEAR", all.x = T, all.y = T)

#subset columns and rename
clo = clos[ , c("Year", "NAO", "AVG")]
names(clo) = c("year", "nao", "soi")

#combine clo with dat.sub
datclo = merge(clo, dat.sub, by.x = "year", by.y = "Year")

#sort by variable and year
datclo = datclo[order(datclo$var, datclo$year), ]

#subset to only include response variables of interest
dat.res = datclo[ , c("sd", "ac", "sk", "kt")]

#subset to only include driver variables of interest
dat.driv = datclo[ , c("soi", "nao")]

#select start and end rows to run correlations for each variable in turn
#by omitting duplicate variables from the beginning and the end of datclo
start.rows = !duplicated(datclo$var)
end.rows = !duplicated(datclo$var, fromLast = TRUE)

#extract original row indexes associated with the start and end of unique state variable
sr.ind = seq_along(datclo$var)[start.rows]
er.ind = seq_along(datclo$var)[end.rows]

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#create containers to hold results of correlation tests 
cat = list()
var = list()
X_sds = list()
Y_sds = list()
pval = list()
tau = list()

#analyze correlations between hypothesized driver and response variables
for (h in seq(1,length(sr.ind))) {
  
  #subset dat.driv
  driv.sub = dat.driv[sr.ind[h]:er.ind[h], ]
  
  #subset dat.res
  res.sub = dat.res[sr.ind[h]:er.ind[h], ]
  
  for(i in 1:ncol(driv.sub)){
    for(j in 1:ncol(res.sub)){
      
      driv = driv.sub[[i]]
      res = res.sub[[j]]
      
      df = data.frame(driv, res)
      
      cat = c(cat, unique(datclo[sr.ind[h]:er.ind[h], ]$cat))
      var = c(var, unique(datclo[sr.ind[h]:er.ind[h], ]$var))
      
      X_sds = append(X_sds, colnames(driv.sub)[i])
      Y_sds = append(Y_sds, colnames(res.sub)[j])
      pval = append(pval, Kendall(df$res, df$driv)[[2]][1])
      tau = append(tau, Kendall(df$res, df$driv)[[1]][1])
      
    }
  }
}

#extract correlation results
cat = unlist(cat)
var = unlist(var)
xvs = unlist(X_sds)
yvs = unlist(Y_sds)
pvs = as.numeric(unlist(pval))
tas = as.numeric(unlist(tau))

#combine results into dataframe  
allclo = data.frame(cbind(pvs, tas))
names(allclo) = c("pval", "tau")
allclo$cat = cat
allclo$var = var
allclo$xvs = xvs
allclo$yvs = yvs

#adjust p-values for multiple comparisons
allclo$padj = p.adjust(allclo$pval, method = "fdr")

#plot the distribution of p-values showing the correlation between climate oscillations and EWS
#and write to folder

#redefine and rename variables for plotting in the order you want
allclo$yvs2 = ifelse(allclo$yvs == "sd", "A",
                   ifelse(allclo$yvs == "ac", "B",
                          ifelse(allclo$yvs == "sk", "C",
                                 ifelse(allclo$yvs == "kt", "D", NA
                                 ))))

allclo$yvs2 = factor(allclo$yvs2, levels = c("A", "B", "C", "D"), 
                   labels = c("SD", "AC", "SK", "KT"))


#call pdf options
pdf.options(width= 4.5, height= 4.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="EWS Summary Figures\\Figure_S5.pdf")

plotpval = ggplot(allclo, aes(x = padj))+
  xlim(c(0,1))+
  geom_density(alpha = 0.5)+
  labs(x = "Adjusted P-Value", y = "Density") +
  ylab("Density")
plotpval + facet_wrap(vars(yvs2))

dev.off()

########################################################################################
#examine relationship between climate and precipitation chemistry drivers and responses#
########################################################################################

#combine dat.sub with clim and prec
dat.sub.c = merge(dat.sub, clim, by.x = "Year", by.y = "WaterYear")
dat.cp = merge(dat.sub.c, ws6prcp, by.x = "Year", by.y = "WaterYear")

#sort by variable and year
dat.cp = dat.cp[order(dat.cp$var, dat.cp$Year), ]

#subset to only include response variables of interest
dat.res = dat.cp[ , c("sd", "ac", "sk", "kt")]

#subset to only include driver variables of interest
dat.driv = dat.cp[ , c("TMEAN", "PRCP", "volwt_pH", "ANC", "NH4_flux",
                   "SO4_flux", "NO3_flux", "totN")]

#select start and end rows to run correlations for each variable in turn
#by omitting duplicate variables from the beginning and the end of dat.cp
start.rows = !duplicated(dat.cp$var)
end.rows = !duplicated(dat.cp$var, fromLast = TRUE)

#extract original row indexes associated with the start and end of unique state variable
sr.ind = seq_along(dat.cp$var)[start.rows]
er.ind = seq_along(dat.cp$var)[end.rows]

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#create containers to hold results of correlation tests 
cat = list()
var = list()
X_sds = list()
Y_sds = list()
pval = list()
tau = list()

#analyze correlations between hypothesized driver and response variables
for (h in seq(1,length(sr.ind))) {
  
  #subset dat.driv
  driv.sub = dat.driv[sr.ind[h]:er.ind[h], ]
  
  #subset dat.res
  res.sub = dat.res[sr.ind[h]:er.ind[h], ]
  

for(i in 1:ncol(driv.sub)){
  for(j in 1:ncol(res.sub)){
    
    driv = driv.sub[[i]]
    res = res.sub[[j]]
    
    df = data.frame(driv, res)
    
    cat = append(cat, unique(dat.cp[sr.ind[h]:er.ind[h], ]$cat))
    var = append(var, unique(dat.cp[sr.ind[h]:er.ind[h], ]$var))
    
    X_sds = append(X_sds, colnames(driv.sub)[i])
    Y_sds = append(Y_sds, colnames(res.sub)[j])
    pval = append(pval, Kendall(df$res, df$driv)[[2]][1])
    tau = append(tau, Kendall(df$res, df$driv)[[1]][1])
    
  }
}
}
  
#extract correlation results
cat = unlist(cat)
var = unlist(var)
xvs = unlist(X_sds)
yvs = unlist(Y_sds)
pvs = as.numeric(unlist(pval))
tas = as.numeric(unlist(tau))

#combine results into dataframe  
allregs = data.frame(cbind(pvs, tas))
names(allregs) = c("pval", "tau")
allregs$cat = cat
allregs$var = var
allregs$xvs = xvs
allregs$yvs = yvs

#adjust p-values for multiple comparisons
allregs$pval = p.adjust(allregs$pval, method = "fdr")

#add columns to dr.cors.1 for formatting final manuscript table
varnames = data.frame(matrix(vector(), nrow = 31, ncol = 2))
names(varnames) = c("vars", "nums")

#create vector with variable names in the order you want them to appear in tabs2
varnames$vars = c("ws1_volwt_Ca", "ws2_volwt_Ca", "ws4_volwt_Ca", "ws5_volwt_Ca", "ws6_volwt_Ca", 
                  "ws1_volwt_NO3", "ws2_volwt_NO3", "ws4_volwt_NO3", "ws5_volwt_NO3", "ws6_volwt_NO3", 
                  "ws1_flux_Ca", "ws2_flux_Ca",  "ws4_flux_Ca", "ws5_flux_Ca", "ws6_flux_Ca",
                  "ws1_flux_NO3", "ws2_flux_NO3",  "ws4_flux_NO3", "ws5_flux_NO3", "ws6_flux_NO3",
                  "BIOC", "BION",
                  "lepyear", "FAGR", "ACSA",
                  "annamre", "annbtbw", "annbtgw", "annlefl", "annoven", "annrevi")


#add index number for ordering response variables as they will appear in supplemental table
varnames$nums = c(1:31)

#merge varnames with allregs
tabs3 = merge(allregs, varnames, by.x = "var", by.y = "vars")

#reconfigure table to make it easier to put in publication format
#make separate dataframes for each driver variable

#make unique ID for each variable by EWS combination for downstream dataframe merging
tabs3$id = paste0(tabs3$var, tabs3$yvs)

regTMEAN = tabs3[tabs3$xvs == "TMEAN", c("id", "var", "nums", "yvs", "pval", "tau")]
names(regTMEAN) = c("id", "var", "nums", "yvs",  "pTMEAN", "tauTMEAN")

regPRCP = tabs3[tabs3$xvs == "PRCP", c("id", "pval", "tau")]
names(regPRCP) = c("id", "pPRCP", "tauPRCP")

regpH = tabs3[tabs3$xvs == "volwt_pH", c("id", "pval", "tau")]
names(regpH) = c("id", "ppH", "taupH")

regANC = tabs3[tabs3$xvs == "ANC", c("id", "pval", "tau")]
names(regANC) = c("id", "pANC", "tauANC")

regSO4 = tabs3[tabs3$xvs == "SO4_flux", c("id", "pval", "tau")]
names(regSO4) = c("id", "pSO4", "tauSO4")

regNO3 = tabs3[tabs3$xvs == "NO3_flux", c("id", "pval", "tau")]
names(regNO3) = c("id", "pNO3", "tauNO3")

regNH4 = tabs3[tabs3$xvs == "NH4_flux", c("id", "pval", "tau")]
names(regNH4) = c("id", "pNH4", "tauNH4")

regtotN = tabs3[tabs3$xvs == "totN", c("id", "pval", "tau")]
names(regtotN) = c("id", "ptotN", "tautotN")

#merge tables together
TMEAN_PRCP = merge(regTMEAN, regPRCP, by.x = "id", by.y = "id")
PRCP_pH = merge(TMEAN_PRCP, regpH, by.x = "id", by.y = "id")
pH_ANC = merge(PRCP_pH, regANC, by.x = "id", by.y = "id")
ANC_SO4 = merge(pH_ANC, regSO4, by.x = "id", by.y = "id")
SO4_NO3 = merge(ANC_SO4, regNO3, by.x = "id", by.y = "id")
NO3_NH4 = merge(SO4_NO3, regNH4, by.x = "id", by.y = "id")
regtab = merge(NO3_NH4, regtotN, by.x = "id", by.y = "id")

#write regtab table to file
write.table(regtab, "EWS Summary Tables\\EWS_regtab.csv", row.names = F, col.names = T, sep = ",")

###################################################
#Make figure of hypothesized drivers for manuscript
#Figure S2
#####################################################

#call pdf options
pdf.options(width= 4.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="EWS Summary Figures\\Figure_S2.pdf")

#set graphical parameters
par(mfrow = c (2,1), mar = c(0.5,0.5,0.5,0.5), oma = c(5,5,5,5))

plot(TMEAN ~ WaterYear, data = clim, type = "l",
     col = "red", ylab = " ", xlab = " ", axes = F)
axis(2, at = c(4.5, 6.5, 8.5, 10.5))
par(new = T)
plot(PRCP ~ WaterYear, data = clim, type = "l",
     axes = F, col = "dodgerblue", ylab = " ", xlab = " ", lty = 2)
axis(4, at = c(600, 900, 1200, 1500, 1800))
box(lty = 1)
mtext(side = 2, expression("Mean Annual Temperature ("*degree~C*")"), line = 2.5)
mtext(side = 4, "Total Annual Precipitation (mm)", line = 2.5)
legend("topleft", bg = "white", bty = "n", cex = 0.8, "A")
legend("top", bg = "white", 
       lty = c(1,2), bty = "n", cex = 0.8,
       col = c("red", "dodgerblue"),
       legend = c("MAT", "PRCP"))

plot(volwt_pH ~ WaterYear, data = ws6prcp, type = "l",
     col = "forestgreen", axes = F, 
     ylab = " ", xlab = "Year", ylim = c(4, 5))
axis(2, at = c(4.0, 4.3, 4.6, 4.9))
par(new = T)
plot(totN ~ WaterYear, data = ws6prcp, type = "l",
     axes = F, col = "gray30", lty = 2,
     xlab = " ", ylab = " ")
axis(1, at = c(1965, 1980, 1995, 2010))
axis(4, at = c(8000, 16000, 24000, 32000), lab = c(8, 16, 24, 32))
box(lty = 1)
mtext(side = 1, "Year", line = 2.5)
mtext(side = 2, "Mean Annual Precipitation pH", line = 2.5)
mtext(side = 4, "Total Annual Nitrogen Flux (kg / ha)", line = 2.5)
legend("topleft", bg = "white", bty = "n", cex = 0.8, "B")
legend("top", bg = "white", 
       lty = c(1,2), bty = "n", cex = 0.8,
       col = c("forestgreen", "gray30"),
       legend = c("pH", "N"))

dev.off()

#####################################################################
#####################################################################
#Make Phase-State Plots Showing Relationship Between Drivers and EWS#
#####################################################################
#####################################################################

####################################################
# Figure 3 Ca Concentration versus Precipitation pH#
####################################################

#subset dat.cp to include just Ca_volwt data 
dat.cp.Ca_volwt = dat.cp[dat.cp$cat == "volwt_Ca", ]

#rename variables for plotting
dat.cp.Ca_volwt$var = factor(dat.cp.Ca_volwt$var, levels = c("ws1_volwt_Ca",
                                                           "ws2_volwt_Ca",
                                                           "ws4_volwt_Ca",
                                                           "ws5_volwt_Ca",
                                                           "ws6_volwt_Ca"), 
                            labels = c("Watershed 1", "Watershed 2", "Watershed 4", "Watershed 5", "Watershed 6"))

Ca.volwt.pH.sd.plot = ggplot(data = dat.cp.Ca_volwt, aes(x = volwt_pH, y = sd, color = Year))+
  geom_point()+
  labs(x = " ", y = "Standard Deviation") +
  theme(axis.title.x=element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5)) +
  facet_grid(. ~ var)
Ca.volwt.pH.sd.plot$labels$colour = "Year"


Ca.volwt.pH.ac.plot = ggplot(data = dat.cp.Ca_volwt, aes(x = volwt_pH, y = ac, color = Year))+
  geom_point()+
  labs(x = " ", y = "Autocorrelation") +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
 facet_grid(. ~ var)
Ca.volwt.pH.ac.plot$labels$colour = "Year"


Ca.volwt.pH.sk.plot = ggplot(data = dat.cp.Ca_volwt, aes(x = volwt_pH, y = sk, color = Year))+
  geom_point()+
  labs(x = " ", y = "Skewness") +
  scale_y_continuous(
    labels = c("0.6", "1.2", "1.8", "2.4"), breaks = c(0.6, 1.2, 1.8, 2.4)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
 facet_grid(. ~ var)
Ca.volwt.pH.sk.plot$labels$colour = "Year"


Ca.volwt.pH.kt.plot = ggplot(data = dat.cp.Ca_volwt, aes(x = volwt_pH, y = kt, color = Year))+
  geom_point()+
  labs(x = " ", y = "Kurtosis") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  facet_grid(. ~ var)
Ca.volwt.pH.kt.plot$labels$colour = "Year" 

#setwd for figure export
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\EWS Summary Figures")

#call pdf options
pdf.options(width= 6.5, height= 7.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="Figure_3.pdf")


#Combine plots
Ca.volwt.pH = plot_grid(Ca.volwt.pH.sd.plot, Ca.volwt.pH.ac.plot, Ca.volwt.pH.sk.plot, Ca.volwt.pH.kt.plot,
              labels = c("a", "b", "c", "d"), axis = 'lb', ncol = 1, rel_heights = c(1,1))

#create common x and y labels
y.grob <- textGrob("EWS of Stream Ca Concentration", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Annual Precipitation pH", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))

#add to plot
grid.arrange(arrangeGrob(Ca.volwt.pH, left = y.grob, bottom = x.grob,
                         vp=viewport(width=1, height=1, clip = TRUE)))

dev.off()

#################################################
# Figure S3 NO3 Flux versus Precipitation N Flux#
#################################################

#subset dat.cp to include just NO3_flux data 
dat.cp.NO3_flux = dat.cp[dat.cp$cat == "flux_NO3", ]

#convert total precipitation N flux to from g / ha to kg / ha 
dat.cp.NO3_flux$totN_kg = dat.cp.NO3_flux$totN / 1000

#rename variables for plotting
dat.cp.NO3_flux$var = factor(dat.cp.NO3_flux$var, levels = c("ws1_flux_NO3",
                                                           "ws2_flux_NO3",
                                                           "ws4_flux_NO3",
                                                           "ws5_flux_NO3",
                                                           "ws6_flux_NO3"), 
                            labels = c("Watershed 1", "Watershed 2", "Watershed 4", "Watershed 5", "Watershed 6"))

NO3.flux.N.sd.plot = ggplot(data = dat.cp.NO3_flux, aes(x = totN_kg, y = sd / 1000, color = Year))+
  geom_point()+
  labs(x = " ", y = "SD / 1000") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  facet_grid(. ~ var)
NO3.flux.N.sd.plot$labels$colour = "Year"


NO3.flux.N.ac.plot = ggplot(data = dat.cp.NO3_flux, aes(x = totN_kg, y = ac, color = Year))+
  geom_point()+
  labs(x = " ", y = "Autocorrelation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  facet_grid(. ~ var)
NO3.flux.N.ac.plot$labels$colour = "Year"


NO3.flux.N.sk.plot = ggplot(data = dat.cp.NO3_flux, aes(x = totN_kg, y = sk, color = Year))+
  geom_point()+
  labs(x = " ", y = "Skewness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  facet_grid(. ~ var)
NO3.flux.N.sk.plot$labels$colour = "Year"


NO3.flux.N.kt.plot = ggplot(data = dat.cp.NO3_flux, aes(x = totN_kg, y = kt, color = Year))+
  geom_point()+
  labs(x = " ", y = "Kurtosis") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)) +
  facet_grid(. ~ var)
NO3.flux.N.kt.plot$labels$colour = "Year"


#setwd for figure export
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience\\EWS Summary Figures")

#call pdf options
pdf.options(width= 6.5, height= 7.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="Figure_S3.pdf")


#Combine plots
NO3.flux.N = plot_grid(NO3.flux.N.sd.plot, NO3.flux.N.ac.plot, NO3.flux.N.sk.plot, NO3.flux.N.kt.plot,
                        labels = c("a", "b", "c", "d"), axis = 'lb', ncol = 1, rel_heights = c(1,1))

#create common x and y labels
y.grob <- textGrob(expression(paste(bold("EWS of Watershed NO"[3]*"Flux"))), 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Annual Precipitation N Flux (kg / ha)", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))

#add to plot
grid.arrange(arrangeGrob(NO3.flux.N, left = y.grob, bottom = x.grob,
                         vp=viewport(width=1, height=1, clip = TRUE)))

dev.off()


#####################
#####################
#multivariate models#
#####################
#####################

#evaluate collinearity between predictor variables
#by calculating the VIF (variance inflation factor)

#fit models with each predictor variable in turn as the 
#response variable
lm1.TMEAN = lm(TMEAN ~ PRCP + ANC + volwt_pH + SO4_flux + NO3_flux + NH4_flux + totN, data = dat.cp)
lm1.PRCP = lm(PRCP ~ TMEAN + ANC + volwt_pH + SO4_flux + NO3_flux + NH4_flux + totN, data = dat.cp)
lm1.ANC = lm(ANC ~ TMEAN + PRCP + volwt_pH + SO4_flux + NO3_flux + NH4_flux + totN, data = dat.cp)
lm1.pH = lm(volwt_pH ~ TMEAN + PRCP + ANC + SO4_flux + NO3_flux + NH4_flux + totN, data = dat.cp)
lm1.SO4 = lm(SO4_flux ~ TMEAN + PRCP + ANC + volwt_pH + NO3_flux + NH4_flux + totN, data = dat.cp)
lm1.NO3 = lm(NO3_flux ~ TMEAN + PRCP + ANC + volwt_pH + SO4_flux + NH4_flux + totN, data = dat.cp)
lm1.NH4 = lm(NH4_flux ~ TMEAN + PRCP + ANC + volwt_pH + SO4_flux + NO3_flux + totN, data = dat.cp)
lm1.totN = lm(totN ~ TMEAN + PRCP + ANC + volwt_pH + SO4_flux + NO3_flux + NH4_flux, data = dat.cp)

#calculate variance inflation statistic.
vif1 = data.frame(cbind("TMEAN" = sqrt(1/(1-(summary(lm1.TMEAN)[[8]]))),
                        "PRCP" = sqrt(1/(1-(summary(lm1.PRCP)[[8]]))),
                        "ANC" = sqrt(1/(1-(summary(lm1.ANC)[[8]]))),
                        "pH" = sqrt(1/(1-(summary(lm1.pH)[[8]]))),
                        "SO4" = sqrt(1/(1-(summary(lm1.SO4)[[8]]))),
                        "NO3" = sqrt(1/(1-(summary(lm1.NO3)[[8]]))),
                        "NH4" = sqrt(1/(1-(summary(lm1.NH4)[[8]]))),
                        "totN" = sqrt(1/(1-(summary(lm1.totN)[[8]])))
))

print(vif1)

#TMEAN     PRCP      ANC       pH       SO4      NO3 NH4 totN
#1.258574  2.069159  6.231187  6.734527 6.184236 Inf Inf  Inf

#run models again with NH4, NO3, ANC, and SO4 removed
lm2.TMEAN = lm(TMEAN ~ PRCP + volwt_pH +  totN, data = dat.cp)
lm2.PRCP = lm(PRCP ~ TMEAN + volwt_pH + totN, data = dat.cp)
lm2.pH = lm(volwt_pH ~ TMEAN + PRCP + totN, data = dat.cp)
lm2.totN = lm(totN ~ TMEAN + PRCP +  volwt_pH, data = dat.cp)

#calculate variance inflation statistic.
vif2 = data.frame(cbind("TMEAN" = sqrt(1/(1-(summary(lm2.TMEAN)[[8]]))),
                        "PRCP" = sqrt(1/(1-(summary(lm2.PRCP)[[8]]))),
                        "pH" = sqrt(1/(1-(summary(lm2.pH)[[8]]))),
                        "totN" = sqrt(1/(1-(summary(lm2.totN)[[8]])))
))

print(vif2)

#TMEAN     PRCP       pH        totN
#1.223364  1.232017   3.247667  2.9133

#all vifs now ~3

####################
####################
#multivariate birds#
####################
####################

#subset bird data
bird.dr = dat.cp[dat.cp$cat == "bird", ]

#subset sd, ac, sk, and kt and standardize data on a 0-1 scale

#sd
bird.sd = data.frame("amre" =  bird.dr[bird.dr$var == "annamre", ]$sd)
bird.sd$btbw = bird.dr[bird.dr$var == "annbtbw", ]$sd
bird.sd$btgw = bird.dr[bird.dr$var == "annbtgw", ]$sd
bird.sd$lefl = bird.dr[bird.dr$var == "annlefl", ]$sd
bird.sd$oven = bird.dr[bird.dr$var == "annoven", ]$sd
bird.sd$revi = bird.dr[bird.dr$var == "annrevi", ]$sd

bird.sd.zs = apply(as.matrix(bird.sd), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#ac
bird.ac = data.frame("amre" =  bird.dr[bird.dr$var == "annamre", ]$ac)
bird.ac$btbw = bird.dr[bird.dr$var == "annbtbw", ]$ac
bird.ac$btgw = bird.dr[bird.dr$var == "annbtgw", ]$ac
bird.ac$lefl = bird.dr[bird.dr$var == "annlefl", ]$ac
bird.ac$oven = bird.dr[bird.dr$var == "annoven", ]$ac
bird.ac$revi = bird.dr[bird.dr$var == "annrevi", ]$ac

bird.ac.zs = apply(as.matrix(bird.ac), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#sk
bird.sk = data.frame("amre" =  bird.dr[bird.dr$var == "annamre", ]$sk)
bird.sk$btbw = bird.dr[bird.dr$var == "annbtbw", ]$sk
bird.sk$btgw = bird.dr[bird.dr$var == "annbtgw", ]$sk
bird.sk$lefl = bird.dr[bird.dr$var == "annlefl", ]$sk
bird.sk$oven = bird.dr[bird.dr$var == "annoven", ]$sk
bird.sk$revi = bird.dr[bird.dr$var == "annrevi", ]$sk

bird.sk.zs = apply(as.matrix(bird.sk), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#kt
bird.kt = data.frame("amre" =  bird.dr[bird.dr$var == "annamre", ]$kt)
bird.kt$btbw = bird.dr[bird.dr$var == "annbtbw", ]$kt
bird.kt$btgw = bird.dr[bird.dr$var == "annbtgw", ]$kt
bird.kt$lefl = bird.dr[bird.dr$var == "annlefl", ]$kt
bird.kt$oven = bird.dr[bird.dr$var == "annoven", ]$kt
bird.kt$revi = bird.dr[bird.dr$var == "annrevi", ]$kt

bird.kt.zs = apply(as.matrix(bird.kt), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))


#combine standardized ews
bird.res.zs = rbind(bird.sd.zs, bird.ac.zs, bird.sk.zs, bird.kt.zs)

#create vector with ews for downstream analysis
ews = rep(c("SD", "AC", "SK", "KT"), 
          each = nrow(bird.sd))

#calculate z scores for driver variables

#first subset driver variables to only include unique rows since drivers repeat for each variable
bird.driv = bird.dr[bird.dr$var == "annamre", ]

bird.dr.zs = data.frame(apply(as.matrix(bird.driv[, c("TMEAN", "PRCP", "volwt_pH", "totN")]), 
        MARGIN = 2, FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T))))

#######################
#run NMDS ordinations#
#######################

#combine rows of predictor variables
dr.zc = rbind(bird.dr.zs, bird.dr.zs, bird.dr.zs, bird.dr.zs)

#determine the best dissimilarity matrix for data
rankindex(grad = dr.zc, veg = bird.res.zs, indices = c("euc", "man", "bray", "gow", "jac", "kul"))
#bray, jac, and kul have the same ranks

#use vegdist function to create distance matrices of bird data
BC_bird = vegdist(bird.res.zs, method = "bray")

#fit NMDS ordinations using 1-5 dimensions
bird.meta1 = metaMDS(BC_bird, k = 1, autotransform = FALSE, trace = FALSE)
bird.meta2 = metaMDS(BC_bird, k = 2, autotransform = FALSE, trace = FALSE)
bird.meta3 = metaMDS(BC_bird, k = 3, autotransform = FALSE, trace = FALSE)
bird.meta4 = metaMDS(BC_bird, k = 4, autotransform = FALSE, trace = FALSE)
bird.meta5 = metaMDS(BC_bird, k = 5, autotransform = FALSE, trace = FALSE)

#extract stress values for each ordination
stress = c(bird.meta1$stress, bird.meta2$stress, bird.meta3$stress, bird.meta4$stress, bird.meta5$stress)
dimen = c(1, 2, 3, 4, 5)
#creates a scree plot for stress
plot(dimen, stress, type = "l")
abline(0.1, 0)
#ordination gets below 0.1 at k = 3 dimensions

#make NMDS object into ordination plots
ordi.bird = ordiplot(bird.meta3)


#make NMDS objects into dataframes, and then add ews information
nmds.bird.df = data.frame(NMDS1 = ordi.bird$sites[,1], NMDS2 =ordi.bird$sites[,2]) #site = rows in the data
nmds.bird.ews = cbind(ews, nmds.bird.df) 
nmds.bird.ews$ews2 = ifelse(nmds.bird.ews$ews == "AC", 2,
                            ifelse(nmds.bird.ews$ews == "KT", 4,
                                   ifelse(nmds.bird.ews$ews == "SD", 1,
                                          ifelse(nmds.bird.ews$ews == "SK", 3, NA))))


#calculate mean values for the axes for plotting results
NMDS1.birds = describeBy(nmds.bird.ews$NMDS1, 
                         group = list(nmds.bird.ews$ews2),
                         na.rm = TRUE, mat = TRUE)
NMDS2.birds = describeBy(nmds.bird.ews$NMDS2, 
                         group = list(nmds.bird.ews$ews2),
                         na.rm = TRUE, mat = TRUE)   

NMDS.birdp= data.frame(NMDS1 = NMDS1.birds$mean, NMDS2 = NMDS2.birds$mean, ews = NMDS1.birds$group1,
                       NMDS1se = NMDS1.birds$se, NMDS2se = NMDS2.birds$se)

######################################################
#envfit to examine gradients in environmental drivers#
######################################################

#run envfit
envfit.bird = envfit(bird.meta3, dr.zc, perms = 999, na.rm = T)

#extract scores 
bird.vecs = as.data.frame(scores(envfit.bird, display = "vectors")) 
bird.vecs$pvals = envfit.bird$vectors$pvals
bird.vecs$r2 = envfit.bird$vectors$r

#scale for plotting. This can be adjusted a lot
bird.vecs$axis1 = bird.vecs$NMDS1 * 0.09
bird.vecs$axis2 = bird.vecs$NMDS2 * 0.09

#create labels for the ordination plot
bird.vecs$label1 = c(" ", " ", " ", "Total N")
bird.vecs$label2 = c(" ", "Prcp", "pH", " ")
bird.vecs$label3 = c("MAT", " ", " ", " ")

#factorize and relevel the ews names so they
#show up in the right order on the plot
NMDS.birdp$ews = factor(NMDS.birdp$ews, levels = c(1, 2, 3, 4), labels =  c("SD", "AC", "SK", "KT"))

#create ordination plot

#note that some of the arguments customizing color have been commented out. The new version of ggplot2 does not
#like them

bird.ord = ggplot(NMDS.birdp, aes(x=NMDS1, y = NMDS2, color = ews, fill = ews,
                                  shape = ews))+
  geom_point(aes(size = ews2), shape = c(21, 22, 23, 24), size = 4)+
  geom_vline(xintercept = 0, color = "red", linetype = 2)+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_x_continuous(
    labels = c("-0.06", "0.00", "0.06", "0.12"), breaks = c(-0.06, 0, 0.06, 0.12)) +
  scale_y_continuous(
    labels = c("-0.06", "0.00", "0.06", "0.12"), breaks = c(-0.06, 0.00, 0.06, 0.12)) +
  #scale_colour_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                               "mediumpurple3", aesthetics = c("colour", "fill")))+
  #scale_fill_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                             "mediumpurple3"))+
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_errorbar(aes(ymax=NMDS2+NMDS2se, ymin=NMDS2-NMDS2se), width = 0)+
  geom_errorbarh(aes(xmax=NMDS1+NMDS1se, xmin=NMDS1-NMDS1se), height = 0)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.text.y = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.y = element_text(family="sans", size = 12, color = "black"))+
  labs(x = "NMDS1 (stress = 0.13)", y = "NMDS2")+
  theme(legend.title = element_blank(), 
        legend.text = element_blank(),
        legend.position = "none") +
  annotate("segment", x = rep(0, 4), y = rep(0,4), xend = bird.vecs$axis1, yend = bird.vecs$axis2,
           arrow = arrow(length=unit(0.025, 'npc')))+ #npc is for changing arrowhead size and the value changes that size
  annotate("text", x = (bird.vecs$axis1 - 0.02), y = (1.5 * bird.vecs$axis2), label = bird.vecs$label1)+
  annotate("text", x = (bird.vecs$axis1 - 0.01), y = (1.1 * bird.vecs$axis2), label = bird.vecs$label2)+
  annotate("text", x = (bird.vecs$axis1 + 0.02), y = (-0.01 + bird.vecs$axis2), label = bird.vecs$label3)+
  annotate("text", x = -0.12, y = 0.14, label = "Bird\nAbundance", fontface = "bold", hjust = 0)

####################
####################
#multivariate trees#
####################
####################

#subset tree data
tree.dr = dat.cp[dat.cp$cat == "bai", ]

#subset sd, ac, sk, and kt and standardize data on a 0-1 scale

#sd
tree.sd = data.frame("ACSA" =  tree.dr[tree.dr$var == "ACSA", ]$sd)
tree.sd$FAGR = tree.dr[tree.dr$var == "FAGR", ]$sd

tree.sd.zs = apply(as.matrix(tree.sd), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#ac
tree.ac = data.frame("ACSA" =  tree.dr[tree.dr$var == "ACSA", ]$ac)
tree.ac$FAGR = tree.dr[tree.dr$var == "FAGR", ]$ac

tree.ac.zs = apply(as.matrix(tree.ac), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#sk
tree.sk = data.frame("ACSA" =  tree.dr[tree.dr$var == "ACSA", ]$sk)
tree.sk$FAGR = tree.dr[tree.dr$var == "FAGR", ]$sk

tree.sk.zs = apply(as.matrix(tree.sk), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#kt
tree.kt = data.frame("ACSA" =  tree.dr[tree.dr$var == "ACSA", ]$kt)
tree.kt$FAGR = tree.dr[tree.dr$var == "FAGR", ]$kt

tree.kt.zs = apply(as.matrix(tree.kt), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#combine standardized ews
tree.res.zs = rbind(tree.sd.zs, tree.ac.zs, tree.sk.zs, tree.kt.zs)

#create vector with ews for downstream analysis
ews = rep(c("SD", "AC", "SK", "KT"), 
          each = nrow(tree.sd))

#calculate z scores for driver variables

#first subset driver variables to only include unique rows since drivers repeat for each variable
tree.driv = tree.dr[tree.dr$var == "ACSA", ]

tree.dr.zs = data.frame(apply(as.matrix(tree.driv[, c("TMEAN", "PRCP", "volwt_pH", "totN")]), 
                              MARGIN = 2, FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T))))

#######################
#run NMDS ordinations#
#######################

#determine the best dissimilarity matrix for data

#combine rows of predictor variables
dr.zc = rbind(tree.dr.zs, tree.dr.zs, tree.dr.zs, tree.dr.zs)

rankindex(grad = dr.zc, veg = tree.res.zs, indices = c("euc", "man", "bray", "gow", "jac", "kul"))
#bray, jac, and kl have the same ranks

#use vegdist function to create distance matrices of basal area increment
BC_tree = vegdist(tree.res.zs, method = "bray")

#fit NMDS ordinations using 1-5 dimensions
tree.meta1 = metaMDS(BC_tree, k = 1, autotransform = FALSE, trace = FALSE)
tree.meta2 = metaMDS(BC_tree, k = 2, autotransform = FALSE, trace = FALSE)
tree.meta3 = metaMDS(BC_tree, k = 3, autotransform = FALSE, trace = FALSE)
tree.meta4 = metaMDS(BC_tree, k = 4, autotransform = FALSE, trace = FALSE)
tree.meta5 = metaMDS(BC_tree, k = 5, autotransform = FALSE, trace = FALSE)

#extract stress values for each ordination
stress = c(tree.meta1$stress, tree.meta2$stress, tree.meta3$stress, tree.meta4$stress, tree.meta5$stress)
dimen = c(1, 2, 3, 4, 5)
#creates a scree plot for stress
plot(dimen, stress, type = "l")
abline(0.1, 0)
#ordination goes below 0.1 with two dimensions until
#there are four dimensions. Use tree.meta2
#in downstream analyses

#make NMDS object into ordination plots
ordi.tree = ordiplot(tree.meta2)

#make NMDS objects into dataframes, and then add ews information
nmds.tree.df = data.frame(NMDS1 = ordi.tree$sites[,1], NMDS2 =ordi.tree$sites[,2]) #site = rows in the data
nmds.tree.ews = cbind(ews, nmds.tree.df) 
nmds.tree.ews$ews2 = ifelse(nmds.tree.ews$ews == "AC", 2,
                            ifelse(nmds.tree.ews$ews == "KT", 4,
                                   ifelse(nmds.tree.ews$ews == "SD", 1,
                                          ifelse(nmds.tree.ews$ews == "SK", 3, NA))))


#calculate mean values for the axes for plotting results
NMDS1.trees = describeBy(nmds.tree.ews$NMDS1, 
                         group = list(nmds.tree.ews$ews2),
                         na.rm = TRUE, mat = TRUE)
NMDS2.trees = describeBy(nmds.tree.ews$NMDS2, 
                         group = list(nmds.tree.ews$ews2),
                         na.rm = TRUE, mat = TRUE)   

NMDS.treep= data.frame(NMDS1 = NMDS1.trees$mean, NMDS2 = NMDS2.trees$mean, ews = NMDS1.trees$group1,
                       NMDS1se = NMDS1.trees$se, NMDS2se = NMDS2.trees$se)

######################################################
#envfit to examine gradients in environmental drivers#
######################################################

#run envfit
envfit.tree = envfit(tree.meta2, dr.zc, perms = 999, na.rm = TRUE)

#extract scores
tree.vecs = as.data.frame(scores(envfit.tree, display = "vectors")) 
tree.vecs$pvals = envfit.tree$vectors$pvals
tree.vecs$r2 = envfit.tree$vectors$r

#scale for plotting. This can be adjusted a lot
tree.vecs$axis1 = tree.vecs$NMDS1 * 0.2
tree.vecs$axis2 = tree.vecs$NMDS2 * 0.2

#create labels for the ordination plot
tree.vecs$label1 = c(" ", " ", "pH", " ")
tree.vecs$label2 = c("MAT", "Prcp", " ", " ")
tree.vecs$label3 = c(" ", " ", " ", "Total N")

#factorize and relevel the ews names so they
#show up in the right order on the plot
NMDS.treep$ews = factor(NMDS.treep$ews, levels = c(1, 2, 3, 4), labels =  c("SD", "AC", "SK", "KT"))

#create ordination plot
tree.ord = ggplot(NMDS.treep, aes(x=NMDS1, y = NMDS2, color = ews, fill = ews,
                                  shape = ews))+
  geom_point(aes(size = ews2), shape = c(21, 22, 23, 24), size = 4)+
  geom_vline(xintercept = 0, color = "red", linetype = 2)+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_y_continuous(breaks = c(-0.09, -0.03, 0.03, 0.09), labels = c("-0.09", "-0.03", "0.03", "0.09"))+
  #scale_x_continuous(breaks = c(-0.07, 0, 0.07, 0.14), labels = c("-0.07", "0.00", "0.07", "0.14"))+
  #scale_colour_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                               "mediumpurple3", aesthetics = c("colour", "fill")))+
  #scale_fill_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                             "mediumpurple3"))+
  #scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = c(21, 22, 23, 24))))+#,
                        #shape = guide_legend(override.aes = list(fill = c("tomato3", "dodgerblue3", "darkolivegreen3", 
                        #                                                  "mediumpurple3")))))+
  geom_errorbar(aes(ymax=NMDS2+NMDS2se, ymin=NMDS2-NMDS2se), width = 0)+
  geom_errorbarh(aes(xmax=NMDS1+NMDS1se, xmin=NMDS1-NMDS1se), height = 0)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.text.y = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.y = element_text(family="sans", size = 12, color = "black"))+
  labs(x = "NMDS1 (stress = 0.09)", y = "NMDS2")+
  theme(legend.position = c(0.65, 0.35), legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 10)) +
  annotate("segment", x = rep(0, 4), y = rep(0,4), xend = tree.vecs$axis1, yend = tree.vecs$axis2,
           arrow = arrow(length=unit(0.025, 'npc')))+ #npc is for changing arrowhead size and the value changes that size
  annotate("text", x = (tree.vecs$axis1 * 1.3), y = (1 * tree.vecs$axis2), label = tree.vecs$label1)+
  annotate("text", x = (tree.vecs$axis1 + 0.02), y = (1 * tree.vecs$axis2), label = tree.vecs$label2)+
  annotate("text", x = (tree.vecs$axis1 - 0.01), y = (0.5 * tree.vecs$axis2), label = tree.vecs$label3)+
  annotate("text", x = -0.08, y = 0.07, label = "Basal Area\nIncrement", fontface = "bold", hjust = 0)

####################
####################
#multivariate concs#
####################
####################

#subset conc data and remove records before 2003 so that all data encompass same period of record
conc.dr = dat.cp[dat.cp$cat == "volwt_Ca" | dat.cp$cat == "volwt_NO3", ]
conc.dr = conc.dr[conc.dr$Year >= 2003 & conc.dr$Year <= 2017, ]


#subset sd, ac, sk, and kt and standardize data on a 0-1 scale

#sd
conc.sd = data.frame("ws1_volwt_Ca" =  conc.dr[conc.dr$var == "ws1_volwt_Ca", ]$sd)
conc.sd$ws2_volwt_Ca = conc.dr[conc.dr$var == "ws2_volwt_Ca", ]$sd
conc.sd$ws4_volwt_Ca = conc.dr[conc.dr$var == "ws4_volwt_Ca", ]$sd
conc.sd$ws5_volwt_Ca = conc.dr[conc.dr$var == "ws5_volwt_Ca", ]$sd
conc.sd$ws6_volwt_Ca = conc.dr[conc.dr$var == "ws6_volwt_Ca", ]$sd

conc.sd.zs = apply(as.matrix(conc.sd), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#ac
conc.ac = data.frame("ws1_volwt_Ca" =  conc.dr[conc.dr$var == "ws1_volwt_Ca", ]$ac)
conc.ac$ws2_volwt_Ca = conc.dr[conc.dr$var == "ws2_volwt_Ca", ]$ac
conc.ac$ws4_volwt_Ca = conc.dr[conc.dr$var == "ws4_volwt_Ca", ]$ac
conc.ac$ws5_volwt_Ca = conc.dr[conc.dr$var == "ws5_volwt_Ca", ]$ac
conc.ac$ws6_volwt_Ca = conc.dr[conc.dr$var == "ws6_volwt_Ca", ]$ac

conc.ac.zs = apply(as.matrix(conc.ac), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))


#sk
conc.sk = data.frame("ws1_volwt_Ca" =  conc.dr[conc.dr$var == "ws1_volwt_Ca", ]$sk)
conc.sk$ws2_volwt_Ca = conc.dr[conc.dr$var == "ws2_volwt_Ca", ]$sk
conc.sk$ws4_volwt_Ca = conc.dr[conc.dr$var == "ws4_volwt_Ca", ]$sk
conc.sk$ws5_volwt_Ca = conc.dr[conc.dr$var == "ws5_volwt_Ca", ]$sk
conc.sk$ws6_volwt_Ca = conc.dr[conc.dr$var == "ws6_volwt_Ca", ]$sk

conc.sk.zs = apply(as.matrix(conc.sk), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#kt
conc.kt = data.frame("ws1_volwt_Ca" =  conc.dr[conc.dr$var == "ws1_volwt_Ca", ]$kt)
conc.kt$ws2_volwt_Ca = conc.dr[conc.dr$var == "ws2_volwt_Ca", ]$kt
conc.kt$ws4_volwt_Ca = conc.dr[conc.dr$var == "ws4_volwt_Ca", ]$kt
conc.kt$ws5_volwt_Ca = conc.dr[conc.dr$var == "ws5_volwt_Ca", ]$kt
conc.kt$ws6_volwt_Ca = conc.dr[conc.dr$var == "ws6_volwt_Ca", ]$kt

conc.kt.zs = apply(as.matrix(conc.kt), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))


#combine standardized ews
conc.res.zs = rbind(conc.sd.zs, conc.ac.zs, conc.sk.zs, conc.kt.zs)

#create vector with ews for downstream analysis
ews = rep(c("SD", "AC", "SK", "KT"), 
          each = nrow(conc.sd))

#calculate z scores for driver variables

#first subset driver variables to only include unique rows since drivers repeat for each variable
conc.driv = conc.dr[conc.dr$var == "ws1_volwt_Ca", ]

conc.dr.zs = data.frame(apply(as.matrix(conc.driv[, c("TMEAN", "PRCP", "volwt_pH", "totN")]), 
                              MARGIN = 2, FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T))))

######################
#run NMDS ordinations#
######################

#determine the best dissimilarity matrix for data

#combine rows of predictor variables
dr.zc = rbind(conc.dr.zs, conc.dr.zs, conc.dr.zs, conc.dr.zs)

rankindex(grad = dr.zc, veg = conc.res.zs, indices = c("euc", "man", "bray", "gow", "jac", "kul"))
#bray, jac, and kl have similar ranks

#use vegdist function to create distance matrices of foliar chemistry data
BC_conc = vegdist(conc.res.zs, method = "bray")

#fit NMDS ordinations using 1-5 dimensions
conc.meta1 = metaMDS(BC_conc, k = 1, autotransform = FALSE, trace = FALSE)
conc.meta2 = metaMDS(BC_conc, k = 2, autotransform = FALSE, trace = FALSE)
conc.meta3 = metaMDS(BC_conc, k = 3, autotransform = FALSE, trace = FALSE)
conc.meta4 = metaMDS(BC_conc, k = 4, autotransform = FALSE, trace = FALSE)
conc.meta5 = metaMDS(BC_conc, k = 5, autotransform = FALSE, trace = FALSE)

#extract stress values for each ordination
stress = c(conc.meta1$stress, conc.meta2$stress, conc.meta3$stress, conc.meta4$stress, conc.meta5$stress)
dimen = c(1, 2, 3, 4, 5)
#creates a scree plot for stress
plot(dimen, stress, type = "l")
abline(0.1, 0)
#ordination goes below 0.1 with three dimensions until
#there are four dimensions. Use conc.meta3
#in downstream analyses

#make NMDS object into ordination plots
ordi.conc = ordiplot(conc.meta3)

#make NMDS objects into dataframes, and then add ews information
nmds.conc.df = data.frame(NMDS1 = ordi.conc$sites[,1], NMDS2 =ordi.conc$sites[,2]) #site = rows in the data
nmds.conc.ews = cbind(ews, nmds.conc.df) 
nmds.conc.ews$ews2 = ifelse(nmds.conc.ews$ews == "AC", 2,
                            ifelse(nmds.conc.ews$ews == "KT", 4,
                                   ifelse(nmds.conc.ews$ews == "SD", 1,
                                          ifelse(nmds.conc.ews$ews == "SK", 3, NA))))


#calculate mean values for the axes for plotting results
NMDS1.concs = describeBy(nmds.conc.ews$NMDS1, 
                         group = list(nmds.conc.ews$ews2),
                         na.rm = TRUE, mat = TRUE)
NMDS2.concs = describeBy(nmds.conc.ews$NMDS2, 
                         group = list(nmds.conc.ews$ews2),
                         na.rm = TRUE, mat = TRUE)   

NMDS.concp= data.frame(NMDS1 = NMDS1.concs$mean, NMDS2 = NMDS2.concs$mean, ews = NMDS1.concs$group1,
                       NMDS1se = NMDS1.concs$se, NMDS2se = NMDS2.concs$se)

######################################################
#envfit to examine gradients in environmental drivers#
######################################################

#combine rows of predictor variables for fitting 
#vectors onto the ordination
dr.zc = rbind(conc.dr.zs, conc.dr.zs, conc.dr.zs, conc.dr.zs)

#run envfit
envfit.conc = envfit(conc.meta3, dr.zc, perms = 999, na.rm = TRUE)

#extract scores
conc.vecs = as.data.frame(scores(envfit.conc, display = "vectors")) 
conc.vecs$pvals = envfit.conc$vectors$pvals
conc.vecs$r2 = envfit.conc$vectors$r

#scale for plotting. This can be adjusted a lot
conc.vecs$axis1 = conc.vecs$NMDS1 * 0.05
conc.vecs$axis2 = conc.vecs$NMDS2 * 0.05

#create labels for the ordination plot
conc.vecs$label1 = c("MAT", " ", "pH", " ")
conc.vecs$label2 = c(" ", "Prcp", " ", " ")
conc.vecs$label3 = c(" ", " ", " ", "Total N")

#factorize and relevel the ews names so they
#show up in the right order on the plot
NMDS.concp$ews = factor(NMDS.concp$ews, levels = c(1, 2, 3, 4), labels =  c("SD", "AC", "SK", "KT"))

#create ordination plot
conc.ord = ggplot(NMDS.concp, aes(x=NMDS1, y = NMDS2, color = ews, fill = ews,
                                  shape = ews))+
  geom_point(aes(size = ews2), shape = c(21, 22, 23, 24), size = 4)+
  geom_vline(xintercept = 0, color = "red", linetype = 2)+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  #scale_y_continuous(breaks = c(-0.10, -0.05, 0.00, 0.05), labels = c("-0.10", "-0.05", "0.00", "0.05"))+
  #scale_x_continuous(breaks = c(-0.10, -0.05, 0.05, 0.10), labels = c("-0.10", "-0.05", "0.05", "0.10"))+
  #scale_colour_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                               "mediumpurple3", aesthetics = c("colour", "fill")))+
  #scale_fill_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                             "mediumpurple3"))+
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  #guides(fill = guide_legend(override.aes = list(shape = c(21, 22, 23, 24)),
  #                           shape = guide_legend(override.aes = list(fill = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                                                                             "mediumpurple3")))))+
  geom_errorbar(aes(ymax=NMDS2+NMDS2se, ymin=NMDS2-NMDS2se), width = 0)+
  geom_errorbarh(aes(xmax=NMDS1+NMDS1se, xmin=NMDS1-NMDS1se), height = 0)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.text.y = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.y = element_text(family="sans", size = 12, color = "black"))+
  labs(x = "NMDS1 (stress = 0.07)", y = "NMDS2")+
  theme(legend.position = "none", legend.title = element_blank(), 
        legend.text = element_blank()) +
  annotate("segment", x = rep(0, 4), y = rep(0,4), xend = conc.vecs$axis1, yend = conc.vecs$axis2,
           arrow = arrow(length=unit(0.025, 'npc')))+ #npc is for changing arrowhead size and the value changes that size
  annotate("text", x = (conc.vecs$axis1 + 0.04), y = (0.8 * conc.vecs$axis2), label = conc.vecs$label1)+
  annotate("text", x = (conc.vecs$axis1 - 0.04), y = (1 * conc.vecs$axis2), label = conc.vecs$label2)+
  annotate("text", x = (conc.vecs$axis1 + 0.03), y = (conc.vecs$axis2 - 0.01), label = conc.vecs$label3)+
  annotate("text", x =  0.02, y = 0.08, label = "Solute\nConcentrations", fontface = "bold", hjust = 0)

#####################
#####################
#multivariate fluxes#
#####################
#####################

#subset flux data and remove records before 2003 so that all data encompass same period of record
flux.dr = dat.cp[dat.cp$cat == "flux_Ca" | dat.cp$cat == "flux_NO3", ]
flux.dr = flux.dr[flux.dr$Year >= 2003 & flux.dr$Year <= 2017, ]

#subset sd, ac, sk, and kt and standardize data on a 0-1 scale

#sd
flux.sd = data.frame("ws1_flux_Ca" =  flux.dr[flux.dr$var == "ws1_flux_Ca", ]$sd)
flux.sd$ws2_flux_Ca = flux.dr[flux.dr$var == "ws2_flux_Ca", ]$sd
flux.sd$ws4_flux_Ca = flux.dr[flux.dr$var == "ws4_flux_Ca", ]$sd
flux.sd$ws5_flux_Ca = flux.dr[flux.dr$var == "ws5_flux_Ca", ]$sd
flux.sd$ws6_flux_Ca = flux.dr[flux.dr$var == "ws6_flux_Ca", ]$sd

flux.sd.zs = apply(as.matrix(flux.sd), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#ac
flux.ac = data.frame("ws1_flux_Ca" =  flux.dr[flux.dr$var == "ws1_flux_Ca", ]$ac)
flux.ac$ws2_flux_Ca = flux.dr[flux.dr$var == "ws2_flux_Ca", ]$ac
flux.ac$ws4_flux_Ca = flux.dr[flux.dr$var == "ws4_flux_Ca", ]$ac
flux.ac$ws5_flux_Ca = flux.dr[flux.dr$var == "ws5_flux_Ca", ]$ac
flux.ac$ws6_flux_Ca = flux.dr[flux.dr$var == "ws6_flux_Ca", ]$ac

flux.ac.zs = apply(as.matrix(flux.ac), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))


#sk
flux.sk = data.frame("ws1_flux_Ca" =  flux.dr[flux.dr$var == "ws1_flux_Ca", ]$sk)
flux.sk$ws2_flux_Ca = flux.dr[flux.dr$var == "ws2_flux_Ca", ]$sk
flux.sk$ws4_flux_Ca = flux.dr[flux.dr$var == "ws4_flux_Ca", ]$sk
flux.sk$ws5_flux_Ca = flux.dr[flux.dr$var == "ws5_flux_Ca", ]$sk
flux.sk$ws6_flux_Ca = flux.dr[flux.dr$var == "ws6_flux_Ca", ]$sk

flux.sk.zs = apply(as.matrix(flux.sk), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))

#kt
flux.kt = data.frame("ws1_flux_Ca" =  flux.dr[flux.dr$var == "ws1_flux_Ca", ]$kt)
flux.kt$ws2_flux_Ca = flux.dr[flux.dr$var == "ws2_flux_Ca", ]$kt
flux.kt$ws4_flux_Ca = flux.dr[flux.dr$var == "ws4_flux_Ca", ]$kt
flux.kt$ws5_flux_Ca = flux.dr[flux.dr$var == "ws5_flux_Ca", ]$kt
flux.kt$ws6_flux_Ca = flux.dr[flux.dr$var == "ws6_flux_Ca", ]$kt

flux.kt.zs = apply(as.matrix(flux.kt), MARGIN = 2, 
                   FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T)))


#combine standardized ews
flux.res.zs = rbind(flux.sd.zs, flux.ac.zs, flux.sk.zs, flux.kt.zs)

#create vector with ews for downstream analysis
ews = rep(c("SD", "AC", "SK", "KT"), 
          each = nrow(flux.sd))

#calculate z scores for driver variables

#first subset driver variables to only include unique rows since drivers repeat for each variable
flux.driv = flux.dr[flux.dr$var == "ws1_flux_Ca", ]

flux.dr.zs = data.frame(apply(as.matrix(flux.driv[, c("TMEAN", "PRCP", "volwt_pH", "totN")]), 
                              MARGIN = 2, FUN = function(X) (X - min(X, na.rm = T))/diff(range(X, na.rm = T))))

######################
#run NMDS ordinations#
######################

#determine the best dissimilarity matrix for data

#combine rows of predictor variables
dr.zc = rbind(flux.dr.zs, flux.dr.zs, flux.dr.zs, flux.dr.zs)

rankindex(grad = dr.zc, veg = flux.res.zs, indices = c("euc", "man", "bray", "gow", "jac", "kul"))
#bray, jac, and kl have the same ranks

#use vegdist function to create distance matrices of foliar chemistry data
BC_flux = vegdist(flux.res.zs, method = "bray")

#fit NMDS ordinations using 1-5 dimensions
flux.meta1 = metaMDS(BC_flux, k = 1, autotransform = FALSE, trace = FALSE)
flux.meta2 = metaMDS(BC_flux, k = 2, autotransform = FALSE, trace = FALSE)
flux.meta3 = metaMDS(BC_flux, k = 3, autotransform = FALSE, trace = FALSE)
flux.meta4 = metaMDS(BC_flux, k = 4, autotransform = FALSE, trace = FALSE)
flux.meta5 = metaMDS(BC_flux, k = 5, autotransform = FALSE, trace = FALSE)

#extract stress values for each ordination
stress = c(flux.meta1$stress, flux.meta2$stress, flux.meta3$stress, flux.meta4$stress, flux.meta5$stress)
dimen = c(1, 2, 3, 4, 5)
#creates a scree plot for stress
plot(dimen, stress, type = "l")
abline(0.1, 0)
#ordination does not goes below 0.1 until
#there are four dimensions. Use flux.meta3
#in downstream analyses

#make NMDS object into ordination plots
ordi.flux = ordiplot(flux.meta3)

#make NMDS objects into dataframes, and then add ews information
nmds.flux.df = data.frame(NMDS1 = ordi.flux$sites[,1], NMDS2 =ordi.flux$sites[,2]) #site = rows in the data
nmds.flux.ews = cbind(ews, nmds.flux.df) 
nmds.flux.ews$ews2 = ifelse(nmds.flux.ews$ews == "AC", 2,
                            ifelse(nmds.flux.ews$ews == "KT", 4,
                                   ifelse(nmds.flux.ews$ews == "SD", 1,
                                          ifelse(nmds.flux.ews$ews == "SK", 3, NA))))


#calculate mean values for the axes for plotting results
NMDS1.fluxs = describeBy(nmds.flux.ews$NMDS1, 
                         group = list(nmds.flux.ews$ews2),
                         na.rm = TRUE, mat = TRUE)
NMDS2.fluxs = describeBy(nmds.flux.ews$NMDS2, 
                         group = list(nmds.flux.ews$ews2),
                         na.rm = TRUE, mat = TRUE)   

NMDS.fluxp= data.frame(NMDS1 = NMDS1.fluxs$mean, NMDS2 = NMDS2.fluxs$mean, ews = NMDS1.fluxs$group1,
                       NMDS1se = NMDS1.fluxs$se, NMDS2se = NMDS2.fluxs$se)

######################################################
#envfit to examine gradients in environmental drivers#
######################################################

#run envfit
envfit.flux = envfit(flux.meta3, dr.zc, perms = 999, na.rm = TRUE)

#extract scores
flux.vecs = as.data.frame(scores(envfit.flux, display = "vectors")) 
flux.vecs$pvals = envfit.flux$vectors$pvals
flux.vecs$r2 = envfit.flux$vectors$r

#scale for plotting. This can be adjusted a lot
flux.vecs$axis1 = flux.vecs$NMDS1 * 0.18
flux.vecs$axis2 = flux.vecs$NMDS2 * 0.18

#create labels for the ordination plot
flux.vecs$label1 = c("MAT", " ", "pH", " ")
flux.vecs$label2 = c(" ", " ", " ", "Total N")
flux.vecs$label3 = c(" ", "Prcp", " ", " ")

#factorize and relevel the ews names so they
#show up in the right order on the plot
NMDS.fluxp$ews = factor(NMDS.fluxp$ews, levels = c(1, 2, 3, 4), labels =  c("SD", "AC", "SK", "KT"))

#create ordination plot
flux.ord = ggplot(NMDS.fluxp, aes(x=NMDS1, y = NMDS2, color = ews, fill = ews,
                                  shape = ews))+
  geom_point(aes(size = ews2), shape = c(21, 22, 23, 24), size = 4)+
  geom_vline(xintercept = 0, color = "red", linetype = 2)+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  #scale_y_continuous(breaks = c(-0.10, -0.05, 0.00, 0.05), labels = c("-0.10", "-0.05", "0.00", "0.05"))+
  #scale_x_continuous(breaks = c(-0.10, -0.05, 0.05, 0.10), labels = c("-0.10", "-0.05", "0.05", "0.10"))+
  #scale_colour_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                               "mediumpurple3", aesthetics = c("colour", "fill")))+
  #scale_fill_manual(values = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                             "mediumpurple3"))+
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  #guides(fill = guide_legend(override.aes = list(shape = c(21, 22, 23, 24)),
  #                           shape = guide_legend(override.aes = list(fill = c("tomato3", "dodgerblue3", "darkolivegreen3", 
  #                                                                             "mediumpurple3")))))+
  geom_errorbar(aes(ymax=NMDS2+NMDS2se, ymin=NMDS2-NMDS2se), width = 0)+
  geom_errorbarh(aes(xmax=NMDS1+NMDS1se, xmin=NMDS1-NMDS1se), height = 0)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.text.y = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.x = element_text(family="sans", size = 12, color = "black"))+
  theme(axis.title.y = element_text(family="sans", size = 12, color = "black"))+
  labs(x = "NMDS1 (stress = 0.09)", y = "NMDS2")+
  theme(legend.position = "none", legend.title = element_blank(), 
        legend.text = element_blank()) +
  annotate("segment", x = rep(0, 4), y = rep(0,4), xend = flux.vecs$axis1, yend = flux.vecs$axis2,
           arrow = arrow(length=unit(0.025, 'npc')))+ #npc is for changing arrowhead size and the value changes that size
  annotate("text", x = (flux.vecs$axis1 - 0.02), y = (0.98 * flux.vecs$axis2), label = flux.vecs$label1)+
  annotate("text", x = (flux.vecs$axis1 - 0.03), y = (1 * flux.vecs$axis2), label = flux.vecs$label2)+
  annotate("text", x = (flux.vecs$axis1 - 0.04), y = (1 * flux.vecs$axis2), label = flux.vecs$label3)+
  annotate("text", x = -0.12, y = 0.09, label = "Solute\nFluxes", fontface = "bold", hjust = 0)

###################################
#Creat output table of all envfits#
###################################

all.vecs = data.frame(cbind("var" = c("TMEAN", "PRCP", "pH", "totN"),
                            conc.p = conc.vecs$pvals, conc.r2 = conc.vecs$r2,
                            flux.p = flux.vecs$pvals, flux.r2 = flux.vecs$r2,
                            bird.p = bird.vecs$pvals, bird.r2 = bird.vecs$r2,
                            tree.p = tree.vecs$pvals, tree.r2 = tree.vecs$r2))
#export table
write.table(all.vecs, "EWS Summary Tables\\allenvfits.csv", col.names = T,
            row.names = F, sep = ",")

###################################################
#make multi-panel figure (Figure 5 for manuscript)#
#of ordination plots#                             #
###################################################

  
#call pdf options
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="EWS Summary Figures\\Figure_5.pdf")

ggarrange(conc.ord, flux.ord, bird.ord, tree.ord, #  +rremove("y.text")+ rremove("ylab"), 
          labels = c("a", "b", "c", "d"), #label.x = c(0.15, 0.075), hjust = c(0.1, 0.5),
          #vjust = 2,
          ncol = 2, nrow = 2)

dev.off()

