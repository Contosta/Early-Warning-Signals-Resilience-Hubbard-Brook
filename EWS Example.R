################################################################################################################
################################################################################################################
################################################################################################################

#This script produces a figure (Figure_S5) showing original NO3 flux data, detrended NO3 flux data, and 
#distributions of detrended NO3 flux data for the biogeochemical reference watershed (Watershed 6)
#and the whole-tree harvested watershed (Watershed 5) at Hubbard Brook as an example of how early warning signals
#such as standard deviation, skewness, and kurtosis can change over time in different directions

#Inputs include original and detrended NO3 flux data in the Detrended Data folder. 

#Outputs include Figure_S5 in the manuscript, which is written to the EWS Summary 
#Figures folder. 

#code developed by A Contosta
#most recent version 7/11/2023

################################################################################################################
################################################################################################################
################################################################################################################

#Initial set-up

#call libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(slider)

#set working directory and import files
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\PES_LTER\\Data\\R Projects\\Hubbard-Brook-Resilience")
ws6_flux =  read.table("Detrended Data\\ws6_flux_dt.csv", head = T, sep = ",")
ws5_flux =  read.table("Detrended Data\\ws5_flux_dt.csv", head = T, sep = ",")

#subset dataframes to only include columns needed for making violin plots of distributions of 
#detrended data
df6 = ws6_flux[ , c("WaterYear", "ksm.NO3_flux.5")]
df5 = ws5_flux[ , c("WaterYear", "ksm.NO3_flux.5")]

#Define the size of windows for generating overlapping distributions
window_size6 = nrow(ws6_flux)/2
window_size5 = round(nrow(ws5_flux)/2)

#Get moving windows of subsets
subsets6 = slide(df6, ~.x, .before = (window_size6 - 1), .complete = TRUE)
subsets5 = slide(df5, ~.x, .before = (window_size5 - 1), .complete = TRUE)

#Combine subsets into a new dataframe in long form
new_df6 = bind_rows(subsets6, .id = "window_id6")
new_df5 = bind_rows(subsets5, .id = "window_id5")

#create new column where window_id is numeric and then sort
new_df6$window6 = as.numeric(new_df6$window_id6)
new_df6 = new_df6[order(new_df6$window6), ]

new_df5$window5 = as.numeric(new_df5$window_id5)
new_df5 = new_df5[order(new_df5$window5), ]

################################################################################################################
################################################################################################################
################################################################################################################

#Make plots

##########################
#plot original NO3 fluxes#
##########################

#W6
orig_NO36 = ggplot(data = ws6_flux, aes(x = WaterYear, y = NO3_flux))+
  geom_point(size = 4)+
  theme_bw()+
  geom_line()+
  theme(legend.position = c(0.5, 0.9))+
  labs(x = "Water Year", 
       y = expression("N Flux (g NO"[3]*~ha^-1*~year^-1*")")
  )

#W5
orig_NO35 = ggplot(data = ws5_flux, aes(x = WaterYear, y = NO3_flux))+
  geom_point(size = 4)+
  theme_bw()+
  geom_line()+
  theme(legend.position = c(0.5, 0.9))+
  labs(x = "Water Year", 
       y = expression("N Flux (g NO"[3]*~ha^-1*~year^-1*")")
  )

###########################
#plot detrended NO3 fluxes#
###########################

#W6
dt_NO36 = ggplot(data = ws6_flux, aes(x = WaterYear, y = ksm.NO3_flux.5))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()+
  theme(legend.position = "blank")+
  labs(x = "Water Year", 
       y = "Detrended N Flux"
  )

#W5
dt_NO35 = ggplot(data = ws5_flux, aes(x = WaterYear, y = ksm.NO3_flux.5))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()+
  theme(legend.position = "blank")+
  labs(x = "Water Year", 
       y = "Detrended N Flux"
  )

#######################################
#make violin plots of the distribution# 
#of detrended NO3 fluxes              #
#######################################

viol_NO36 = ggplot(data = new_df6, aes(x = window6, y = ksm.NO3_flux.5, fill = window_id6))+
  geom_violin(alpha = 1)+
  theme_bw()+
  theme(legend.position = "blank")+
  labs(x = "Window", 
       y = "Distribution of N Flux",
       fill = "Window")

viol_NO35 = ggplot(data = new_df5, aes(x = window5, y = ksm.NO3_flux.5, fill = window_id5))+
  geom_violin(alpha = 1)+
  theme_bw()+
  theme(legend.position = "blank")+
  labs(x = "Window", 
       y = "Distribution of N Flux",
       fill = "Window")


####################################
#arrange panels and write to folder#
####################################

#call pdf options
pdf.options(width= 6.5, height= 8.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="EWS Summary Figures\\Figure_S5.pdf")

ggarrange(orig_NO36, orig_NO35,
          dt_NO36, dt_NO35,
          viol_NO36, viol_NO35,
          labels = c("a", "d", "b", "e", "c", "f"),  ncol = 2, nrow = 3,
          common.legend = F)

dev.off()

