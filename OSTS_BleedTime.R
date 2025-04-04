#########
#Packages
#########
library(car)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(emmeans)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)

setwd("/Users/ben/Desktop/OSTS/OSTS 2.0/Data Analysis")
####################
#Loading in the data
####################

datum1_short = read.csv("Experiment1_Short.csv") #this reads the data from the csv
head(datum1_short)  #This shows you your data to double check it is the right dataset
datum1_short$PopulationF = as.factor(datum1_short$Pop2) #as.factor() makes sure R reads this as a categorical variable when graphing or analyzing
datum1_short$TempF = as.factor(datum1_short$Temp)

datum2_short = read.csv("Experiment2_Short.csv")
head(datum2_short)
datum2_short$PopulationF = as.factor(datum2_short$Population) #as.factor() makes sure R reads this as a categorical variable when graphing or analyzing
datum2_short$TreatmentF = as.factor(datum2_short$Treatment)
datum2_short$BaseBleedTime = datum2_short$Baseline.Bleed.Duration..mm.ss.
datum2_short$TreatBleedTime = datum2_short$TreatmentTime

#########
#Analysis
#########

#Convert the time format from excel into seconds
datum2_short$TreatBleedTime = sapply(strsplit(as.character(datum2_short$TreatmentTime), ":"), 
                                     function(x) as.numeric(x[1]) * 60 + as.numeric(x[2]))

#Divide by 60 to get to minutes
datum2_short$TreatBleedTime <- datum2_short$TreatBleedTime / 60

#The actual t-test
t.test(TreatBleedTime ~ Treatment, data = datum2_short, var.equal = TRUE)

#averages and sd for the groups
datum2_short %>%
  group_by(Treatment) %>%  # Group by Cold/Warmup
  summarise(
    Mean_BleedTime = mean(TreatBleedTime, na.rm = TRUE),  # Calculate mean
    SD_BleedTime = sd(TreatBleedTime, na.rm = TRUE),      # Calculate standard deviation
    Count = n())  # Optional: Count number of observations per group
    







