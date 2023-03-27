# Data Cleanup and Manipulation #####
# 2021-06-14
# Author: Monserrat Garcia 


#### Data Cleaning 2018
#Loading Packages needed
require(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
# Data 

meta18<-read.csv("C:/Users/monse/OneDrive/Oyster16S_Data23/Data/metadata_de18 - Copy.csv")
genetics_data2018 <- read.csv("Data/DE2018_alldata - Copy.csv")

# Adding Na's in blank columns
genetics_data2018[genetics_data2018 == "" | genetics_data2018 == " "] <- NA

# Dropping NA rows in "Bucket" column
genetics_data2018 <- genetics_data2018 %>%
  drop_na(Bucket,Color.Number)

#using ifelse for the new species abbreviations

meta18$Species2<- ifelse(meta18$Species== "MB", "LP", 
                         ifelse(meta18$Species== "IR", "IR", 
                                ifelse(meta18$Species=="CV","CV", "AM")))

genetics_data2018$Species2<- ifelse(genetics_data2018$Species== "MB", "LP", 
                                    ifelse(genetics_data2018$Species== "IR", "IR", 
                                           ifelse(genetics_data2018$Species=="CV","CV", "AM")))

#Using ifelse for creating "High/low, poly/mono" ID's 
genetics_data2018$Bucket2<- ifelse(genetics_data2018$Bucket== "HM1", "HIGH_MONO", 
                                   ifelse(genetics_data2018$Bucket=="HM2", "HIGH_MONO", 
                                          ifelse(genetics_data2018$Bucket== "HM3", "HIGH_MONO", 
                                                 ifelse(genetics_data2018$Bucket=="HM4", "HIGH_MONO", 
                                                        ifelse(genetics_data2018$Bucket=="HM5", "HIGH_MONO", 
                                                               ifelse(genetics_data2018$Bucket=="HM6", "HIGH_MONO", 
                                                                      ifelse(genetics_data2018$Bucket=="HM7", "HIGH_MONO", 
                                                                             ifelse(genetics_data2018$Bucket=="HM8", "HIGH_MONO", 
                                                                                    ifelse(genetics_data2018$Bucket=="HM9", "HIGH_MONO", ifelse(genetics_data2018$Bucket== "LM1", "LOW_MONO", 
                                                                                                                                                ifelse(genetics_data2018$Bucket=="LM2", "LOW_MONO", 
                                                                                                                                                       ifelse(genetics_data2018$Bucket== "LM3", "LOW_MONO", 
                                                                                                                                                              ifelse(genetics_data2018$Bucket=="LM4", "LOW_MONO", 
                                                                                                                                                                     ifelse(genetics_data2018$Bucket=="LM5", "LOW_MONO", 
                                                                                                                                                                            ifelse(genetics_data2018$Bucket=="LM6", "LOW_MONO", 
                                                                                                                                                                                   ifelse(genetics_data2018$Bucket=="LM7", "LOW_MONO", 
                                                                                                                                                                                          ifelse(genetics_data2018$Bucket=="LM8", "LOW_MONO", 
                                                                                                                                                                                                 ifelse(genetics_data2018$Bucket=="LM9", "LOW_MONO", ifelse(genetics_data2018$Bucket== "LP1", "LOW_POLY", 
                                                                                                                                                                                                                                                            ifelse(genetics_data2018$Bucket=="LP2", "LOW_POLY", 
                                                                                                                                                                                                                                                                   ifelse(genetics_data2018$Bucket== "LP3", "LOW_POLY", 
                                                                                                                                                                                                                                                                          ifelse(genetics_data2018$Bucket=="LP4", "LOW_POLY", 
                                                                                                                                                                                                                                                                                 ifelse(genetics_data2018$Bucket=="LP5", "LOW_POLY", 
                                                                                                                                                                                                                                                                                        ifelse(genetics_data2018$Bucket=="LP6", "LOW_POLY", 
                                                                                                                                                                                                                                                                                               ifelse(genetics_data2018$Bucket=="LP7", "LOW_POLY", 
                                                                                                                                                                                                                                                                                                      ifelse(genetics_data2018$Bucket=="LP8", "LOW_POLY", 
                                                                                                                                                                                                                                                                                                             ifelse(genetics_data2018$Bucket=="LP9", "LOW_POLY", "HIGH_POLY")))))))))))))))))))))))))))

# Combining Bucket and Color.Number columns

genetics_data2018$Bucket_colnum <- paste0(genetics_data2018$Bucket, genetics_data2018$Color.Number, sep = "")


# Creating UniqueID's for genetics_data2018, Example: 2018__HIGH_POLY_HP1W5_CV

genetics_data2018$UniqueID <- paste("2018", genetics_data2018$Bucket2, genetics_data2018$Bucket_colnum, genetics_data2018$Species, sep = "_")


# Merging the datasets 
meta_gen18_data <- merge(meta18, genetics_data2018, by= "UniqueID", all.x=TRUE)

#Taking out unecessary columns in merging data
meta_gen18_data <- subset(meta_gen18_data, select = -c(X.x, Site, V1, Phase_1_DO, Phase_2_DO, Phase_1_temp, Phase_2_Temp, Overall_treatment, Date_initial_measure, 
                                                       Mortality_Date, Date_FinalMeasurement, RFTM_Date, Parasites, X.y, Species.x, Species.y))

#Saving the new data
write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")

#uploading the asvtable_18
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")


# 2017 

meta17<-read.csv("C:/Users/monse/OneDrive/Oyster16S_Data23/Data/metadata_de17 - Copy.csv")
data<-read.csv("Data/DE_DATA_ForGenetics - Copy.csv")


#Adding new Column for treatment
data$Treatment2<-ifelse(data$Treatment=="HH", "HIGH_POLY", ifelse(data$Treatment=="HL", "HIGH_MONO", ifelse(data$Treatment== "LL", "LOW_MONO","LOW_POLY")))

data$Colornumber<- paste0(data$Color, data$Number)


#Creating unique IDs, Example: 2017_NW_HIGH_MONO_B11_CV

data$UniqueID <- paste("2017", data$Site, data$Treatment2, data$Colornumber,data$Species, sep = "_")

meta17_data<- merge(meta17, data, by= "UniqueID", all.x= TRUE) # if you want to keep all the rows in meta17 all.x=TRUE, but for data you do all.y=TRUE

#Taking out unnecessary Columns
meta17_data2 <- subset(meta17_data, select = -c(X,V1,Phase_1_DO,Phase_1_temp,Phase_2_DO,Phase_2_Temp,
                                                Overall_treatment,Notes_pre,Notes_post,Date_post,
                                                Dry_weight_shell,POST_DEAD_ALIVE,Dry_Weight_plate,Dry_weight_final,Genetics_Weight))

