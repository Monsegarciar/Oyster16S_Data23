#Final Document
#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

#Libraries ####
library("tidyverse")
library("data.table")
library("dplyr")
library("plyr")
library("tidyr")
library("stringr")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("GenomicRanges")
library("GenomeInfoDb") 
library("genefilter")
library("vegan")
library("taxa")
library("ggpubr")
library("RColorBrewer")
library("ggh4x")
library("readr")
library("metacoder")
display.brewer.all()
theme_set(theme_bw())

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
#Data Cleanup + Phyloseq object ####
#Loading the data (Original Data is called DE_DATA_ForGenetics_17.csv)
##de_data17 <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/DE_DATA_ForGenetics_17.csv")
de_data17 <- read.csv("Data/DE2018_alldata - Copy.csv")

#Loading the data (Original Data is called metadata_de17.csv)
##meta17 <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/metadata_de17.csv")
meta17 <- read.csv("Data/metadata_de17 - Copy.csv")

#Loading the data (Original Data Name = asvtable_de17.csv)
##asv17 <- fread("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/asvtable_de17.csv")
asv17 <- fread("Data/asvtable_de17 - Copy.csv")

#Renaming the Treatment Names
de_data17$Treatment
de_data17$Treatment2 <- ifelse(de_data17$Treatment =="HH", "HIGH_POLY",
                               ifelse(de_data17$Treatment == "HL", "HIGH_MONO", 
                                      ifelse(de_data17$Treatment == "LL", "LOW_MONO", "LOW_POLY")))

#Creating a new column names Colornumber 
de_data17$Colornumber <- paste0(de_data17$Color, de_data17$Number) 

#Creating the UniqueIDs in de_data
de_data17$UniqueID <- paste("2017", de_data17$Site, de_data17$Treatment2, de_data17$Colornumber, de_data17$Species, sep = "_")

#Using Merge to combine the two data frames
meta17data <- merge(meta17, de_data17, by = "UniqueID", all.x = TRUE) #Matching by column UniqueID, all.x referrers to Meta17 because it was on the X place

#Deleting columns in the new data frame
meta17data <- select(meta17data, 
                     -"X",
                     -"V1", 
                     -"Phase_1_DO",
                     -"Phase_1_temp", 
                     -"Phase_2_DO", 
                     -"Phase_2_Temp",
                     -"Overall_treatment",
                     -"Date_pre", 
                     -"Date_post",
                     -"Site.y", 
                     -"Number.y", 
                     -"Species.y", 
                     -"RFTM_score.y", 
                     -"peacrabs.y",
                     -"Notes_pre", 
                     -"POST_DEAD_ALIVE",
                     -"Dry_Weight_plate", 
                     -"Dry_weight_final", 
                     -"Dry_weight_shell", 
                     -"Notes_post", 
                     -"Genetics_Weight")

#Making Unique IDs the new row names for Phyloseq
rownames(meta17data) = meta17data$UniqueID
rownames(meta17data) #UniqueID
meta17data$UniqueID=NULL
write.csv(meta17data, file = "Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/cleanmetadata17")


## PHYLOSEQ ####
c_meta17data <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/cleanmetadata17")
asvtable17 <- fread("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/asvtable_de17.csv")
run23 <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Oyster_data_raw/Run123_taxa_complete.csv")

## CHANGING ROW NAMES FOR EACH DATA SET
rownames(c_meta17data) = c_meta17data$X
c_meta17data$X=NULL
rownames(c_meta17data)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable17) = asvtable17$V1
rownames(asvtable17)
asvtable17$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$Row.names
run23$Row.names = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 

## CONVERTING TO MATRICIES
meta17_matrix <- as.data.frame(c_meta17data, rownames("X"))
rownames(meta17_matrix)
#STILL UNIQUE ID

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 

rownames(run23)
taxmat_matrix <- as.matrix(run23) 
colnames(taxmat_matrix) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(taxmat_matrix)
#STILL SEQUENCE 


## SETTING OTU, TAX, SAMP
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
rownames(OTU) #UniqueID

TAX <- tax_table(taxmat_matrix)
rownames(TAX) #Sequence

SAMP <- sample_data(meta17_matrix)
rownames(SAMP)#UniqueID

## INSPECTING SAMPLE NAMES
sample_names(SAMP) #UniqueID
sample_names(OTU) #UniqueID
sample_names(TAX) #NULL


## EVENING OUT THE DATA
OTU=transform_sample_counts(OTU, function(OTU) 1E6 * OTU/sum(OTU))

## READING THROUGH PHYLOSEQ
physeq_class = phyloseq(OTU, TAX, SAMP) 

#Final Object (Original)
physeq_class

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 16383 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 24 sample variables ]
#tax_table()   Taxonomy Table:    [ 16383 taxa by 9 taxonomic ranks ]

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Identifying Intercept Variables ####  

#Make peacrabs into a factor
SAMP$peacrabs.f <- factor(SAMP$peacrabs.x)
typeof(SAMP$peacrabs.f)#Remains "integer"
is.factor(SAMP$peacrabs.f) #True

#Make RFTM_score into a factor
SAMP$RFTM_score.f <- factor(SAMP$RFTM_score.x)
typeof(SAMP$RFTM_score.f) #Remains "integer"
is.factor(SAMP$RFTM_score.f)

typeof(SAMP$peacrabs.x) #Integer
typeof(SAMP$RFTM_score.x) #Double

#make a simpler rftm factor without 0.5 and combining 4 and 5
#changing 0.5 to 1 and 5 to 4, we are left with 1,2,3,4
SAMP$RFTM_simp <- SAMP$RFTM_score.f
SAMP$RFTM_simp <- sub("0.5", "1", SAMP$RFTM_simp)
SAMP$RFTM_simp <- sub("5", "4", SAMP$RFTM_simp)

#make RFTM presence absence
#Anything equal to O is absence, everything else=(1,2,3,4) is presence
SAMP$RFTM_pa <- ifelse(SAMP$RFTM_score.x=="0", 0, 1)
SAMP$RFTM_pa <- factor(SAMP$RFTM_pa)

physeq = phyloseq(OTU, TAX, SAMP) 

dds_rftmxpeasite <- phyloseq_to_deseq2(physeq, ~ RFTM_pa * peacrabs.f + Site.x)
dds_rftmxpeasite <- DESeq(dds_rftmxpeasite, test="Wald", fitType="parametric")
resultsNames(dds_rftmxpeasite)
#RFTM_pa_1_vs_0   use this intercept
#peacrabs.f_1_vs_0   use this intercept
#Site.x_OY_vs_NW
#Site.x_SW_vs_NW
#RFTM_pa1.peacrabs.f1

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ## RFTM_PA
  resPA_rftm_3 <- results(dds_rftmxpeasite, name="RFTM_pa_1_vs_0")
sigPA_rftm_3 <- resPA_rftm_3[which(resPA_rftm_3$padj < 0.05), ]
sigPA_rftm_3 <- cbind(as(sigPA_rftm_3, "data.frame"), as(tax_table(physeq)[rownames(sigPA_rftm_3), ], "matrix"))
sigPA_rftm_3 <- select(sigPA_rftm_3, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_rftm_3)[names(sigPA_rftm_3) == "Genus.x"] <- "Genus"
view(sigPA_rftm_3)
write.table(sigPA_rftm_3, file="Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_RFTMpa_compLog.csv", quote=FALSE,sep = ",", col.names=NA)

sigPA_RFTM_tax3 <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_RFTMpa_compLog.csv")
rownames(sigPA_RFTM_tax3) <- sigPA_RFTM_tax3$X
sigPA_RFTM_tax3$X = NULL
rownames(sigPA_RFTM_tax3)
#Rownames are sequences

sigPA_RFTM_tax3 <- as.matrix(sigPA_RFTM_tax3)
colnames(sigPA_RFTM_tax3) <- c("log2FoldChange", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(sigPA_RFTM_tax3)
#Still sequences 

OTU_R <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
sigPA_RFTM_tax3 <- tax_table(sigPA_RFTM_tax3)
SAMP
view(SAMP)

sample_names(SAMP) #Unique ID
sample_names(OTU_R) #Unique ID
sample_names(sigPA_RFTM_tax3) #NULL

#RFTMpa Complete
#With Default dimensions
F3_RFTMPA_wlog <- phyloseq(OTU_R, sigPA_RFTM_tax3, SAMP) 
F3_RFTMPA_wlog
#otu_table()   OTU Table:         [ 584 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 28 sample variables ]
#tax_table()   Taxonomy Table:    [ 584 taxa by 8 taxonomic ranks ]


##PeaCrab_PA
#peacrabs.f_1_vs_0
resPA_pea_F3 <- results(dds_rftmxpeasite, name="peacrabs.f_1_vs_0")
sigPA_pea_F3 <- resPA_pea_F3[which(resPA_pea_F3$padj < 0.05), ]
sigPA_pea_F3 <- cbind(as(sigPA_pea_F3, "data.frame"), as(tax_table(physeq)[rownames(sigPA_pea_F3), ], "matrix"))
sigPA_pea_F3 <- select(sigPA_pea_F3, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_pea_F3)[names(sigPA_pea_F3) == "Genus.x"] <- "Genus"
view(sigPA_pea_F3)
write.table(sigPA_pea_F3, file="Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_PEApa_compLog.csv", quote=FALSE,sep = ",", col.names=NA)

sigPA_PEA_tax3 <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_PEApa_compLog.csv")
rownames(sigPA_PEA_tax3) <- sigPA_PEA_tax3$X
sigPA_PEA_tax3$X = NULL
rownames(sigPA_PEA_tax3)
#Rownames are sequences

sigPA_PEA_tax3 <- as.matrix(sigPA_PEA_tax3)
colnames(sigPA_PEA_tax3) <- c("log2FoldChange", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(sigPA_PEA_tax3)
#Still sequences 

OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
sigPA_PEA_tax3 <- tax_table(sigPA_PEA_tax3)
SAMP

sample_names(SAMP) #Unique ID
sample_names(OTU) #Unique ID
sample_names(sigPA_PEA_tax3) #NULL

##PEApa Complete
#With Default dimensions
F3_PEAPA_wlog <- phyloseq(OTU, sigPA_PEA_tax3, SAMP) 
F3_PEAPA_wlog
#otu_table()   OTU Table:         [ 326 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 28 sample variables ]
#tax_table()   Taxonomy Table:    [ 326 taxa by 8 taxonomic ranks ]

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #33.3 percent filtering ####  

##RFTMpa Pruned
#Remove OTUs that do not appear more than once in more than one third of the samples
F3_RFTMPA_wlog_filter <- genefilter_sample(F3_RFTMPA_wlog, filterfun_sample(function(x) x > 1), A=0.33*nsamples(F3_RFTMPA_wlog))
F3_RFTMPA_wlog_prune <- prune_taxa(F3_RFTMPA_wlog_filter, F3_RFTMPA_wlog)

F3_RFTMPA_wlog_prune
#otu_table()   OTU Table:         [ 13 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 28 sample variables ]
#tax_table()   Taxonomy Table:    [ 13 taxa by 8 taxonomic ranks ]

F3_RFTMPA_wlog_prune_tax <- tax_table(F3_RFTMPA_wlog_prune)
write.table(F3_RFTMPA_wlog_prune_tax, file="Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_RFTMpa_pruned_log.csv", quote=FALSE,sep = ",", col.names=NA)


##PEApa Pruned
#Remove OTUs that do not appear more than once in more than one third of the samples
F3_PEAPA_wlog_filter <- genefilter_sample(F3_PEAPA_wlog, filterfun_sample(function(x) x > 1), A=0.33*nsamples(F3_PEAPA_wlog))
F3_PEAPA_wlog_prune <- prune_taxa(F3_PEAPA_wlog_filter, F3_PEAPA_wlog)

F3_PEAPA_wlog_prune
#otu_table()   OTU Table:         [ 7 taxa and 112 samples ]
#sample_data() Sample Data:       [ 112 samples by 28 sample variables ]
#tax_table()   Taxonomy Table:    [ 7 taxa by 8 taxonomic ranks ]

F3_PEAPA_wlog_prune_tax <- tax_table(F3_PEAPA_wlog_prune)
write.table(F3_PEAPA_wlog_prune_tax, file="Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_PEApa_pruned_log.csv", quote=FALSE,sep = ",", col.names=NA)

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Filling in Taxonomy #### 

##RFTM
RFTM_log <- read.csv("/Users/dianaportugal/Desktop/16S/Condensed Materials/Oyster-16S-22/Log2Fold/3F_RFTMpa_pruned_log.csv")
#Transposing the data set 
RFTM_log <- t(RFTM_log)
RFTM_log <- as.data.frame(RFTM_log)
RFTM_log <- tibble::rownames_to_column(RFTM_log, "Rank") #Make 1 column for row names
RFTM_log[1, 1] <- " " #remove X title
RFTM_log[2, 1] <- " " #remove Log2fold title
View(RFTM_log)
#For loop to get previous name 
for (i in colnames(RFTM_log)[2:NCOL(RFTM_log)]){
  RFTM_log[[i]] <- (str_c(RFTM_log$Rank, "_", RFTM_log[[i]]))}
RFTM_log <- fill(RFTM_log, names(RFTM_log)) #Fill func
RFTM_log <- t(RFTM_log) #Transpose
RFTM_log[1, 1] <- "OTU" #Adding OTU label
RFTM_log[1, 2] <- "L2F" #Adding L2F label
colnames(RFTM_log) <- RFTM_log[1,] #Setting column names
RFTM_log <- RFTM_log[-1,] #Removing column names template column 
rownames(RFTM_log) = NULL #Removing rownames 
RFTM_log <- as.data.frame(RFTM_log)
RFTM_log$OTU <- gsub("_","",as.character(RFTM_log$OTU)) #Removing the extra "_" in the begining of the column
RFTM_log$L2F <- gsub("_","",as.character(RFTM_log$L2F)) #Removing the extra "_" in the begining of the column
RFTM_log$L2F <- as.numeric(as.character(RFTM_log$L2F)) 
RFTM_log$L2F <- signif(RFTM_log$L2F, 4) 
View(RFTM_log)


## PEACRABS
PEA_log <- read.csv("/Users/dianaportugal/Desktop/16S/Condensed Materials/Oyster-16S-22/Log2Fold/3F_PEApa_pruned_log.csv")
#Transposing the data set 
PEA_log <- t(PEA_log)
PEA_log <- as.data.frame(PEA_log)
PEA_log <- tibble::rownames_to_column(PEA_log, "Rank") #Make 1 column for row names
PEA_log[1, 1] <- " " #remove X title
PEA_log[2, 1] <- " " #remove Log2fold title
for (i in colnames(PEA_log)[2:NCOL(PEA_log)]){
  PEA_log[[i]] <- (str_c(PEA_log$Rank, "_", PEA_log[[i]]))}
PEA_log <- fill(PEA_log, names(PEA_log)) #Fill func
PEA_log <- t(PEA_log) #Transpose
PEA_log[1, 1] <- "OTU" #Adding OTU label
PEA_log[1, 2] <- "L2F" #Adding L2F label
colnames(PEA_log) <- PEA_log[1,] #Setting column names
PEA_log <- PEA_log[-1,] #Removing column names template column 
rownames(PEA_log) = NULL #Removing rownames 
PEA_log <- as.data.frame(PEA_log)
PEA_log$OTU <- gsub("_","",as.character(PEA_log$OTU)) #Removing the extra "_" in the begining of the column
PEA_log$L2F <- gsub("_","",as.character(PEA_log$L2F)) #Removing the extra "_" in the begining of the column
PEA_log$L2F <- as.numeric(as.character(PEA_log$L2F)) 
PEA_log$L2F <- signif(PEA_log$L2F, 4) 
View(PEA_log)

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # Log2Fold Taxa ####

RFTMCC <- length(unique(RFTM_log$Genus))
RFTMgp <- colorRampPalette(brewer.pal(8, "Set1"))

ggplot(RFTM_log,aes(x = OTU, y=L2F, fill = Genus)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=L2F), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values = RFTMgp(RFTMCC))+
  scale_y_continuous(limits = c(-5, 5))+
  geom_hline(yintercept = 0)+
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=8),
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(face = "bold", size = 15, margin=margin(15, 15, 15, 15)),
        axis.title.y = element_text(face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        panel.border = element_rect(color = "black",fill = NA,size = 1))+
  labs(x = "OTU",
       y = "Log2Fold Change")
#ggsave(filename = "RFTM_Log2Fold_barplot.jpeg", plot=last_plot(), path="/Users/dianaportugal/Desktop", width = 11, height = 8)  


PEACC <- length(unique(PEA_log$Genus))
PEAgp <- colorRampPalette(brewer.pal(8, "Set1"))

ggplot(PEA_log,aes(x = OTU, y=L2F, fill = Genus)) +
  geom_bar(stat="identity")+
  #scale_fill_brewer(palette="Dark2")+
  geom_text(aes(label=L2F), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values = PEAgp(PEACC))+
  scale_y_continuous(limits = c(-31, 31))+
  geom_hline(yintercept = 0)+
  theme(legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=8),
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(face = "bold", size = 15, margin=margin(15, 15, 15, 15)),
        axis.title.y = element_text(face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        panel.border = element_rect(color = "black",fill = NA,size = 1))+
  labs(x = "OTU",
       y = "Log2Fold Change")
#ggsave(filename = "RFTM_Log2Fold_barplot.jpeg", plot=last_plot(), path="/Users/dianaportugal/Desktop", width = 11, height = 8)  

----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Comp_RFTM_log <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_RFTMpa_compLog.csv")
Comp_PEA_log <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_PEApa_compLog.csv")

onethird_PEA_log <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_PEApa_pruned_log.csv")
onethird_RFTM_log <- read.csv("Documents/Smithsonian_2021/Oyster_16S/Log2Fold/3F_RFTMpa_pruned_log.csv")

F3_RFTMPA_wlog_prune_tax
F3_PEAPA_wlog_prune_tax


#Graphing - Heat Tree####
## RFTM
#Complete with log
F3_RFTMPA_wlog <- phyloseq(OTU, F3_RFTMPA_wlog_prune_tax, SAMP) 
F3_RFTMPA_wlog

#Same thing without log2fold (removing the first column)
sigPA_RFTM_tax3_r <- F3_RFTMPA_wlog_prune_tax[,-1]
view(sigPA_RFTM_tax3_r)

F3_RFTMpa_pruned <- phyloseq(OTU, sigPA_RFTM_tax3_r, SAMP) 
F3_RFTMpa_pruned

taxmapRFTMPA <- parse_phyloseq(F3_RFTMpa_pruned)
tax_table(F3_RFTMpa_pruned)
view(tax_table(F3_RFTMpa_pruned))
F3_RFTMpa_pruned

taxmapRFTMPA %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapRFTMPA),
            node_color = n_obs(taxmapRFTMPA),
            layout = "davidson-harel", initial_layout = "reingold-tilford", 
            node_legend_title	= "Taxon Observation", 
            node_color_axis_label ="Taxon Observation",
            margin_size = c(0.1, 0.1, 0.1, 0.1))+
  theme(legend.position = "none", 
        legend.text=element_text(size=8),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  labs(title = "OTU Taxonomy Abundance for Oysters with RFTM", 
       subtitle = "Presence Absence",
       caption = "Data source: Oyster 16s 2017")
#ggsave(filename = "RFTM_Comp_taxa.jpeg", plot=last_plot(), path="/Users/dianaportugal/Desktop", width = 11, height = 8)  


## Peacrabs
#Complete with log
F3_PEAPA_wlog <- phyloseq(OTU, F3_PEAPA_wlog_prune_tax, SAMP) 
F3_PEAPA_wlog

#Same thing without log2fold (removing the first column)
sigPA_PEA_tax3_r <- F3_PEAPA_wlog_prune_tax[,-1]
view(sigPA_PEA_tax3_r)

F3_PEApa_pruned <- phyloseq(OTU, sigPA_PEA_tax3_r, SAMP) 
F3_PEApa_pruned

taxmapPEAPA <- parse_phyloseq(F3_PEApa_pruned)
tax_table(F3_PEApa_pruned)
view(tax_table(F3_PEApa_pruned))
F3_PEApa_pruned

taxmapPEAPA %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapPEAPA),
            node_color = n_obs(taxmapPEAPA),
            layout = "davidson-harel", initial_layout = "reingold-tilford", 
            node_legend_title	= "Taxon Observation", 
            node_color_axis_label ="Taxon Observation",
            margin_size = c(0.1, 0.1, 0.1, 0.1))+
  theme(legend.position = "none", 
        legend.text=element_text(size=8),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  labs(title = "OTU Taxonomy Abundance for Oysters with Peacrabs", 
       subtitle = "Presence Absence",
       caption = "Data source: Oyster 16s 2017")
#ggsave(filename = "Peacrabs_Comp_taxa.jpeg", plot=last_plot(), path="/Users/dianaportugal/Desktop", width = 11, height = 8)  
