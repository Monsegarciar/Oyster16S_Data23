#### 2017 Analysis
# 2023-26-03
# Author: Monserrat Garcia 

BiocManager::install("ggtree")
# Required Packages ####
require(phyloseq)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggtree)
require(metacoder)

#Loading Data 

meta17_data <- read.csv("Data/meta17_data_update.csv")


asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")
Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")


#Changing row names in "Run23_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))


#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]
otumat17=asvtable_17

#Converting to matrix
otu_matrix17= as.matrix(otumat17, rownames = rownames(asvtable_17))

tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")
meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= otu_table(otu_matrix17, taxa_are_rows = FALSE)

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)


OTU_count17=transform_sample_counts(OTU17, function(x) 1E6 * x/sum(x))


physeq_class17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_class17

physeq_count17 = phyloseq(OTU_count17, TAX17, SAMP17)
physeq_count17

