#### 2018 Analysis
# 2023-26-03
# Author: Monserrat Garcia 

# Phyloseq Analysis ####

#Loading Data

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")
Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run123_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
taxmat18=Run123_taxa
taxmat18=Run123_taxa[-c(1)]
otumat18=asvtable_18

#Converting to matrix
otu_matrix18= as.matrix(otumat18, rownames = "V1")

tax_matrix18=as.matrix(taxmat18, rownames = "V2")
colnames(tax_matrix18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18= otu_table(otu_matrix18, taxa_are_rows = FALSE)

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

OTU_count18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))

physeq_class18 = phyloseq(OTU18, TAX18, SAMP18)
physeq_class18

physeq_count18 = phyloseq(OTU_count18, TAX18, SAMP18)
physeq_count18



