#### Taxonomy Associated with Growth ####
# 2022-09-22
# Author: Monserrat Garcia

# Packages ####
require(microViz)
require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(tidyr)
require(DESeq2)
require(data.table)
library(genefilter)
require(ggplot2)
require(tidyverse)
require(phyloseq)
library(metacoder)

# Filtering out OTU tables 

physeq_count18 <- readRDS("Data/physeq_count18.rds")
physeq_count18

physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

#2017: Site and Scaled Volume
factor_site <- as.factors(sample_data(physeq_count17)$Site.x)

physeq_count17 = subset_samples(physeq_count17, Volume_scale != "NA")
des_count17 <- phyloseq_to_deseq2(physeq_count17, ~ Volume_scale+Site.x)
des_count17 <- DESeq(des_count17, test="Wald", fitType = "parametric")

results17 <- results(des_count17, name = "Volume_scale")
results17
significant17 <- results17[which(results17$padj <0.05), ]
sigtab17 = cbind(as(significant17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_17_vol <- subset(sigtab17, select = -c(baseMean, lfcSE, stat, pvalue, padj))


###### Log2fold Change  ####

#####Bar Plot#####
mycolors1= brewer.pal(11, "Paired")

plot_bar(physeq_v17, x= "Volume_scale", fill= "Order") +
  geom_bar(aes(color=Order, fill = Order), stat = "identity", position = "stack") +
  scale_fill_manual(values = mycolors1) +
  scale_color_manual(values = mycolors1) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_blank(), 
        text = element_text(size=10))

# Turning into a phyloseq ####
physeq_count17
taxa_volume <- as.matrix(sigtab_17_vol)
taxtable_volume <- tax_table(taxa_volume)
view(taxtable_volume)
physeq17 = subset_taxa(prune_taxa(rownames(taxtable_volume), physeq_count17))

## Checking Differences ##
physeq_count17
physeq17_vol

# with site = 3871, and taxa is 10 

physeq_si17_vol = tax_filter(physeq17_vol, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_si17_vol

genefilter_17vol <- genefilter_sample(physeq17_vol, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq17_vol))
genefilter_17vol
# 2193, 6

physeq17_newvol = prune_taxa(genefilter_17vol, physeq17_vol)
physeq17_newvol
# 8 taxa 

# Adding 1 to OTU Tables due to Error
otu <- otu_table(physeq17_newvol)
otu17 <- otu +1
View(otu17)

taxatable_df <- as.data.frame(tax_table(physeq17_v))
view(taxatable_df)
sigvol_17 <- merge(taxatable_df, sigtab_17_vol, by ='row.names', all = TRUE)
view(sigvol_17)
newsigvol_17 <- subset.data.frame(sigvol_17, Kingdom.x != "NA")
view(newsigvol_17)
rownames(newsigvol_17)= sss17$Row.names
sigvolume_17 <- subset(newsigvol_17, select = -c(Row.names,Kingdom.y, Phylum.y, Class.y, Order.y, Genus.y.y, Species.y, Family.y, Genus.x.y))

log2_17vol <- merge(sigvolume_17, sigtab_17_vol,by ="row.names")
log2fold_volume17 <- subset(log2_17vol, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Genus.y.x, Species.x, Family.x, Genus.x.x, log2FoldChange.y))

write.csv(log2fold_volume17, file = "Data/log2fold2017.csv")


###### Turning into Phyloseq Object #####

meta17_data <- read.csv("Data/meta17_data_update.csv")

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

#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")

#Setting OTU, TAX, and SAMP
OTU17= otu17

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)

physeq_v17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_v17

saveRDS(physeq_v17, "Data/physeq_sitexvol17.rds")

taxtable_17 = as.data.frame(tax_table(physeq17_v))

taxtable17 = as.data.frame(tax_table(physeq_v17))

write.csv(taxtable_17, file = "Data/taxtable_sitexvol17.csv")

##### Heatmap #####

plot_heatmap(physeq_v17, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", taxa.label = "Phylum", sample.label = "Volume_scale", sample.order = "Volume_scale")

sub_significant17 <- subset_taxa(prune_taxa(rownames(significant17), physeq_count17))

#### Positive OTUs

significant17_pos <- significant17[significant17$log2FoldChange>0,]
dim(significant17_pos)
dd <- as.data.frame(significant17_pos)
# 25, 6

sub_significant17_pos <- subset_taxa(prune_taxa(rownames(significant17_pos), physeq_count17))
sub_significant17_pos

heatmap17 = parse_phyloseq(physeq17_v)

heatmap17 %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")
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


#2018 Volume and Bucket ####
physeq_count18 = subset_samples(physeq_count18, Volume_scale != "NA")
des_count18 <- phyloseq_to_deseq2(physeq_count18, ~ Volume_scale + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18 = estimateSizeFactors(des_count18, geoMeans=geoMeans, locfunc=shorth)
deseq18_v = DESeq(deseq18, test="Wald", fitType="parametric")

results18 <- results(deseq18_v, name = "Volume_scale")
significant18 <- results18[which(results18$padj <0.05), ]
sigtab18 = cbind(as(significant18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant18), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_18_vol <- subset(sigtab18, select = -c(baseMean, 
                                              lfcSE, stat, pvalue, padj))


# Turning it into a phyloseq ####
physeq_count18
taxa18 <- as.matrix(sigtab_18_vol)
taxa_vol18 <- tax_table(taxa18)
physeq18 = subset_taxa(prune_taxa(rownames(taxa_vol18), physeq_count18)) 

#Checking the Difference 
physeq_count18
physeq18

# With site = 2471, 5 taxa  
physeq_si18 = tax_filter(physeq18, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_si18

genefiltervol_18 <- genefilter_sample(physeq18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq18))
# 1277, 6

physeq18_v = prune_taxa(genefiltervol_18, physeq18)
physeq18_v
# 5 taxa 

volume_df18 <- as.data.frame(tax_table(physeq18_v))
sigvol18 <- merge(volume_df18, sigtab_18_vol, by ='row.names', all = TRUE)
newsigvol18 <- subset.data.frame(sigvol18, Kingdom.x != "NA")
sigvolume_18 <- subset(newsigvol18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(sigvolume_18, file = "Data/log2fold2018_volxbuck.csv")

# Adding 1 to OTU Table due to Error
otu18 <- otu_table(physeq_si18)
otu_volume18 <- otu18 +1
View(otu_volume18)

# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run123_taxa"
Run123_taxa$V1=NULL
rownames(table18)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Setting taxmat and otumat
taxmat18=table18
taxmat18=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix18=as.matrix(taxmat18, rownames = "V2")
colnames(tax_matrix18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")

#Setting OTU, TAX, and SAMP
OTU18= otu_v18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

physeq_sig18 = phyloseq(OTU18, TAX18, SAMP18)
physeq_sig18

saveRDS(physeq_sig18, "Data/physeq_volxbuck18.rds")
taxtable_18 <- as.data.frame(tax_table(physeq_si18))

taxtable18 <- as.data.frame(tax_table(physeq_sig18))

write.csv(taxtable18, file = "Data/taxtable_volxbuck18.csv")

# Graphs ####
plot_heatmap(physeq_sig18, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", taxa.label = "Family", sample.label = "Volume_scale", sample.order = "Volume_scale")

ggplot(sigtab18, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2018")

plot_richness(physeq_sig18, x= "Volume_scale", color = "Bucket2", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for Treatment and Species 2017")

# Taking out positive OTU's
neg_otus_weight <- sigtab18 %>% 
  filter(log2FoldChange < 0)
pos_otus_vol <- sigtab18 %>% 
  filter(log2FoldChange > 0)

pos_otus <- subset(pos_otus_vol, select = -c(baseMean, log2FoldChange, 
                                             lfcSE, stat, pvalue, padj))

write.csv(pos_otus, file = "Data/po_vol18.csv")
pos_otus17 <- read.csv("Data/pos_otus17.csv")



Phy.ord <- ordinate(physeq_sig18, "NMDS", "bray")
plot_ordination(physeq_sig18, Phy.ord, type = "biplot", color = "Volume_scale", shape = "Bucket2", title = "biplot")

physeq18_v <- readRDS("Data/physeq_volxbuck18.rds")

heatmap18_vol = parse_phyloseq(physeq18_v)

heatmap18_vol %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")

taxtable18 = as.data.frame(tax_table(physeq18_v))


# Weight Analysis#####

physeq_weight = subset_samples(physeq_count17, Weight_delta != "NA")
des_weight17 <- phyloseq_to_deseq2(physeq_weight, ~ Weight_delta + Site.x)
des_weight17 <- DESeq(des_weight17, test="Wald", fitType = "parametric")

results_weight17 <- results(des_weight17, name = "Weight_delta")
significant_weight17 <- results_weight17[which(results_weight17$padj <0.05), ]
sigtab_weight17 = cbind(as(significant_weight17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant_weight17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_weigh17 <- subset(sigtab_weight17, select = -c(baseMean, log2FoldChange, 
                                                      lfcSE, stat, pvalue, padj))

physeq_count17
taxa_weigh <- as.matrix(sigtab_weigh17)
taxa_weight <- tax_table(taxa_weigh)
physeq_weight17 = subset_taxa(prune_taxa(rownames(taxa_weight), physeq_count17))
physeq_count17
physeq_weight17

physeq_filter_weight17 = tax_filter(physeq_weight17, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_filter_weight17
sigfigweight_17 <- genefilter_sample(physeq_weight17, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_weight17))
# 820, 6

phy_weight17 = prune_taxa(sigfigweight_17, physeq_weight17)
# 8 taxa 

#### *Note: All the taxa were filtered out when filtering for 1/3 of samples and at least appeared once ####

# 2018 Weight

physeq_weight18 = subset_samples(physeq_count18, Weight_delta !="NA")

des_weight18 <- phyloseq_to_deseq2(physeq_weight18, ~ Weight_delta + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_weight18 = estimateSizeFactors(des_weight18, geoMeans=geoMeans, locfunc=shorth)
des_weight18 <- DESeq(deseq_weight18, test="Wald", fitType = "parametric")

results_weight18 <- results(des_weight18, name = "Weight_delta")
significant_weight18 <- results_weight18[which(results_weight18$padj <0.05), ]
sigtab_weight18 = cbind(as(significant_weight18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_weight18), ], "matrix"))

sigtab_18_weight <- subset(sigtab_weight18, select = -c(baseMean, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_w18 <- as.matrix(sigtab_18_weight)
taxa_weight18 <- tax_table(taxa_w18)
physeq_weight18 = subset_taxa(prune_taxa(rownames(taxa_weight18), physeq_count18))
physeq_count18
physeq_weight18

weight18_df <- as.data.frame(tax_table(phy_weight18))
sigweigh_18 <- merge(weight18_df, sigtab_18_weight, by ='row.names', all = TRUE)
newsigweigh_18 <- subset.data.frame(sigweigh_18, Kingdom.x != "NA")
sigweight_18 <- subset(newsigweigh_18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(sigweight_18, file = "Data/log2fold2018_weight.csv")

# Weight and bucket analysis gave 2730 and 6 taxa 
phys_weight18 = tax_filter(physeq_weight18, min_prevalence = 0.3, min_sample_abundance = 1)
phys_weight18

sigfigweight_18 <- genefilter_sample(physeq_weight18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_weight18))
# 1302, 6

phy_weight18 = prune_taxa(sigfigweight_18, physeq_weight18)
phy_weight18
#7, 8

# Adding 1 to OTU Table due to Error
otu_weight18 <- otu_table(phy_weight18)
otu_w18 <- otu_weight18 +1
View(otu_w18)

#Setting OTU, TAX, and SAMP with new OTU table
OTU18= otu_w18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

phy_weight_sig18 = phyloseq(OTU18, TAX18, SAMP18)
phy_weight_sig18

saveRDS(phy_weight_sig18, "Data/physeq_buckxweight18.rds")

physeq_weight_sig18

taxtable_weight18 <- as.data.frame(tax_table(phy_weight_sig18))

taxtable_BuckxWeigh <- as.data.frame(tax_table(phy_weight_sig18))

write.csv(taxtable_BuckxWeigh, file = "Data/taxtable_Bucket&Weight.csv")

# Plot Analysis ####

physeq_buckxweight18 <- readRDS("Data/physeq_buckxweight18.rds")
physeq_buckxweight18

heatmap = parse_phyloseq(physeq_buckxweight18)

heatmap %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")


# Merging Tax Tables with significant OTUs for Weight and Volume ####

#### Loading Tax Tables ####
taxtable_BuckxWeigh <- read.csv("Data/taxtable_Bucket&Weight.csv")
taxtable_sitexvol17 <- read.csv("Data/taxtable_sitexvol17.csv")
tabtable_volxbuck18 <- read.csv("Data/taxtable_volxbuck18.csv")

#### Combining Tax Tables ####
sig_OTU_combined <- rbind(taxtable_BuckxWeigh, taxtable_sitexvol17, tabtable_volxbuck18)

write.csv(sig_OTU_combined, file = "Data/sig_OTU_combined.csv")



# Non-scaled Volume ####

physeq_count17 = subset_samples(physeq_count17, Volume_delta != "NA")
des_volume17 <- phyloseq_to_deseq2(physeq_count17, ~ Volume_delta+Site.x)
des_volume17 <- DESeq(des_volume17, test="Wald", fitType = "parametric")

results_volume17 <- results(des_volume17, name = "Volume_delta")
results_volume17
significant_volume17 <- results_volume17[which(results_volume17$padj <0.05), ]
sigtab_volume17 = cbind(as(significant_volume17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant_volume17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_Volume17 <- subset(sigtab_volume17, select = -c(baseMean, log2FoldChange, 
                                                       lfcSE, stat, pvalue, padj))

# Turning it into a phyloseq ####
physeq_count17
taxa_volume17 <- as.matrix(sigtab_Volume17)
taxa_vol17 <- tax_table(taxa_volume17)
physeq_volume17 = subset_taxa(prune_taxa(rownames(taxa_vol17), physeq_count17))

#Checkign Difference
physeq_count17
physeq_volume17

#Look for which sites are more prominent or appear more in the significant OTUs#

physeq_sigvol17 = tax_filter(physeq_volume17, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_sigvol17

sf_vol17 <- genefilter_sample(physeq_volume17, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_volume17))
# 2193, 6

physeq17_volume = prune_taxa(sf_vol17, physeq_volume17)
physeq17_volume

## *Note: All 2017 filtered out with the filter settings above ##

physeq_volume18 = subset_samples(physeq_count18, Volume_delta !="NA")

des_volume18 <- phyloseq_to_deseq2(physeq_volume18, ~ Volume_delta + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_volume18 = estimateSizeFactors(des_volume18, geoMeans=geoMeans, locfunc=shorth)
des_volume18 <- DESeq(deseq_volume18, test="Wald", fitType = "parametric")

results_volume18 <- results(des_volume18, name = "Volume_delta")
significant_volume18 <- results_volume18[which(results_volume18$padj <0.05), ]
sigtab_volume18 = cbind(as(significant_volume18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_volume18), ], "matrix"))

sigtab_18_volume <- subset(sigtab_volume18, select = -c(baseMean, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_v18 <- as.matrix(sigtab_18_volume)
taxa_volume18 <- tax_table(taxa_v18)
physeq_volume18 = subset_taxa(prune_taxa(rownames(taxa_volume18), physeq_count18))

# Checking Difference
physeq_count18
physeq_volume18


deltavol_df18 <- as.data.frame(tax_table(physeq18_volume))
delvol_18 <- merge(deltavoldf18, sigtab_volume18, by ='row.names', all = TRUE)
newdelvol_18 <- subset.data.frame(delvol_18, Kingdom.x != "NA")
deltavolume_18 <- subset(newdelvol_18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.y.x,Genus.x.x, Genus.y.y,Species.x,baseMean,lfcSE, stat,pvalue,padj))
write.csv(deltavolume_18, file = "Data/log2fold2018_voldelta.csv")



physeq_sigvol18 = tax_filter(physeq_volume18, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_sigvol18

sigdel_vol18 <- genefilter_sample(physeq_volume18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_volume18))
# 3 taxa

physeq18_volume = prune_taxa(sigdel_vol18, physeq_volume18)
physeq18_volume

taxtable_voldelxbuck18 <- as.data.frame(tax_table(physeq_sigvol18))


heatmap18 = parse_phyloseq(physeq_sigvol18)

heatmap18 = parse_phyloseq(taxtable_sitexvol17)

heatmap18 %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")
write.csv(taxtable_voldelxbuck18, file = "Data/taxtable_voldelxbuck18.csv")

taxtable_voldelxbuck18 <- read.csv("Data/taxtable_voldelxbuck18.csv")


# Sites in the Significant OTUs ####

physeq_buckxweight18 <- readRDS("Data/physeq_buckxweight18.rds")
site_data18 <- sample_data(physeq_buckxweight18)


physeq_sitexvol17 <- readRDS("Data/physeq_sitexvol17.rds")
site_data17 <- sample_data(physeq_sitexvol17)

physeq_volxbuck18 <- readRDS("Data/physeq_volxbuck18.rds")

##### Log2Fold Change Graphs ####
log2fold18_vol <- read.csv("Data/log2fold2018_volxbuck.csv")
log2fold18_vol

log2fold18_weight <- read.csv("Data/log2fold2018_weight.csv")
log2fold2017 <- read.csv("Data/log2fold2017.csv")


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

##### 2018 Volume and Bucket ####
x = tapply(log2fold18_vol$log2FoldChange, log2fold18_vol$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold18_vol$Genus = factor(as.character(log2fold18_vol$Genus), levels=names(x))

ggplot(log2fold18_vol, aes(x=Genus, y=log2FoldChange, color=Genus)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ labs(title = "Log2FoldChange Volume 2018")+ ylab("Log2FoldChange")+geom_bar(stat="identity", aes(fill = Genus))

##### 2018 Weight ####

x = tapply(log2fold18_weight$log2FoldChange, log2fold18_weight$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold18_weight$Genus = factor(as.character(log2fold18_weight$Genus), levels=names(x))

ggplot(log2fold18_weight, aes(x=Genus, y=log2FoldChange, color=Genus)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ labs(title = "Log2FoldChange Weight 2018")+ylab("Log2FoldChange")+geom_bar(stat="identity", aes(fill = Genus))

##### 2017 ####

x = tapply(log2fold2017$log2FoldChange.x, log2fold2017$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold2017$Genus = factor(as.character(log2fold2017$Genus), levels=names(x))

ggplot(log2fold2017, aes(x=Genus, y=log2FoldChange.x, color=Genus)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Log2FoldChange Volume 2017") +  ylab("Log2FoldChange")+geom_bar(stat="identity", aes(fill = Genus))
