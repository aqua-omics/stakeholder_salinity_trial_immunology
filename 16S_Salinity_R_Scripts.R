###Please see Figure_Scripts_FiMS.R for actual scripts used for figures, those seen here are precursors that underwent adjustments for publication

###PREPPING YOUR R ENVIRONMENT###

##Set your working directory##

#You are gonna be exporting a lot of stuff from R (tables, pictures, graphs, etc...), it is best to be somewhere where you are happy to place all those items

#Use getwd() to find out where you currently are
getwd()

#Use setwd() to set your working directory to whereever you want
setwd("C:/Users/dbrad/Documents/Bioinformatic_Analysis/LAB_salinity_trial/16S_R")

#Rerun getwd() to check it worked
getwd()

##Load up SummarySE function##

#To run this put cursor in front of summarySE below and click Run above. This is a function created by someone else that we borrow, never seen a citation for it.#

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Install packages##

#Some packages require BiocManager to be installed before being able to actually install the package you want, think of it as like a non-default conda channel#
install.packages("BiocManager")

#phyloseq is the main package you will use to do most of the analysis scripts. Like QIIME2 it is a wrapper of a lot of other tools. You will also get vegan (statistics), ggplot2(graphics), and plyr(data manipulation) by installing this#
BiocManager::install("phyloseq")

#Install other useful packages#
install.packages(c("tidyverse", "FSA", "dplyr", "reshape", "rcompanion"))
BiocManager::install("DESeq2")
BiocManager::install("metagenomeSeq")
devtools::install_github("yanlinlin82/ggvenn")

##Load packages##
library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(tidyverse)
library(FSA)
library(dplyr)
library(reshape)
library(DESeq2)
library(metagenomeSeq)
library(ggvenn)
library(rcompanion)

###GETTING YOUR DATA INTO R###

##Import files##
BIOM <- import_biom(file.choose()) #phyloseq.biom
TREE =  read_tree(file.choose()) #tree/tree.nwk
META <- import_qiime_sample_data(file.choose()) #16S_metadata.txt

#Upload metada in human-readable format in R#
META2  = read.delim(file.choose(), row.names=1) #16S_metadata.txt

#check that sample names are the same between all files, if they do not match, then fix it in the mapping file, honestly just because it is easier# 
sample_names(META)
sample_names(BIOM)

#Merge three items into one phyloseq object#
data_w_LSE1DP2a <- merge_phyloseq (BIOM,TREE,META)
data_w_LSE1DP2a

###DATA MANIPULATION AND SUMMARIZING SEQUENCES### 

##Change taxa names if needed##

#See what present names are#
colnames(tax_table(data_w_LSE1DP2a))

#Switch them to classic format#
colnames(tax_table(data_w_LSE1DP2a))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

#Check that it worked#
colnames(tax_table(data_w_LSE1DP2a))

##Visualize bad sample##

#Check number of taxa and samples#
ntaxa(data_w_LSE1DP2a) #30337#
nsamples(data_w_LSE1DP2a) #195#

#Filter Eukaryotes
data_w_LSE1DP2a <- subset_taxa(data_w_LSE1DP2a, Kingdom!="d__Eukaryota")
data_w_LSE1DP2a
ntaxa(data_w_LSE1DP2a) #30255#
nsamples(data_w_LSE1DP2a) #195#

#Filter low abudance sequences
Fdata_w_LSE1DP2a <- filter_taxa(data_w_LSE1DP2a, function(x) sum(x) >10, TRUE)
ntaxa(Fdata_w_LSE1DP2a) #15688#

#Subset samples to explore dataset for contamination
W_Fdata_w_LSE1DP2a <- subset_samples(Fdata_w_LSE1DP2a, Sample_Type =="Water")
nsamples(W_Fdata_w_LSE1DP2a) #126#

W_1D_Fdata_w_LSE1DP2a = subset_samples(W_Fdata_w_LSE1DP2a, Days_Post_Hatch=="D1")
nsamples(W_1D_Fdata_w_LSE1DP2a)#27

W_1D_w_LSE1DP2a_phy <- transform_sample_counts(W_1D_Fdata_w_LSE1DP2a, function(x) 100*x/sum(x))

##Create table ready for making stacked bar graph for Genuses <1%##
# agglomerate taxa
glom <- tax_glom(W_1D_w_LSE1DP2a_phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Sample_ID, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Sample_ID","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Sample_ID, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Sample_ID, Abundance, Genus))
#combine with original table
Genus_W_1D_w_LSE1DP2a_Sample_ID <- rbind(dat, Abundance)

spatial_plot_Genus_W_1D_w_LSE1DP2a_data_Sample_ID <- ggplot(data=Genus_W_1D_w_LSE1DP2a_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "rosybrown"))+
  ggtitle("1 DPH Water Genera by Sample ID with LSE1DP2a") +
  xlab("Sample ID") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) 

spatial_plot_Genus_W_1D_w_LSE1DP2a_data_Sample_ID

tiff('1DPH Water Genera by Sample ID with LSE1DP2a.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Genus_W_1D_w_LSE1DP2a_data_Sample_ID
dev.off()

#Make a distance matrix#
W_1D_w_LSE1DP2a_Bray<-phyloseq::distance(W_1D_Fdata_w_LSE1DP2a, "bray")

W_1D_w_LSE1DP2a_Bray_df <- as.data.frame(as.matrix(W_1D_w_LSE1DP2a_Bray))

mean(W_1D_w_LSE1DP2a_Bray_df$LSE1DP2a) #0.961217
mean(W_1D_w_LSE1DP2a_Bray_df$LSE1DP2b) #0.5428415
mean(W_1D_w_LSE1DP2a_Bray_df$LSE1DP2c) #0.5642056

#Remove a sample that has been proven to be weird
data <- subset_samples(vdata_w_LSE1DP2a, Sample_ID != "LSE1DP2a")
data <- filter_taxa(data, function(x) sum(x) >0, TRUE)

ntaxa(data) #30040#
nsamples(data) #194#


#Filter taxa that have occurred less than 10 times across all samples, this is an alternative to Deblur's default setting of removing any taxa that did not show up in at least 10 samples. For my studies, I had some sample types that has less than 10 occurances and I did not want to lose anything that could have been exclusive to those samples. 
Fdata = filter_taxa(data, function(x) sum(x) >10, TRUE)
ntaxa(Fdata) #15642#

write.csv(tax_table(Fdata), 'Fdata_taxa.csv')

#Check number of samples
nsamples(Fdata) #Should be 194#

##Summarize sequences by ASVs (filtered and unfiltered)## 

#Unfiltered version first#

#export table of sequence sums#
data_seqs_per_ESV <- as.data.frame(taxa_sums(data))

#Change colnames to Sequences#
colnames(data_seqs_per_ESV) <- c("Sequences")
sum(data_seqs_per_ESV$Sequences) #18857782#

#Use cbind to get taxonomy added in#
data_seqs_per_ESV = cbind(as(data_seqs_per_ESV, "data.frame"), as(tax_table(data)[rownames(data_seqs_per_ESV), ], "matrix"))

#export it as a csv, and tell it that the the rownames should be named OTUID#
write.csv(data.frame("OTUID" =rownames(data_seqs_per_ESV), data_seqs_per_ESV) , "seqs_per_ESV.csv", row.names=FALSE)

#Filtered version second#

Fdata_seqs_per_ESV <- as.data.frame(taxa_sums(Fdata))
colnames(Fdata_seqs_per_ESV) <- c("Sequences")
sum(Fdata_seqs_per_ESV$Sequences) #18794113#
Fdata_seqs_per_ESV = cbind(as(Fdata_seqs_per_ESV, "data.frame"), as(tax_table(Fdata)[rownames(Fdata_seqs_per_ESV), ], "matrix"))
write.csv(data.frame("OTUID" =rownames(Fdata_seqs_per_ESV), Fdata_seqs_per_ESV) , "filtered_seqs_per_ESV.csv", row.names=FALSE)

##Summarize sequences by samples (filtered and unfiltered)##

#Unfiltered first
data_seqs_per_sample <- as.data.frame(sample_sums(data))
colnames(data_seqs_per_sample) <- c("Full_Sequences")
sum(data_seqs_per_sample$Full_Sequences) #18857782#
Fdata_seqs_per_sample <- as.data.frame(sample_sums(Fdata))
colnames(Fdata_seqs_per_sample) <- c("Trimmed_Sequences")
sum(Fdata_seqs_per_sample$Trimmed_Sequences) #18798301#

#Sum of sequences should be same as above. Instead of exporting two different tables you can instead combine them with cbind and export that#
Comdata_seqs_per_sample = cbind(as(data_seqs_per_sample, "data.frame"), as(Fdata_seqs_per_sample, "data.frame"))
write.csv(data.frame("OTUID" =rownames(Comdata_seqs_per_sample), Comdata_seqs_per_sample) , "sequences_per_sample.csv", row.names=FALSE)

##Determine numbers of various taxonomic levels for all samples

##INFORMATION USED IN SUPPLEMENTARY TABLE 2

Fdatataxa <- as.data.frame(tax_table(Fdata))

Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("d__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("p__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("c__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("o__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("f__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("g__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("s__", "", x)}))

Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #38
Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #37

Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #95
Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #84

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #245
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #199

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #468
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #325

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #1048
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #680

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #1939
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #280

#Split samples between Sample Typles
W_Fdata = subset_samples(Fdata, Sample_Type=="Water")
W_Fdata #15642  taxa and 125 samples
W_Fdata = filter_taxa(W_Fdata, function(x) sum(x) >0, TRUE)
W_Fdata #8369 taxa and 125 samples

T_Fdata = subset_samples(Fdata, Sample_Type=="Tissue")
T_Fdata #15642  taxa and 69 samples
T_Fdata = filter_taxa(T_Fdata, function(x) sum(x) >0, TRUE)
T_Fdata #10386 taxa and 69 samples

#Extract out and copy otu/taxonomy tables for analysis in Excel
Uniq_Fdata_taxa_species_level <- distinct(Fdatataxa, Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE) #1946

write.csv(Uniq_Fdata_taxa_species_level, "Dada2Uniq_Fdata_taxa_species_level.csv")

W_Fdata_asv_table <- as.data.frame(otu_table(W_Fdata))

W_Fdata_asv_tax_table = cbind(as(W_Fdata_asv_table, "matrix"), as(tax_table(W_Fdata)[rownames(W_Fdata_asv_table), ], "matrix"))
clipr::write_clip(W_Fdata_asv_tax_table)

T_Fdata_asv_table <- as.data.frame(otu_table(T_Fdata))

T_Fdata_asv_tax_table = cbind(as(T_Fdata_asv_table, "matrix"), as(tax_table(T_Fdata)[rownames(T_Fdata_asv_table), ], "matrix"))
clipr::write_clip(T_Fdata_asv_tax_table)


##Determine numbers of various taxonomic levels for Water samples

##INFORMATION USED IN SUPPLEMENTARY TABLE 2

W_Fdatataxa <- as.data.frame(tax_table(W_Fdata))

W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("d__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("p__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("c__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("o__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("f__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("g__", "", x)}))
W_Fdatataxa <- data.frame(lapply(W_Fdatataxa, function(x) {gsub("s__", "", x)}))

W_Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #33
W_Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #32

W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #64
W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #55

W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #157
W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #135

W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #286
W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #215

W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #583
W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #394

W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #971
W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #139


##Determine numbers of various taxonomic levels for Tissue samples
T_Fdatataxa <- as.data.frame(tax_table(T_Fdata))

T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("d__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("p__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("c__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("o__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("f__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("g__", "", x)}))
T_Fdatataxa <- data.frame(lapply(T_Fdatataxa, function(x) {gsub("s__", "", x)}))

T_Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #37
T_Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #36

T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #93
T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #82

T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #230
T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #186

T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #438
T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #308

T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #946
T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #605

T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #1696
T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #244

###Venn Diagrams###

#Note: Although some of these venn diagrams were not used in the manuscript the information gained from them was useful

##Generate Sample Type based Venn Diagrams
#Extract out ASVs as a list from the asv table stored in phyloseq objects
water_asvs <- as.list(rownames(otu_table(W_Fdata)))
tissue_asvs <- as.list(rownames(otu_table(T_Fdata)))

#Create a table agglomerated to genus level
W_Fdatataxa_genus <- W_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
W_Fdatataxa_genus$Taxonomy <- paste(W_Fdatataxa_genus$Kingdom,W_Fdatataxa_genus$Phylum,W_Fdatataxa_genus$Class,W_Fdatataxa_genus$Order,W_Fdatataxa_genus$Family,W_Fdatataxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
water_genera_list <- as.list(W_Fdatataxa_genus$Taxonomy)

#Do it for the tissue
T_Fdatataxa_genus <- T_Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
T_Fdatataxa_genus$Taxonomy <- paste(T_Fdatataxa_genus$Kingdom,T_Fdatataxa_genus$Phylum,T_Fdatataxa_genus$Class,T_Fdatataxa_genus$Order,T_Fdatataxa_genus$Family,T_Fdatataxa_genus$Genus,sep="-")
tissue_genera_list <- as.list(T_Fdatataxa_genus$Taxonomy)

#Combine the list into a list of lists
ASV_Genera_Lists <- list('Water (8369)' = water_asvs,
                         'Tissue (10386)' = tissue_asvs,
                         'Water (583)' = water_genera_list,
                         'Tissue (946)' = tissue_genera_list)

#Create a Venn Diagram comparing ASVs and save it
ggvenn(ASV_Genera_Lists, c("Water (8369)", "Tissue (10386)"), fill_color = c("lightblue", "gray"), show_percentage = FALSE)+
  ggtitle("ASVs by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('ASV level venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

#Create a Venn Diagram comparing Genera and save it
ggvenn(ASV_Genera_Lists, c("Water (583)", "Tissue (946)"), fill_color = c("lightblue", "gray"), show_percentage = FALSE)+
  ggtitle("Genera by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Genus level venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

##Generate Salinity based Venn Diagrams
#ASV Level Salinity Venn Diagrams

S10_W_Fdata <- subset_samples(W_Fdata, Salinity=="S10")
S10_W_Fdata <- filter_taxa(S10_W_Fdata, function(x) sum(x) >0, TRUE)
water_S10_asvs <- as.list(rownames(otu_table(S10_W_Fdata)))

S20_W_Fdata <- subset_samples(W_Fdata, Salinity=="S20")
S20_W_Fdata <- filter_taxa(S20_W_Fdata, function(x) sum(x) >0, TRUE)
water_S20_asvs <- as.list(rownames(otu_table(S20_W_Fdata)))

S30_W_Fdata <- subset_samples(W_Fdata, Salinity=="S30")
S30_W_Fdata <- subset_samples(S30_W_Fdata, Days_Post_Hatch!="D1")
S30_W_Fdata <- filter_taxa(S30_W_Fdata, function(x) sum(x) >0, TRUE)
water_S30_asvs <- as.list(rownames(otu_table(S30_W_Fdata)))

W_Salinities_Lists <- list('S10 (2694)' = water_S10_asvs,
                           'S20 (3215)' = water_S20_asvs,
                           'S30 (3491)' = water_S30_asvs)

ggvenn(W_Salinities_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  ggtitle("Water ASVs by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Water ASV level Salinity venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

S10_T_Fdata <- subset_samples(T_Fdata, Salinity=="S10")
S10_T_Fdata <- filter_taxa(S10_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S10_asvs <- as.list(rownames(otu_table(S10_T_Fdata)))

S20_T_Fdata <- subset_samples(T_Fdata, Salinity=="S20")
S20_T_Fdata <- filter_taxa(S20_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S20_asvs <- as.list(rownames(otu_table(S20_T_Fdata)))

S30_T_Fdata <- subset_samples(T_Fdata, Salinity=="S30")
S30_T_Fdata <- subset_samples(S30_T_Fdata, Tank_Number!="P0")
S30_T_Fdata <- filter_taxa(S30_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S30_asvs <- as.list(rownames(otu_table(S30_T_Fdata)))

T_Salinities_Lists <- list('S10 (4234)' = tissue_S10_asvs,
                           'S20 (4174)' = tissue_S20_asvs,
                           'S30 (4899)' = tissue_S30_asvs)

ggvenn(T_Salinities_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  ggtitle("Tissue ASVs by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Tissue ASV level Salinity venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()


##Genus level Salinity Venn Diagrams

#PRECURSOR TO SUPPLEMENTARY FIGURE 3
#Water Genus Level Salinity Venn Diagram
Genus_W_Fdata <- tax_glom(W_Fdata, taxrank = "Genus", NArm = FALSE)

S10_Genus_W_Fdata <- subset_samples(Genus_W_Fdata, Salinity=="S10")
S10_Genus_W_Fdata <- filter_taxa(S10_Genus_W_Fdata, function(x) sum(x) >0, TRUE)
water_S10_genera <- as.list(rownames(otu_table(S10_Genus_W_Fdata)))

S20_Genus_W_Fdata <- subset_samples(Genus_W_Fdata, Salinity=="S20")
S20_Genus_W_Fdata <- filter_taxa(S20_Genus_W_Fdata, function(x) sum(x) >0, TRUE)
water_S20_genera <- as.list(rownames(otu_table(S20_Genus_W_Fdata)))

S30_Genus_W_Fdata <- subset_samples(Genus_W_Fdata, Salinity=="S30")
S30_Genus_W_Fdata <- subset_samples(S30_Genus_W_Fdata, Days_Post_Hatch!="D1")
S30_Genus_W_Fdata <- filter_taxa(S30_Genus_W_Fdata, function(x) sum(x) >0, TRUE)
water_S30_genera <- as.list(rownames(otu_table(S30_Genus_W_Fdata)))


W_Salinities_Genera_Lists <- list('10 ppt (352)' = water_S10_genera,
                                  '20 ppt (365)' = water_S20_genera,
                                  '30 ppt (331)' = water_S30_genera)

W_Genera_Salinity_Venn <- ggvenn(W_Salinities_Genera_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  #ggtitle("Water Genera by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))
W_Genera_Salinity_Venn

tiff('Water Genera by Salinity Venn Diagram.tiff', units="in", width=6, height=6, res=300)
W_Genera_Salinity_Venn
dev.off()

#PRECURSOR TO SUPPLEMENTARY FIGURE 6
#Larvae Genus Level Salinity Venn Diagram
Genus_T_Fdata <- tax_glom(T_Fdata, taxrank = "Genus", NArm = FALSE)

S10_Genus_T_Fdata <- subset_samples(Genus_T_Fdata, Salinity=="S10")
S10_Genus_T_Fdata <- filter_taxa(S10_Genus_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S10_genera <- as.list(rownames(otu_table(S10_Genus_T_Fdata)))

S20_Genus_T_Fdata <- subset_samples(Genus_T_Fdata, Salinity=="S20")
S20_Genus_T_Fdata <- filter_taxa(S20_Genus_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S20_genera <- as.list(rownames(otu_table(S20_Genus_T_Fdata)))

S30_Genus_T_Fdata <- subset_samples(Genus_T_Fdata, Salinity=="S30")
S30_Genus_T_Fdata <- subset_samples(S30_Genus_T_Fdata, Tank_Number!="P0")
S30_Genus_T_Fdata <- filter_taxa(S30_Genus_T_Fdata, function(x) sum(x) >0, TRUE)
tissue_S30_genera <- as.list(rownames(otu_table(S30_Genus_T_Fdata)))

T_Salinities_Genera_Lists <- list('10 ppt (668)' = tissue_S10_genera,
                                  '20 ppt (606)' = tissue_S20_genera,
                                  '30 ppt (675)' = tissue_S30_genera)

T_Genera_Salinity_Venn <- ggvenn(T_Salinities_Genera_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  #ggtitle("Tissue Genera by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))
T_Genera_Salinity_Venn

tiff('Larvae Genera by Salinity Venn Diagram.tiff', units="in", width=6, height=6, res=300)
T_Genera_Salinity_Venn
dev.off()

###ALPHA DIVERSITY###

##Calculate alpha diversity##

#Estimate Richness (7 alpha diversity metrics chosen by phyloseq, we will only focus on four because we filtered our data and thus cannot use Chao1 or ACE) and export it#
Falphadiv = estimate_richness(Fdata, split = TRUE)

#Sometimes phyloseq notices you do not have any singletons and yells at you for doing that#

write.csv(Falphadiv, file='Falphadiv.csv')

#cbind in the metadata from the filtered phyloseq#
Falphadiv_metadata = cbind(as(Falphadiv, "data.frame"), as(sample_data(Fdata)[rownames(Falphadiv), ], "data.frame"))

##Test for normality##

#Null hypothesis of Shapiro-Wilk test is that the data is normally distributed. IF the p-value is LESS THAN 0.05 then the null hypothesis is rejected and the data is likely NOT normally distributed.#
shapiro.test(Falphadiv_metadata$Shannon)
shapiro.test(Falphadiv_metadata$Observed)
shapiro.test(Falphadiv_metadata$Fisher)
shapiro.test(Falphadiv_metadata$Simpson)

##Correlation Tests##

#Test to see how well the other alpha diversity metrics are correlated to Shannon Diversity. Use the non-parametric spearman rank test due to all of the alpha diversity values not being normally distributed (p-value < 0.05). Pearson is the equivalent parametric test#
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Observed, method = "spearman", exact = FALSE)
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Fisher, method = "spearman", exact = FALSE)
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Simpson, method = "spearman", exact = FALSE)

##Multiple test adjustment##

#Make a list of all the p vales from the correlation tests, make it a dataframe, add column of adjusted p values#
corr_alpha_div_p.value=c(2.2e-16,2.2e-16,2.2e-16)
corr_alpha_div_p.value=data.frame(corr_alpha_div_p.value)
corr_alpha_div_p.value$padj <- p.adjust(corr_alpha_div_p.value$corr_alpha_div_p.value, method = "BH")
corr_alpha_div_p.value

# corr_alpha_div_p.value    padj
# 1                2.2e-16 2.2e-16
# 2                2.2e-16 2.2e-16
# 3                2.2e-16 2.2e-16

##Basic visualization of Shannon statistics##

#You have determined that Shannon is a good representation of the alpha diversity since it is correlated with the other metrics. Now can determine summary of its statistics in the four Site.Survey categories using ddply to get the mean, standard deviation, median, and interquartile range#

#Note: Samples associated with days post hat 0 and 1 are removed from the datasets because they are not relevant to the analysis since fish were nether in the tanks yet when water samples were taken, nor were they at different salinties. 

Falphadiv_no0D_metadata <- filter(Falphadiv_metadata, Days_Post_Hatch!="D0")
Falphadiv_no0D1D_metadata <- filter(Falphadiv_no0D_metadata, Days_Post_Hatch!="D1")


W_Falphadiv_metadata <- filter(Falphadiv_metadata, Sample_Type=="Water")
T_Falphadiv_metadata <- filter(Falphadiv_metadata, Sample_Type=="Tissue")

W_Falphadiv_no1D_metadata <- filter(W_Falphadiv_metadata, Days_Post_Hatch!="D1")
T_Falphadiv_no0D1D_metadata <- filter(T_Falphadiv_metadata, Days_Post_Hatch!="D0")
T_Falphadiv_no0D1D_metadata <- filter(T_Falphadiv_no0D1D_metadata, Days_Post_Hatch!="D1")

ddply(Falphadiv_metadata, .(Sample_Type), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

# Sample_Type     mean        sd   median      IQR
# 1      Tissue 3.798622 0.8593483 3.833657 1.351252
# 2       Water 2.830504 0.7100213 2.773253 1.166025

#Make a simple boxplot to get a visual summary of these statistics#
#Bars denote largest and smallest values within 1.5 times the interquartile range, middle line is the median, ends of boxes are the first and third quartiles (IQR length)#
ggplot(Falphadiv_metadata, aes(x=Sample_Type, y=Shannon, color=Sample_Type)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

ggplot(W_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("legend", values=c("deeppink", "midnightblue", "blue"))+
  ggtitle("Water Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24"))


ggplot(T_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+ scale_color_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))


ggplot(Falphadiv_no0D1D_metadata, aes(x=Sample_Type, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Shannon Diversity by Salinity") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))


###TESTS FOR TOTAL AND PAIRWISE SIGNIFICANCE###
##Kruskal-Wallis: Non-parametric test for overall significance##
#Use IF not normally distributed per Shapiro#

##INFORMATION USED IN SUPPLEMENTARY TABLE 3

kruskal.test(Shannon ~ Sample_Type, data = Falphadiv_no0D1D_metadata)

kruskal.test(Shannon ~ Salinity, data = W_Falphadiv_no1D_metadata)
kruskal.test(Shannon ~ Days_Post_Hatch, data = W_Falphadiv_no1D_metadata)
kruskal.test(Shannon ~ Tank_Number, data = W_Falphadiv_no1D_metadata)

kruskal.test(Shannon ~ Salinity, data = T_Falphadiv_no0D1D_metadata)
kruskal.test(Shannon ~ Days_Post_Hatch, data = T_Falphadiv_no0D1D_metadata)
kruskal.test(Shannon ~ Tank_Number, data = T_Falphadiv_no0D1D_metadata)


kruskal.test(Shannon ~ Salinity, data = filter(W_Falphadiv_no1D_metadata, Days_Post_Hatch == "D6"))
kruskal.test(Shannon ~ Salinity, data = filter(W_Falphadiv_no1D_metadata, Days_Post_Hatch == "D12"))
kruskal.test(Shannon ~ Salinity, data = filter(W_Falphadiv_no1D_metadata, Days_Post_Hatch == "D18"))
kruskal.test(Shannon ~ Salinity, data = filter(W_Falphadiv_no1D_metadata, Days_Post_Hatch == "D24"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D3"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D6"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D9"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D12"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D15"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D18"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D20"))
kruskal.test(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D24"))

alpha_div_test = c("Sample_Type", "W_Salinity", "W_DPH", "W_Tank#", "T_Salinity", "T_DPH", "T_Tank#", "W_6DPH_Salinity", "W_12DPH_Salinity","W_18DPH_Salinity","W_24DPH_Salinity", "T_3DPH_Salinity", "T_6DPH_Salinity","T_9DPH_Salinity","T_12DPH_Salinity", "T_15DPH_Salinity", "T_18DPH_Salinity","T_20DPH_Salinity","T_24DPH_Salinity")
alpha_div_test=data.frame(alpha_div_test)
alpha_div_test$p.value <- c(6.195e-15,0.05744,1.144e-14, 0.604, 0.9319, 4.343e-08, 0.9301, 0.03307, 0.0003193, 4.007e-05, 0.0002556, 0.6703, 0.4054, 0.06081, 0.5611, 0.4298, 0.3889, 0.5547, 0.1534)

alpha_div_test$padj <- p.adjust(alpha_div_test$p.value, method = "BH")
alpha_div_test
# > alpha_div_test
# alpha_div_test   p.value         padj
# 1       Sample_Type 6.195e-15 1.086800e-13
# 2        W_Salinity 5.744e-02 1.283767e-01
# 3             W_DPH 1.144e-14 1.086800e-13
# 4           W_Tank# 6.040e-01 7.172500e-01
# 5        T_Salinity 9.319e-01 9.319000e-01
# 6             T_DPH 4.343e-08 2.750567e-07
# 7           T_Tank# 9.301e-01 9.319000e-01
# 8   W_6DPH_Salinity 3.307e-02 8.976143e-02
# 9  W_12DPH_Salinity 3.193e-04 1.011117e-03
# 10 W_18DPH_Salinity 4.007e-05 1.903325e-04
# 11 W_24DPH_Salinity 2.556e-04 9.712800e-04
# 12  T_3DPH_Salinity 6.703e-01 7.491588e-01
# 13  T_6DPH_Salinity 4.054e-01 6.281692e-01
# 14  T_9DPH_Salinity 6.081e-02 1.283767e-01
# 15 T_12DPH_Salinity 5.611e-01 7.107267e-01
# 16 T_15DPH_Salinity 4.298e-01 6.281692e-01
# 17 T_18DPH_Salinity 3.889e-01 6.281692e-01
# 18 T_20DPH_Salinity 5.547e-01 7.107267e-01
# 19 T_24DPH_Salinity 1.534e-01 2.914600e-01
##Dunn: Non-parametric test for pair-wise significance##

dunnTest(Shannon ~ Salinity, data = W_Falphadiv_no1D_metadata, method="bh")
dunnTest(Shannon ~ Days_Post_Hatch, data = W_Falphadiv_no1D_metadata, method="bh")
Dunn <- dunnTest(Shannon ~ Tank_Number, data = W_Falphadiv_no1D_metadata, method="bh")
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
#No Significant differences



dunnTest(Shannon ~ Salinity, data = T_Falphadiv_no0D1D_metadata, method="bh")
Dunn <- dunnTest(Shannon ~ Days_Post_Hatch, data = T_Falphadiv_no0D1D_metadata, method="bh")

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 6    D3     cd         cd
# 7    D6      d          d
# 8    D9      d          d
# 1   D12    abc       abc 
# 2   D15      a       a   
# 3   D18     ab       ab  
# 4    D2     ab       ab  
# 5   D24    bcd        bcd

Dunn <- dunnTest(Shannon ~ Tank_Number, data = T_Falphadiv_no0D1D_metadata, method="bh")
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
#No Significant differences


dunnTest(Shannon ~ Salinity, data = filter(W_Falphadiv_metadata, Days_Post_Hatch == "D6"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(W_Falphadiv_metadata, Days_Post_Hatch == "D12"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(W_Falphadiv_metadata, Days_Post_Hatch == "D18"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(W_Falphadiv_metadata, Days_Post_Hatch == "D24"), method="bh")

dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D3"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D6"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D9"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D12"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D15"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D18"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D20"), method="bh")
dunnTest(Shannon ~ Salinity, data = filter(T_Falphadiv_metadata, Days_Post_Hatch == "D24"), method="bh")


###COMBINING STATISTICAL AND ALPHA DIVERSITY TESTING###

##PRECURSOR TO FIGURE 2

#Water by DPH
W_Shannon_DPH_no_stats_boxplot <- ggplot(W_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  #scale_color_manual("Days Post Hatch", values=c("deeppink", "midnightblue", "blue", "red", "orange", "purple"), limits=c("D1", "D6", "D12","D18", "D24"))+
  # ggtitle("Water Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24")) +
  annotate("text", x = c(2.0, 3, 4, 5) , y = c(3.1, 2.9, 4.0, 3.5), label = c("a", "b", "c", "c"), size=5) + 
  annotate("segment", x=0.5, xend=3.0, y=4.25, yend=4.25, size=2, color = "forestgreen")+ #0-12 dph
  annotate("text", x=1.75, y=4.375, label="Green Water", size=5, color="forestgreen")+
  annotate("segment", x=1, xend=2.75, y=4.5, yend=4.5, size=2, color = "red")+ #1-10 dph
  annotate("text", x=1.875, y=4.625, label="Enriched Rotifers", size=5, color ="red") +  
  annotate("segment", x=2.5, xend=3.5, y=4.75, yend=4.75, size=2, color = "blue")+ #9-15 dph
  annotate("text", x=3.0, y=4.875, label="Artemia nauplii", size=5, color ="blue")+
  annotate("segment", x=3.0, xend=4.5, y=5, yend=5, size=2, color = "orange")+ #12-21 dph
  annotate("text", x=3.75, y=5.125, label="Enriched Artemia", size=5, color ="orange")+
  annotate("segment", x=3.25, xend=5, y=5.25, yend=5.25, size=2, color = "purple")+ # 13-24 dph
  annotate("text", x=4.125, y=5.375, label="Microfeeds", size=5, color ="purple")+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) 
W_Shannon_DPH_no_stats_boxplot

tiff('Water Shannon Diversity by Days Post Hatch wo stats.tiff', units="in", width=10, height=6, res=300)
W_Shannon_DPH_no_stats_boxplot
dev.off()


##PRECURSOR TO FIGURE 5

#Tissue by DPH
T_Shannon_DPH_boxplot_no_stats_title <- ggplot(T_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon)) + 
  geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1))+ 
  #scale_color_manual("Days Post Hatch", values=c("lightgray", "deeppink", "lightblue", "midnightblue", "maroon", "blue", "rosybrown", "red", "orange", "darkgray", "purple"), limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
  #ggtitle("Tissue Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  annotate("text", x = c(3, 4, 5, 6, 7, 8, 9, 10) , y = c(5.3, 5.3, 5.35, 4.65, 3.3, 4.5, 3.55, 4.8), label = c("cd", "d", "d", "abc", "a", "ab", "ab", "bcd"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  annotate("segment", x=1, xend=6.0, y=5.5, yend=5.5, size=2, color = "forestgreen")+ #0-12 dph
  annotate("text", x=3.5, y=5.625, label="Green Water", size=5, color="forestgreen")+
  annotate("segment", x=2, xend=5.5, y=5.75, yend=5.75, size=2, color = "red") + #1-10 dph
  annotate("text", x=3.75, y=5.875, label="Enriched Rotifers", size=5, color ="red") +  
  annotate("segment", x=5, xend=7, y=6, yend=6, size=2, color = "blue") + #9-15 dph
  annotate("text", x=6, y=6.125, label="Artemia nauplii", size=5, color ="blue")+
  annotate("segment", x=6, xend=9.5, y=6.25, yend=6.25, size=2, color = "orange")+ #12-21 dph
  annotate("text", x=7.75, y=6.375, label="Enriched Artemia", size=5, color ="orange")+
  annotate("segment", x=6.25, xend=10, y=6.5, yend=6.5, size=2, color = "purple")+ #13-24 dph
  annotate("text", x=8.125, y=6.625, label="Microfeeds", size=5, color ="purple")
T_Shannon_DPH_boxplot_no_stats_title

tiff('Tissue Shannon Diversity by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
T_Shannon_DPH_boxplot_no_stats_title
dev.off()

###BETA DIVERSITY ANALYSIS##

##Determining best transformation##

#Transform data using square root
SR_Fdata <- transform_sample_counts(Fdata, function(x){x^(1/2)})


#Alternative:
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/229

#Install and load DESeq2 so you can use their rlog based transfomation
library("DESeq2")

#Load alternative gm_mean function to handle zeros based upon github issue shown above#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
dds = phyloseq_to_deseq2(Fdata, ~1)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
Fdata
DS_Fdata = Fdata
DS_Fdata

#Switch the asv table with the DESeq2 transformed data
otu_table(DS_Fdata) <- otu_table(getVarianceStabilizedData(dds), taxa_are_rows = TRUE)

#Check to see if your basic phyloseq information was not changed
DS_Fdata

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445

#Make any negative values equal to 0 since the more negative they are the more likely they were to be zero over very small and unlikely to affect your results

#Make a copy of DS_Fdata to manipulate
ZDS_Fdata <- DS_Fdata
ZDS_Fdata

#extract out asv table
DESeq2_otu_table <- as.data.frame(otu_table(ZDS_Fdata))

#Change all negatives to zero
DESeq2_otu_table[DESeq2_otu_table < 0.0] <- 0.0

#Switch out the asv table in phyloseq object
otu_table(ZDS_Fdata) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

#Check to make sure basic phyloseq info did not change
ZDS_Fdata

#Show how amount of positive numbers changed throughout the transformation 
z <- otu_table(Fdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# # 0.98250415 0.01749585 

z <- otu_table(DS_Fdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.98250481 0.01749519    

z <- otu_table(ZDS_Fdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.98250481 0.01749519 

#Standardize with CSS:
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Install and load
library(metagenomeSeq)

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(Fdata) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGS_Fdata = Fdata
MGS_Fdata

#Switch out old ASV table for transformed one
otu_table(MGS_Fdata) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGS_Fdata

#Plot each of the different methods to see the differences#
plot_ordination(SR_Fdata, ordinate(SR_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Sample_Type") 
plot_ordination(ZDS_Fdata, ordinate(ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Sample_Type")
plot_ordination(MGS_Fdata, ordinate(MGS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Sample_Type")

#Check your library size differences by dividing the sequences sum of the sample with the most sequences by the one with the least sequences typically want it to be <10x difference

max(sample_sums(Fdata))/min(sample_sums(Fdata))#61.90364
max(sample_sums(SR_Fdata))/min(sample_sums(SR_Fdata)) #30.43625
max(sample_sums(ZDS_Fdata))/min(sample_sums(ZDS_Fdata)) #16.08665
max(sample_sums(MGS_Fdata))/min(sample_sums(MGS_Fdata)) #36.37826

#Check to see if they are statistically normal, if >0.05 then it is likely normal
shapiro.test(sample_sums(Fdata)) #8.6e-06
shapiro.test(sample_sums(SR_Fdata)) #1.992e-06
shapiro.test(sample_sums(ZDS_Fdata)) #0.0002333
shapiro.test(sample_sums(MGS_Fdata)) #1.883e-07

#Visual representation of the sample sums
hist(sample_sums(Fdata))
hist(sample_sums(SR_Fdata))
hist(sample_sums(ZDS_Fdata))
hist(sample_sums(MGS_Fdata))

#Zeroed Variance Stabilizing Transformed data performed thye best, thus it is used for further analysis

##Export ASV and  table for PRIMER7##

##INFORMATION USED IN SUPPLEMENTARY TABLE 4 AND 5

# Extract abundance matrix from the phyloseq object#
LAB_Salinity_ZDS_Fdata_Bio = as(otu_table(ZDS_Fdata), "matrix")

# Coerce to data.frame#
LAB_Salinity_ZDS_Fdata_Bio = as.data.frame(LAB_Salinity_ZDS_Fdata_Bio)

#Make a csv file#
write.csv(LAB_Salinity_ZDS_Fdata_Bio, file = 'LAB_Salinity_ZDS_Fdata_Bio.csv')

##PCoA Analysis##

W_ZDS_Fdata = subset_samples(ZDS_Fdata, Sample_Type=="Water")
W_ZDS_Fdata #15642 taxa and 125 samples
W_ZDS_Fdata = filter_taxa(W_ZDS_Fdata, function(x) sum(x) >0, TRUE)
W_ZDS_Fdata #8369 taxa and 125 samples

T_ZDS_Fdata = subset_samples(ZDS_Fdata, Sample_Type=="Tissue")
T_ZDS_Fdata #15642 taxa and 69 samples
T_ZDS_Fdata = filter_taxa(T_ZDS_Fdata, function(x) sum(x) >0, TRUE)
T_ZDS_Fdata #10385 taxa and 69 samples

#Subset water samples based upon Days post hatch
W_1D_ZDS_Fdata = subset_samples(W_ZDS_Fdata, Days_Post_Hatch=="D1")
nsamples(W_1D_ZDS_Fdata)#26

W_6D_ZDS_Fdata = subset_samples(W_ZDS_Fdata, Days_Post_Hatch=="D6")
nsamples(W_6D_ZDS_Fdata)#27

W_12D_ZDS_Fdata = subset_samples(W_ZDS_Fdata, Days_Post_Hatch=="D12")
nsamples(W_12D_ZDS_Fdata)#27

W_18D_ZDS_Fdata = subset_samples(W_ZDS_Fdata, Days_Post_Hatch=="D18")
nsamples(W_18D_ZDS_Fdata)#24

W_24D_ZDS_Fdata = subset_samples(W_ZDS_Fdata, Days_Post_Hatch=="D24")
nsamples(W_24D_ZDS_Fdata)#21


#Subset tissue samples based upon Days post hatch
T_0D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D0")
nsamples(T_0D_ZDS_Fdata)#1

T_1D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D1")
nsamples(T_1D_ZDS_Fdata)#2

T_3D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D3")
nsamples(T_3D_ZDS_Fdata)#9

T_6D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D6")
nsamples(T_6D_ZDS_Fdata)#8

T_9D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D9")
nsamples(T_9D_ZDS_Fdata)#9

T_12D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D12")
nsamples(T_12D_ZDS_Fdata)#9

T_15D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D15")
nsamples(T_15D_ZDS_Fdata)#9

T_18D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D18")
nsamples(T_18D_ZDS_Fdata)#8

T_20D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D20")
nsamples(T_20D_ZDS_Fdata)#7

T_24D_ZDS_Fdata = subset_samples(T_ZDS_Fdata, Days_Post_Hatch=="D24")
nsamples(T_24D_ZDS_Fdata)#7

#Create days post hatch phyloseqs
TW_1D_Fdata = subset_samples(ZDS_Fdata, Days_Post_Hatch=="D1") #29
nsamples(TW_1D_Fdata) #28

TW_6D_Fdata = subset_samples(ZDS_Fdata, Days_Post_Hatch=="D6")
nsamples(TW_6D_Fdata) #35

TW_12D_Fdata = subset_samples(ZDS_Fdata, Days_Post_Hatch=="D12")
nsamples(TW_12D_Fdata) #36

TW_18D_Fdata = subset_samples(ZDS_Fdata, Days_Post_Hatch=="D18")
nsamples(TW_18D_Fdata) #32

TW_24D_Fdata = subset_samples(ZDS_Fdata, Days_Post_Hatch=="D24")
nsamples(TW_24D_Fdata) #28

##PRECURSOR TO SUPPLEMENTARY FIGURE 2

#Overall Sample Type PCoA
plot_ordination(ZDS_Fdata, ordinate(ZDS_Fdata, "PCoA", "bray"), color = "Sample_Type")+ 
  scale_color_manual(name="Sample Type", values = c("gray", "lightblue"))+
  ggtitle("Overall PCoA by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))

tiff('Overall PCoA by Sample_Type.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

#Subset Water data without D1 control
nsamples(W_ZDS_Fdata)#125
W_noD1_ZDS_Fdata <- subset_samples(W_ZDS_Fdata, Days_Post_Hatch != "D1")
nsamples(W_noD1_ZDS_Fdata) #99

##PRECURSOR TO FIGURE 3

#Water PCoA without 1 DPH samples
plot_ordination(W_noD1_ZDS_Fdata, ordinate(W_noD1_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D6", "D12", "D18", "D24"), values=c("red", "forestgreen", "purple", "orange")) +
  ggtitle("Water PCoA by DPH & Salinity") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 

tiff('Water PCoA by DPH and Salinity wo D1.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

##TISSUE

#Subset Tissue data without D0 & D1 control
nsamples(T_ZDS_Fdata)#69
T_noD01_ZDS_Fdata <- subset_samples(T_ZDS_Fdata, Days_Post_Hatch != "D0")
T_noD01_ZDS_Fdata <- subset_samples(T_noD01_ZDS_Fdata, Days_Post_Hatch != "D1")
nsamples(T_noD01_ZDS_Fdata) #66

#PRECURSOR TO FIGURE 6

# Tissue PCoA without 0 & 1 DPH samples
plot_ordination(T_noD01_ZDS_Fdata, ordinate(T_noD01_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red",  "black", "forestgreen", "brown", "purple", "dark gray", "orange"))+
  ggtitle("Tissue PCoA by DPH and Salinity") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) +
  annotate("text", x = c(0.15) , y = c(-0.2), label = "Salinity PERMANOVA", size=4)+
  annotate("text", x = c(0.15) , y = c(-0.225), label = "Pseudo-F = 1.52", size=4)+
  annotate("text", x = c(0.15) , y = c(-0.25), label = "P(MC) = 0.0021", size=4) +
  annotate("text", x = c(0) , y = c(0.25), label = "DPH PERMANOVA", size=4)+
  annotate("text", x = c(0) , y = c(0.225), label = "Pseudo-F = 2.09", size=4)+
  annotate("text", x = c(0) , y = c(0.2), label = "P(MC) = 0.0001", size=4)

tiff('Tissue PCoA by DPH and Salinity wo D0-1.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()


###TAXONOMIC BAR PLOTS###

#Use filtered, untransformed data for this

#Both sets of scripts for the Phylum and Genus graphs use the same names for intermediate steps and only has unique name for final data table. You do NOT need to keep the intermediates for each of the steps for each of the graphs you intend to make. You only care about keeping the last data table for graphing, no need to waste space in R.#

# get abundance in relative percentage#
phy <- transform_sample_counts(Fdata, function(x) 100*x/sum(x))
nsamples(phy)#194

noD1_phy <- subset_samples(phy, Days_Post_Hatch != "D0")
noD01_phy <- subset_samples(noD1_phy, Days_Post_Hatch != "D1")
nsamples(noD01_phy)#165


Wphy <- transform_sample_counts(W_Fdata, function(x) 100*x/sum(x))
nsamples(Wphy)#125

W_noD1_Fdataphy <- subset_samples(W_Fdataphy, Days_Post_Hatch != "D1")
nsamples(W_noD1_Fdataphy) #99


Tphy <- transform_sample_counts(T_Fdata, function(x) 100*x/sum(x))
nsamples(Tphy) #69

T_noD1_Fdataphy <- subset_samples(T_Fdataphy, Days_Post_Hatch != "D0")
T_noD01_Fdataphy <- subset_samples(T_noD1_Fdataphy, Days_Post_Hatch != "D1")
nsamples(T_noD01_Fdataphy)#66

##Make Sample Type Graph

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(noD01_phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Sample_Type, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Sample_Type","Genus"), na.rm = FALSE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Sample_Type, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Sample_Type, Abundance, Genus))
#combine with original table
Genus_noD01_Fdata_Sample_Type <- rbind(dat, Abundance)

##PRECURSOR TO FIGURE 1

spatial_plot_Genus_noD01_Fdata_Sample_Type <- ggplot(data=Genus_noD01_Fdata_Sample_Type, aes(x=Sample_Type, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("red","midnightblue","orange","maroon","purple","lightblue","firebrick3","darkseagreen3","seagreen","turquoise","goldenrod","wheat","black","rosybrown"))+
  ggtitle("Genera by Sample Type") +
  xlab("Sample Type") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  annotate("text", x = c(1:2 ) , y = 105, label = c("a", "b"), size=5)

spatial_plot_Genus_noD01_Fdata_Sample_Type

tiff('Genus by Sample Type no D1.tiff', units="in", width=6, height=6, res=300)
spatial_plot_Genus_noD01_Fdata_Sample_Type
dev.off()


## Make Water by Salinity graph

# agglomerate taxa
glom <- tax_glom(W_noD1_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Salinity, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columnsA
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#combine with original table
Genus_W_noD1_Fdata_Salinity <- rbind(dat, Abundance)

##PRECURSOR TO SUPPLEMENTARY FIGURE 4

#Water by Salinity
spatial_plot_Genus_NoD1_W_Fdata_Salinity <- ggplot(data=Genus_W_noD1_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("yellow","red","midnightblue","orange","green","maroon","purple","lightblue","orangered","firebrick3","darkseagreen3","seagreen","white","turquoise","wheat","goldenrod","black","rosybrown"))+
  ggtitle("Water Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "c"), size=5)
spatial_plot_Genus_NoD1_W_Fdata_Salinity

tiff('Water Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_NoD1_W_Fdata_Salinity
dev.off()


##Make Tissue by Salinity Graph

# agglomerate taxa
glom <- tax_glom(T_noD01_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Salinity, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#combine with original table
Genus_T_noD01_Fdata_Salinity <- rbind(dat, Abundance)

##PRECURSOR TO SUPPLEMTNARY FIGURE 7

#Tissue by Salinity
spatial_plot_Genus_T_noD01_Fdata_Salinity <- ggplot(data=Genus_T_noD01_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","purple","lightblue","firebrick3","darkseagreen3","seagreen","wheat","goldenrod","lightgray","black","darkgray","rosybrown"))+
  ggtitle("Tissue Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "b"), size=5)
spatial_plot_Genus_T_noD01_Fdata_Salinity

tiff('Tissue Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_T_noD01_Fdata_Salinity
dev.off()

###GENERA TRENDS

##Sample Type
#Switch the structure of dataframe from long to wide
#i.e. summarize by expanding a column (here Sample_Type) into its subcategories (here Water adn Tissue) and places abundances in each column
#This makes a dataframe that is 28 rows by 3 columns to on that is 14 rows by 3 columns
#Typically this creates a dataframe that is longer because there are more than two subcategories (hence wide)
wide_Genus_noD01_Fdata_Sample_Type <- spread(Genus_noD01_Fdata_Sample_Type, Sample_Type, Abundance) 

#Add a new column that is the Tissue abundance of a particular genera divided by the Water abundance of that same genera
wide_Genus_noD01_Fdata_Sample_Type$TW <- wide_Genus_noD01_Fdata_Sample_Type$Tissue / wide_Genus_noD01_Fdata_Sample_Type$Water 

#Same thing except that its Water divided by Tissue
wide_Genus_noD01_Fdata_Sample_Type$WT <- wide_Genus_noD01_Fdata_Sample_Type$Water / wide_Genus_noD01_Fdata_Sample_Type$Tissue


##Water by Salinity
#Switch the structure of dataframe from long to wide
wide_Genus_W_noD1_Fdata_Salinity <- spread(Genus_W_noD1_Fdata_Salinity, Salinity, Abundance) 

#Add a new column that is the 10 ppt abundance of a particular genera divided by the 20 ppt abundance of that same genera
wide_Genus_W_noD1_Fdata_Salinity$S10S20 <- wide_Genus_W_noD1_Fdata_Salinity$S10 / wide_Genus_W_noD1_Fdata_Salinity$S20 

#Same thing except that its S20 divided by S10
wide_Genus_W_noD1_Fdata_Salinity$S20S10 <- wide_Genus_W_noD1_Fdata_Salinity$S20 / wide_Genus_W_noD1_Fdata_Salinity$S10

#Add a new column that is the 10 ppt abundance of a particular genera divided by the 30 ppt abundance of that same genera
wide_Genus_W_noD1_Fdata_Salinity$S10S30 <- wide_Genus_W_noD1_Fdata_Salinity$S10 / wide_Genus_W_noD1_Fdata_Salinity$S30 

#Same thing except that its S30 divided by S10
wide_Genus_W_noD1_Fdata_Salinity$S30S10 <- wide_Genus_W_noD1_Fdata_Salinity$S30 / wide_Genus_W_noD1_Fdata_Salinity$S10

#Add a new column that is the 20 ppt abundance of a particular genera divided by the 30 ppt abundance of that same genera
wide_Genus_W_noD1_Fdata_Salinity$S20S30 <- wide_Genus_W_noD1_Fdata_Salinity$S20 / wide_Genus_W_noD1_Fdata_Salinity$S30 

#Same thing except that its S30 divided by S20
wide_Genus_W_noD1_Fdata_Salinity$S30S20 <- wide_Genus_W_noD1_Fdata_Salinity$S30 / wide_Genus_W_noD1_Fdata_Salinity$S20


##Tissue by Salinity
#Switch the structure of dataframe from long to wide
wide_Genus_T_noD01_Fdata_Salinity <- spread(Full_Genus_T_noD1_Salinity_dat, Salinity, Abundance) 

#Add a new column that is the 10 ppt abundance of a particular genera divided by the 20 ppt abundance of that same genera
wide_Genus_T_noD01_Fdata_Salinity$S10S20 <- wide_Genus_T_noD01_Fdata_Salinity$S10 / wide_Genus_T_noD01_Fdata_Salinity$S20 

#Same thing except that its S20 divided by S10
wide_Genus_T_noD01_Fdata_Salinity$S20S10 <- wide_Genus_T_noD01_Fdata_Salinity$S20 / wide_Genus_T_noD01_Fdata_Salinity$S10

#Add a new column that is the 10 ppt abundance of a particular genera divided by the 30 ppt abundance of that same genera
wide_Genus_T_noD01_Fdata_Salinity$S10S30 <- wide_Genus_T_noD01_Fdata_Salinity$S10 / wide_Genus_T_noD01_Fdata_Salinity$S30 

#Same thing except that its S30 divided by S10
wide_Genus_T_noD01_Fdata_Salinity$S30S10 <- wide_Genus_T_noD01_Fdata_Salinity$S30 / wide_Genus_T_noD01_Fdata_Salinity$S10

#Add a new column that is the 20 ppt abundance of a particular genera divided by the 30 ppt abundance of that same genera
wide_Genus_T_noD01_Fdata_Salinity$S20S30 <- wide_Genus_T_noD01_Fdata_Salinity$S20 / wide_Genus_T_noD01_Fdata_Salinity$S30 

#Same thing except that its S30 divided by S20
wide_Genus_T_noD01_Fdata_Salinity$S30S20 <- wide_Genus_T_noD01_Fdata_Salinity$S30 / wide_Genus_T_noD01_Fdata_Salinity$S20

### POTENIALLY PATHOGENIC GENERA

#Potentially pathogenic genera for Trachinotus species found in the literature include: Vibrio, Aeromonas, Photobacterium, Pseudomonas, Flexibacter, Mycobacterium, Nocardia, Streptococcus


###Make a new set of tables focusing only on these organisms and plot them as a bar plot

##Water Salinity PPG

#Make a table with all genera in Water (no 1 DPH) samples by Salinity
W_noD1_Fdataphy <- subset_samples(W_Fdataphy, Days_Post_Hatch != "D1")
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(W_noD1_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#combine with original table
Full_Genus_W_noD1_Salinity_dat <- rbind(dat)

#Keep only PPG
PPFG_Full_Genus_W_noD1_Salinity_dat <- Full_Genus_W_noD1_Salinity_dat[Full_Genus_W_noD1_Salinity_dat$Genus %in% PPFG,]

##PRECURSOR TO FIGURE 4B

spatial_plot_PPFG_Full_Genus_W_noD1_Salinity <- ggplot(data=PPFG_Full_Genus_W_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink","midnightblue","chocolate","lightblue","olivedrab","orchid","black"
  ))+
  # scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D20", "D24"))+
  xlab("Salinity") +
  ggtitle("Water Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity

tiff('Water Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity
dev.off()

#Water DPH PPG

#Make a table with all genera in Water samples by DPH

W_Fdataphy <- transform_sample_counts(W_Fdata, function(x) 100*x/sum(x))
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(W_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#combine with original table
Full_Genus_W_Days_Post_Hatch_dat <- rbind(dat)

#Keep only the PPG 
PPFG_Full_Genus_W_Days_Post_Hatch_dat <- Full_Genus_W_Days_Post_Hatch_dat[Full_Genus_W_Days_Post_Hatch_dat$Genus %in% PPFG,]

##PRECURSOR TO FIGURE 4D

spatial_plot_PPFG_Full_Genus_W_DPH <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","chocolate","lightblue","olivedrab","orchid","black"
))+
   scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D24"))+
  xlab("Days Post Hatch") +
  ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
    theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_W_DPH

tiff('Water Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=7, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH
dev.off()

##Tissue Salinity PPG

#Make a table with all genera in Tissue (no 0 & 1 DPH) samples by Salinity

T_noD1_Fdataphy <- subset_samples(T_Fdataphy, Days_Post_Hatch != "D0")
T_noD01_Fdataphy <- subset_samples(T_noD1_Fdataphy, Days_Post_Hatch != "D1")
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(T_noD1_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Salinity, Abundance, Genus))
#combine with original table
Full_Genus_T_noD1_Salinity_dat <- rbind(dat)

#Keep only PPG
PPFG_Full_Genus_T_noD1_Salinity_dat <- Full_Genus_T_noD1_Salinity_dat[Full_Genus_T_noD1_Salinity_dat$Genus %in% PPFG,]

##PRECURSOR TO FIGURE 7B

spatial_plot_PPFG_Full_Genus_T_noD1_Salinity <- ggplot(data=PPFG_Full_Genus_T_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  # scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D20", "D24"))+
  xlab("Salinity") +
  ggtitle("Tissue Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity

tiff('Tissue Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity
dev.off()

##Tissue DPH PPG

#Make a table with all genera in Tissue by DPH

T_Fdataphy <- transform_sample_counts(T_Fdata, function(x) 100*x/sum(x))
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(T_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#combine with original table
Full_Genus_T_Days_Post_Hatch_dat <- rbind(dat)

#Keep only PPG
PPFG_Full_Genus_T_Days_Post_Hatch_dat <- Full_Genus_T_Days_Post_Hatch_dat[Full_Genus_T_Days_Post_Hatch_dat$Genus %in% PPFG,]

##PRECURSOR TO FIGURE 7D

spatial_plot_PPFG_Full_Genus_T_DPH <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
  xlab("Days Post Hatch") +
  ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
    theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_DPH

tiff('Tissue Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH
dev.off()

###Do statistical analysis on PPGs and create boxplots of sums

##Create series of tables to use for testing

##INFORMATION USED IN SUPPLEMENTARY TABLE 3

# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#Make a list of potentially pathogenic genera
PPFG <- c("f__Vibrionaceae;g__Vibrio", 
          "f__Aeromonadaceae;g__Aeromonas", 
          "f__Vibrionaceae;g__Photobacterium",
          "f__Pseudomonadaceae;g__Pseudomonas", 
          "f__Mycobacteriaceae;g__Mycobacterium",
          "f__Nocardiaceae;g__Nocardia", 
          "f__Streptococcaceae;g__Streptococcus",
          "f__Flavobacteriaceae;g__Flavobacterium")
#chan name of potentially pathogenic genera
dat[dat$Genus %in% PPFG,]$Genus <- 'Genera with Potentially Pathogenic Species'
#keep only PPG
dat <- filter(dat, Genus=="Genera with Potentially Pathogenic Species")
#Summarize based upon target parameter
PPG_dat <- dat %>% group_by(Sample, Sample_Type, Days_Post_Hatch, Tank_Number, Salinity, Sample_Type_Days_Post_Hatch, Sample_Type_Tank_Number,Sample_Type_Salinity,Days_Post_Hatch_Tank_Number,Days_Post_Hatch_Salinity,Tank_Number_Salinity) %>%  summarise(sum_PPG = sum(Abundance))

noD0_PPG_dat <- filter(PPG_dat, Days_Post_Hatch!="D0")
noD01_PPG_dat <- filter(noD0_PPG_dat, Days_Post_Hatch!="D1")


W_PPG_dat <- filter(PPG_dat, Sample_Type=="Water")
W_noD1_PPG_dat <- filter(W_PPG_dat, Days_Post_Hatch!="D1")

T_PPG_dat <- filter(PPG_dat, Sample_Type=="Tissue")
T_noD1_PPG_dat <- filter(T_PPG_dat, Days_Post_Hatch!="D0")
T_noD01_PPG_dat <- filter(T_noD1_PPG_dat, Days_Post_Hatch!="D1")

#Water vs Tissue PPG
kruskal.test(sum_PPG ~ Sample_Type, data = noD01_PPG_dat)
# Kruskal-Wallis rank sum test
# 
# data:  sum_PPG by Sample_Type
# Kruskal-Wallis chi-squared = 53.938, df = 1, p-value = 2.07e-13

#Water PPG by Salinity
kruskal.test(sum_PPG ~ Salinity, data = W_noD1_PPG_dat)
# Kruskal-Wallis chi-squared = 45.794, df = 2, p-value = 1.137e-10

dunnTest(sum_PPG ~ Salinity, data = W_noD1_PPG_dat, method="bh") 
# Comparison        Z     P.unadj      P.adj
# 1  S10 - S20 5.077628 3.821763e-07 5.732645e-07
# 2  S10 - S30 6.494876 8.310199e-11 2.493060e-10
# 3  S20 - S30 1.347227 1.779071e-01 1.779071e-01

##PRECURSOR TO FIGURE 4A

ggplot(W_noD1_PPG_dat, aes(x=Salinity, y=sum_PPG, color=Salinity)) +   geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("Salinity", values=c("deeppink", "midnightblue", "blue"))+
  ggtitle("Water PPG% by Salinity") +
  xlab("Salinity") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(73, 38, 8), label = c("a", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

tiff('Water Potentially Pathogenic Genera by Salinity BP.tiff', units="in", width=6, height=6, res=300)
#plot
dev.off()

#Tissue PPG by Salinity
kruskal.test(sum_PPG ~ Salinity, data = T_noD01_PPG_dat)
# Kruskal-Wallis chi-squared = 0.70324, df = 2, p-value = 0.7035

dunnTest(sum_PPG ~ Salinity, data = T_noD01_PPG_dat, method="bh") 
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 -0.8279304 0.4077099 1.0000000
# 2  S10 - S30 -0.3133920 0.7539829 0.7539829
# 3  S20 - S30  0.5416908 0.5880315 0.8820473

##PRECURSOR TO FIGURE 7A

ggplot(T_noD01_PPG_dat, aes(x=Salinity, y=sum_PPG, color=Salinity)) +   geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("Salinity", values=c("deeppink", "midnightblue", "blue"))+
  ggtitle("Tissue PPG% by Salinity") +
  xlab("Salinity") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(97, 97, 99), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

tiff('Tissue Potentially Pathogenic Genera by Salinity BP.tiff', units="in", width=6, height=6, res=300)
#plot
dev.off()


#Water PPG by DPH
kruskal.test(sum_PPG ~ Days_Post_Hatch, data = W_noD1_PPG_dat)
# Kruskal-Wallis chi-squared = 33.045, df = 3, p-value = 3.151e-07

Dunn <- dunnTest(sum_PPG ~ Days_Post_Hatch, data = W_noD1_PPG_dat, method="bh") 
# Comparison          Z      P.unadj        P.adj
# 1  D12 - D18 -1.1203567 2.625618e-01 3.150742e-01
# 2  D12 - D24  0.9895577 3.223904e-01 3.223904e-01
# 3  D18 - D24  2.0154330 4.385930e-02 6.578895e-02
# 4   D12 - D6 -4.5388064 5.657356e-06 1.697207e-05
# 5   D18 - D6 -3.2829323 1.027333e-03 2.054667e-03
# 6   D24 - D6 -5.2352223 1.647864e-07 9.887187e-07

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a         a 
# 2   D18      a         a 
# 3   D24      a         a 
# 4    D6      b          b

##PRECURSOR TO FIGURE 4C

ggplot(W_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG, color=Days_Post_Hatch)) +   geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("Days Post Hatch", values=c("black", "red", "forestgreen", "purple", "orange"), limits=c("D1", "D6", "D12","D18", "D24"))+
  ggtitle("Water PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24")) +
  annotate("text", x = c(2, 3, 4, 5) , y = c(73, 6, 9, 36), label = c("a", "b", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

tiff('Water Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
#plot
dev.off()

#Tissue PPG by DPH
kruskal.test(sum_PPG ~ Days_Post_Hatch, data = T_noD01_PPG_dat)
# Kruskal-Wallis chi-squared = 46.747, df = 7, p-value = 6.252e-08

Dunn <-dunnTest(sum_PPG ~ Days_Post_Hatch, data = T_noD01_PPG_dat, method="bh") 
# Comparison          Z      P.unadj        P.adj
# 1   D12 - D15 -3.1801292 1.472094e-03 5.152329e-03
# 2   D12 - D18 -1.2194793 2.226623e-01 2.968831e-01
# 3   D15 - D18  1.8656992 6.208347e-02 1.158891e-01
# 4   D12 - D20 -1.6834427 9.228942e-02 1.520061e-01
# 5   D15 - D20  1.2912957 1.966011e-01 2.752416e-01
# 6   D18 - D20 -0.4942814 6.211074e-01 6.688849e-01
# 7   D12 - D24  1.8901813 5.873371e-02 1.174674e-01
# 8   D15 - D24  4.8649198 1.145032e-06 1.603045e-05
# 9   D18 - D24  2.9854599 2.831523e-03 8.809181e-03
# 10  D20 - D24  3.3692451 7.537438e-04 3.517471e-03
# 11   D12 - D3  1.0068363 3.140135e-01 3.822773e-01
# 12   D15 - D3  4.1869655 2.827088e-05 2.638615e-04
# 13   D18 - D3  2.1962540 2.807377e-02 6.550546e-02
# 14   D20 - D3  2.6252518 8.658490e-03 2.424377e-02
# 15   D24 - D3 -0.9483722 3.429400e-01 4.000967e-01
# 16   D12 - D6  0.5494357 5.827064e-01 6.526312e-01
# 17   D15 - D6  3.6346142 2.783969e-04 1.559022e-03
# 18   D18 - D6  1.7190765 8.560045e-02 1.498008e-01
# 19   D20 - D6  2.1550671 3.115658e-02 6.710649e-02
# 20   D24 - D6 -1.3246743 1.852792e-01 2.730431e-01
# 21    D3 - D6 -0.4273389 6.691325e-01 6.939152e-01
# 22   D12 - D9  2.2101284 2.709625e-02 6.897228e-02
# 23   D15 - D9  5.3902576 7.035675e-08 1.969989e-06
# 24   D18 - D9  3.3636188 7.692771e-04 3.077108e-03
# 25   D20 - D9  3.7508286 1.762512e-04 1.233758e-03
# 26   D24 - D9  0.1772045 8.593478e-01 8.593478e-01
# 27    D3 - D9  1.2032921 2.288633e-01 2.912806e-01
# 28    D6 - D9  1.5947038 1.107785e-01 1.723221e-01

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

#   Group Letter MonoLetter
# 1   D12    abc       abc 
# 2   D15      d          d
# 3   D18    abd       ab d
# 4    D2     ad       a  d
# 5   D24      c         c 
# 6    D3     bc        bc 
# 7    D6    abc       abc 
# 8    D9      c         c 

##PRECURSOR TO FIGURE 7C

ggplot(T_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG, color=Days_Post_Hatch)) +   geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("wheat", "black", "rosybrown", "red",  "blue", "forestgreen", "brown", "purple", "dark gray", "orange"))+
  ggtitle("Tissue PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24")) +
  annotate("text", x = c(3, 4, 5, 6, 7, 8, 9, 10) , y = c(38, 47, 18, 79, 99, 87, 92, 32), label = c("bc", "abc", "c", "abc", "d", "abd", "ad", "c"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

tiff('Tissue Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
#plot
dev.off()                 

##Conduct multiple testing adjustments
ppg_test = c("Sample_Type", "W_Salinity", "W_DPH", "T_Salinity", "T_DPH")
ppg_test=data.frame(ppg_test)
ppg_test$p.value <- c(2.07e-13, 1.137e-10, 3.151e-07,0.7035, 6.252e-08 )

ppg_test$padj <- p.adjust(ppg_test$p.value, method = "BH")
ppg_test

###Create Bar Plots with every Salinity/DPH shown

##Water by DPH by Salinity

W_Fdataphy <- transform_sample_counts(W_Fdata, function(x) 100*x/sum(x))
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(W_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch_Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch_Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch_Salinity, Abundance, Genus))
#combine with original table
Full_Genus_W_Days_Post_Hatch_Salinity_dat <- rbind(dat)

#Keep only PPG
PPFG_Full_Genus_W_Days_Post_Hatch_Salinity_dat <- Full_Genus_W_Days_Post_Hatch_Salinity_dat[Full_Genus_W_Days_Post_Hatch_Salinity_dat$Genus %in% PPFG,]

##PRECURSOR TO SUPPLEMENTARY FIGURE 5

spatial_plot_PPFG_Full_Genus_W_DPH_Salinity <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink","midnightblue", "chocolate","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c(
    "D1_S30", 
    "D6_S10", "D6_S20", "D6_S30",
    "D12_S10", "D12_S20", "D12_S30",
    "D18_S10", "D18_S20", "D18_S30",
    "D24_S10", "D24_S20", "D24_S30"))+
  xlab("Days Post Hatch and Salinity") +
  ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
  spatial_plot_PPFG_Full_Genus_W_DPH_Salinity

tiff('Water Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity
dev.off()

##Tissue by DPH by Salinity

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(T_Fdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch_Salinity, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch_Salinity","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch_Salinity, Abundance, Genus))
#combine with original table
Full_Genus_T_Days_Post_Hatch_Salinity_dat <- rbind(dat)

#Keep only PPG
PPFG_Full_Genus_T_Days_Post_Hatch_Salinity_dat <- Full_Genus_T_Days_Post_Hatch_Salinity_dat[Full_Genus_T_Days_Post_Hatch_Salinity_dat$Genus %in% PPFG,]

##PRECURSOR TO SUPPLEMENTARY FIGURE 8

spatial_plot_PPFG_Full_Genus_T_DPH_Salinity <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c(
    "D0_S30", 
    "D1_S30", 
    "D3_S10", "D3_S20", "D3_S30",
    "D6_S10", "D6_S20", "D6_S30",
    "D9_S10", "D9_S20", "D9_S30",
    "D12_S10", "D12_S20", "D12_S30",
    "D15_S10", "D15_S20", "D15_S30",
    "D18_S10", "D18_S20", "D18_S30",
    "D20_S10", "D20_S20", "D20_S30",
    "D24_S10", "D24_S20", "D24_S30"))+
  xlab("Days Post Hatch and Salinity") +
  ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
    theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity

tiff('Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity
dev.off()