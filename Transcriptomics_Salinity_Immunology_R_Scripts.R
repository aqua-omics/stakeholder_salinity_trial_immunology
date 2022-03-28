###Prepatory scripts used in previous publication###

#Importing information into R was done using the following tximport vignette
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

#Analysis using DESeq2 was done using the following vignette
#http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

#Setting the working directory
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/R")


#Install packages#
BiocManager::install("tximportData")
BiocManager::install("BUSpaRse")
BiocManager::install("tximport")
BiocManager::install("edgeR")
BiocManager::install("apeglm")
install.packages("tidyverse")
install.packages("readxl")
BiocManager::install("topGO")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("org.Mn.eg.db")
BiocManager::install("org.Hs.eg.db")
install.packages("dendextend")

BiocManager::install("RnaSeqGeneEdgeRQL")


#Load libraries#
library(RnaSeqGeneEdgeRQL)
library(topGO)
library(phyloseq)
library(DESeq2)
library(tximportData)
library(tximport)
library(edgeR)
library(apeglm)
library(tidyverse)
library(readxl)
library(phyloseq)
library(tximeta)
library(gtools)
library(pheatmap)
library(dendextend)
library(FSA)
library(biomaRt)
library(GO.db)
library(clusterProfiler)
library(rcompanion)
library(ggpubr)
library(pathview)
library(ggvenn)

###PREPARE FOR IMPORTING OF CLC GENOMICS DATA AND PHYLOSEQ###

##Prepare your working environment

#In your working directory you should have:
#CLC Genomics Excel file - each sheet should be a different sample
#sample.txt - tab deliminated file that has your metadata for your samples including a column named "sample" that lists the names that you want your samples to have in the ultimate abundance table
#tx2gene.txt - tab deliminated file with the first column being transcript names (TXNAME) and the second column being gene names (GENEID), these are based upon the CLC Genomics columns of "Name" and "Gene name" or if you want to get technical "Transcript ID" and "GeneID" 

#Set your working directory
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis")
getwd()

#Upload your metadata file
samples <- read.csv("salinity_transcriptomics_metadata.csv", row.names = 1)

#Create a new folder  in that work directory
dir.create("CLC_Genomics_Files")

#Set your working directory to that new folder for this next section
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/CLC_Genomics_Files")
getwd()

##Prep abundance files##

#Load function that will write every sheet in an Excel file (designated path here) to a series of tsv files#
read_then_tsv <- function(sheet, path) {
  pathbase <- path %>%
    basename() %>%
    tools::file_path_sans_ext()
  path %>%
    read_excel(sheet = sheet) %>%
    write_tsv(paste0(pathbase, "-", sheet, ".tsv"))
}

#Set path to equal the Excel with all your sample sheets
path <- "C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/Salinity_TE_all.xlsx"

#Make a tsv file for every sheet in the Excel file
path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_then_tsv, path = path)

##Rename all the files using your metadata table

#list.files() - just lists files in your present working directory (pwd) similar to ls in Linux
myFiles <- list.files()

#See what this list looks like
myFiles

#mixed sort allows numbers to be sorted in correct order instead of treating them like a factor, rev reverses the order
myFiles <- rev(mixedsort(myFiles))

#Check that it worked
myFiles

#use my files to change the name of all files in the pwd
file.rename(myFiles, paste0(samples$sample,".tsv"))

#Check to see that it worked
list.files()

#Set your working directory back to the other folder
getwd()
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis")
getwd()

##Make a vector pointing to quantification files

#Make a value that points to folder with your new set of tsv files for each sample 
dir <- "C:/Users/bradshawd/Documents/Bioinformatic_Analysis/LAB_salinity_trial/Transcriptomics_Analysis/CLC_Genomics_Files"
list.files(dir)

#Create a dataframe of the file names to use later
names <- as.data.frame(list.files(dir))
names$names <- names$`list.files(dir)`
names

#Create a list of files and their associated paths
files <- file.path(dir, names$names)
files
all(file.exists(files))

#Use the metadata table to add labels to each of these files
names(files) <- paste0(samples$sample)
names(files)

#Load in your transcripts to genes file
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)

#Use tximport to import all the files into a list of data frames that have transcript information

#You are using the type = "none" because tximport does not have a built in way to deal with CLC Genomics outputs. That means you must define the following information by telling the tool the names of the columns in the tables that have this information: column with transcript name that must match whats in your tx2gene table (txIdCol), column with TPM or FPKM abundances (abundanceCol), column with estimate counts (countsCol), and column with length information (lengthCol)

#You must also define the importer used to get your tables into R, use read_tsv here since it is opposite of the write_tsv used earlier

#countsFromAbundance - allows you to get counts scaled to library size (scaledTPM), or scaled using average transcript length over samples and then by library size(lengthScaledTPM), or scaled using the median transcript length among isoforms fo the gene and then by library size (dtuScaledTPM) that is generally used with txOut=TRUE and if you want to do differential transcript usage analysis

#txOut = TRUE - means that you are just generating transcript level outputs

#Make a set of tables with transcript information
txi.CLC.tx <- tximport(files, type="none", txIdCol = "Name", tx2gene=tx2gene, abundanceCol = "TPM", countsCol = "Total transcript reads", lengthCol = "Transcript length", importer = read_tsv, txOut=TRUE)

names(txi.CLC.tx)
#"abundance" "counts" "length" "countsFromAbundance"

#Extract tables of each of the elements of this tximport object and save them
head(txi.CLC.tx$counts)
tx_counts_table <- as.data.frame(txi.CLC.tx$counts)
write.csv(tx_counts_table, "salinity_tx_counts.csv")

head(txi.CLC.tx$abundance)
tx_abundance_table <- as.data.frame(txi.CLC.tx$abundance)
write.csv(tx_abundance_table, "salinity_tx_abundance.csv")

head(txi.CLC.tx$length)
tx_length_table <- as.data.frame(txi.CLC.tx$length)
write.csv(tx_length_table, "salinity_tx_length.csv")

#Make a gene-level object
txi.CLC.gn <- summarizeToGene(txi.CLC.tx, tx2gene)

names(txi.CLC.gn)
#"abundance" "counts" "length" "countsFromAbundance"

#Make a set of tables with gene information
head(txi.CLC.gn$counts)
gn_counts_table <- as.data.frame(txi.CLC.gn$counts)
write.csv(gn_counts_table, "salinity_gn_counts.csv")

head(txi.CLC.gn$abundance)
gn_abundance_table <- as.data.frame(txi.CLC.gn$abundance)
write.csv(gn_abundance_table, "salinity_gn_abundance.csv")

#Take just one column of the length table to be used later and relate it length
head(txi.CLC.gn$length)
gn_length_table <- as.data.frame(txi.CLC.gn$length)
gn_length_table <- gn_length_table[,1, drop=FALSE]
colnames(gn_length_table)=c("length")
head(gn_length_table)
write.csv(gn_length_table, "salinity_gn_length.csv")

##Get sequence information for transcripts and genes
tx_Totals <- as.data.frame(colSums(tx_counts_table))
colnames(tx_Totals) <- "Sum"
max(tx_Totals$Sum)/min(tx_Totals$Sum) #3.57964
sum(tx_Totals$Sum) #553,366,960
dim(tx_counts_table) #34893    57

gene_Totals <- as.data.frame(colSums(gn_counts_table))
colnames(gene_Totals) <- "Sum"
max(gene_Totals$Sum)/min(gene_Totals$Sum) #3.579674
sum(gene_Totals$Sum) #553366960
dim(gn_counts_table) #23800    57

write.csv(tx_Totals, "transcrip_totals.csv")

##Prepare transcript-level data for phyloseq
#Convert counts and metadata into phyloseq friendly versions#
META <- import_qiime_sample_data("salinity_transcriptomics_metadata.txt")

#Make a table of the transcript counts into phyloseq format
TXCTS = otu_table(tx_counts_table, taxa_are_rows = TRUE)

#Create a phyloseq object
TX_dataA <- merge_phyloseq(TXCTS, META)

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
nofactor_tx_dds <- DESeqDataSetFromTximport(txi.CLC.tx, samples, ~1)

#Conduct DESEQ2 test#
nofactor_tx_dds = DESeq(nofactor_dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
TX_dataA
VST_TX_dataA = TX_dataA
VST_TX_dataA

#Switch the asv table with the DESeq2 transformed data
otu_table(VST_TX_dataA) <- otu_table(getVarianceStabilizedData(nofactor_tx_dds), taxa_are_rows = TRUE)

#Extract out VST transcript table
VST_TX_Table <- as.data.frame(otu_table(VST_TX_dataA))

#Print the VST transcript table
write.csv(otu_table(VST_TX_dataA), "VST_TX_Table.csv")

#Make a PCoA
TX_PCoA_Ord <- ordinate(VST_TX_dataA,"PCoA")
TX_PCoA_Plot = plot_ordination(VST_TX_dataA, TX_PCoA_Ord, color="Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red", "forestgreen", "brown", "purple", "dark gray", "orange"), labels=c("3", "6", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PCoA by DPH and Salinity") +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 
TX_PCoA_Plot


##Prepare gene-level data for phyloseq
GNCTS = otu_table(gn_counts_table, taxa_are_rows = TRUE)
GN_dataA <- merge_phyloseq(GNCTS, META)

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
nofactor_gn_dds <- DESeqDataSetFromTximport(txi.CLC.gn, samples, ~1)

#Conduct DESEQ2 test#
nofactor_gn_dds = DESeq(nofactor_gn_dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
GN_dataA
VST_GN_dataA = GN_dataA
VST_GN_dataA

#Switch the asv table with the DESeq2 transformed data
otu_table(VST_GN_dataA) <- otu_table(getVarianceStabilizedData(nofactor_gn_dds), taxa_are_rows = TRUE)

#Extract out VST gene table
VST_GN_Table <- as.data.frame(otu_table(VST_GN_dataA))

#Print the VST gene table for PRIMER7/PERMANOVA+ analysis (TABLE T4)
write.csv(VST_GN_Table, "VST_GN_Table.csv")

#Extract out the VST abundance data from the phyloseq object
data <- as.data.frame(otu_table(VST_GN_dataA))

#Check how different the sample sums are between samples
max(rowSums(t(data)))/min(rowSums(t(data))) #1.046653

#Visualize the samples sums
hist(rowSums(t(data)))


###Gene Set Enrichment Analysis
#Because there was such a variation in the number of DEGs at the dph level and there were little DEGs at the pooled level, this lead to similar varied results in the GO and KEGG analysis (See Other Transcriptome Exploration Scripts Not In Manuscript.R). Thus it was decided to test all genes at once instead of just the DEGs using GSEA and jus to focus on the pooled trends.

##S10 vs S30
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric
full_S10_S30 <- read.csv("Total All DPH 10 vs. 30 ppt.csv", stringsAsFactors=FALSE)
full_S10_S30 <- full_S10_S30[which(full_S10_S30$Fold.change!="#N/A"),]
full_S10_S30$Fold.change <- as.numeric(full_S10_S30$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S10_S30[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S10_S30[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_10vs30 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_10vs30)

df_GSEA_10vs30 <- as.data.frame(GSEA_10vs30)

##S10 vs S20
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric
full_S10_S20 <- read.csv("Total All DPH 10 vs. 20 ppt.csv", stringsAsFactors=FALSE)
full_S10_S20 <- full_S10_S20[which(full_S10_S20$Fold.change!="#N/A"),]
full_S10_S20$Fold.change <- as.numeric(full_S10_S20$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S10_S20[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S10_S20[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_10vs20 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_10vs20)

df_GSEA_10vs20 <- as.data.frame(GSEA_10vs20)

##S20 vs S30
#Upload the full list of genes post CLC-Genomics DEGs analysis, get rid of lines with #N/A as the Fold Change, and make the fold change numeric

full_S20_S30 <- read.csv("Total All DPH 20 vs. 30 ppt.csv", stringsAsFactors=FALSE)
full_S20_S30 <- full_S20_S30[which(full_S20_S30$Fold.change!="#N/A"),]
full_S20_S30$Fold.change <- as.numeric(full_S20_S30$Fold.change)

#Prep for GSEA
# feature 1: create a numeric vector using the column with Fold Change values
geneList = full_S20_S30[,6]
# feature 2: name the values in this vector with the Gene Id column
names(geneList) = as.character(full_S20_S30[,10])
# feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#Run the GSEA with normal parameters make a dataframe from the results
GSEA_20vs30 <- gseKEGG(geneList     = geneList,
                       organism     = 'sdu',
                       minGSSize    = 120,
                       pvalueCutoff = 0.05)
head(GSEA_20vs30)

df_GSEA_20vs30 <- as.data.frame(GSEA_20vs30)

##Plot the results

#S10 vs S30

#Create a new column denoting activated/suppressed status
df_GSEA_10vs30$Direction <- "Activated"

df_GSEA_10vs30$Direction[which(df_GSEA_10vs30$NES < 0)] <- "Suppressed"

##FIGURE 1

#Plot and publish the results
S10vsS30_GSEA_point_plot <- ggplot(df_GSEA_10vs30)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_10vs30[which(df_GSEA_10vs30$p.adjust < 0.05),], n = 37))+
  theme_bw()+
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S10vsS30_GSEA_point_plot

tiff('Fig 1 Gene Sets Enriched between 10 ppt and 30 ppt.tiff', units="in", width=6, height=6, res=300)
S10vsS30_GSEA_point_plot
dev.off()

#S20 vs S30

#Create a new column denoting activated/suppressed status
df_GSEA_20vs30$Direction <- "Activated"

df_GSEA_20vs30$Direction[which(df_GSEA_20vs30$NES < 0)] <- "Suppressed"

##SUPP FIGURE 2

#Plot and publish the results
S20vsS30_GSEA_point_plot <- ggplot(df_GSEA_20vs30)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_20vs30[which(df_GSEA_20vs30$p.adjust < 0.05),], n = 37))+
  theme_bw() +
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S20vsS30_GSEA_point_plot

tiff('Supp Fig 2 Gene Sets Enriched between 20 ppt and 30 ppt.tiff', units="in", width=6, height=6, res=300)
S20vsS30_GSEA_point_plot
dev.off()

#S10 vs S20

#Create a new column denoting activated/suppressed status
df_GSEA_10vs20$Direction <- "Activated"

df_GSEA_10vs20$Direction[which(df_GSEA_10vs20$NES < 0)] <- "Suppressed"

##SUPP FIGURE 3

#Plot and publish the results
S10vsS20_GSEA_point_plot <- ggplot(df_GSEA_10vs20)+
  geom_point(mapping = aes(y= reorder(Description, NES), x=(rank), color = p.adjust, size=abs(NES), shape=Direction, fill = p.adjust),
             data = head(df_GSEA_10vs20[which(df_GSEA_10vs20$p.adjust < 0.05),], n = 37))+
  theme_bw()+
  scale_shape_manual(values=c(24,25))+
  scale_fill_continuous(low="red", high="blue") +
  scale_color_continuous(low="red", high="blue", guide="none") +
  labs(y = "Enriched Gene Sets", x = "rank", color = "p.adjust", size = "NES", shape = "Expression") +
  theme(axis.text=element_text(size=8)) +
  #ggtitle('24 DPH 10ppt vs 30 ppt') +
  theme(legend.key.size = unit(1.2, "line"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
S10vsS20_GSEA_point_plot

tiff('Supp Fig 3 Gene Sets Enriched between 10 ppt and 20 ppt.tiff', units="in", width=6, height=6, res=300)
S10vsS20_GSEA_point_plot
dev.off()


###GENE TREND STATISTICS AND FIGURES

#Copy data to manipulate it, change rownames to from gene ids to gene names, Transpose table, and add in metadata
VST_cts_boxplot <- data
rownames(VST_cts_boxplot) <- name2id$Gene.name[match(rownames(VST_cts_boxplot), name2id$GeneID)]
VST_cts_boxplot <- as.data.frame(t(VST_cts_boxplot))
VST_cts_boxplot <- cbind(VST_cts_boxplot, samples)




###Scripts unique to this publication###

##Immune Regulation

#cat
#Salinity
kruskal.test(cat ~ Salinity, data=VST_cts_boxplot)
dunnTest(cat ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(cat ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(cat ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a        a  
# 2   D15    abc        abc
# 3   D18     bc         bc
# 4    D2      b         b 
# 5   D24    abc        abc
# 6    D3     ac        a c
# 7    D6    abc        abc

#"ac", "abc", "a", "abc", "bc", "b", "abc"

##FIGURE T3A

cat <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=cat, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("cat - catalase") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.25, 9.85, 7.55, 9.1, 9.1, 9.3, 8.85), label = c("ac", "abc", "a", "abc", "bc", "b", "abc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
cat

tiff('cat VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
cat
dev.off()

#lyz
#Salinity
kruskal.test(lyz ~ Salinity, data=VST_cts_boxplot)
dunnTest(lyz ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(lyz ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(lyz ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a        a  
# 2   D15     ab        ab 
# 3   D18     bc         bc
# 4    D2    abc        abc
# 5   D24    abc        abc
# 6    D3      c          c
# 7    D6    abc        abc

#"a", "abc", "c", "bc", "ab", "abc", "abc"

##FIGURE T3B

lyz <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=lyz, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("lyz - lysozyme activity") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(6.4, 4.85, 4.2, 4.8, 5.15, 7.05, 5.8), label = c("a", "abc", "c", "bc", "ab", "abc", "abc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
lyz

tiff('lyz VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
lyz
dev.off()


#nod1
#Salinity
kruskal.test(nod1 ~ Salinity, data=VST_cts_boxplot)
dunnTest(nod1 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(nod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(nod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab        ab 
# 2   D15    abc        abc
# 3   D18     ac        a c
# 4    D2      c          c
# 5   D24     ac        a c
# 6    D3      b         b 
# 7    D6     ab        ab 

#"a", "ab", "ab", "abc", "bc", "c", "bc"

##FIGURE T3C

nod1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=nod1, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("nod1 - regulation of apoptotic process") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(3.5, 4.3, 3.9, 4.8, 5.05, 5.1, 4.7), label = c("a", "ab", "ab", "abc", "bc", "c", "bc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
nod1

tiff('nod1 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
nod1
dev.off()

#sod1
#Salinity
kruskal.test(sod1 ~ Salinity, data=VST_cts_boxplot)
dunnTest(sod1 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(sod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(sod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab         ab
# 2   D15      a         a 
# 3   D18      b          b
# 4    D2     ab         ab
# 5   D24     ab         ab
# 6    D3      b          b
# 7    D6     ab         ab

#"a", "ab", "ab", "b", "a", "ab", "ab"

##FIGURE T3D

sod1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=sod1, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("sod1 - superoxide dismutase activity") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(10.45, 11, 11.15, 11.35, 10.55, 10.9, 10.95), label = c("a", "ab", "ab", "b", "a", "ab", "ab"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
sod1

tiff('sod1 VST Gene Counts by Days Post Hatch.tiff', units="in", width=6, height=6, res=300)
sod1
dev.off()

#tlr3
#Salinity
kruskal.test(tlr3 ~ Salinity, data=VST_cts_boxplot)
dunnTest(tlr3 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(tlr3 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(tlr3 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12      a         a 
# 2   D15      a         a 
# 3   D18      b          b
# 4    D2      b          b
# 5   D24     ab         ab
# 6    D3     ab         ab
# 7    D6     ab         ab

#"ab", "ab", "a", "a", "b", "b", "ab"

##FIGURE T3E

tlr3 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=tlr3, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("tlr3 - positive regulation of inflammatory response") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(5.2, 5.2, 4.9, 5.1, 5.7, 5.2, 6), label = c("ab", "ab", "a", "a", "b", "b", "ab"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
tlr3

tiff('tlr3 VST Gene Counts by Days Post Hatch.tiff', units="in", width=6, height=6, res=300)
tlr3
dev.off()

#tnfaip3
#Salinity
kruskal.test(tnfaip3 ~ Salinity, data=VST_cts_boxplot)
dunnTest(tnfaip3 ~ Salinity, data=VST_cts_boxplot)

#DPH
kruskal.test(tnfaip3 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn <- dunnTest(tnfaip3 ~ Days_Post_Hatch, data=VST_cts_boxplot)

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1   D12     ab         ab
# 2   D15      a         a 
# 3   D18      a         a 
# 4    D2      a         a 
# 5   D24     ab         ab
# 6    D3      b          b
# 7    D6     ab         ab

#"a", "ab", "ab", "b", "b", "b", "ab"

##FIGURE T3F

tnfaip3 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=tnfaip3, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("green", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw()+
  ggtitle("tnfaip3 - inhibit apoptosis") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.35, 8.05, 7.53, 6.95, 6.8, 6.7, 7.28), label = c("a", "ab", "ab", "b", "b", "b", "ab"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
tnfaip3

tiff('tnfaip3 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
tnfaip3
dev.off()

##FIGURE T3

#Make a multi-part figure, removing the axis labels from each separate figure
target_immunity_gene_figure <- ggarrange(cat + rremove("ylab") + rremove("xlab"), lyz + rremove("ylab") + rremove("xlab"), nod1 + rremove("ylab") + rremove("xlab"), sod1 + rremove("ylab") + rremove("xlab"), tlr3 + rremove("ylab") + rremove("xlab"), tnfaip3 + rremove("ylab") + rremove("xlab"),
                                         labels = c("A", "B", "C", "D", "E", "F"),
                                         ncol = 2, nrow = 3,
                                         common.legend = TRUE)
target_immunity_gene_figure

#Add annotations to the above figure
target_immunity_gene_figure <- annotate_figure(target_immunity_gene_figure, 
                                               bottom = text_grob("Days Post Hatch", size = 16),
                                               left = text_grob("VST Gene Count", rot = 90, size =16))
target_immunity_gene_figure

tiff('Figure 3 Target Immunity Genes Figure.tiff', units="in", width=10, height=10, res=300)
target_immunity_gene_figure
dev.off()
