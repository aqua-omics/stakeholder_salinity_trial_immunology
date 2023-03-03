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

##FIGURE 5 of Bradshaw et al., 2023 in Aquaculture

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

tiff('Fig 5 Gene Sets Enriched between 10 ppt and 30 ppt.tiff', units="in", width=6, height=6, res=300)
S10vsS30_GSEA_point_plot
dev.off()

#S20 vs S30

#Create a new column denoting activated/suppressed status
df_GSEA_20vs30$Direction <- "Activated"

df_GSEA_20vs30$Direction[which(df_GSEA_20vs30$NES < 0)] <- "Suppressed"

##SUPP FIGURE 2 of Bradshaw et al., 2023 in Aquaculture

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

##SUPP FIGURE 3 of Bradshaw et al., 2023 in Aquaculture

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

#Create dataframe to store Kruskal-Wallis Tests results
#INFORMATION USED IN SUPPLEMENTARY TABLE 6

gene_stats <- data.frame()

##INFECTION RELATED GENES##

#nod1
#Salinity
kruskal.test(nod1 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 8.511, df = 2, p-value = 0.01419
gene_stats = rbind(gene_stats, c("nod1 - Salinity", 8.511, 0.01419))

Dunn <- dunnTest(nod1 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z     P.unadj      P.adj
# 1  S10 - S20 -0.8334241 0.404605603 0.40460560
# 2  S10 - S30 -2.8224131 0.004766373 0.01429912
# 3  S20 - S30 -1.9575281 0.050285413 0.10057083

# a ab b

#Create data frame to save Dunn results
#INFORMATION USED IN SUPPLEMENTARY TABLE 6

Dunn_Full <- data.frame()

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("nod1 - Sal", "nod1 - Sal", "nod1 - Sal", "nod1 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(nod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 31.668, df = 6, p-value = 1.889e-05
gene_stats = rbind(gene_stats, c("nod1 - DPH", 31.668, 1.889e-05))

#dunnTest - FSA package
#cldList - rcompanion package
Dunn <- dunnTest(nod1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn
# Comparison          Z      P.unadj        P.adj
# 1   D12 - D15  1.2780430 2.012343e-01 1.0000000000
# 2   D12 - D18 -2.3230601 2.017593e-02 0.2824629625
# 3   D15 - D18 -3.5629439 3.667191e-04 0.0066009441
# 4   D12 - D20 -1.7040615 8.836959e-02 0.8836959499
# 5   D15 - D20 -2.8995612 3.736853e-03 0.0597896531
# 6   D18 - D20  0.5217630 6.018354e-01 1.0000000000
# 7   D12 - D24 -0.5768761 5.640232e-01 1.0000000000
# 8   D15 - D24 -1.7723758 7.633220e-02 0.8396541801
# 9   D18 - D24  1.6193362 1.053750e-01 0.9483745962
# 10  D20 - D24  1.0627206 2.879087e-01 1.0000000000
# 11   D12 - D3  2.1016707 3.558214e-02 0.4625677613
# 12   D15 - D3  0.8236277 4.101511e-01 1.0000000000
# 13   D18 - D3  4.3619801 1.288906e-05 0.0002706703
# 14   D20 - D3  3.6699944 2.425558e-04 0.0046085598
# 15   D24 - D3  2.5428089 1.099653e-02 0.1649480078
# 16   D12 - D6  1.4895826 1.363340e-01 1.0000000000
# 17   D15 - D6  0.2496988 8.028203e-01 0.8028202700
# 18   D18 - D6  3.7052227 2.112051e-04 0.0042241019
# 19   D20 - D6  3.0578222 2.229518e-03 0.0379018101
# 20   D24 - D6  1.9602490 4.996669e-02 0.5996002805
# 21    D3 - D6 -0.5493374 5.827739e-01 1.0000000000
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

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("nod1 - DPH", "nod1 - DPH", "nod1 - DPH", "nod1 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8A
nod1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=nod1, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("nod1 - pathogen recognition") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(3.5, 4.3, 3.9, 4.8, 5.05, 5.1, 4.7), label = c("a", "ab", "ab", "abc", "bc", "c", "bc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
nod1

tiff('Figure 8A nod1 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
nod1
dev.off()

##SUPPLEMENTARY FIGURE 9
nod1_Sal <- ggplot(VST_cts_boxplot, aes(x=Salinity, y=nod1)) +   
  geom_boxplot() + 
  theme_bw() +
  ggtitle("nod1 - pathogen recognition") +
  xlab("Salinity") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("10", "20", "30")) +
  annotate("text", x = 1:3 , y = c(5.1, 4.85, 5), label = c("a", "ab", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
nod1_Sal

tiff('Supplementary Figure 9 nod1 VST Gene Counts by Salinity.tiff', units="in", width=6, height=6, res=300)
nod1_Sal
dev.off()

#ripk2
#Salinity
kruskal.test(ripk2 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 0.79731, df = 2, p-value = 0.6712
gene_stats = rbind(gene_stats, c("ripk2 - Salinity", 0.79731, 0.6712))

Dunn <- dunnTest(ripk2 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 0.75309404 0.4513934 0.9027868
# 2  S10 - S30 0.80087461 0.4232042 1.0000000
# 3  S20 - S30 0.01935199 0.9845603 0.9845603

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("ripk2 - Sal", "ripk2 - Sal", "ripk2 - Sal", "ripk2 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(ripk2 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 31.287, df = 6, p-value = 2.234e-05
gene_stats = rbind(gene_stats, c("ripk2 - DPH", 31.287, 2.234e-05))


Dunn <- dunnTest(ripk2 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3      c          c
# 7    D6    abc        abc
# 1   D12     ab        ab 
# 2   D15      a        a  
# 3   D18     ab        ab 
# 4   D20     ab        ab 
# 5   D24     bc         bc


#a, abc, cb, c, cb, cb, ab

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("ripk2 - DPH", "ripk2 - DPH", "ripk2 - DPH", "ripk2 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8B
ripk2 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=ripk2, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("ripk2 - nod1 signal mediator") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(6.7, 6.7, 6.05, 5.6, 5.95, 5.6, 6.35), label = c("a", "abc", "bc", "c", "bc", "bc", "ab"), size=5) +   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
ripk2

tiff('Figure 8B ripk2 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
ripk2
dev.off()

#IL-6 Like
#Salinity
kruskal.test(LOC111223243 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 3.4898, df = 2, p-value = 0.1747
gene_stats = rbind(gene_stats, c("il6l - Salinity", 3.4898, 0.1747))

Dunn <- dunnTest(LOC111223243 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison         Z   P.unadj     P.adj
# 1  S10 - S20  0.000000 1.0000000 1.0000000
# 2  S10 - S30 -1.596946 0.1102778 0.3308335
# 3  S20 - S30 -1.596946 0.1102778 0.2205556

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("il6l - Sal", "il6l - Sal", "il6l - Sal", "il6l - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(LOC111223243 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 4.4127, df = 6, p-value = 0.621
gene_stats = rbind(gene_stats, c("il6l - DPH", 4.4127, 0.621))

Dunn <- dunnTest(LOC111223243 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 1   D12      a          a
# 2   D15      a          a
# 3   D18      a          a
# 4    D2      a          a
# 5   D24      a          a
# 6    D3      a          a
# 7    D6      a          a

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("il6l - DPH", "il6l - DPH", "il6l - DPH", "il6l - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8C
il6l <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=LOC111223243, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("il6-like - proinflammatory cytokine") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(1.95, 1.15, 1.95, 1.15, 1.15, 1.15, 1.15), label = c("a", "a", "a", "a", "a", "a", "a"), size=5) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
il6l

tiff('Figure 8C il6l VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
il6l
dev.off()

#IL-8 LIKE
#Salinity
kruskal.test(LOC111231462 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 1.9649, df = 2, p-value = 0.3744
gene_stats = rbind(gene_stats, c("il8l - Salinity", 1.9649, 0.3744))

Dunn <- dunnTest(LOC111231462 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20  0.6828053 0.4947299 0.4947299
# 2  S10 - S30 -0.6922058 0.4888081 0.9776163
# 3  S20 - S30 -1.4007863 0.1612780 0.4838340

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("il8l - Sal", "il8l - Sal", "il8l - Sal", "il8l - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(LOC111231462 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 10.066, df = 6, p-value = 0.1219
gene_stats = rbind(gene_stats, c("il8l - DPH", 10.066, 0.1219))

Dunn <- dunnTest(LOC111231462 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3     ab         ab
# 7    D6     ab         ab
# 1   D12      a         a 
# 2   D15     ab         ab
# 3   D18     ab         ab
# 4   D20     ab         ab
# 5   D24      b          b

#ab, ab, a, ab, ab, ab, b

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("il8l - DPH", "il8l - DPH", "il8l - DPH", "il8l - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8D
il8l <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=LOC111231462, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("il8-like - proinflammatory cytokine") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(3.55, 3.5, 3, 3.75, 4.1, 3.4, 4.05), label = c("ab", "ab", "a", "ab", "ab", "ab", "b"), size=5) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
il8l

tiff('Figure 8D il8l VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
il8l
dev.off()

#tlr5
#Salinity
kruskal.test(tlr5 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 1.2464, df = 2, p-value = 0.5362
gene_stats = rbind(gene_stats, c("tlr5 - Salinity", 1.2464, 0.5362))

Dunn <- dunnTest(tlr5 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison           Z   P.unadj     P.adj
# 1  S10 - S20 -0.96427293 0.3349091 0.6698182
# 2  S10 - S30 -0.98131505 0.3264374 0.9793122
# 3  S20 - S30  0.01935826 0.9845553 0.9845553

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("tlr5 - Sal", "tlr5 - Sal", "tlr5 - Sal", "tlr5 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(tlr5 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 19.15, df = 6, p-value = 0.003917
gene_stats = rbind(gene_stats, c("tlr5 - DPH", 19.15, 0.003917))

Dunn <- dunnTest(tlr5 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3      b          b
# 7    D6     ab         ab
# 1   D12     ab         ab
# 2   D15     ab         ab
# 3   D18     ab         ab
# 4   D20      a         a 
# 5   D24     ab         ab


#a, ab, ab, ab, ab, b, ab

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("tlr5 - DPH", "tlr5 - DPH", "tlr5 - DPH", "tlr5 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8E
tlr5 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=tlr5, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("tlr5 - pathogen recognition") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(4.05, 3.9, 3.75, 4.15, 4.05, 4.45, 4.35), label = c("a", "ab", "ab", "ab", "ab", "b", "ab"), size=5) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
tlr5

tiff('Figure 8E tlr5 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
tlr5
dev.off()

#tab1
#Salinity
kruskal.test(tab1 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 0.12641, df = 2, p-value = 0.9388
gene_stats = rbind(gene_stats, c("tab1 - Salinity", 0.12641, 0.9388))

Dunn <- dunnTest(tab1 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 -0.3414026 0.7328005 1.0000000
# 2  S10 - S30 -0.0922941 0.9264644 0.9264644
# 3  S20 - S30  0.2619962 0.7933244 1.0000000

# a a a

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("tab1 - Sal", "tab1 - Sal", "tab1 - Sal", "tab1 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(tab1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 29.539, df = 6, p-value = 4.809e-05
gene_stats = rbind(gene_stats, c("tab1 - DPH", 29.539, 4.809e-05))

#dunnTest - FSA package
#cldList - rcompanion package
Dunn <- dunnTest(tab1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn
# Comparison          Z      P.unadj        P.adj
# 1   D12 - D15  1.0224344 3.065754e-01 1.0000000000
# 2   D12 - D18  2.2593438 2.386201e-02 0.3340680901
# 3   D15 - D18  1.2674368 2.049992e-01 1.0000000000
# 4   D12 - D20  1.8710520 6.133788e-02 0.6747166789
# 5   D15 - D20  0.9146522 3.603743e-01 1.0000000000
# 6   D18 - D20 -0.2993381 7.646821e-01 0.7646820607
# 7   D12 - D24 -0.9810688 3.265588e-01 1.0000000000
# 8   D15 - D24 -1.9374686 5.268809e-02 0.6849451332
# 9   D18 - D24 -3.0765308 2.094246e-03 0.0376964367
# 10  D20 - D24 -2.6890053 7.166529e-03 0.1218309930
# 11   D12 - D3 -2.3146778 2.063057e-02 0.3094585567
# 12   D15 - D3 -3.3371122 8.465377e-04 0.0160842159
# 13   D18 - D3 -4.5049111 6.640072e-06 0.0001394415
# 14   D20 - D3 -4.0362348 5.431586e-05 0.0010863172
# 15   D24 - D3 -1.1841140 2.363679e-01 1.0000000000
# 16   D12 - D6  0.3065268 7.592036e-01 1.0000000000
# 17   D15 - D6 -0.6853802 4.931040e-01 1.0000000000
# 18   D18 - D6 -1.8977970 5.772283e-02 0.6926739676
# 19   D20 - D6 -1.5341079 1.250031e-01 1.0000000000
# 20   D24 - D6  1.2430848 2.138366e-01 1.0000000000
# 21    D3 - D6  2.5520942 1.070776e-02 0.1713241234

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# 6    D3      c          c
# 7    D6    abc        abc 
# 1   D12    abc        abc
# 2   D15     ab        ab 
# 3   D18      a        a  
# 4    D2     ab        ab 
# 5   D24     bc         bc


#"a", "abc", "abc", "bc", "c", "bc", "ab"

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("tab1 - DPH", "tab1 - DPH", "tab1 - DPH", "tab1 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 8F
tab1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=tab1, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("tab1 - tlr5 signal mediator") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(7.9, 7.3, 7.4, 7.25, 7.1, 7, 7.65), label = c("a", "abc", "abc", "bc", "c", "bc", "ab"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
tab1

tiff('Figure 8F tab1 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
tab1
dev.off()

##FIGURE 8

#Make a multi-part figure, removing the axis labels from each separate figure
target_infection_immunity_gene_figure <- ggarrange(nod1 + rremove("ylab") + rremove("xlab"), ripk2 + rremove("ylab") + rremove("xlab"), il6l + rremove("ylab") + rremove("xlab"), il8l + rremove("ylab") + rremove("xlab"), tlr5 + rremove("ylab") + rremove("xlab"), tab1 + rremove("ylab") + rremove("xlab"),
                                                   labels = c("A", "B", "C", "D", "E", "F"),
                                                   ncol = 2, nrow = 3,
                                                   common.legend = TRUE)
target_infection_immunity_gene_figure

#Add annotations to the above figure
target_infection_immunity_gene_figure <- annotate_figure(target_infection_immunity_gene_figure, 
                                                         bottom = text_grob("Days Post Hatch", size = 16),
                                                         left = text_grob("VST Gene Count", rot = 90, size =16))
target_infection_immunity_gene_figure

tiff('Figure 8 Target Infection Immunity Genes Figure.tiff', units="in", width=10, height=10, res=300)
target_infection_immunity_gene_figure
dev.off()

##SALINITY RELATED GENES##

#tlr2
#Salinity
kruskal.test(tlr2 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 0.37271, df = 2, p-value = 0.83
gene_stats = rbind(gene_stats, c("tlr2 - Salinity", 0.37271, 0.83))

Dunn <- dunnTest(tlr2 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 -0.3464719 0.7289881 1.0000000
# 2  S10 - S30  0.2499224 0.8026474 0.8026474
# 3  S20 - S30  0.6094733 0.5422108 1.0000000

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("tlr2 - Sal", "tlr2 - Sal", "tlr2 - Sal", "tlr2 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(tlr2 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 12.578, df = 6, p-value = 0.05024
gene_stats = rbind(gene_stats, c("tlr2 - DPH", 12.578, 0.05024))

Dunn <- dunnTest(tlr2 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 1   D12      a          a
# 2   D15      a          a
# 3   D18      a          a
# 4    D2      a          a
# 5   D24      a          a
# 6    D3      a          a
# 7    D6      a          a

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("tlr2 - DPH", "tlr2 - DPH", "tlr2 - DPH", "tlr2 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9A
tlr2 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=tlr2, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("tlr2 - pathogen recognition") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(2.75, 2.8, 2.65, 3.05, 2.8, 2.6, 2.55), label = c("a", "a", "a", "a", "a", "a", "a"), size=5) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
tlr2

tiff('Figure 9A tlr2 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
tlr2
dev.off()

#lyz
#Salinity
kruskal.test(lyz ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 0.066315, df = 2, p-value = 0.9674
gene_stats = rbind(gene_stats, c("lyz - Salinity", 0.066315, 0.9674))

Dunn <- dunnTest(lyz ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison            Z   P.unadj     P.adj
# 1  S10 - S20 -0.220907585 0.8251644 1.0000000
# 2  S10 - S30 -0.227758021 0.8198344 1.0000000
# 3  S20 - S30  0.001488615 0.9988123 0.9988123

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("lyz - Sal", "lyz - Sal", "lyz - Sal", "lyz - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(lyz ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 26.98, df = 6, p-value = 0.000146
gene_stats = rbind(gene_stats, c("lyz - DPH", 26.98, 0.000146))

Dunn <- dunnTest(lyz ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# 6    D3      c          c
# 7    D6    abc        abc
# 1   D12      a        a  
# 2   D15     ab        ab 
# 3   D18     bc         bc
# 4   D20    abc        abc
# 5   D24    abc        abc


#a, abc, c, bc, ab, abc, abc

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("lyz - DPH", "lyz - DPH", "lyz - DPH", "lyz - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9B
lyz <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=lyz, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("lyz - antimicrobial activity") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(6.4, 4.85, 4.2, 4.8, 5.15, 7.05, 5.8), label = c("a", "abc", "c", "bc", "ab", "abc", "abc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
lyz

tiff('Figure 9B lyz VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
lyz
dev.off()

#cat
#Salinity
kruskal.test(cat ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 1.6632, df = 2, p-value = 0.4354
gene_stats = rbind(gene_stats, c("cat - Salinity", 1.6632, 0.4354))

Dunn <- dunnTest(cat ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 -0.8233828 0.4102903 0.8205807
# 2  S10 - S30 -1.2757426 0.2020465 0.6061395
# 3  S20 - S30 -0.4212779 0.6735522 0.6735522

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("cat - Sal", "cat - Sal", "cat - Sal", "cat - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(cat ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 28.425, df = 6, p-value = 7.815e-05
gene_stats = rbind(gene_stats, c("cat - DPH", 28.425, 7.815e-05))

Dunn <- dunnTest(cat ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

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

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("cat - DPH", "cat - DPH", "cat - DPH", "cat - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9C
cat <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=cat, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("cat - antioxidant defense") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.25, 9.85, 7.55, 9.1, 9.1, 9.3, 8.85), label = c("ac", "abc", "a", "abc", "bc", "b", "abc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
cat

tiff('Figure 9C cat VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
cat
dev.off()

#gsto1
#Salinity
kruskal.test(gsto1 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 3.8221, df = 2, p-value = 0.1479

gene_stats = rbind(gene_stats, c("gsto1 - Salinity", 3.8221, 0.1479))

Dunn <- dunnTest(gsto1 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison         Z    P.unadj     P.adj
# 1  S10 - S20 -1.184868 0.23606966 0.4721393
# 2  S10 - S30 -1.944131 0.05187971 0.1556391
# 3  S20 - S30 -0.714535 0.47489644 0.4748964

# a a a

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("gsto1 - Sal", "gsto1 - Sal", "gsto1 - Sal", "gsto1 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(gsto1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 38.697, df = 6, p-value = 8.207e-07
gene_stats = rbind(gene_stats, c("gsto1 - DPH", 38.697, 8.207e-07))

#dunnTest - FSA package
#cldList - rcompanion package
Dunn <- dunnTest(gsto1 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn
# Comparison           Z      P.unadj        P.adj
# 1   D12 - D15 -1.37744631 1.683743e-01 1.0000000000
# 2   D12 - D18 -3.07043447 2.137476e-03 0.0299246622
# 3   D15 - D18 -1.73411526 8.289763e-02 0.7460787095
# 4   D12 - D20 -2.77052318 5.596632e-03 0.0727562114
# 5   D15 - D20 -1.48204014 1.383296e-01 1.0000000000
# 6   D18 - D20  0.18500760 8.532231e-01 1.0000000000
# 7   D12 - D24 -3.48782302 4.869703e-04 0.0077915249
# 8   D15 - D24 -2.19933998 2.785376e-02 0.3342450859
# 9   D18 - D24 -0.51344805 6.076379e-01 1.0000000000
# 10  D20 - D24 -0.67627677 4.988649e-01 1.0000000000
# 11   D12 - D3  0.78102626 4.347871e-01 1.0000000000
# 12   D15 - D3  2.15847257 3.089111e-02 0.3398021978
# 13   D18 - D3  3.82814123 1.291147e-04 0.0024531786
# 14   D20 - D3  3.50110635 4.633309e-04 0.0078766255
# 15   D24 - D3  4.21840618 2.460353e-05 0.0005166741
# 16   D12 - D6  0.66471548 5.062325e-01 1.0000000000
# 17   D15 - D6  2.00103469 4.538865e-02 0.4538865166
# 18   D18 - D6  3.62991328 2.835165e-04 0.0051032962
# 19   D20 - D6  3.32182180 8.943179e-04 0.0134147687
# 20   D24 - D6  4.02027744 5.812964e-05 0.0011625928
# 21    D3 - D6 -0.09299129 9.259105e-01 0.9259104847


#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3      b         b 
# 7    D6      b         b 
# 1   D12     ab        ab 
# 2   D15    abc        abc
# 3   D18      c          c
# 4    D2     ac        a c
# 5   D24      c          c


#"a", "a", "ab", "abc", "c", "bc", "c"

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("gsto1 - DPH", "gsto1 - DPH", "gsto1 - DPH", "gsto1 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9D
gsto1 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=gsto1, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("gsto1 - glutathione transferase activity") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(8.5, 8.3, 8.7, 9.65, 9.9, 10.05, 10.3), label = c("a", "a", "ab", "abc", "c", "bc", "c"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
gsto1

tiff('Figure 9D gsto1 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
gsto1
dev.off()


#gss
#Salinity
kruskal.test(gss ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 3.2834, df = 2, p-value = 0.1936

gene_stats = rbind(gene_stats, c("gss - Salinity", 3.2834, 0.1936))

Dunn <- dunnTest(gss ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z    P.unadj     P.adj
# 1  S10 - S20 -1.2551567 0.20942187 0.4188437
# 2  S10 - S30 -1.7684740 0.07698169 0.2309451
# 3  S20 - S30 -0.4659363 0.64126107 0.6412611

# a a a

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("gss - Sal", "gss - Sal", "gss - Sal", "gss - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(gss ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 33.163, df = 6, p-value = 9.753e-06
gene_stats = rbind(gene_stats, c("gss - DPH", 33.163, 9.753e-06))

#dunnTest - FSA package
#cldList - rcompanion package
Dunn <- dunnTest(gss ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn
# Comparison            Z      P.unadj        P.adj
# 1   D12 - D15 -0.766825781 4.431851e-01 1.000000e+00
# 2   D12 - D18 -2.455658750 1.406266e-02 1.968773e-01
# 3   D15 - D18 -1.711728469 8.694672e-02 7.825205e-01
# 4   D12 - D20 -3.243030216 1.182657e-03 2.247048e-02
# 5   D15 - D20 -2.525730379 1.154581e-02 1.847329e-01
# 6   D18 - D20 -0.852282188 3.940575e-01 1.000000e+00
# 7   D12 - D24 -0.715402218 4.743605e-01 1.000000e+00
# 8   D15 - D24  0.001897619 9.984859e-01 9.984859e-01
# 9   D18 - D24  1.608942472 1.076289e-01 8.610314e-01
# 10  D20 - D24  2.383070530 1.716890e-02 2.231958e-01
# 11   D12 - D3  1.917064453 5.522974e-02 5.522974e-01
# 12   D15 - D3  2.683890234 7.277099e-03 1.237107e-01
# 13   D18 - D3  4.315484451 1.592533e-05 3.185065e-04
# 14   D20 - D3  5.036279808 4.746666e-07 9.967998e-06
# 15   D24 - D3  2.508651811 1.211929e-02 1.817893e-01
# 16   D12 - D6 -0.115378076 9.081455e-01 1.000000e+00
# 17   D15 - D6  0.628552205 5.296423e-01 1.000000e+00
# 18   D18 - D6  2.274344005 2.294531e-02 2.753437e-01
# 19   D20 - D6  3.049507243 2.292171e-03 4.125908e-02
# 20   D24 - D6  0.588282583 5.563426e-01 1.000000e+00
# 21    D3 - D6 -1.975203777 4.824503e-02 5.306953e-01

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3      b         b 
# 7    D6     ab        ab
# 1   D12     ab        ab 
# 2   D15    abc        abc
# 3   D18     ac        a c
# 4    D2      c          c
# 5   D24    abc        abc


#"a", "ab", "ab", "abc", "bc", "c", "abc"

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("gss - DPH", "gss - DPH", "gss - DPH", "gss - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9E
gss <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=gss, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("gss - glutathione synthase activity") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(5.6, 7.1, 6.7, 6.85, 6.95, 7.4, 7.15), label = c("a", "ab", "ab", "abc", "bc", "c", "abc"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
gss

tiff('Figure 9E gss VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
gss
dev.off()

#scara3
#Salinity
kruskal.test(scara3 ~ Salinity, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 3.0639, df = 2, p-value = 0.2161
gene_stats = rbind(gene_stats, c("scara3 - Salinity", 3.0639, 0.2161))

Dunn <- dunnTest(scara3 ~ Salinity, data=VST_cts_boxplot)
Dunn
# Comparison          Z   P.unadj     P.adj
# 1  S10 - S20 -1.0643729 0.2871598 0.5743196
# 2  S10 - S30 -1.7401904 0.0818256 0.2454768
# 3  S20 - S30 -0.6356384 0.5250121 0.5250121

#Save Dunn results as a data frame, create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn = as.data.frame(Dunn$res)
Dunn_Full = rbind(Dunn_Full, c("scara3 - Sal", "scara3 - Sal", "scara3 - Sal", "scara3 - Sal"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

#DPH
kruskal.test(scara3 ~ Days_Post_Hatch, data=VST_cts_boxplot)
#Kruskal-Wallis chi-squared = 32.076, df = 6, p-value = 1.578e-05
gene_stats = rbind(gene_stats, c("scara3 - DPH", 32.076, 1.578e-05))

Dunn <- dunnTest(scara3 ~ Days_Post_Hatch, data=VST_cts_boxplot)
Dunn

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)
# Group Letter MonoLetter
# 6    D3      b         b 
# 7    D6      b         b
# 1   D12    abc        abc
# 2   D15     ab        ab 
# 3   D18      c          c
# 4   D20     ac        a c
# 5   D24    abc        abc


#a, a, abc, ab, c, bc, abc

#Create spacing/testing line in full Dunn data frame, match column names between Dunn and full Dunn data frames, and then add Dunn results to full Dunn data frame
Dunn_Full = rbind(Dunn_Full, c("scara3 - DPH", "scara3 - DPH", "scara3 - DPH", "scara3 - DPH"))
colnames(Dunn_Full) = colnames(Dunn)
Dunn_Full = rbind(Dunn_Full, Dunn)

##FIGURE 9F
scara3 <- ggplot(VST_cts_boxplot, aes(x=Days_Post_Hatch, y=scara3, fill=Feeding)) +   
  geom_boxplot() + 
  scale_fill_manual(name="Feeding Type", limits=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "Eart + MF", "MF"), values=c("forestgreen", "brown", "blue", "orange", "purple"),labels=c("GW + Rot", "GW + Art + EArt", "Art + EArt + MF", "EArt + MF", "MF"))+
  theme_bw() +
  ggtitle("scara3 - deplete reactive oxidative species") +
  xlab("Days Post Hatch") +
  ylab("VST Gene Count") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D3", "D6", "D12", "D15", "D18", "D20", "D24"), labels=c("3", "6", "12", "15", "18", "20", "24")) +
  annotate("text", x = 1:7 , y = c(5.1, 6, 7.35, 6.45, 7.35, 7.4, 7.25), label = c("a", "a", "abc", "ab", "c", "bc", "abc"), size=5) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
scara3

tiff('Figure 9F scara3 VST Gene Counts by DPH.tiff', units="in", width=6, height=6, res=300)
scara3
dev.off()

##FIGURE 9

#Make a multi-part figure, removing the axis labels from each separate figure
target_salinity_immunity_gene_figure <- ggarrange(tlr2 + rremove("ylab") + rremove("xlab"), lyz + rremove("ylab") + rremove("xlab"), cat + rremove("ylab") + rremove("xlab"), gsto1 + rremove("ylab") + rremove("xlab"), gss + rremove("ylab") + rremove("xlab"), scara3 + rremove("ylab") + rremove("xlab"),
                                                  labels = c("A", "B", "C", "D", "E", "F"),
                                                  ncol = 2, nrow = 3,
                                                  common.legend = TRUE)
target_salinity_immunity_gene_figure

#Add annotations to the above figure
target_salinity_immunity_gene_figure <- annotate_figure(target_salinity_immunity_gene_figure, 
                                                        bottom = text_grob("Days Post Hatch", size = 16),
                                                        left = text_grob("VST Gene Count", rot = 90, size =16))
target_salinity_immunity_gene_figure

tiff('Figure 9 Target Salinity Immunity Genes Figure.tiff', units="in", width=10, height=10, res=300)
target_salinity_immunity_gene_figure
dev.off()

colnames(gene_stats) <- c("Category", "Statistic", "p.value")
gene_stats

gene_stats$padj <- p.adjust(gene_stats$p.value, method = "BH")
gene_stats

clipr::write_clip(gene_stats)
clipr::write_clip(Dunn_Full)
