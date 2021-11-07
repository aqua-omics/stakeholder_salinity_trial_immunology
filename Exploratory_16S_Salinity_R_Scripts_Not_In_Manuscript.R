

##Extra Shannon Diversity Plots

#Sample Type
Shannon_x_Sample_Type_boxplot <- ggplot(Falphadiv_no0D1D_metadata, aes(x=Sample_Type, y=Shannon, color=Sample_Type)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Sample Type", values=c("gray", "lightblue"))+
  ggtitle("Shannon Diversity by Sample Type") +
  xlab("Sample Type") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1.0, 2.0), y = c(5.35, 4.0), label = c("a", "b"), size=5)+ 
  annotate("text", x = c(2), y = c(5.3), label = "Kruskal-Wallis", size=3)+ 
  annotate("text", x = c(2), y = c(5), label = "chi2 = 60.839", size=3)+ 
  annotate("text", x = c(2), y = c(4.7), label = "BH adj. p-value = 9.72E-14", size=3)+ 
  
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
Shannon_x_Sample_Type_boxplot

tiff('Shannon Diversity by Sample Type.tiff', units="in", width=6, height=6, res=300)
Shannon_x_Sample_Type_boxplot
dev.off()

#Water by Salinity
Water_Shannon_Salinity_boxplot <- ggplot(W_Falphadiv_no1D_metadata, aes(x=Salinity, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  ggtitle("Water Shannon Diversity by Salinity") +
  xlab("Salinity") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1.0, 2.0, 3.0), y = c(3.4, 3.5, 4), label = c("a", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  annotate("text", x = c(2), y = c(4.2), label = "Kruskal-Wallis", size=3)+ 
  annotate("text", x = c(2), y = c(4.0), label = "chi^2 = 5.71", size=3)+ 
  annotate("text", x = c(2), y = c(3.8), label = "BH adj. p-value = 0.115", size=3)
Water_Shannon_Salinity_boxplot

tiff('Water Shannon Diversity by Salinity.tiff', units="in", width=10, height=6, res=300)
Water_Shannon_Salinity_boxplot
dev.off()

#Tissue by Salinity
Tissue_Shannon_Salinity_boxplot <- ggplot(T_Falphadiv_no0D1D_metadata, aes(x=Salinity, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  ggtitle("Tissue Shannon Diversity by Salinity") +
  xlab("Salinity") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1.0, 2.0, 3.0), y = c(5.35, 5.3, 5.3), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  annotate("text", x = c(2), y = c(5.9), label = "Kruskal-Wallis", size=3)+ 
  annotate("text", x = c(2), y = c(5.7), label = "chi^2 = 0.141", size=3)+ 
  annotate("text", x = c(2), y = c(5.5), label = "BH adj. p-value = 0.932", size=3)
Tissue_Shannon_Salinity_boxplot

tiff('Tisue Shannon Diversity by Salinity.tiff', units="in", width=10, height=6, res=300)
Tissue_Shannon_Salinity_boxplot
dev.off()


#Water DPH by Salinity
W_Shannon_DPH_Salinity_boxplot <- ggplot(W_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"))+
  ggtitle("Water Shannon Diversity by Days Post Hatch and Salinity") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24")) +
  annotate("text", x = c(1.75, 2.0, 2.25, 2.75, 3.0, 3.25, 3.75, 4.0, 4.25, 4.75, 5.0, 5.25 ) , y = c(2.85, 3.1, 2.55, 2.15, 2.85, 2.3, 2.65, 3.4, 3.95, 3.35, 3.05, 3.5), label = c("a", "ab", "b", "a", "b", "b", "a", "b", "c", "a", "a", "b"), size=5)+ 
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
W_Shannon_DPH_Salinity_boxplot

tiff('Water Shannon Diversity by Days Post Hatch by Salinity.tiff', units="in", width=10, height=6, res=300)
W_Shannon_DPH_Salinity_boxplot
dev.off()

#Tissue DPH by Salinity
T_Shannon_DPH_Salinity_boxplot <- ggplot(T_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon, color=Salinity)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+ 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  ggtitle("Tissue Shannon Diversity by Days Post Hatch and Salinity") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
  annotate("text", x = c(2.75, 3.0, 3.25, 3.75, 4.0, 4.25, 4.75, 5.0, 5.25, 5.75, 6.0, 6.25, 6.75, 7.0, 7.25, 7.75, 8.0, 8.25, 8.75, 9.0, 9.25, 9.75, 10.0, 10.25 ) , y = c(4.85, 5.3, 4.85, 5.1, 4.75, 5.25, 5.35, 4.65, 5.25, 4.6, 4.5, 4.25, 2.85, 3.2, 2.75, 3.8, 3.65, 4.4, 3.4, 3.45, 3.5, 4.5, 4.8, 4.3), label = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"), size=5)+
  annotate("segment", x=1, xend=6.0, y=5.5, yend=5.5, size=2, color = "forestgreen")+ #0-12 dph
  annotate("text", x=3.5, y=5.625, label="Green Water", size=5, color="forestgreen")+
  annotate("segment", x=2, xend=5.5, y=5.75, yend=5.75, size=2, color = "red") + #1-10 dph
  annotate("text", x=3.75, y=5.875, label="Enriched Rotifers", size=5, color ="red") +  
  annotate("segment", x=5, xend=7, y=6, yend=6, size=2, color = "blue") + #9-15 dph
  annotate("text", x=6, y=6.125, label="Artemia nauplii", size=5, color ="blue")+
  annotate("segment", x=6, xend=9.5, y=6.25, yend=6.25, size=2, color = "orange")+ #12-21 dph
  annotate("text", x=7.75, y=6.375, label="Enriched Artemia", size=5, color ="orange")+
  annotate("segment", x=6.25, xend=10, y=6.5, yend=6.5, size=2, color = "purple")+ #13-24 dph
  annotate("text", x=8.125, y=6.625, label="Microfeeds", size=5, color ="purple")+ 
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
T_Shannon_DPH_Salinity_boxplot

tiff('Tissue Shannon Diversity by Days Post Hatch by Salinity.tiff', units="in", width=10, height=6, res=300)
T_Shannon_DPH_Salinity_boxplot
dev.off()

###Other PCoAs/NMDSs

#Note: In general NMDSs had large stress values and thus were not used over the PCoA. This was tested because other data was visualized with NMDSs such as the fatty acid information.

#Sample Type NMDS 
ZDS_Fdata_NMDS_Bray <- ordinate(ZDS_Fdata, "NMDS", "bray")
ZDS_Fdata_NMDS_Bray #stress = 0.2147626 

plot_ordination(ZDS_Fdata, ZDS_Fdata_NMDS_Bray, color = "Sample_Type")+ 
  scale_color_manual(name="Sample Type", values = c("gray", "lightblue"))+
  ggtitle("Overall PCoA by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Overall PCoA by Sample_Type.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

#Water NMDS with 1 DPH samples
W_ZDS_Fdata_NMDS_Bray <- ordinate(W_ZDS_Fdata, "NMDS", "bray")
W_ZDS_Fdata_NMDS_Bray #stress = 0.1817207  

plot_ordination(W_ZDS_Fdata, W_ZDS_Fdata_NMDS_Bray, color = "Days_Post_Hatch") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D1", "D6", "D12", "D18", "D24"), values=c("midnightblue", "hotpink", "red", "yellow", "orange")) +
  ggtitle("Water PCoA by Days Post Hatch") +
  theme(plot.title = element_text(hjust = 0.5))


tiff('Water PCoA by DPH & Salinity.tiff', units="in", width=10, height=6, res=300)
W_PCoA_DPH_Salinity
dev.off()

#Water PCoA with 1 DPH samples
W_PCoA_DPH_Salinity <- plot_ordination(W_ZDS_Fdata, ordinate(W_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D1", "D6", "D12", "D18", "D24"), values=c("midnightblue", "red", "forestgreen", "purple", "orange")) +
  ggtitle("Water PCoA by Days Post Hatch") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) +
  annotate("text", x = c(-0.3) , y = c(0), label = "Salinity PERMANOVA", size=4)+
  annotate("text", x = c(-0.3) , y = c(-0.025), label = "Pseudo-F = 4.96", size=4)+
  annotate("text", x = c(-0.3) , y = c(-0.05), label = "P(MC) = 0.0001", size=4) +
  annotate("text", x = c(-0.3) , y = c(-0.1), label = "DPH PERMANOVA", size=4)+
  annotate("text", x = c(-0.3) , y = c(-0.125), label = "Pseudo-F = 9.34", size=4)+
  annotate("text", x = c(-0.3) , y = c(-0.15), label = "P(MC) = 0.0001", size=4)
W_PCoA_DPH_Salinity

#Water NMDS without 1 DPH samples
W_noD1_ZDS_Fdata_NMDS_Bray <- ordinate(W_noD1_ZDS_Fdata, "NMDS", "bray")
W_noD1_ZDS_Fdata_NMDS_Bray #stress = 0.1738464   

plot_ordination(W_noD1_ZDS_Fdata, W_noD1_ZDS_Fdata_NMDS_Bray, color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D6", "D12", "D18", "D24"), values=c("red", "forestgreen", "purple", "orange")) +
  ggtitle("Water PCoA by DPH & Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Water PCoA by DPH & Salinity wo D1.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

# Tissue PCoA with 0 & 1 DPH samples
plot_ordination(T_ZDS_Fdata, ordinate(T_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("blue", "midnightblue", "rosybrown", "hotpink",  "green", "red", "purple", "yellow", "dark gray", "orange"))+
  ggtitle("Tissue PCoA by Days Post Hatch") +
  theme(plot.title = element_text(hjust = 0.5))

# Tissue NMDS with 0 & 1 DPH samples
T_ZDS_Fdata_NMDS_Bray <- ordinate(T_ZDS_Fdata, "NMDS", "bray") 
T_ZDS_Fdata_NMDS_Bray #stress = 0.2427299 
plot_ordination(T_ZDS_Fdata, T_ZDS_Fdata_NMDS_Bray, color = "Days_Post_Hatch") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("blue", "midnightblue", "rosybrown", "hotpink",  "green", "red", "purple", "yellow", "dark gray", "orange"))+
  ggtitle("Tissue PCoA by Days Post Hatch") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Tissue PCoA by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

# Tissue NMDS without 0 & 1 DPH samples
T_noD01_ZDS_Fdata_NMDS_Bray <- ordinate(T_noD01_ZDS_Fdata, "NMDS", "bray")
T_noD01_ZDS_Fdata_NMDS_Bray #stress = 0.2400142   

plot_ordination(T_noD01_ZDS_Fdata, T_noD01_ZDS_Fdata_NMDS_Bray, color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red",  "black", "forestgreen", "brown", "purple", "dark gray", "orange"))+
  ggtitle("Tissue PCoA by DPH and Salinity") +
  theme(plot.title = element_text(hjust = 0.5))

#Run through and make Water and Tissue PCoAs by DPH
plot_ordination(W_1D_ZDS_Fdata, ordinate(W_1D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Replicate") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "blue", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))

plot_ordination(W_6D_ZDS_Fdata, ordinate(W_6D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "blue", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))

plot_ordination(T_6D_ZDS_Fdata, ordinate(T_6D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("tomato", "blue", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))#Missing Tank 1

plot_ordination(W_12D_ZDS_Fdata, ordinate(W_12D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "blue", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))

plot_ordination(T_12D_ZDS_Fdata, ordinate(T_12D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "blue", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))

plot_ordination(W_18D_ZDS_Fdata, ordinate(W_18D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))#Missing Tank P3

plot_ordination(T_18D_ZDS_Fdata, ordinate(T_18D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "tomato", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))#Missing Tank P3

plot_ordination(W_24D_ZDS_Fdata, ordinate(W_24D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab"))#Missing Tanks P2 and P3

plot_ordination(T_24D_ZDS_Fdata, ordinate(T_24D_ZDS_Fdata, "PCoA", "bray"), color = "Tank_Number", shape = "Salinity") + 
  scale_color_manual(name="Tank Number", values=c("red", "darkgreen", "steelblue", "skyblue", "indianred", "springgreen", "olivedrab")) #Missing Tanks P2 and P3


##Other Taxonomic Bar Plots
T_1D_phy <- transform_sample_counts(T_1D_Fdata, function(x) 100*x/sum(x))
W_1D_phy <- transform_sample_counts(W_1D_Fdata, function(x) 100*x/sum(x))

T_6D_phy <- transform_sample_counts(T_6D_Fdata, function(x) 100*x/sum(x))
W_6D_phy <- transform_sample_counts(W_6D_Fdata, function(x) 100*x/sum(x))
T_12D_phy <- transform_sample_counts(T_12D_Fdata, function(x) 100*x/sum(x))
W_12D_phy <- transform_sample_counts(W_12D_Fdata, function(x) 100*x/sum(x))
T_18D_phy <- transform_sample_counts(T_18D_Fdata, function(x) 100*x/sum(x))
W_18D_phy <- transform_sample_counts(W_18D_Fdata, function(x) 100*x/sum(x))
T_24D_phy <- transform_sample_counts(T_24D_Fdata, function(x) 100*x/sum(x))
W_24D_phy <- transform_sample_counts(W_24D_Fdata, function(x) 100*x/sum(x))


##Taxonomic Bar Plots


#Make Sample_Type Phylum graph#
##Create table ready for making stacked bar graph for Phyla >1% by Site-Survey combo##

# agglomerate taxa such that the ASVs with the same upper taxonomic ranks are merged at whatever rank you choose#
glom <- tax_glom(phy, taxrank = 'Phylum')

# Create dataframe from phyloseq object, just melts the phyloseq object to a full table#
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R#
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate mean rel. abundance#
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%#
Other <- means[means$mean <= 1,]$Phylum

# change their name to "Other Prokaryotes"#
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'

#remove all Phylums labeled Other Prokaryotes#
dat <-dat[!dat$Phylum=="Other Prokaryotes",]

#remove unnessary columns#
dat <- subset(dat, select=c(Sample_Type, Abundance, Phylum))

#Summarize based upon target parameter#
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Sample_Type","Phylum"), na.rm = TRUE)

#Arrange by Abundance#
dat <- arrange(dat, Abundance)

#Create a table that is the leftover Prokaryotes#
Abundance <- ddply(dat, ~Sample_Type, function(x) c(Abundance=100-sum(x$Abundance)))

#Add a column labeling the leftover Prokaryotes#
Abundance$Phylum<- "Other Prokaryotes"

#remove unnessary columns#
dat <- subset(dat, select=c(Sample_Type, Abundance, Phylum))

#combine with original table#
Phylum_Fdata_Sample_Type <- rbind(dat, Abundance)


#Plot a bar graph by chosen factor (here Sample.Survey) and abundance with fill (taxonomic level) as whatever you chose#
spatial_plot_Phylum_Fdata_Sample_Type <- ggplot(data=Phylum_Fdata_Sample_Type, aes(x=Sample_Type, y=Abundance, fill=Phylum)) +
  
  #Make it a stacked bar graph using the data you calculate (stat="identity" )#
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  
  #Make x axis labels at 45 degrees to make easier to read#
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  
  #Name the fill(legend) whatever you like ("Phylum" here) and uses the values parameter to pick your colors to make it pretty#
  scale_fill_manual("Phylum", values=c("rosybrown", "purple","forestgreen", "orange", "blue"))+
  
  #rename x and y axes#
  xlab("Sample Type") +
  ylab("Percentage") +
  
  #Give it a title#
  ggtitle("Phyla by Sample Type")+
  
  #Center the title#
  theme(plot.title = element_text(hjust = 0.5)) 

#Call the stored plot#
spatial_plot_Phylum_Fdata_Sample_Type

##Save Graph for Publication##

tiff('Phylum by Sample Type.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Phylum_Fdata_Sample_Type
dev.off()

## Make Water by DPH graph
# agglomerate taxa
glom <- tax_glom(Wphy, taxrank = 'Genus')
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
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#combine with original table
Genus_W_Fdata_Days_Post_Hatch <- rbind(dat, Abundance)

spatial_plot_Genus_W_Fdata_Days_Post_Hatch <- ggplot(data=Genus_W_Fdata_Days_Post_Hatch, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "red", "yellow", "orange", "blue", "green", "wheat", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "seagreen", "turquoise", "black", "rosybrown"))+
  ggtitle("Water Genera by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12", "D18", "D24")) +
  annotate("text", x = c(1:5 ) , y = 105, label = c("a", "b", "c", "d", "e"), size=5)
spatial_plot_Genus_W_Fdata_Days_Post_Hatch

tiff('Water Genera by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Genus_W_Fdata_Days_Post_Hatch
dev.off()

#Make Tissue by DPH Graph
# agglomerate taxa
glom <- tax_glom(Tphy, taxrank = 'Genus')
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
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Genus))
#combine with original table
Genus_T_Fdata_Days_Post_Hatch <- rbind(dat, Abundance)

spatial_plot_Genus_T_Fdata_Days_Post_Hatch <- ggplot(data=Genus_T_Fdata_Days_Post_Hatch, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "red", "green", "maroon", "purple", "lightblue", "orangered", "seagreen", "white",  "gray", "black", "darkviolet", "rosybrown"))+
  ggtitle("Tissue Genera by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24")) +
  annotate("text", x = c(3:10 ) , y = 105, label = c("a", "ac", "c", "c", "d", "b", "b", "c"), size=5)

spatial_plot_Genus_T_Fdata_Days_Post_Hatch

tiff('Tissue Genera by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Genus_T_Fdata_Days_Post_Hatch
dev.off()

## Make 6DPH Water by Salinity graph
# agglomerate taxa
glom <- tax_glom(W_6D_phy, taxrank = 'Genus')
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
Genus_W_6D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_W_6D_Fdata_Salinity <- ggplot(data=Genus_W_6D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Genera 6 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_W_6D_Fdata_Salinity

#Make 6 DPH Tissue by Salinity Graph
# agglomerate taxa
glom <- tax_glom(T_6D_phy, taxrank = 'Genus')
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
Genus_T_6D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_T_6D_Fdata_Salinity <- ggplot(data=Genus_T_6D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Tissue Genera 6 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_T_6D_Fdata_Salinity

## Make 12DPH Water by Salinity graph
# agglomerate taxa
glom <- tax_glom(W_12D_phy, taxrank = 'Genus')
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
Genus_W_12D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_W_12D_Fdata_Salinity <- ggplot(data=Genus_W_12D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Genera 12 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_W_12D_Fdata_Salinity

#Make 12 DPH Tissue by Salinity Graph
# agglomerate taxa
glom <- tax_glom(T_12D_phy, taxrank = 'Genus')
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
Genus_T_12D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_T_12D_Fdata_Salinity <- ggplot(data=Genus_T_12D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Tissue Genera 12 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_T_12D_Fdata_Salinity

## Make 18DPH Water by Salinity graph
# agglomerate taxa
glom <- tax_glom(W_18D_phy, taxrank = 'Genus')
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
Genus_W_18D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_W_18D_Fdata_Salinity <- ggplot(data=Genus_W_18D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Genera 18 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_W_18D_Fdata_Salinity

#Make 18 DPH Tissue by Salinity Graph
# agglomerate taxa
glom <- tax_glom(T_18D_phy, taxrank = 'Genus')
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
Genus_T_18D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_T_18D_Fdata_Salinity <- ggplot(data=Genus_T_18D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Tissue Genera 18 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_T_18D_Fdata_Salinity

## Make 24DPH Water by Salinity graph
# agglomerate taxa
glom <- tax_glom(W_24D_phy, taxrank = 'Genus')
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
Genus_W_24D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_W_24D_Fdata_Salinity <- ggplot(data=Genus_W_24D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Water Genera 24 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_W_24D_Fdata_Salinity

#Make 24 DPH Tissue by Salinity Graph
# agglomerate taxa
glom <- tax_glom(T_24D_phy, taxrank = 'Genus')
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
Genus_T_24D_Fdata_Salinity <- rbind(dat, Abundance)

spatial_plot_Genus_T_24D_Fdata_Salinity <- ggplot(data=Genus_T_24D_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("Tissue Genera 24 DPH by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))
spatial_plot_Genus_T_24D_Fdata_Salinity


##Sample-level analysis##

#Tissue 1 DPH

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(T_1D_phy, taxrank = 'Genus')
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
Genus_T_1D_data_Sample_ID <- rbind(dat, Abundance)

spatial_plot_Genus_T_1D_data_Sample_ID <- ggplot(data=Genus_T_1D_data_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("1 DPH Tissue Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("LST1DP0a", "LST1DP0b"))

spatial_plot_Genus_T_1D_data_Sample_ID

#Water 1 DPH 

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(W_1D_phy, taxrank = 'Genus')
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
Genus_W_1D_data_Sample_ID <- rbind(dat, Abundance)

spatial_plot_Genus_W_1D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("1 DPH Water Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c(
    "LSE1DP3a", "LSE1DP3b", "LSE1DP3c", 
    "LSE1DP5a", "LSE1DP5b", "LSE1DP5c",
    "LSE1DP6a", "LSE1DP6b", "LSE1DP6c",
    "LSE1DP1a", "LSE1DP1b", "LSE1DP1c", 
    "LSE1DP2a", "LSE1DP2b", "LSE1DP2c",
    "LSE1DP7a", "LSE1DP7b", "LSE1DP7c",
    "LSE1DP4a", "LSE1DP4b", "LSE1DP4c",
    "LSE1DP8a", "LSE1DP8b", "LSE1DP8c",
    "LSE1DP9a", "LSE1DP9b", "LSE1DP9c"))+
  annotate("segment", x=0.5, xend=9.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=9.5, xend=18.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=18.5, xend=27.5, y=102, yend=102, color = "blue") +
  annotate("text", x=5, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=14, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=23, y=105, label="30 ppt", size=3, color ="blue")
spatial_plot_Genus_W_1D_data_Sample_ID

#Sample ID Table

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
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
Genus_Fdata_Sample_ID <- rbind(dat, Abundance)

#Plot each set of samples by sample type and DPH using  scale_x_discrete to just select samples you are interested in

#6 DPH Tissue by Sample ID
spatial_plot_Genus_T_6D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("6 DPH Tissue Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("LST6DP3a", "LST6DP5a","LST6DP6a","LST6DP1a","LST6DP2a","LST6DP7a","LST6DP4a","LST6DP8a","LST6DP9a"))+
  annotate("segment", x=0.5, xend=3.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=3.5, xend=6.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=6.5, xend=9.5, y=102, yend=102, color = "blue") +
  annotate("text", x=2, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=5, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=8, y=105, label="30 ppt", size=3, color ="blue")

spatial_plot_Genus_T_6D_data_Sample_ID

#6 DPH Water by Sample ID
spatial_plot_Genus_W_6D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("6 DPH Water Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c(
    "LSE6DP3a", "LSE6DP3b", "LSE6DP3c", 
    "LSE6DP5a", "LSE6DP5b", "LSE6DP5c",
    "LSE6DP6a", "LSE6DP6b", "LSE6DP6c",
    "LSE6DP1a", "LSE6DP1b", "LSE6DP1c", 
    "LSE6DP2a", "LSE6DP2b", "LSE6DP2c",
    "LSE6DP7a", "LSE6DP7b", "LSE6DP7c",
    "LSE6DP4a", "LSE6DP4b", "LSE6DP4c",
    "LSE6DP8a", "LSE6DP8b", "LSE6DP8c",
    "LSE6DP9a", "LSE6DP9b", "LSE6DP9c"))+
  annotate("segment", x=0.5, xend=9.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=9.5, xend=18.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=18.5, xend=27.5, y=102, yend=102, color = "blue") +
  annotate("text", x=5, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=14, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=23, y=105, label="30 ppt", size=3, color ="blue")
spatial_plot_Genus_W_6D_data_Sample_ID




#12 DPH Tissue by Sample ID
spatial_plot_Genus_T_12D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("12 DPH Tissue Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("LST12DP3a", "LST12DP5a","LST12DP6a","LST12DP1a","LST12DP2a","LST12DP7a","LST12DP4a","LST12DP8a","LST12DP9a"))+
  annotate("segment", x=0.5, xend=3.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=3.5, xend=6.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=6.5, xend=9.5, y=102, yend=102, color = "blue") +
  annotate("text", x=2, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=5, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=8, y=105, label="30 ppt", size=3, color ="blue")

spatial_plot_Genus_T_12D_data_Sample_ID

#12 DPH Water by Sample ID
spatial_plot_Genus_W_12D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("12 DPH Water Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c(
    "LSE12DP3a", "LSE12DP3b", "LSE12DP3c", 
    "LSE12DP5a", "LSE12DP5b", "LSE12DP5c",
    "LSE12DP6a", "LSE12DP6b", "LSE12DP6c",
    "LSE12DP1a", "LSE12DP1b", "LSE12DP1c", 
    "LSE12DP2a", "LSE12DP2b", "LSE12DP2c",
    "LSE12DP7a", "LSE12DP7b", "LSE12DP7c",
    "LSE12DP4a", "LSE12DP4b", "LSE12DP4c",
    "LSE12DP8a", "LSE12DP8b", "LSE12DP8c",
    "LSE12DP9a", "LSE12DP9b", "LSE12DP9c"))+
  annotate("segment", x=0.5, xend=9.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=9.5, xend=18.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=18.5, xend=27.5, y=102, yend=102, color = "blue") +
  annotate("text", x=5, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=14, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=23, y=105, label="30 ppt", size=3, color ="blue")
spatial_plot_Genus_W_12D_data_Sample_ID


#18 DPH Tissue by Sample ID
spatial_plot_Genus_T_18D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("18 DPH Tissue Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("LST18DP3a", "LST18DP5a","LST18DP6a","LST18DP1a","LST18DP2a","LST18DP7a","LST18DP4a","LST18DP8a","LST18DP9a"))+
  annotate("segment", x=0.5, xend=3.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=3.5, xend=6.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=6.5, xend=9.5, y=102, yend=102, color = "blue") +
  annotate("text", x=2, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=5, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=8, y=105, label="30 ppt", size=3, color ="blue")

spatial_plot_Genus_T_18D_data_Sample_ID

#18 DPH Water by Sample ID
spatial_plot_Genus_W_18D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("18 DPH Water Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c(
    "LSE18DP3a", "LSE18DP3b", "LSE18DP3c", 
    "LSE18DP5a", "LSE18DP5b", "LSE18DP5c",
    "LSE18DP6a", "LSE18DP6b", "LSE18DP6c",
    "LSE18DP1a", "LSE18DP1b", "LSE18DP1c", 
    "LSE18DP2a", "LSE18DP2b", "LSE18DP2c",
    "LSE18DP7a", "LSE18DP7b", "LSE18DP7c",
    "LSE18DP4a", "LSE18DP4b", "LSE18DP4c",
    "LSE18DP8a", "LSE18DP8b", "LSE18DP8c",
    "LSE18DP9a", "LSE18DP9b", "LSE18DP9c"))+
  annotate("segment", x=0.5, xend=9.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=9.5, xend=18.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=18.5, xend=27.5, y=102, yend=102, color = "blue") +
  annotate("text", x=5, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=14, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=23, y=105, label="30 ppt", size=3, color ="blue")
spatial_plot_Genus_W_18D_data_Sample_ID


#24 DPH Tissue by Sample ID
spatial_plot_Genus_T_24D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("24 DPH Tissue Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("LST24DP3a", "LST24DP5a","LST24DP6a","LST24DP1a","LST24DP2a","LST24DP7a","LST24DP4a","LST24DP8a","LST24DP9a"))+
  annotate("segment", x=0.5, xend=3.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=3.5, xend=6.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=6.5, xend=9.5, y=102, yend=102, color = "blue") +
  annotate("text", x=2, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=5, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=8, y=105, label="30 ppt", size=3, color ="blue")

spatial_plot_Genus_T_24D_data_Sample_ID

#24 DPH Water by Sample ID
spatial_plot_Genus_W_24D_data_Sample_ID <- ggplot(data=Genus_Fdata_Sample_ID, aes(x=Sample_ID, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  ggtitle("24 DPH Water Genera by Sample ID") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c(
    "LSE24DP3a", "LSE24DP3b", "LSE24DP3c", 
    "LSE24DP5a", "LSE24DP5b", "LSE24DP5c",
    "LSE24DP6a", "LSE24DP6b", "LSE24DP6c",
    "LSE24DP1a", "LSE24DP1b", "LSE24DP1c", 
    "LSE24DP2a", "LSE24DP2b", "LSE24DP2c",
    "LSE24DP7a", "LSE24DP7b", "LSE24DP7c",
    "LSE24DP4a", "LSE24DP4b", "LSE24DP4c",
    "LSE24DP8a", "LSE24DP8b", "LSE24DP8c",
    "LSE24DP9a", "LSE24DP9b", "LSE24DP9c"))+
  annotate("segment", x=0.5, xend=9.5, y=102, yend=102, color = "hotpink")+
  annotate("segment", x=9.5, xend=18.5, y=102, yend=102, color = "midnightblue")+
  annotate("segment", x=18.5, xend=27.5, y=102, yend=102, color = "blue") +
  annotate("text", x=5, y=105, label="10 ppt", size=3, color="hotpink")+
  annotate("text", x=14, y=105, label="20 ppt", size=3, color ="midnightblue") +
  annotate("text", x=23, y=105, label="30 ppt", size=3, color ="blue")
spatial_plot_Genus_W_24D_data_Sample_ID




###PPG Focused Species Level Graphs

##Water Pseudomonas

#Keep only Pseudomonas
ntaxa(W_Fdata) #166944#
Pseudomonas_W_Fdata <- subset_taxa(W_Fdata, Genus=="g__Pseudomonas")
ntaxa(Pseudomonas_W_Fdata)#111

##Create a table of the Pseudomonas Species by DPH
Pseudomonas_W_Fdataphy <- transform_sample_counts(Pseudomonas_W_Fdata, function(x) 100*x/sum(x))
nsamples(Pseudomonas_W_Fdataphy)#125
Pseudomonas_W_Fdataphy2 <- prune_samples(sample_sums(Pseudomonas_W_Fdataphy)>0, Pseudomonas_W_Fdataphy)
nsamples(Pseudomonas_W_Fdataphy2)

# agglomerate taxa
glom <- tax_glom(Pseudomonas_W_Fdataphy2, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 0.01,]$Species
# change their name to "Other Pseudomonas"
dat[dat$Species %in% Other,]$Species <- 'Other Pseudomonas'
#remove all Specieses labeled Other Pseudomonas
dat <-dat[!dat$Species=="Other Pseudomonas",]
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Pseudomonas"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Pseudomonas_Species_W_DPH_dat <- rbind(dat, Abundance)

##Create a table of the All Species by DPH

W_Fdataphy <- transform_sample_counts(W_Fdata, function(x) 100*x/sum(x))
# agglomerate taxa
glom <- tax_glom(W_Fdataphy, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Full_Species_W_DPH_dat <- rbind(dat)


#Summarize the Pseudomonas Species Table by DPH in wide format
Wide_Pseudomonas_Species_W_DPH_dat<- spread(Pseudomonas_Species_W_DPH_dat, Species, Abundance)

#Save the Pseudomonas species names as a list
pseudomonas_species <- colnames(Wide_Pseudomonas_Species_W_DPH_dat)

#Use this list to create a filtered version of the Full Species list
Psudomonas_Full_Species_W_DPH_dat <- Full_Species_W_DPH_dat[Full_Species_W_DPH_dat$Species %in% pseudomonas_species,]

#Use this table to create a bar graph of just Pseudomonas Species
spatial_plot_Psudomonas_Full_Species_W_DPH <- ggplot(data=Psudomonas_Full_Species_W_DPH_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Species)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "rosybrown", "black", "seagreen", "turquoise", "orchid", "gold3", "darkolivegreen", "hotpink"))+
  #scale_x_discrete(limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018")) + 
  xlab("Estuary by Sampling Period") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Psudomonas_Full_Species_W_DPH

tiff('Water Pseudomonas Species by DPH.tiff', units="in", width=6, height=6, res=300)
spatial_plot_Psudomonas_Full_Species_W_DPH
dev.off()


## Water Vibrio

##Create a table of the Vibrio Species by DPH
Vibrio_W_Fdataphy <- transform_sample_counts(Vibrio_W_Fdata, function(x) 100*x/sum(x))
nsamples(Vibrio_W_Fdataphy)
Vibrio_W_Fdataphy2 <- prune_samples(sample_sums(Vibrio_W_Fdataphy)>0, Vibrio_W_Fdataphy)
nsamples(Vibrio_W_Fdataphy2)#124

# agglomerate taxa
glom <- tax_glom(Vibrio_W_Fdataphy2, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 0.01,]$Species
# change their name to "Other Vibrio"
dat[dat$Species %in% Other,]$Species <- 'Other Vibrio'
#remove all Specieses labeled Other Vibrio
dat <-dat[!dat$Species=="Other Vibrio",]
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Vibrio"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Vibrio_Species_W_DPH_dat <- rbind(dat, Abundance)

#Summarize the Vibrio Species Table by DPH in wide format
Wide_Vibrio_Species_W_DPH_dat<- spread(Vibrio_Species_W_DPH_dat, Species, Abundance)

#Save the Vibrio species names as a list
Vibrio_species <- colnames(Wide_Vibrio_Species_W_DPH_dat)

#Use this list to create a filtered version of the Full Species list
Vibrio_Full_Species_W_DPH_dat <- Full_Species_W_DPH_dat[Full_Species_W_DPH_dat$Species %in% Vibrio_species,]

#Use this table to create a bar graph of just Vibrio Species
spatial_plot_Vibrio_Full_Species_W_DPH <- ggplot(data=Vibrio_Full_Species_W_DPH_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Species)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Species", values=c("purple","lightblue","orangered","firebrick3","darkseagreen3","seagreen","white"), labels=c("s__uncultured_bacterium", "s__Vibrio_cidicii", "s__Vibrio_fluvialis", "s__Vibrio_fortis", "s__Vibrio_furnissii", "s__Vibrio_mediterranei", "s__Vibrio_vulnificus"))+
  scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D24"))+
  xlab("Days Post Hatch") +
  ggtitle("Water Vibrio Species by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Vibrio_Full_Species_W_DPH

tiff('Water Vibrio Species by DPH.tiff', units="in", width=6, height=6, res=300)
spatial_plot_Vibrio_Full_Species_W_DPH
dev.off()

##Tissue Pseudomonas

##Create a table of the Pseudomonas Species by DPH
Pseudomonas_T_Fdataphy <- transform_sample_counts(Pseudomonas_T_Fdata, function(x) 100*x/sum(x))
nsamples(Pseudomonas_T_Fdataphy)
Pseudomonas_T_Fdataphy2 <- prune_samples(sample_sums(Pseudomonas_T_Fdataphy)>0, Pseudomonas_T_Fdataphy)
nsamples(Pseudomonas_T_Fdataphy2)

# agglomerate taxa
glom <- tax_glom(Pseudomonas_T_Fdataphy2, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 0.01,]$Species
# change their name to "Other Pseudomonas"
dat[dat$Species %in% Other,]$Species <- 'Other Pseudomonas'
#remove all Specieses labeled Other Pseudomonas
dat <-dat[!dat$Species=="Other Pseudomonas",]
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Pseudomonas"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Pseudomonas_Species_T_DPH_dat <- rbind(dat, Abundance)

##Create a table of the All Species by DPH

T_Fdataphy <- transform_sample_counts(T_Fdata, function(x) 100*x/sum(x))
##Create table ready for making stacked bar graph for Specieses <0.5%##
# agglomerate taxa
glom <- tax_glom(T_Fdataphy, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Full_Species_T_DPH_dat <- rbind(dat)

#Summarize the Pseudomonas Species Table by DPH in wide format
Wide_Pseudomonas_Species_T_DPH_dat<- spread(Pseudomonas_Species_T_DPH_dat, Species, Abundance)

#Save the Pseudomonas species names as a list
pseudomonas_species <- colnames(Wide_Pseudomonas_Species_T_DPH_dat)

#Use this list to create a filtered version of the Full Species list
Psudomonas_Full_Species_T_DPH_dat <- Full_Species_T_DPH_dat[Full_Species_T_DPH_dat$Species %in% pseudomonas_species,]

#Use this table to create a bar graph of just Pseudomonas Species
spatial_plot_Psudomonas_Full_Species_T_DPH <- ggplot(data=Psudomonas_Full_Species_T_DPH_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Species)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus:Species", values=c("deeppink","green","maroon","black","midnightblue","rosybrown","yellow","olivedrab","orchid","khaki"), labels=c("s__Pseudomonas_alcaligenes", "s__Pseudomonas_balearica", "s__Pseudomonas_indica", "s__Pseudomonas_matsuisoli", "s__Pseudomonas_oleovorans", "s__Pseudomonas_otitidis", "s__Pseudomonas_pachastrellae", "s__Pseudomonas_peli", "s__Pseudomonas_psychrotolerans", "s__Spumella-like_flagellate"))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
  xlab("Days Post Hatch") +
  ggtitle("Tissue Pseudomonas Species by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Psudomonas_Full_Species_T_DPH

tiff('Tissue Pseudomonas Species by DPH.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Psudomonas_Full_Species_T_DPH
dev.off()

##Tissue Vibrio

##Create a table of the Vibrio Species by DPH
Vibrio_T_Fdataphy <- transform_sample_counts(Vibrio_T_Fdata, function(x) 100*x/sum(x))
nsamples(Vibrio_T_Fdataphy)
Vibrio_T_Fdataphy2 <- prune_samples(sample_sums(Vibrio_T_Fdataphy)>0, Vibrio_T_Fdataphy)
nsamples(Vibrio_T_Fdataphy2)#124

##Create table ready for making stacked bar graph for Specieses <0.5%##
# agglomerate taxa
glom <- tax_glom(Vibrio_T_Fdataphy2, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Genus and Species columns#
dat<- unite(dat, Species, Genus:Species, sep='; ')
#Combine Family and Genus/Species columns#
dat<- unite(dat, Species, Family:Species, sep='; ')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 0.01,]$Species
# change their name to "Other Vibrio"
dat[dat$Species %in% Other,]$Species <- 'Other Vibrio'
#remove all Specieses labeled Other Vibrio
dat <-dat[!dat$Species=="Other Vibrio",]
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Days_Post_Hatch","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Days_Post_Hatch, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Vibrio"
#remove unnessary columns
dat <- subset(dat, select=c(Days_Post_Hatch, Abundance, Species))
#combine with original table
Vibrio_Species_T_DPH_dat <- rbind(dat, Abundance)

#Summarize the Pseudomonas Species Table by DPH in wide format
Wide_Vibrio_Species_T_DPH_dat<- spread(Vibrio_Species_T_DPH_dat, Species, Abundance)

#Save the Vibrio species names as a list
Vibrio_species <- colnames(Wide_Vibrio_Species_T_DPH_dat)

#Use this list to create a filtered version of the Full Species list
Vibrio_Full_Species_T_DPH_dat <- Full_Species_T_DPH_dat[Full_Species_T_DPH_dat$Species %in% Vibrio_species,]

#Use this table to create a bar graph of just Vibiro Species
spatial_plot_Vibrio_Full_Species_T_DPH <- ggplot(data=Vibrio_Full_Species_T_DPH_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Species)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Species", values=c("lightblue","orangered","firebrick3","darkseagreen3","seagreen","turquoise","white"), labels=c("s__Vibrio_cidicii", "s__Vibrio_fluvialis", "s__Vibrio_fortis", "s__Vibrio_furnissii", "s__Vibrio_mediterranei", "s__Vibrio_ponticus", "s__Vibrio_vulnificus"))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
  xlab("Days Post Hatch") +
  ggtitle("Tissue Vibrio Species by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Vibrio_Full_Species_T_DPH

tiff('Tissue Vibrio Species by DPH.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Vibrio_Full_Species_T_DPH
dev.off()