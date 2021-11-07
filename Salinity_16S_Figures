#MICROBIOME FIGURE 1
PCoA_Sample_Type_no_stats <- plot_ordination(ZDS_Fdata, ordinate(ZDS_Fdata, "PCoA", "bray"), color = "Sample_Type")+ 
  scale_color_manual(name="Sample Type", values = c("gray", "lightblue"))+
  # ggtitle("Overall PCoA by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))

PCoA_Sample_Type_no_stats

tiff('Figure 1 - Overall PCoA by Sample Type wo stats.tiff', units="in", width=10, height=6, res=300)
PCoA_Sample_Type_no_stats
dev.off()

#MICROBIOME FIGURE 2
spatial_plot_Genus_noD01_Fdata_Sample_Type <- ggplot(data=Genus_noD01_Fdata_Sample_Type, aes(x=Sample_Type, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("red","midnightblue","orange","maroon","purple","lightblue","firebrick3","darkseagreen3","seagreen","turquoise","goldenrod","wheat","black","rosybrown"))+
  #ggtitle("Genera by Sample Type") +
  xlab("Sample Type") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  annotate("text", x = c(1:2 ) , y = 105, label = c("a", "b"), size=5)

spatial_plot_Genus_noD01_Fdata_Sample_Type

tiff('Figure 2 - Genus by Sample Type.tiff', units="in", width=6, height=6, res=300)
spatial_plot_Genus_noD01_Fdata_Sample_Type
dev.off()

#MICROBIOME FIGURE 3
W_Shannon_DPH_no_stats_boxplot <- ggplot(W_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon)) + 
  geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  #scale_color_manual("Days Post Hatch", values=c("deeppink", "midnightblue", "blue", "red", "orange", "purple"), limits=c("D1", "D6", "D12","D18", "D24"))+
  # ggtitle("Water Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24")) +
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

tiff('Figure 3 - Water Shannon Diversity by Days Post Hatch wo stats.tiff', units="in", width=10, height=6, res=300)
W_Shannon_DPH_no_stats_boxplot
dev.off()

jpeg('Figure 3 - Water Shannon Diversity by Days Post Hatch wo stats.jpeg', units="in", width=10, height=6, res=300)
W_Shannon_DPH_no_stats_boxplot
dev.off()

#MICROBIOME FIGURE 4
W_PCoA_DPH_no_stats_Salinity <- plot_ordination(W_noD1_ZDS_Fdata, ordinate(W_noD1_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D6", "D12", "D18", "D24"), values=c("red", "forestgreen", "purple", "orange"), labels=c("6", "12","18", "24")) +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  # ggtitle("Water PCoA by DPH & Salinity") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
W_PCoA_DPH_no_stats_Salinity

tiff('Figure 4 - Water PCoA by DPH & Salinity wo stats.tiff', units="in", width=10, height=6, res=300)
W_PCoA_DPH_no_stats_Salinity
dev.off()

jpeg('Figure 4 - Water PCoA by DPH & Salinity wo stats.jpeg', units="in", width=10, height=6, res=300)
W_PCoA_DPH_no_stats_Salinity
dev.off()

#MICROBIOME FIGURE 5
#Water by Salinity
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title <- ggplot(data=Genus_W_noD1_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("yellow","red","midnightblue","orange","green","maroon","purple","lightblue","orangered","firebrick3","darkseagreen3","seagreen","white","turquoise","wheat","goldenrod","black","rosybrown"))+
  #ggtitle("Water Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "c"), size=5)
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title

tiff('Figure 5 - Water Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title
dev.off()

#MICROBIOME FIGURE 6a
W_PPG_Salinity_BP_no_stats <- ggplot(W_noD1_PPG_dat, aes(x=Salinity, y=sum_PPG, color=Salinity)) +   geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  #ggtitle("Water PPG% by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(73, 38, 8), label = c("a", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
W_PPG_Salinity_BP_no_stats

tiff('Figure 6a - Water Potentially Pathogenic Genera by Salinity BP no stats.tiff', units="in", width=6, height=6, res=300)
W_PPG_Salinity_BP_no_stats
dev.off()

jpeg('Figure 6a - Water Potentially Pathogenic Genera by Salinity BP no stats.jpeg', units="in", width=10, height=6, res=300)
W_PPG_Salinity_BP_no_stats
dev.off()

#MICROBIOME FIGURE 6b
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_W_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","chocolate","lightblue","olivedrab","orchid","black"
  ))+
  scale_x_discrete(labels=c("10","20", "30"))+
  xlab("Salinity (ppt)") +
  #ggtitle("Water Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
  
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title

tiff('Figure 6b - Water Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title
dev.off()

#MICROBIOME FIGURE 6c
W_PPG_DPH_BP_no_stats <- ggplot(W_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG, color=Days_Post_Hatch)) +   geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Days Post Hatch", values=c("black", "red", "forestgreen", "purple", "orange"), limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24"))+
  #ggtitle("Water PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24")) +
  annotate("text", x = c(2, 3, 4, 5) , y = c(73, 6, 9, 36), label = c("a", "b", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))

W_PPG_DPH_BP_no_stats

tiff('Figure 6c - Water Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
W_PPG_DPH_BP_no_stats
dev.off()

jpeg('Figure 6c - Water Potentially Pathogenic Genera by DPH BP.jpeg', units="in", width=10, height=6, res=300)
W_PPG_DPH_BP_no_stats
dev.off()

#MICROBIOME FIGURE 6d
spatial_plot_PPFG_Full_Genus_W_DPH_no_title <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","chocolate","lightblue","olivedrab","orchid","black"
  ))+
  scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D24"), labels=c("1", "6", "12","18", "24"))+
  xlab("Days Post Hatch") +
  #ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_W_DPH_no_title

tiff('Figure 6d - Water Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=7, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_no_title
dev.off()

#MICROBIOME FIGURE 7
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c(
    "D1_S30", 
    "D6_S10", "D6_S20", "D6_S30",
    "D12_S10", "D12_S20", "D12_S30",
    "D18_S10", "D18_S20", "D18_S30",
    "D24_S10", "D24_S20", "D24_S30"), labels=c(
      "1 dph; 30 ppt", 
      "6 dph; 10 ppt", "6 dph; 20 ppt", "6 dph; 30 ppt",
      "12 dph; 10 ppt", "12 dph; 20 ppt", "12 dph; 30 ppt",
      "18 dph; 10 ppt", "18 dph; 20 ppt", "18 dph; 30 ppt",
      "24 dph; 10 ppt", "24 dph; 20 ppt", "24 dph; 30 ppt"))+
  xlab("Days Post Hatch and Salinity") +
  #ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title

tiff('Figure 7 - Water Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title
dev.off()

jpeg('Figure 7 - Water Potentially Pathogenic Genera by Days Post Hatch and Salinity.jpeg', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title
dev.off()

#MICROBIOME FIGURE 8
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

tiff('Figure 8 - Tissue Shannon Diversity by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
T_Shannon_DPH_boxplot_no_stats_title
dev.off()

jpeg('Figure 8 - Tissue Shannon Diversity by Days Post Hatch.jpeg', units="in", width=10, height=6, res=300)
T_Shannon_DPH_boxplot_no_stats_title
dev.off()

#MICROBIOME FIGURE 9
T_PCoA_DPH_no_stats_Salinity <- plot_ordination(T_noD01_ZDS_Fdata, ordinate(T_noD01_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("rosybrown", "red",  "black", "forestgreen", "brown", "purple", "dark gray", "orange"), labels=c("3", "6", "9", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PCoA by DPH and Salinity") +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 
T_PCoA_DPH_no_stats_Salinity

tiff('Figure 9 - Tissue PCoA by DPH and Salinity wo D0-1.tiff', units="in", width=10, height=6, res=300)
T_PCoA_DPH_no_stats_Salinity
dev.off()

jpeg('Figure 9 - Tissue PCoA by DPH and Salinity wo D0-1h.jpeg', units="in", width=10, height=6, res=300)
T_PCoA_DPH_no_stats_Salinity
dev.off()

#MICROBIOME FIGURE 10
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title <- ggplot(data=Genus_T_noD01_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","purple","lightblue","firebrick3","darkseagreen3","seagreen","wheat","goldenrod","lightgray","black","darkgray","rosybrown"))+
  #ggtitle("Tissue Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "b"), size=5)
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title


tiff('Figure 10 - Tissue Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title
dev.off()

#MICROBIOME FIGURE 11A
T_PPG_Salinity_BP_no_stats <- ggplot(T_noD01_PPG_dat, aes(x=Salinity, y=sum_PPG, color=Salinity)) +   geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels = c("10", "20", "30"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  #ggtitle("Tissue PPG% by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(97, 97, 99), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
T_PPG_Salinity_BP_no_stats

tiff('Figure 11A - Tissue Potentially Pathogenic Genera by Salinity BP.tiff', units="in", width=6, height=6, res=300)
T_PPG_Salinity_BP_no_stats
dev.off()

jpeg('Figure 11A - Tissue Potentially Pathogenic Genera by Salinity BPh.jpeg', units="in", width=10, height=6, res=300)
T_PPG_Salinity_BP_no_stats
dev.off()

#MICROBIOME FIGURE 11B
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_T_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  xlab("Salinity(ppt)") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title

tiff('Figure 11B - Tissue Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title
dev.off()

#MICROBIOME FIGURE 11C
T_PPG_DPH_BP_no_stats <- ggplot(T_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG, color=Days_Post_Hatch)) +   geom_boxplot() +
 # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("wheat", "black", "rosybrown", "red",  "blue", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24")) +
  annotate("text", x = c(3, 4, 5, 6, 7, 8, 9, 10) , y = c(38, 47, 18, 79, 99, 87, 92, 32), label = c("bc", "abc", "c", "abc", "d", "abd", "ad", "c"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
T_PPG_DPH_BP_no_stats

tiff('Figure 11C - Tissue Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
T_PPG_DPH_BP_no_stats
dev.off()

jpeg('Figure 11C - Tissue Potentially Pathogenic Genera by DPH BP.jpeg', units="in", width=10, height=6, res=300)
T_PPG_DPH_BP_no_stats
dev.off()

#MICROBIOME FIGURE 11D
spatial_plot_PPFG_Full_Genus_T_DPH_no_title <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"),labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  xlab("Days Post Hatch") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_DPH_no_title 

tiff('Figure 11D - Tissue Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_no_title
dev.off()

#MICROBIOME FIGURE 12
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
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
    "D24_S10", "D24_S20", "D24_S30"),
    labels=c("0 dph; 30 ppt", 
      "1 dph; 30 ppt", 
      "3 dph; 10 ppt", "3 dph; 20 ppt", "3 dph; 30 ppt",
      "6 dph; 10 ppt", "6 dph; 20 ppt", "6 dph; 30 ppt",
      "9 dph; 10 ppt", "9 dph; 20 ppt", "9 dph; 30 ppt",
      "12 dph; 10 ppt", "12 dph; 20 ppt", "12 dph; 30 ppt",
      "15 dph; 10 ppt", "15 dph; 20 ppt", "15 dph; 30 ppt",
      "18 dph; 10 ppt", "18 dph; 20 ppt", "18 dph; 30 ppt",
      "20 dph; 10 ppt", "20 dph; 20 ppt", "20 dph; 30 ppt",
      "24 dph; 10 ppt", "24 dph; 20 ppt", "24 dph; 30 ppt"))+
  xlab("Days Post Hatch and Salinity") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title

tiff('Figure 12 - Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title
dev.off()

jpeg('Figure 12 - Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity.jpeg', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title
dev.off()











spatial_plot_PPFG_Full_Genus_T_DPH_no_title_feeding <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue", "chocolate","khaki","lightblue","olivedrab","orchid","black"))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"),labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  xlab("Days Post Hatch") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("segment", x=1, xend=6.0, y=100, yend=100, size=2, color = "forestgreen")+ #0-12 dph
  annotate("text", x=3.5, y=105, label="Green Water", size=5, color="forestgreen")+
  annotate("segment", x=2, xend=5.5, y=110, yend=110, size=2, color = "red") + #1-10 dph
  annotate("text", x=3.75, y=115, label="Enriched Rotifers", size=5, color ="red") +  
  annotate("segment", x=5, xend=7, y=120, yend=120, size=2, color = "blue") + #9-15 dph
  annotate("text", x=6, y=125, label="Artemia nauplii", size=5, color ="blue")+
  annotate("segment", x=6, xend=9.5, y=130, yend=130, size=2, color = "orange")+ #12-21 dph
  annotate("text", x=7.75, y=135, label="Enriched Artemia", size=5, color ="orange")+
  annotate("segment", x=6.25, xend=10, y=140, yend=140, size=2, color = "purple")+ #13-24 dph
  annotate("text", x=8.125, y=145, label="Microfeeds", size=5, color ="purple")

spatial_plot_PPFG_Full_Genus_T_DPH_no_title_feeding


tiff('Tissue Potentially Pathogenic Genera by Days Post Hatch with Feeding.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_no_title_feeding
dev.off()
