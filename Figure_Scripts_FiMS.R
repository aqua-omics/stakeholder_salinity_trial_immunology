#SUPPLEMENTARY FIGURE 2
PCoA_Sample_Type_no_stats <- plot_ordination(ZDS_Fdata, ordinate(ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Sample_Type")+ 
  scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("purple", "deeppink",  "gold", "blue",  "black", "forestgreen", "brown", "lightblue", "darkgray", "orange"), labels=c("0", "1", "3", "6", "9", "12", "15", "18", "20", "24"))+
  # ggtitle("Overall PCoA by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  geom_point(size=6)+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))

PCoA_Sample_Type_no_stats

tiff('Supplementary Figure 2 - Sample Type PCoA.tiff', units="in", width=10, height=6, res=300)
PCoA_Sample_Type_no_stats
dev.off()

#FIGURE 1
spatial_plot_Genus_noD01_Fdata_Sample_Type <- ggplot(data=Genus_noD01_Fdata_Sample_Type, aes(x=Sample_Type, y=Abundance, fill=Genus)) +
  theme_bw()+
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("red","midnightblue","orange","maroon","purple","lightblue","blue","darkseagreen3","seagreen","turquoise","goldenrod","wheat","black","rosybrown"))+
  scale_x_discrete(labels=c("Larvae", "Water"))+
  #ggtitle("Genera by Sample Type") +
  xlab("Sample Type") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))  + 
  annotate("text", x = c(1:2 ) , y = 105, label = c("a", "b"), size=5)+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_Genus_noD01_Fdata_Sample_Type

tiff('Figure 1 - Genus by Sample Type.tiff', units="in", width=6, height=6, res=300)
spatial_plot_Genus_noD01_Fdata_Sample_Type
dev.off()

#SUPPLEMENTARY FIGURE 3
W_Genera_Salinity_Venn <- ggvenn(W_Salinities_Genera_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  #ggtitle("Water Genera by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))
W_Genera_Salinity_Venn

tiff('Supplementary Figure 3 - Water Genera by Salinity Venn Diagram.tiff', units="in", width=6, height=6, res=300)
W_Genera_Salinity_Venn
dev.off()


#FIGURE 2
W_Shannon_DPH_no_stats_boxplot <- ggplot(W_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon)) + 
  geom_boxplot(fill="gray90") +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Days Post Hatch", values=c("deeppink", "midnightblue", "blue", "red", "orange", "purple"), limits=c("D1", "D6", "D12","D18", "D24"))+
  # ggtitle("Water Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24")) +
  annotate("text", x = c(2.0, 3, 4, 5) , y = c(3.1, 2.9, 4.0, 3.5), label = c("a", "b", "c", "c"), size=5) + 
  annotate("segment", x=0.5, xend=3.0, y=4.25, yend=4.25, size=2)+ #0-12 DPH
  annotate("text", x=1.75, y=4.375, label="Green Water", size=5)+
  annotate("segment", x=1, xend=2.75, y=4.5, yend=4.5, size=2)+ #1-10 DPH
  annotate("text", x=1.875, y=4.625, label="Enriched Rotifers", size=5) +  
  annotate("segment", x=2.5, xend=3.5, y=4.75, yend=4.75, size=2)+ #9-15 DPH
  annotate("text", x=3.0, y=4.875, label="Artemia nauplii", size=5)+
  annotate("segment", x=3.0, xend=4.5, y=5, yend=5, size=2)+ #12-21 DPH
  annotate("text", x=3.75, y=5.125, label="Enriched Artemia", size=5)+
  annotate("segment", x=3.25, xend=5, y=5.25, yend=5.25, size=2)+ # 13-24 DPH
  annotate("text", x=4.125, y=5.375, label="Microfeeds", size=5)+  theme_bw() +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
W_Shannon_DPH_no_stats_boxplot

tiff('Figure 2 - Water Shannon Diversity by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
W_Shannon_DPH_no_stats_boxplot
dev.off()

jpeg('Figure 2 - Water Shannon Diversity by Days Post Hatch.jpeg', units="in", width=10, height=6, res=300)
W_Shannon_DPH_no_stats_boxplot
dev.off()

#MICROBIOME FIGURE 3
W_PCoA_DPH_no_stats_Salinity <- plot_ordination(W_noD1_ZDS_Fdata, ordinate(W_noD1_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D6", "D12", "D18", "D24"), values=c("blue", "forestgreen", "lightblue", "orange"), labels=c("6", "12","18", "24")) +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  # ggtitle("Water PCoA by DPH & Salinity") +
theme_bw() +
  geom_point(size=6)+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
W_PCoA_DPH_no_stats_Salinity

tiff('Figure 3 - Water PCoA by DPH & Salinity.tiff', units="in", width=10, height=6, res=300)
W_PCoA_DPH_no_stats_Salinity
dev.off()

jpeg('Figure 3 - Water PCoA by DPH & Salinity.jpeg', units="in", width=10, height=6, res=300)
W_PCoA_DPH_no_stats_Salinity
dev.off()

#SUPPLEMENTARY FIGURE 4
#Water by Salinity
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title <- ggplot(data=Genus_W_noD1_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  scale_fill_manual("Family;Genus", values=c("yellow","red","midnightblue","orange","green","maroon","purple","lightblue","orangered","blue","darkseagreen3","seagreen","white","turquoise","wheat","goldenrod","black","rosybrown"))+
  #ggtitle("Water Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
theme_bw()+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "c"), size=5)+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title

tiff('Supplementary Figure 4 - Water Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title
dev.off()

jpeg('Supplementary Figure 4 - Water Genera by Salinity.jpeg', units="in", width=7, height=6, res=300)
spatial_plot_Genus_NoD1_W_Fdata_Salinity_no_title
dev.off()

#FIGURE 4

#FIGURE 4a
W_PPG_Salinity_BP_no_stats <- ggplot(W_noD1_PPG_dat, aes(x=Salinity, y=sum_PPG)) +   geom_boxplot(fill="gray90") +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  #scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels=c("10", "20", "30"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  #ggtitle("Water PPG% by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Potentially Pathogenic Genera Percentage") +
theme_bw()+
  annotate("text", x = c(1, 2, 3) , y = c(73, 38, 8), label = c("a", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
W_PPG_Salinity_BP_no_stats

tiff('Figure 4a - Water Potentially Pathogenic Genera by Salinity BP no stats.tiff', units="in", width=6, height=6, res=300)
W_PPG_Salinity_BP_no_stats
dev.off()

jpeg('Figure 4a - Water Potentially Pathogenic Genera by Salinity BP no stats.jpeg', units="in", width=10, height=6, res=300)
W_PPG_Salinity_BP_no_stats
dev.off()

#FIGURE 4b
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_W_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","darkseagreen3","plum","darkorange","black"))+
  scale_x_discrete(labels=c("10","20", "30"))+
  xlab("Salinity (ppt)") +
  #ggtitle("Water Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme_bw() +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
  
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title

tiff('Figure 4b - Water Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title
dev.off()

jpeg('Figure 4b - Water Potentially Pathogenic Genera by Salinity.jpeg', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title
dev.off()

#FIGURE 4c
W_PPG_DPH_BP_no_stats <- ggplot(W_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG)) +   geom_boxplot(fill="gray90") +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  #scale_color_manual("Days Post Hatch", values=c("black", "red", "forestgreen", "purple", "orange"), limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24"))+
  #ggtitle("Water PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme_bw() +
  scale_x_discrete(limits=c("D1", "D6", "D12","D18", "D24"), labels=c("1", "6", "12","18", "24")) +
  annotate("text", x = c(2, 3, 4, 5) , y = c(73, 6, 9, 36), label = c("a", "b", "b", "b"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())

W_PPG_DPH_BP_no_stats

tiff('Figure 4c - Water Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
W_PPG_DPH_BP_no_stats
dev.off()

jpeg('Figure 4c - Water Potentially Pathogenic Genera by DPH BP.jpeg', units="in", width=10, height=6, res=300)
W_PPG_DPH_BP_no_stats
dev.off()

#FIGURE 4d
spatial_plot_PPFG_Full_Genus_W_DPH_no_title <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","darkseagreen3","plum","darkorange","black"  ))+
  scale_x_discrete(limits=c("D1","D6", "D12", "D18", "D24"), labels=c("1", "6", "12","18", "24"))+
  xlab("Days Post Hatch") +
  #ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme_bw() +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_PPFG_Full_Genus_W_DPH_no_title

tiff('Figure 4d - Water Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=7, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_no_title
dev.off()

#FIGURE 4
#Make a multi-part figure, removing the axis labels from each separate figure
PPFG_W_Figures <- ggarrange(W_PPG_Salinity_BP_no_stats, spatial_plot_PPFG_Full_Genus_W_noD1_Salinity_no_title, W_PPG_DPH_BP_no_stats, spatial_plot_PPFG_Full_Genus_W_DPH_no_title,
                                         labels = c("A", "B", "C", "D"),
                                         ncol = 2, nrow = 2,
                                         common.legend = FALSE,
                            widths=c(1,2))
PPFG_W_Figures

tiff('Figure 4 - Potentially Pathogenic Genera Water Figures.tiff', units="in", width=12, height=12, res=300)
PPFG_W_Figures
dev.off()

jpeg('Figure 4 - Potentially Pathogenic Genera Water Figures.jpeg', units="in", width=12, height=12, res=300)
PPFG_W_Figures
dev.off()


#SUPPLEMENTARY FIGURE 5
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_W_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  theme_bw() + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","darkseagreen3","plum","darkorange","black"))+
  scale_x_discrete(limits=c(
    "D1_S30", 
    "D6_S10", "D6_S20", "D6_S30",
    "D12_S10", "D12_S20", "D12_S30",
    "D18_S10", "D18_S20", "D18_S30",
    "D24_S10", "D24_S20", "D24_S30"), labels=c(
      "1 DPH; 30 ppt", 
      "6 DPH; 10 ppt", "6 DPH; 20 ppt", "6 DPH; 30 ppt",
      "12 DPH; 10 ppt", "12 DPH; 20 ppt", "12 DPH; 30 ppt",
      "18 DPH; 10 ppt", "18 DPH; 20 ppt", "18 DPH; 30 ppt",
      "24 DPH; 10 ppt", "24 DPH; 20 ppt", "24 DPH; 30 ppt"))+
  xlab("Days Post Hatch and Salinity") +
  #ggtitle("Water Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title

tiff('Supplementary Figure 5 - Water Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title
dev.off()

jpeg('Supplementary Figure 5 - Water Potentially Pathogenic Genera by Days Post Hatch and Salinity.jpeg', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_W_DPH_Salinity_no_title
dev.off()

#SUPPLEMENTARY FIGURE 6
T_Genera_Salinity_Venn <- ggvenn(T_Salinities_Genera_Lists, fill_color = c("deeppink", "midnightblue", "blue"), show_percentage = FALSE)+
  #ggtitle("Tissue Genera by Salinity") +
  theme(plot.title = element_text(hjust = 0.5))
T_Genera_Salinity_Venn

tiff('Supplementary Figure 6 - Larvae Genera by Salinity Venn Diagram.tiff', units="in", width=6, height=6, res=300)
T_Genera_Salinity_Venn
dev.off()


#FIGURE 5
T_Shannon_DPH_boxplot_no_stats_title <- ggplot(T_Falphadiv_metadata, aes(x=Days_Post_Hatch, y=Shannon)) + 
  geom_boxplot(fill="gray90") +
  theme_bw() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1))+ 
  #scale_color_manual("Days Post Hatch", values=c("lightgray", "deeppink", "lightblue", "midnightblue", "maroon", "blue", "rosybrown", "red", "orange", "darkgray", "purple"), limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"))+
#ggtitle("Tissue Shannon Diversity by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Shannon") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1","D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  annotate("text", x = c(3, 4, 5, 6, 7, 8, 9, 10) , y = c(5.3, 5.3, 5.35, 4.65, 3.3, 4.5, 3.55, 4.8), label = c("cd", "d", "d", "abc", "a", "ab", "ab", "bcd"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  annotate("segment", x=1, xend=6.0, y=5.5, yend=5.5, size=2)+ #0-12 DPH
  annotate("text", x=3.5, y=5.625, label="Green Water", size=5)+
  annotate("segment", x=2, xend=5.5, y=5.75, yend=5.75, size=2) + #1-10 DPH
  annotate("text", x=3.75, y=5.875, label="Enriched Rotifers", size=5) +  
  annotate("segment", x=5, xend=7, y=6, yend=6, size=2) + #9-15 DPH
  annotate("text", x=6, y=6.125, label="Artemia nauplii", size=5)+
  annotate("segment", x=6, xend=9.5, y=6.25, yend=6.25, size=2)+ #12-21 DPH
  annotate("text", x=7.75, y=6.375, label="Enriched Artemia", size=5)+
  annotate("segment", x=6.25, xend=10, y=6.5, yend=6.5, size=2)+ #13-24 DPH
  annotate("text", x=8.125, y=6.625, label="Microfeeds", size=5)+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
T_Shannon_DPH_boxplot_no_stats_title

tiff('Figure 5 - Tissue Shannon Diversity by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
T_Shannon_DPH_boxplot_no_stats_title
dev.off()

jpeg('Figure 5 - Tissue Shannon Diversity by Days Post Hatch.jpeg', units="in", width=10, height=6, res=300)
T_Shannon_DPH_boxplot_no_stats_title
dev.off()

#FIGURE 6
T_PCoA_DPH_no_stats_Salinity <- plot_ordination(T_noD01_ZDS_Fdata, ordinate(T_noD01_ZDS_Fdata, "PCoA", "bray"), color = "Days_Post_Hatch", shape = "Salinity") + 
  scale_color_manual(name="Days Post Hatch", limits=c("D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("gold", "blue",  "black", "forestgreen", "brown", "lightblue", "darkgray", "orange"), labels=c("3", "6", "9", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PCoA by DPH and Salinity") +
  scale_shape_discrete(name="Salinity (ppt)", labels=c("10", "20", "30"))+
  theme_bw()+
  geom_point(size=6)+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13)) 
T_PCoA_DPH_no_stats_Salinity

tiff('Figure 6 - Larvae PCoA by DPH and Salinity wo D0-1.tiff', units="in", width=10, height=6, res=300)
T_PCoA_DPH_no_stats_Salinity
dev.off()

jpeg('Figure 6 - Larvae PCoA by DPH and Salinity wo D0-1h.jpeg', units="in", width=10, height=6, res=300)
T_PCoA_DPH_no_stats_Salinity
dev.off()

#MICROBIOME SUPPLEMENTARY FIGURE 7
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title <- ggplot(data=Genus_T_noD01_Fdata_Salinity, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
theme_bw()+
  scale_fill_manual("Family;Genus", values=c("deeppink","midnightblue","purple","lightblue","blue","darkseagreen3","seagreen","wheat","goldenrod","lightgray","black","darkgray","rosybrown"))+
  #ggtitle("Tissue Genera by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  scale_x_discrete(limits=c("S10", "S20", "S30"), labels=c("10", "20", "30"))+
  annotate("text", x = c(1:3 ) , y = 105, label = c("a", "b", "b"), size=5)+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title


tiff('Supplementary Figure 7 - Tissue Genera by Salinity.tiff', units="in", width=7, height=6, res=300)
spatial_plot_Genus_T_noD01_Fdata_Salinity_no_title
dev.off()

#FIGURE 7

#FIGURE 7A
T_PPG_Salinity_BP_no_stats <- ggplot(T_noD01_PPG_dat, aes(x=Salinity, y=sum_PPG)) +   geom_boxplot(fill="gray90") +
  theme_bw() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_manual("Salinity (ppt)", values=c("deeppink", "midnightblue", "blue"), labels = c("10", "20", "30"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  #ggtitle("Tissue PPG% by Salinity") +
  xlab("Salinity (ppt)") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x = c(1, 2, 3) , y = c(97, 97, 99), label = c("a", "a", "a"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
T_PPG_Salinity_BP_no_stats

tiff('Figure 7A - Larvae Potentially Pathogenic Genera by Salinity BP.tiff', units="in", width=6, height=6, res=300)
T_PPG_Salinity_BP_no_stats
dev.off()

jpeg('Figure 7A - Larvae Potentially Pathogenic Genera by Salinity BPh.jpeg', units="in", width=10, height=6, res=300)
T_PPG_Salinity_BP_no_stats
dev.off()

#FIGURE 7B
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_T_noD1_Salinity_dat, aes(x=Salinity, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme_bw() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","saddlebrown","darkseagreen3","plum","darkorange","black"))+
  scale_x_discrete(labels=c("10", "20", "30"))+
  xlab("Salinity(ppt)") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title

tiff('Figure 7B - Tissue Potentially Pathogenic Genera by Salinity.tiff', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title
dev.off()

jpeg('Figure 7B - Tissue Potentially Pathogenic Genera by Salinity.jpeg', units="in", width=6, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title
dev.off()

#FIGURE 7C
T_PPG_DPH_BP_no_stats <- ggplot(T_PPG_dat, aes(x=Days_Post_Hatch, y=sum_PPG)) +   geom_boxplot(fill="gray90") +
  theme_bw() +
 # theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  #scale_color_manual(name="Days Post Hatch", limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), values=c("wheat", "black", "rosybrown", "red",  "blue", "forestgreen", "brown", "purple", "dark gray", "orange"),labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  #ggtitle("Tissue PPG% by Days Post Hatch") +
  xlab("Days Post Hatch") +
  ylab("Potentially Pathogenic Genera Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"), labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24")) +
  annotate("text", x = c(3, 4, 5, 6, 7, 8, 9, 10) , y = c(38, 47, 18, 79, 99, 87, 92, 32), label = c("bc", "abc", "c", "abc", "d", "abd", "ad", "c"), size=5)+   
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
T_PPG_DPH_BP_no_stats

tiff('Figure 7C - Tissue Potentially Pathogenic Genera by DPH BP.tiff', units="in", width=6, height=6, res=300)
T_PPG_DPH_BP_no_stats
dev.off()

jpeg('Figure 7C - Tissue Potentially Pathogenic Genera by DPH BP.jpeg', units="in", width=10, height=6, res=300)
T_PPG_DPH_BP_no_stats
dev.off()

#FIGURE 7D
spatial_plot_PPFG_Full_Genus_T_DPH_no_title <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_dat, aes(x=Days_Post_Hatch, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme_bw() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","saddlebrown","darkseagreen3","plum","darkorange","black"))+
  scale_x_discrete(limits=c("D0", "D1", "D3", "D6", "D9", "D12", "D15", "D18", "D20", "D24"),labels=c("0", "1","3", "6", "9", "12", "15", "18", "20", "24"))+
  xlab("Days Post Hatch") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_PPFG_Full_Genus_T_DPH_no_title 

tiff('Figure 7D - Tissue Potentially Pathogenic Genera by Days Post Hatch.tiff', units="in", width=10, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_no_title
dev.off()

#FIGURE 7
#Make a multi-part figure, removing the axis labels from each separate figure
PPFG_T_Figures <- ggarrange(T_PPG_Salinity_BP_no_stats, spatial_plot_PPFG_Full_Genus_T_noD1_Salinity_no_title, T_PPG_DPH_BP_no_stats, spatial_plot_PPFG_Full_Genus_T_DPH_no_title,
                            labels = c("A", "B", "C", "D"),
                            ncol = 2, nrow = 2,
                            common.legend = FALSE,
                            widths=c(1,2))
PPFG_T_Figures

tiff('Figure 7 - Potentially Pathogenic Genera Larvae Figures.tiff', units="in", width=12, height=12, res=300)
PPFG_T_Figures
dev.off()

jpeg('Figure 7 - Potentially Pathogenic Genera Larvae Figures.jpeg', units="in", width=12, height=12, res=300)
PPFG_T_Figures
dev.off()


#SUPPLEMENTARY FIGURE 8
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title <- ggplot(data=PPFG_Full_Genus_T_Days_Post_Hatch_Salinity_dat, aes(x=Days_Post_Hatch_Salinity, y=Abundance, fill=Genus)) +
  theme_bw() +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("gold","red","cyan","saddlebrown","darkseagreen3","plum","darkorange","black"))+
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
    labels=c("0 DPH; 30 ppt", 
      "1 DPH; 30 ppt", 
      "3 DPH; 10 ppt", "3 DPH; 20 ppt", "3 DPH; 30 ppt",
      "6 DPH; 10 ppt", "6 DPH; 20 ppt", "6 DPH; 30 ppt",
      "9 DPH; 10 ppt", "9 DPH; 20 ppt", "9 DPH; 30 ppt",
      "12 DPH; 10 ppt", "12 DPH; 20 ppt", "12 DPH; 30 ppt",
      "15 DPH; 10 ppt", "15 DPH; 20 ppt", "15 DPH; 30 ppt",
      "18 DPH; 10 ppt", "18 DPH; 20 ppt", "18 DPH; 30 ppt",
      "20 DPH; 10 ppt", "20 DPH; 20 ppt", "20 DPH; 30 ppt",
      "24 DPH; 10 ppt", "24 DPH; 20 ppt", "24 DPH; 30 ppt"))+
  xlab("Days Post Hatch and Salinity") +
  #ggtitle("Tissue Potentially Pathogenic Genera by Days Post Hatch and Salinity") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title

tiff('Supplementary Figure 8 - Larvae Potentially Pathogenic Genera by Days Post Hatch and Salinity.tiff', units="in", width=14, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title
dev.off()

jpeg('Supplementary Figure 8 - Larvae Potentially Pathogenic Genera by Days Post Hatch and Salinity.jpeg', units="in", width=14, height=6, res=300)
spatial_plot_PPFG_Full_Genus_T_DPH_Salinity_no_title
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

#FIGURE 8
#FIGURE 8A
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

#FIGURE 8B
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

#FIGURE 8C
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

#FIGURE 8D
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

#FIGURE 8E
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

#FIGURE 8F
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

#FIGURE 8

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

#FIGURE 9
#FIGURE 9A
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

#FIGURE 9B
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

#FIGURE 9C
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

#FIGURE 9D
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


#FIGURE 9E
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

#FIGURE 9F
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

#FIGURE 9

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