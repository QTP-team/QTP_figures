library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
#install.packages('ggforce')
library(ggforce)
library(ggpubr)
library(stringr)

All_SGBs_BGCs <- read.table("0.data/BGC.Network_Annotations_Full.QTP.tsv", sep = "\t", header = T)
#All_SGBs_BGCs <- All_SGBs_BGCs[,c(3,5)]
All_SGBs_BGCs$SGBs <- str_split_fixed(All_SGBs_BGCs$Description, "_k141", 2)[,1]
#colnames(All_SGBs_BGCs) <- c("Contigs", "BGCs", "SGBs")
#write.csv(SGBs_profile_sample, file = "SGBs_profile_sample.csv")
All_SGBs_taxonomy <- read.csv("0.data/QTP_SGBs_14062_GTDB_r214.csv")
All_SGBs_taxonomy <- separate(data = All_SGBs_taxonomy, col = taxonomy, 
                     into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
All_SGBs_taxonomy <- All_SGBs_taxonomy[,c(1,3)]


### all
# All 饼图
all_SGB_BGC_taxo <- merge(All_SGBs_BGCs, All_SGBs_taxonomy, by.x = "SGBs", by.y = "SGB_ID")
write.csv(all_SGB_BGC_taxo, file = "all_SGB_BGC_taxo.csv")
all_BGC_class <- as.data.frame(table(all_SGB_BGC_taxo$BiG.SCAPE.class))
colnames(all_BGC_class) <- c("BGC_class", "all")
all_BGC_class <- all_BGC_class[order(all_BGC_class$all),]
p_all <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("") + ylab('') +
  scale_fill_manual(values = c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                               "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                               "Saccharides" = "#b3de69", "Terpene" = "#fccde5"))+
  geom_arc_bar(data=all_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = all, fill = BGC_class))+
  theme(legend.position = "bottom")
# four phylum 柱状图
all_BGC_data_annotation_4phylum <- all_SGB_BGC_taxo[which((all_SGB_BGC_taxo$Phylum=="p__Bacillota_A")|(all_SGB_BGC_taxo$Phylum=="p__Bacteroidota")|(all_SGB_BGC_taxo$Phylum=="p__Verrucomicrobiota")|(all_SGB_BGC_taxo$Phylum=="p__Spirochaetota")),]
all_BGC_data_annotation_4phylum$Phylum <- factor(all_BGC_data_annotation_4phylum$Phylum, levels = c("p__Bacillota_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
all_BGC_data_annotation_4phylum$BiG.SCAPE.class <- factor(all_BGC_data_annotation_4phylum$BiG.SCAPE.class, levels = c( "Saccharides", "PKS-NRP_Hybrids","PKSI","PKSother","Terpene","NRPS", "Others","RiPPs"))
p_all_BGC_class <- ggplot(all_BGC_data_annotation_4phylum, aes(x = Phylum, fill = BiG.SCAPE.class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "",
       y = "BGC superfamily percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_all_BGC <- ggarrange(p_all, p_all_BGC_class, widths = c(1,1), nrow = 1)
ggsave(p_all_BGC, filename = "all_BGC.pdf", width = 10, height = 6)

