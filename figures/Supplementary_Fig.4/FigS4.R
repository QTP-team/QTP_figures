# Supplementary Fig. 4. Predicted biosynthetic genes from the SGBs catalogs of the six host species.
library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
#install.packages('ggforce')
library(ggforce)
library(ggpubr)

TAN_BGC_data <- read.csv("0.data/TAN_BGC.csv")
TS_BGC_data <- read.csv("0.data/TS_BGC.csv")
Yak_BGC_data <- read.csv("0.data/Yak_BGC.csv")
TC_BGC_data <- read.csv("0.data/TC_BGC.csv")
TA_BGC_data <- read.csv("0.data/TA_BGC.csv")
TH_BGC_data <- read.csv("0.data/TH_BGC.csv")
All_SGBs_annotation <- read.csv("0.data/TableS3_SGBs_summary.csv")
All_SGBs_annotation <- All_SGBs_annotation[,c(9,20)]

## TAN
# All 饼图
TAN_BGC_class <- as.data.frame(table(TAN_BGC_data$BiG.SCAPE_class))
colnames(TAN_BGC_class) <- c("BGC_class", "TAN")
TAN_BGC_class <- TAN_BGC_class[order(TAN_BGC_class$TAN),]
p_TAN_all <- ggplot()+
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
  geom_arc_bar(data=TAN_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = TAN, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
TAN_BGC_data_new <- separate(TAN_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
TAN_BGC_data_new$BGC_ID <- paste(TAN_BGC_data_new$BGC_ID, '.fa', sep="") 
TAN_BGC_data_annotation <- merge(TAN_BGC_data_new, All_SGBs_annotation, 
                                 by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(TAN_BGC_data_annotation, file = "TAN_BGC_data_annotation.csv")
TAN_BGC_data_annotation_4phylum <- TAN_BGC_data_annotation[which((TAN_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(TAN_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(TAN_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(TAN_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
TAN_BGC_data_annotation_4phylum$GTDB_phylum <- factor(TAN_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
TAN_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(TAN_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c( "Saccharides","PKS-NRP_Hybrids","PKSI", "PKSother","Terpene","NRPS","Others","RiPPs"))
p_TAN_BGC_class <- ggplot(TAN_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "TAN",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_TAN_BGC <- ggarrange(p_TAN_all, p_TAN_BGC_class, widths = c(1,1), nrow = 2)
#ggsave(p_TAN_BGC, filename = "TAN_BGC.pdf", width = 6, height = 8)

## TS
# All 饼图
TS_BGC_class <- as.data.frame(table(TS_BGC_data$BiG.SCAPE_class))
colnames(TS_BGC_class) <- c("BGC_class", "TS")
TS_BGC_class <- TS_BGC_class[order(TS_BGC_class$TS),]
p_TS_all <- ggplot()+
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
  geom_arc_bar(data=TS_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = TS, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
TS_BGC_data_new <- separate(TS_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
TS_BGC_data_new$BGC_ID <- paste(TS_BGC_data_new$BGC_ID, '.fa', sep="") 
TS_BGC_data_annotation <- merge(TS_BGC_data_new, All_SGBs_annotation, 
                                by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(TS_BGC_data_annotation, file = "TS_BGC_data_annotation.csv")
TS_BGC_data_annotation_4phylum <- TS_BGC_data_annotation[which((TS_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(TS_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(TS_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(TS_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
TS_BGC_data_annotation_4phylum$GTDB_phylum <- factor(TS_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
TS_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(TS_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c( "Saccharides","PKS-NRP_Hybrids","PKSI", "PKSother","Terpene","NRPS","Others","RiPPs"))
p_TS_BGC_class <- ggplot(TS_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "TS",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_TS_BGC <- ggarrange(p_TS_all, p_TS_BGC_class, widths = c(1,1), nrow = 2)


## TC
# All 饼图
TC_BGC_class <- as.data.frame(table(TC_BGC_data$BiG.SCAPE_class))
colnames(TC_BGC_class) <- c("BGC_class", "TC")
TC_BGC_class <- TC_BGC_class[order(TC_BGC_class$TC),]
p_TC_all <- ggplot()+
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
  geom_arc_bar(data=TC_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = TC, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
TC_BGC_data_new <- separate(TC_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
TC_BGC_data_new$BGC_ID <- paste(TC_BGC_data_new$BGC_ID, '.fa', sep="") 
TC_BGC_data_annotation <- merge(TC_BGC_data_new, All_SGBs_annotation, 
                                by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(TC_BGC_data_annotation, file = "TC_BGC_data_annotation.csv")
TC_BGC_data_annotation_4phylum <- TC_BGC_data_annotation[which((TC_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(TC_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(TC_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(TC_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
TC_BGC_data_annotation_4phylum$GTDB_phylum <- factor(TC_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
TC_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(TC_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c("PKS-NRP_Hybrids","PKSother","PKSI", "Terpene","NRPS","Others","RiPPs"))
p_TC_BGC_class <- ggplot(TC_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "TC",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_TC_BGC <- ggarrange(p_TC_all, p_TC_BGC_class, widths = c(1,1), nrow = 2)


## TA
# All 饼图
TA_BGC_class <- as.data.frame(table(TA_BGC_data$BiG.SCAPE_class))
colnames(TA_BGC_class) <- c("BGC_class", "TA")
TA_BGC_class <- TA_BGC_class[order(TA_BGC_class$TA),]
p_TA_all <- ggplot()+
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
  geom_arc_bar(data=TA_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = TA, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
TA_BGC_data_new <- separate(TA_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
TA_BGC_data_new$BGC_ID <- paste(TA_BGC_data_new$BGC_ID, '.fa', sep="") 
TA_BGC_data_annotation <- merge(TA_BGC_data_new, All_SGBs_annotation, 
                                by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(TA_BGC_data_annotation, file = "TA_BGC_data_annotation.csv")
TA_BGC_data_annotation_4phylum <- TA_BGC_data_annotation[which((TA_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(TA_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(TA_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(TA_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
TA_BGC_data_annotation_4phylum$GTDB_phylum <- factor(TA_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
TA_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(TA_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c("PKS-NRP_Hybrids","PKSI", "PKSother","Terpene","Others","NRPS","RiPPs"))
p_TA_BGC_class <- ggplot(TA_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "TA",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_TA_BGC <- ggarrange(p_TA_all, p_TA_BGC_class, widths = c(1,1), nrow = 2)

## Yak
# All 饼图
Yak_BGC_class <- as.data.frame(table(Yak_BGC_data$BiG.SCAPE_class))
colnames(Yak_BGC_class) <- c("BGC_class", "Yak")
Yak_BGC_class <- Yak_BGC_class[order(Yak_BGC_class$Yak),]
p_Yak_all <- ggplot()+
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
  geom_arc_bar(data=Yak_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = Yak, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
Yak_BGC_data_new <- separate(Yak_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
Yak_BGC_data_new$BGC_ID <- paste(Yak_BGC_data_new$BGC_ID, '.fa', sep="") 
Yak_BGC_data_annotation <- merge(Yak_BGC_data_new, All_SGBs_annotation, 
                                 by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(Yak_BGC_data_annotation, file = "Yak_BGC_data_annotation.csv")
Yak_BGC_data_annotation_4phylum <- Yak_BGC_data_annotation[which((Yak_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(Yak_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(Yak_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(Yak_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
Yak_BGC_data_annotation_4phylum$GTDB_phylum <- factor(Yak_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
Yak_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(Yak_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c( "Saccharides","PKS-NRP_Hybrids","PKSI", "PKSother","Terpene","NRPS","Others","RiPPs"))
p_Yak_BGC_class <- ggplot(Yak_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "Yak",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_Yak_BGC <- ggarrange(p_Yak_all, p_Yak_BGC_class, widths = c(1,1), nrow = 2)


## TH
# All 饼图
TH_BGC_class <- as.data.frame(table(TH_BGC_data$BiG.SCAPE_class))
colnames(TH_BGC_class) <- c("BGC_class", "TH")
TH_BGC_class <- TH_BGC_class[order(TH_BGC_class$TH),]
p_TH_all <- ggplot()+
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
  geom_arc_bar(data=TH_BGC_class,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount = TH, fill = BGC_class))+
  theme(legend.position = "none")
# four phylum 柱状图
TH_BGC_data_new <- separate(TH_BGC_data, BGC_ID, into = c("BGC_ID","tmp"), sep = ".fa")
TH_BGC_data_new$BGC_ID <- paste(TH_BGC_data_new$BGC_ID, '.fa', sep="") 
TH_BGC_data_annotation <- merge(TH_BGC_data_new, All_SGBs_annotation, 
                                by.x = "BGC_ID", by.y = "Representative_genome")
#write.csv(TH_BGC_data_annotation, file = "TH_BGC_data_annotation.csv")
TH_BGC_data_annotation_4phylum <- TH_BGC_data_annotation[which((TH_BGC_data_annotation$GTDB_phylum=="p__Firmicutes_A")|(TH_BGC_data_annotation$GTDB_phylum=="p__Bacteroidota")|(TH_BGC_data_annotation$GTDB_phylum=="p__Spirochaetota")|(TH_BGC_data_annotation$GTDB_phylum=="p__Verrucomicrobiota")),]
TH_BGC_data_annotation_4phylum$GTDB_phylum <- factor(TH_BGC_data_annotation_4phylum$GTDB_phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota", "p__Spirochaetota", "p__Verrucomicrobiota"))
TH_BGC_data_annotation_4phylum$BiG.SCAPE_class <- factor(TH_BGC_data_annotation_4phylum$BiG.SCAPE_class, levels = c("PKS-NRP_Hybrids","PKSI", "PKSother","Terpene","Others","NRPS","RiPPs"))
p_TH_BGC_class <- ggplot(TH_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count", width = 0.5) + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) +   
  scale_y_continuous(labels = percent) +   
  labs(x = "TH",
       y = "BGC class percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = "none")

p_TH_BGC <- ggarrange(p_TH_all, p_TH_BGC_class, widths = c(1,1), nrow = 2)

p <- ggplot(TS_BGC_data_annotation_4phylum, aes(x = GTDB_phylum, fill = BiG.SCAPE_class)) +
  geom_bar(position="fill", stat = "count") + 
  scale_fill_manual(values=c("NRPS" = "#ffffb3", "Others" = "#fdb462", "PKS-NRP_Hybrids" = "#bebada", 
                             "PKSI" = "#fb8072", "PKSother" = "#8dd3c7", "RiPPs" = "#80b1d3",
                             "Saccharides" = "#b3de69", "Terpene" = "#fccde5")) 
legend <- get_legend(p)

p_All_BGC <- ggarrange(p_TA_BGC, p_TAN_BGC, p_Yak_BGC, 
                       p_TH_BGC, p_TS_BGC, p_TC_BGC, 
                       widths = c(1,1,1,1,1,1), nrow = 2, ncol = 3)

p_All_BGC_leged <- ggarrange(p_All_BGC, legend,  
                             widths = c(6,1), nrow = 1)

ggsave(p_All_BGC_leged, filename = "p_All_BGC_new.pdf", width = 9, height = 10.3)
