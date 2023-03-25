# Figure 6. Functional divergence of representative SGBs with high frequent co-speciation or host-swap events across host species.
# BCDEFG
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(reshape2)

## Sample_data
sample_data <- read.table("sample.txt", header = T)
sample_data$Group1 <- sample_data$Host
sample_data$Group1[sample_data$Group1=="TA"] <- "Perissodactyla"
sample_data$Group1[sample_data$Group1=="TH"] <- "Perissodactyla"
sample_data$Group1[sample_data$Group1=="TAN"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="TS"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="Yak"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="TC"] <- "Artiodactyla"
sample_data$Group2 <- sample_data$Host
sample_data$Group2[sample_data$Group2=="TAN"] <- "Caprinae"
sample_data$Group2[sample_data$Group2=="TS"] <- "Caprinae"
sample_data$Group2[sample_data$Group2=="Yak"] <- "Bovinae"
sample_data$Group2[sample_data$Group2=="TC"] <- "Bovinae"

## HGM05376
# data
HGM05376_data <- read.csv("HGM05376/HGM05376_instrain_result.tsv", header = F)
colnames(HGM05376_data) <- c("Sample", "scaffold", "gene", "coverage",	"breadth", "pNpS_variants")
HGM05376_data_merge <- merge(HGM05376_data, sample_data, by.x = "Sample", by.y = "Sample_ID")
HGM05376_data_merge$Host <- factor(HGM05376_data_merge$Host, levels = c("TA", "TH", "TAN", "TS","Yak", "TC"))
function_HGM05376 <- read.table("HGM05376/function_HGM05376_SGB.txt", sep = "\t")
colnames(function_HGM05376) <- c("Gene", "evalue", "score", "COG_category", "GOs", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "CAZy")
function_HGM05376_merge <- merge(HGM05376_data_merge, function_HGM05376,
                                 by.x = "gene", by.y = "Gene", all.x = T)

# violin gene_num
HGM05376_data_merge_new1 <- HGM05376_data_merge[which(HGM05376_data_merge$pNpS_variants>1),]
HGM05376_data_merge_new1 <- HGM05376_data_merge_new1[HGM05376_data_merge_new1$pNpS_variants!="1",]
HGM05376_data_gene_num <- as.data.frame(table(HGM05376_data_merge_new1$Sample))
colnames(HGM05376_data_gene_num) <- c("Sample", "Gene_num")
HGM05376_data_gene_num_merge <- merge(sample_data, HGM05376_data_gene_num,
                                      by.x = "Sample_ID", by.y = "Sample")

HGM05376_data_gene_num_merge$Group2 <- factor(HGM05376_data_gene_num_merge$Group2, levels = c("Caprinae", "Bovinae"))
p_HGM05376_group2_gene_num_pN <- ggviolin(HGM05376_data_gene_num_merge, x = "Group2", y = "Gene_num", fill = "Group2") +
  geom_boxplot(aes(fill = Group2),width=0.1, outlier.size = 0.5) +
  scale_y_continuous(trans='log10') +
  geom_signif(comparisons = list(c("Caprinae", "Bovinae")), y_position = c(log10(200)),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("HGM05376 gene number (pN/pS > 1)")

HGM05376_data_merge_new1 <- HGM05376_data_merge[HGM05376_data_merge$pNpS_variants<1,]
HGM05376_data_gene_num <- as.data.frame(table(HGM05376_data_merge_new1$Sample))
colnames(HGM05376_data_gene_num) <- c("Sample", "Gene_num")
HGM05376_data_gene_num_merge <- merge(sample_data, HGM05376_data_gene_num,
                                      by.x = "Sample_ID", by.y = "Sample")

HGM05376_data_gene_num_merge$Group2 <- factor(HGM05376_data_gene_num_merge$Group2, levels = c("Caprinae", "Bovinae"))
p_HGM05376_group2_gene_num_pS <- ggviolin(HGM05376_data_gene_num_merge, x = "Group2", y = "Gene_num", fill = "Group2") +
  geom_boxplot(aes(fill = Group2),width=0.1, outlier.size = 0.5) +
  scale_y_continuous(trans='log10') +
  geom_signif(comparisons = list(c("Caprinae", "Bovinae")), y_position = c(log10(5000)),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("HGM05376 gene number (pN/pS < 1)")

## RGIG4079
# data

RGIG4079_data <- read.csv("RGIG4079/RGIG4079_instrain_result.tsv", header = F)
colnames(RGIG4079_data) <- c("Sample", "scaffold", "gene", "coverage",	"breadth", "pNpS_variants")
RGIG4079_data_merge <- merge(RGIG4079_data, sample_data, by.x = "Sample", by.y = "Sample_ID")
RGIG4079_data_merge$Host <- factor(RGIG4079_data_merge$Host, levels = c("TAN", "TS","Yak", "TC"))
function_RGIG4079 <- read.table("RGIG4079/function_RGIG4079_SGB.txt", sep = "\t")
colnames(function_RGIG4079) <- c("Gene", "evalue", "score", "COG_category", "GOs", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "CAZy")
function_RGIG4079_merge <- merge(RGIG4079_data_merge, function_RGIG4079,
                                 by.x = "gene", by.y = "Gene", all.x = T)

# violin gene_num
RGIG4079_data_merge_new1 <- RGIG4079_data_merge[RGIG4079_data_merge$pNpS_variants>1,]
RGIG4079_data_gene_num <- as.data.frame(table(RGIG4079_data_merge_new1$Sample))
colnames(RGIG4079_data_gene_num) <- c("Sample", "Gene_num")
RGIG4079_data_gene_num_merge <- merge(sample_data, RGIG4079_data_gene_num,
                                      by.x = "Sample_ID", by.y = "Sample")

RGIG4079_data_gene_num_merge$Group2 <- factor(RGIG4079_data_gene_num_merge$Group2, levels = c("Caprinae", "Bovinae"))
p_RGIG4079_group2_gene_num_pN <- ggviolin(RGIG4079_data_gene_num_merge, x = "Group2", y = "Gene_num", fill = "Group2") +
  geom_boxplot(aes(fill = Group2),width=0.1, outlier.size = 0.5) +
  scale_y_continuous(trans='log10') +
  geom_signif(comparisons = list(c("Caprinae", "Bovinae")), y_position = c(log10(100)),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("RGIG4079 gene number (pN/pS > 1)")

RGIG4079_data_merge_new1 <- RGIG4079_data_merge[RGIG4079_data_merge$pNpS_variants<1,]
RGIG4079_data_gene_num <- as.data.frame(table(RGIG4079_data_merge_new1$Sample))
colnames(RGIG4079_data_gene_num) <- c("Sample", "Gene_num")
RGIG4079_data_gene_num_merge <- merge(sample_data, RGIG4079_data_gene_num,
                                      by.x = "Sample_ID", by.y = "Sample")

RGIG4079_data_gene_num_merge$Group2 <- factor(RGIG4079_data_gene_num_merge$Group2, levels = c("Caprinae", "Bovinae"))
p_RGIG4079_group2_gene_num_pS <- ggviolin(RGIG4079_data_gene_num_merge, x = "Group2", y = "Gene_num", fill = "Group2") +
  geom_boxplot(aes(fill = Group2),width=0.1, outlier.size = 0.5) +
  scale_y_continuous(trans='log10') +
  geom_signif(comparisons = list(c("Caprinae", "Bovinae")), y_position = c(log10(1000)),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("RGIG4079 gene number (pN/pS < 1)")

gene_num_merge_two_group <- plot_grid(p_HGM05376_group2_gene_num_pN, p_HGM05376_group2_gene_num_pS, 
                                      p_RGIG4079_group2_gene_num_pN, p_RGIG4079_group2_gene_num_pS, 
                                      labels = c("A", "B", "C", "D"))
ggsave(gene_num_merge_two_group, filename = "gene_num_merge_two_group.pdf", width = 10, height = 8)


## Bovinae:Caprinae
HGM05376_pS_function <- read.table("HGM05376/HGM05376_pS_function.txt", header = T, sep = "\t")
HGM05376_pS_function_long <-gather(HGM05376_pS_function, Group, ko_num, Bovinae:Caprinae)
HGM05376_pS_function_long$Pathway <- factor(HGM05376_pS_function_long$Pathway, levels = c("Bacterial secretion system", "Two-component system", "ABC transporters", 
                                                                                          "Ribosome", "Protein export", "Metabolic pathways", "Biosynthesis of secondary metabolites", 
                                                                                          "Biosynthesis of amino acids", "Microbial metabolism in diverse environments", "Carbon metabolism", "Glycolysis / Gluconeogenesis",
                                                                                          "Biosynthesis of cofactors"))
HGM05376_pS_function$Pathway <- factor(HGM05376_pS_function$Pathway, levels = c("Bacterial secretion system", "Two-component system", "ABC transporters", 
                                                                                "Ribosome", "Protein export", "Metabolic pathways", "Biosynthesis of secondary metabolites", 
                                                                                "Biosynthesis of amino acids", "Microbial metabolism in diverse environments", "Carbon metabolism", "Glycolysis / Gluconeogenesis",
                                                                                "Biosynthesis of cofactors"))
p_HGM05376 <- ggplot(HGM05376_pS_function_long) + 
  geom_col(aes(ko_num, Pathway, fill = Group), width = 0.5) +
  facet_grid(.~Group) +
  labs(x = 'ko number', y = '') +
  theme_classic(base_family = "Helvetica") 
ggsave(p_HGM05376, filename = "p_HGM05376.pdf", height = 8, width = 8)

HGM05376_pS_function_1 <- HGM05376_pS_function_long[which(HGM05376_pS_function_long$Function=="Metabolism"),]
p_HGM05376_1 <- ggplot(HGM05376_pS_function_1, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity') +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
  labs(x = '', y = '')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") 
HGM05376_pS_function_2 <- HGM05376_pS_function_long[which(HGM05376_pS_function_long$Function=="Environmental Information Processing"),]
p_HGM05376_2 <- ggplot(HGM05376_pS_function_2, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity') +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
  labs(x = '', y = '')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") 
HGM05376_pS_function_3 <- HGM05376_pS_function_long[which(HGM05376_pS_function_long$Function=="Genetic Information Processing"),]
p_HGM05376_3 <- ggplot(HGM05376_pS_function_3, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity') +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-10, 10), breaks = seq(-10, 10, 5)) +
  labs(x = '', y = '')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom") 
p_HGM05376_3_new <- p_HGM05376_3 + theme(legend.position = "none")
legend <- get_legend(p_HGM05376_3)
p_HGM05376_merge <- plot_grid(p_HGM05376_1, p_HGM05376_2, p_HGM05376_3_new, nrow = 3, align = "v", rel_heights = c(2,2.5,2))
p_merge <- plot_grid(p_HGM05376_merge, legend, nrow = 2, rel_heights = c(10,1))
ggsave(p_merge, filename = "p_HGM05376_merge.pdf", height = 4, width = 6)


RGIG4079_pS_function <- read.table("RGIG4079/RGIG4079_pS_function.txt", header = T, sep = "\t")
RGIG4079_pS_function_long <-gather(RGIG4079_pS_function, Group, ko_num, Bovinae:Caprinae)
RGIG4079_pS_function_long$Pathway <- factor(RGIG4079_pS_function_long$Pathway, levels = c("Two-component system", "Ribosome", "Metabolic pathways", "Biosynthesis of secondary metabolites"))
p_RGIG4079 <- ggplot(RGIG4079_pS_function_long) + 
  geom_col(aes(ko_num, Pathway, fill = Group), width = 0.5) +
  facet_grid(.~Group) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  labs(x = 'ko number', y = '') +
  theme_bw() +
  theme_classic(base_family = "Helvetica") 
ggsave(p_RGIG4079, filename = "p_RGIG4079.pdf", height = 5, width = 8)

p_RGIG4079 <- ggplot(RGIG4079_pS_function_long, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity', width = 0.8) +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-12, 12), breaks = seq(-12, 12, 4)) +
  labs(x = '', y = 'KOs Number')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(p_RGIG4079, filename = "p_RGIG4079_merge.pdf", height = 3, width = 6)

RGIG4079_pS_function_2 <- RGIG4079_pS_function_long[which(RGIG4079_pS_function_long$Function=="Environmental Information Processing"),]
p_RGIG4079_2 <- ggplot(RGIG4079_pS_function_2, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity') +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-12, 12), breaks = seq(-12, 12, 4)) +
  labs(x = '', y = '')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") 
RGIG4079_pS_function_3 <- RGIG4079_pS_function_long[which(RGIG4079_pS_function_long$Function=="Genetic Information Processing"),]
p_RGIG4079_3 <- ggplot(RGIG4079_pS_function_3, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
  geom_bar(stat = 'identity') +  
  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
  coord_flip() + 
  scale_y_continuous(labels = abs, limits = c(-12, 12), breaks = seq(-12, 12, 4)) +
  labs(x = '', y = '')  +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom") 
p_RGIG4079_3_new <- p_RGIG4079_3 + theme(legend.position = "none")
legend <- get_legend(p_RGIG4079_3)
p_RGIG4079_merge <- plot_grid(p_RGIG4079_3_new, p_RGIG4079_2, nrow = 2, align = "v", rel_heights = c(0.5,0.5))
p_merge <- plot_grid(p_RGIG4079_merge, legend, nrow = 2, rel_heights = c(10,1))
ggsave(p_merge, filename = "p_RGIG4079_merge.pdf", height = 2, width = 6)



