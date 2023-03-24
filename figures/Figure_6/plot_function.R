setwd("~/Desktop/BGI/高原项目/pN_pS/inStrain")
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(reshape2)

## Bovinae:Caprinae
HGM05376_pS_function <- read.table("HGM05376_pS_function.txt", header = T, sep = "\t")
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
 
 
 
RGIG4079_pS_function <- read.table("RGIG4079_pS_function.txt", header = T, sep = "\t")
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

#RGIG4079_pS_function_1 <- RGIG4079_pS_function_long[which(RGIG4079_pS_function_long$Function=="Metabolism"),]
#p_RGIG4079_1 <- ggplot(RGIG4079_pS_function_1, aes(x = Pathway, y = ifelse(Group == "Bovinae", ko_num, -ko_num), fill = Group)) +
#  geom_bar(stat = 'identity') +  
#  scale_fill_manual(values=c("#66c2a5", "#8da0cb")) +
#  coord_flip() + 
#  scale_y_continuous(labels = abs, limits = c(-12, 12), breaks = seq(-12, 12, 4)) +
#  labs(x = '', y = '')  +
#  theme_bw() +
#  theme(panel.grid = element_blank(), legend.position = "none") 
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
 
