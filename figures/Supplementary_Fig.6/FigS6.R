library(tidyverse)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(ggtext)
library(reshape2)
library(UpSetR)
#install.packages("ggtext")
profile <- read.csv("0.data/all_sample_profile.csv", check.names=FALSE)
all_taxo <- read.csv("0.data/QTP_19251_SGBs_GTDB_r207.tsv", check.names=FALSE)
all_taxo <- separate(data = all_taxo, col = Classification, 
                     into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
all_taxo_profile <- merge(all_taxo, profile, by.x = "SGBs", by.y = "SGB_ID")
Sample <- read.csv("0.data/TableS0_samples_info.csv", check.names=FALSE)
Sample <- Sample[, c(2, 3)]

### phylum
phylum_profile <- all_taxo_profile[, -c(1, 2, 4, 5, 6, 7, 8, 9)]
phylum_profile_sum <- aggregate(phylum_profile[, -1], by=list(phylum_profile$Phylum), sum)
phylum_profile_sum_str <- as.data.frame(t(phylum_profile_sum))
colnames(phylum_profile_sum_str) <- phylum_profile_sum_str[1,]
phylum_profile_sum_str <- phylum_profile_sum_str[-1,]
phylum_profile_sum_str$Sample <- rownames(phylum_profile_sum_str)
phylum_profile_sample <- merge(phylum_profile_sum_str, Sample, by.x = "Sample", by.y = "Sample_ID")
#write.csv(phylum_profile_sample, file = "phylum_profile_sample.csv")
phylum_profile_sample <- phylum_profile_sample[, -1]
phylum_profile_sample_mean <- aggregate(as.data.frame(lapply(phylum_profile_sample[, -30],as.numeric)), by=list(type=phylum_profile_sample$Host), mean)
#write.csv(phylum_profile_sample_mean, file = "phylum_profile_sample_mean.csv")
phylum_final <- data.frame(matrix(ncol = 0, nrow = 1412))
phylum_final$Host <- phylum_profile_sample$Host
phylum_final$p__Firmicutes_A <- phylum_profile_sample$p__Firmicutes_A
phylum_final$p__Verrucomicobiota <- phylum_profile_sample$p__Verrucomicrobiota
phylum_final$p__Bacteroidota <- phylum_profile_sample$p__Bacteroidota
phylum_final$p__Firmicutes <- phylum_profile_sample$p__Firmicutes
phylum_final$p__Proteobacteria <- phylum_profile_sample$p__Proteobacteria
phylum_final$p__Spirochaetota <- phylum_profile_sample$p__Spirochaetota
#write.csv(phylum_final, file = "phylum_final.csv")
phylum_final <- read.csv("phylum_final.csv")
phylum_final <- phylum_final[,-1]
phylum_final_new <- gather(phylum_final, phylum, profile, p__Firmicutes_A:p__Spirochaetota)
phylum_final_new$Group1 <- phylum_final_new$Host
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Ass"] <- "Perissodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Horse"] <- "Perissodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Antelope"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Sheep"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Yak"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Tibetan_Cattle"] <- "Artiodactyla"
phylum_final_new$Group2 <- phylum_final_new$Host
phylum_final_new$Group2[phylum_final_new$Group2=="Tibetan_Antelope"] <- "Caprinae"
phylum_final_new$Group2[phylum_final_new$Group2=="Tibetan_Sheep"] <- "Caprinae"
phylum_final_new$Group2[phylum_final_new$Group2=="Tibetan_Yak"] <- "Bovinae"
phylum_final_new$Group2[phylum_final_new$Group2=="Tibetan_Cattle"] <- "Bovinae"
phylum_final_new$Host <- factor(phylum_final_new$Host, levels = c("Tibetan_Ass", "Tibetan_Horse",
                                                                  "Tibetan_Antelope", "Tibetan_Sheep",
                                                                  "Tibetan_Yak", "Tibetan_Cattle"))
aggregate(phylum_final_new$profile, by=list(type=phylum_final_new$phylum),mean)
phylum_final_new$phylum <- factor(phylum_final_new$phylum, levels = c("p__Firmicutes_A", "p__Bacteroidota",
                                                                  "p__Firmicutes", "p__Verrucomicobiota",
                                                                  "p__Proteobacteria", "p__Spirochaetota"))

df_p_val1 <- phylum_final_new %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Host) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val1_new1 <- df_p_val1[(df_p_val1$group1=="Tibetan_Ass"&df_p_val1$group2=="Tibetan_Horse"),]
df_p_val1_new2 <- df_p_val1[(df_p_val1$group1=="Tibetan_Antelope"&df_p_val1$group2=="Tibetan_Sheep"),]
df_p_val1_new3 <- df_p_val1[(df_p_val1$group1=="Tibetan_Yak"&df_p_val1$group2=="Tibetan_Cattle"),]
df_p_val1_new <- rbind(df_p_val1_new1,df_p_val1_new2,df_p_val1_new3)
df_p_val1_new$y.position <- 95

df_p_val2 <- phylum_final_new %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Group1) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val2$y.position <- 100

phylum_final3 <- phylum_final_new[phylum_final_new$Group2!="Tibetan_Ass",]
phylum_final3 <- phylum_final3[phylum_final3$Group2!="Tibetan_Horse",]
df_p_val3 <- phylum_final3 %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Group2) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val3$y.position <- 105

plot_phylum <- ggplot(phylum_final_new, aes(phylum, profile)) + 
  stat_boxplot(aes(color = Host), geom = "errorbar",width=0.6,
               position=position_dodge(0.8)) +
  labs(x = '', y = 'Relativate abudance (%)')+
  geom_boxplot(aes(fill = Host),position=position_dodge(0.8),
               outlier.size = 0.1,width=0.6)+
  stat_pvalue_manual(df_p_val1_new,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val2,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val3,label = '{p.signif}',tip.length = 0.01) +
  coord_flip() +
  theme_classic(base_family = "Helvetica") 

#ggsave(plot_phylum, filename = "phylum_6.pdf", height = 10, width = 8)
ggsave(plot_phylum, filename = "phylum_6_test.pdf", height = 8, width = 8)

### genus
genus_profile <- all_taxo_profile[, -c(1, 2, 3, 4, 5, 6, 8, 9)]
genus_profile_sum <- aggregate(genus_profile[, -1], by=list(genus_profile$Genus), sum)
### 0,1
genus_profile_01 <- genus_profile_sum
rownames(genus_profile_01) <- genus_profile_01$Group.1
genus_profile_01 <- genus_profile_01[,-1]
genus_profile_01[genus_profile_01>0] <- 1
genus_profile_01_str <- as.data.frame(t(genus_profile_01))
genus_profile_01_str$Sample <- rownames(genus_profile_01_str)
Sample$Host <- gsub("Tibetan_Ass", "TA", Sample$Host)
Sample$Host <- gsub("Tibetan_Antelope", "TAN", Sample$Host)
Sample$Host <- gsub("Tibetan_Cattle", "TC", Sample$Host)
Sample$Host <- gsub("Tibetan_Horse", "TH", Sample$Host)
Sample$Host <- gsub("Tibetan_Sheep", "TS", Sample$Host)
Sample$Host <- gsub("Tibetan_Yak", "Yak", Sample$Host)
genus_profile_sample <- merge(genus_profile_01_str, Sample, by.x = "Sample", by.y = "Sample_ID")
genus_profile_sample <- genus_profile_sample[,-1]
genus_profile_sample_sum <- aggregate(genus_profile_sample[, -972], by=list(genus_profile_sample$Host), sum)
rownames(genus_profile_sample_sum) <- genus_profile_sample_sum$Group.1
genus_profile_sample_sum <- genus_profile_sample_sum[,-1]
genus_profile_sample_sum_str <- as.data.frame(t(genus_profile_sample_sum))
genus_profile_sample_sum_str <- genus_profile_sample_sum_str[-1,]
#write.csv(genus_profile_sample_sum_str, file = "genus_profile_sample_sum_str.csv")
genus_profile_core50 <- genus_profile_sample_sum_str[which(genus_profile_sample_sum_str$TA>=24 | genus_profile_sample_sum_str$TAN>=127.5 | genus_profile_sample_sum_str$TC>=98
                                                           | genus_profile_sample_sum_str$TH>=39.5 | genus_profile_sample_sum_str$TS>=223 | genus_profile_sample_sum_str$Yak>=194),]
genus_profile_core50[genus_profile_core50>0] <- 1
sample_lst = c("TC","Yak","TS","TAN","TH","TA")
plot <- upset(genus_profile_core50, sets = sample_lst, order.by = c("freq"),
              keep.order = TRUE,
              nintersects = 10,
              queries = list(list(query = intersects, 
                                  params = list(colnames(genus_profile_sample_sum_str)), 
                                  active = T)))
###
genus_profile_sum_str <- as.data.frame(t(genus_profile_sum))
colnames(genus_profile_sum_str) <- genus_profile_sum_str[1,]
genus_profile_sum_str <- genus_profile_sum_str[-1,]
genus_profile_sum_str$Sample <- rownames(genus_profile_sum_str)
genus_profile_sample <- merge(genus_profile_sum_str, Sample, by.x = "Sample", by.y = "Sample_ID")
#write.csv(genus_profile_sample, file = "genus_profile_sample.csv")
genus_profile_sample <- genus_profile_sample[, -1]
genus_profile_sample_mean <- aggregate(as.data.frame(lapply(genus_profile_sample[, -972],as.numeric)), by=list(type=genus_profile_sample$Host), mean)
#write.csv(genus_profile_sample_mean, file = "genus_profile_sample_mean.csv")
genus_final <- data.frame(matrix(ncol = 0, nrow = 1412))
genus_final$Host <- genus_profile_sample$Host
genus_final$g__Faecousia <- genus_profile_sample$g__Faecousia
genus_final$g__Alistipes <- genus_profile_sample$g__Alistipes
genus_final$g__Cryptobacteroides <- genus_profile_sample$g__Cryptobacteroides
genus_final$g__Phocaeicola <- genus_profile_sample$g__Phocaeicola
genus_final$g__Treponema_D <- genus_profile_sample$g__Treponema_D
genus_final$g__UBA3663 <- genus_profile_sample$g__UBA3663
genus_final$g__Prevotella <- genus_profile_sample$g__Prevotella
genus_final$g__RF16 <- genus_profile_sample$g__RF16
genus_final$g__HGM04593 <- genus_profile_sample$g__HGM04593
genus_final$g__UBA1258 <- genus_profile_sample$g__UBA1258
#write.csv(genus_final, file = "genus_final.csv")
genus_final <- read.csv("genus_final.csv")
genus_final <- genus_final[,-1]
genus_final_new <- gather(genus_final, genus, profile, g__Faecousia:g__UBA1258)
genus_final_new$Group1 <- genus_final_new$Host
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Ass"] <- "Perissodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Horse"] <- "Perissodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Antelope"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Sheep"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Yak"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Tibetan_Cattle"] <- "Artiodactyla"
genus_final_new$Group2 <- genus_final_new$Host
genus_final_new$Group2[genus_final_new$Group2=="Tibetan_Antelope"] <- "Caprinae"
genus_final_new$Group2[genus_final_new$Group2=="Tibetan_Sheep"] <- "Caprinae"
genus_final_new$Group2[genus_final_new$Group2=="Tibetan_Yak"] <- "Bovinae"
genus_final_new$Group2[genus_final_new$Group2=="Tibetan_Cattle"] <- "Bovinae"
genus_final_new$Host <- factor(genus_final_new$Host, levels = c("Tibetan_Ass", "Tibetan_Horse",
                                                                "Tibetan_Antelope", "Tibetan_Sheep",
                                                                "Tibetan_Yak", "Tibetan_Cattle"))
aggregate(genus_final_new$profile, by=list(type=genus_final_new$genus),mean)
genus_final_new$genus <- factor(genus_final_new$genus, levels = c("g__Cryptobacteroides", "g__Faecousia",
                                                                  "g__Alistipes", "g__HGM04593",
                                                                  "g__Phocaeicola", "g__RF16",
                                                                  "g__Treponema_D", "g__Prevotella",
                                                                  "g__UBA3663", "g__UBA1258"))

genus_final_new <- genus_final_new[genus_final_new$profile>0,]
df_p_val1 <- genus_final_new %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Host) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')

df_p_val1_new1 <- df_p_val1[(df_p_val1$group1=="Tibetan_Ass"&df_p_val1$group2=="Tibetan_Horse"),]
df_p_val1_new2 <- df_p_val1[(df_p_val1$group1=="Tibetan_Antelope"&df_p_val1$group2=="Tibetan_Sheep"),]
df_p_val1_new3 <- df_p_val1[(df_p_val1$group1=="Tibetan_Yak"&df_p_val1$group2=="Tibetan_Cattle"),]
df_p_val1_new <- rbind(df_p_val1_new1,df_p_val1_new2,df_p_val1_new3)
df_p_val1_new$y.position <- 40

genus_final2 <- genus_final_new[genus_final_new$genus!="g__HGM04593",]
df_p_val2 <- genus_final2 %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Group1) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')
df_p_val2$y.position <- 43

genus_final3 <- genus_final_new[genus_final_new$Group2!="Tibetan_Ass",]
genus_final3 <- genus_final3[genus_final3$Group2!="Tibetan_Horse",]
df_p_val3 <- genus_final3 %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Group2) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')
df_p_val3$y.position <- 46

plot_genus <- ggplot(genus_final_new, aes(genus, profile)) + 
  stat_boxplot(aes(color = Host), geom = "errorbar",width=0.6,
               position=position_dodge(0.8)) +
  geom_boxplot(aes(fill = Host), width=0.6,
               position=position_dodge(0.8), outlier.size = 0.1)+
  labs(x = '', y = 'Relativate abudance (%)')+
  stat_pvalue_manual(df_p_val1_new,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val2,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val3,label = '{p.signif}',tip.length = 0.01) +
  coord_flip() +
  theme_classic(base_family = "Helvetica") 

ggsave(plot_genus, filename = "genus_10_test.pdf", height = 10, width = 8)


## a. Alpha diversity
cal_shannon_richness <- function(infile){
  #infile <- "../00.tables/profile.csv"
  profile.dset <- read.csv(infile, header = T, row.names = 1, stringsAsFactors=FALSE, check.names = F)
  shannon.res <- diversity(profile.dset, index = "shannon")
  simpson.res <- diversity(profile.dset, index = "simpson")
  
  res.dset <- tibble(Sample_ID = rownames(profile.dset), Richness = Reduce(`+`, as.data.frame(profile.dset > 0))) %>%
    left_join(tibble(Sample_ID = names(shannon.res), Shannon_index = shannon.res), by = "Sample_ID") %>%
    left_join(tibble(Sample_ID = names(simpson.res), Simpson_index = simpson.res), by = "Sample_ID")
  return(res.dset) 
}
rownames(genus_profile_sum) <- genus_profile_sum[,1]
genus_profile_sum <- genus_profile_sum[,-1]
genus_profile_sum_str <- as.data.frame(t(genus_profile_sum))
write.csv(genus_profile_sum_str, file = "genus_profile_sum_str.csv")
d.res <- cal_shannon_richness("genus_profile_sum_str.csv")
sample.dset <- Sample
meta.alpha <- left_join(sample.dset, d.res, by = "Sample_ID") %>%
  pivot_longer(Richness : Simpson_index,
               names_to = "Diversity",
               values_to = "value")
# Host
meta.alpha$Host <- gsub("Tibetan_Ass", "TA", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Antelope", "TAN", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Cattle", "TC", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Horse", "TH", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Sheep", "TS", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Yak", "Yak", meta.alpha$Host)
meta.alpha$Host <- factor(meta.alpha$Host, levels = c("TA", "TH", "TAN", "TS", "Yak", "TC"))
Host_comparisons <- list(c("TA", "TH"), c("TAN", "TS"), c("Yak", "TC"))
richness.plot <- ggviolin(filter(meta.alpha, Diversity == "Richness"), x = "Host", y = "value", fill = "Host") +
  geom_boxplot(aes(fill = Host),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = Host_comparisons,
              y_position = c(550),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("Richness") 
ggsave("genus.Richness.pdf", richness.plot, width = 5.5, height = 4)
Shannon.plot <- ggviolin(filter(meta.alpha, Diversity == "Shannon_index"), x = "Host", y = "value", fill = "Host") +
  geom_boxplot(aes(fill = Host),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = Host_comparisons, 
              y_position = c(5.3),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("Shannon index") 
ggsave("genus.Shannon.pdf", Shannon.plot, width = 5.5, height = 4)

# order
meta.alpha$Host <- gsub("Tibetan_Ass", "Perissodactyla", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Antelope", "Artiodactyla", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Cattle", "Artiodactyla", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Horse", "Perissodactyla", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Sheep", "Artiodactyla", meta.alpha$Host)
meta.alpha$Host <- gsub("Tibetan_Yak", "Artiodactyla", meta.alpha$Host)
richness.plot <- ggviolin(filter(meta.alpha, Diversity == "Richness"), x = "Host", y = "value", fill = "Host") +
  geom_boxplot(aes(fill = Host),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("Perissodactyla","Artiodactyla")),
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("Richness") 
Shannon.plot <- ggviolin(filter(meta.alpha, Diversity == "Shannon_index"), x = "Host", y = "value", fill = "Host") +
  geom_boxplot(aes(fill = Host),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("Perissodactyla","Artiodactyla")), 
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  xlab("") +
  ylab("Shannon index") 
