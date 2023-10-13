setwd("~/Desktop/BGI/高原项目/QTP_redo/Figure S6.taxo_plot")
library(tidyverse)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(ggtext)
library(reshape2)
library(UpSetR)
#install.packages("ggtext")
profile <- read.csv("0.data/SGBs_14062_profile.csv", check.names=FALSE)
rownames(profile) <- profile$Sample
profile <- profile[,-1]
profile <- as.data.frame(t(profile))
profile$SGB_ID <- rownames(profile)
all_taxo <- read.csv("0.data/QTP_SGBs_14062_GTDB_r214.csv", check.names=FALSE)
all_taxo <- separate(data = all_taxo, col = taxonomy, 
                     into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
all_taxo_profile <- merge(all_taxo, profile, by = "SGB_ID")
Sample <- read.csv("0.data/TableS0_samples_info.csv", check.names=FALSE)
Sample <- Sample[, c(2, 3)]
Sample$Host <- gsub("Tibetan_Ass", "TA", Sample$Host)
Sample$Host <- gsub("Tibetan_Antelope", "TAN", Sample$Host)
Sample$Host <- gsub("Tibetan_Cattle", "TC", Sample$Host)
Sample$Host <- gsub("Tibetan_Horse", "TH", Sample$Host)
Sample$Host <- gsub("Tibetan_Sheep", "TS", Sample$Host)
Sample$Host <- gsub("Tibetan_Yak", "Yak", Sample$Host)

### phylum
phylum_profile <- all_taxo_profile[, -c(1, 2, 4, 5, 6, 7, 8)]
phylum_profile_sum <- aggregate(phylum_profile[, -1], by=list(phylum_profile$Phylum), sum)
phylum_profile_sum_str <- as.data.frame(t(phylum_profile_sum))
colnames(phylum_profile_sum_str) <- phylum_profile_sum_str[1,]
phylum_profile_sum_str <- phylum_profile_sum_str[-1,]
phylum_profile_sum_str$Sample <- rownames(phylum_profile_sum_str)
phylum_profile_sample <- merge(phylum_profile_sum_str, Sample, by.x = "Sample", by.y = "Sample_ID")
#write.csv(phylum_profile_sample, file = "phylum_profile_sample.csv")
phylum_profile_sample <- phylum_profile_sample[, -1]
phylum_profile_sample_mean <- aggregate(as.data.frame(lapply(phylum_profile_sample[, -29],as.numeric)), by=list(type=phylum_profile_sample$Host), mean)
#write.csv(phylum_profile_sample_mean, file = "phylum_profile_sample_mean.csv")
phylum_final <- data.frame(matrix(ncol = 0, nrow = 1412))
phylum_final$Host <- phylum_profile_sample$Host
phylum_final$p__Bacillota_A <- phylum_profile_sample$p__Bacillota_A
phylum_final$p__Bacteroidota <- phylum_profile_sample$p__Bacteroidota
phylum_final$p__Verrucomicrobiota <- phylum_profile_sample$p__Verrucomicrobiota
phylum_final$p__Bacillota <- phylum_profile_sample$p__Bacillota
phylum_final$p__Pseudomonadota <- phylum_profile_sample$p__Pseudomonadota
phylum_final$p__Spirochaetota <- phylum_profile_sample$p__Spirochaetota
phylum_final$p__Fibrobacterota <- phylum_profile_sample$p__Fibrobacterota
#write.csv(phylum_final, file = "phylum_final.csv")
#phylum_final <- read.csv("phylum_final.csv")
#phylum_final <- phylum_final[,-1]
phylum_final_new <- gather(phylum_final, phylum, profile, p__Bacillota_A:p__Fibrobacterota)
phylum_final_new$Group1 <- phylum_final_new$Host
phylum_final_new$Group1[phylum_final_new$Group1=="TA"] <- "Perissodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="TH"] <- "Perissodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="TAN"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="TS"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="Yak"] <- "Artiodactyla"
phylum_final_new$Group1[phylum_final_new$Group1=="TC"] <- "Artiodactyla"
phylum_final_new$Group2 <- phylum_final_new$Host
phylum_final_new$Group2[phylum_final_new$Group2=="TAN"] <- "Caprinae"
phylum_final_new$Group2[phylum_final_new$Group2=="TS"] <- "Caprinae"
phylum_final_new$Group2[phylum_final_new$Group2=="Yak"] <- "Bovinae"
phylum_final_new$Group2[phylum_final_new$Group2=="TC"] <- "Bovinae"
phylum_final_new$Host <- factor(phylum_final_new$Host, levels = c("TA", "TH",
                                                                  "TAN", "TS",
                                                                  "Yak", "TC"))
phylum_final_new$profile <- as.numeric(phylum_final_new$profile)
aggregate(phylum_final_new$profile, by=list(type=phylum_final_new$phylum),mean)
phylum_final_new$phylum <- factor(phylum_final_new$phylum, levels = c("p__Bacillota_A", "p__Bacteroidota",
                                                                  "p__Verrucomicrobiota", "p__Bacillota",
                                                                  "p__Pseudomonadota", "p__Spirochaetota", "p__Fibrobacterota"))

df_p_val1 <- phylum_final_new %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Host) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val1_new1 <- df_p_val1[(df_p_val1$group1=="TA"&df_p_val1$group2=="TH"),]
df_p_val1_new2 <- df_p_val1[(df_p_val1$group1=="TAN"&df_p_val1$group2=="TS"),]
df_p_val1_new3 <- df_p_val1[(df_p_val1$group1=="Yak"&df_p_val1$group2=="TC"),]
df_p_val1_new <- rbind(df_p_val1_new1,df_p_val1_new2,df_p_val1_new3)
df_p_val1_new$y.position <- 95
phylum_hosts_p_value <- as.data.frame(df_p_val1_new)
phylum_hosts_p_value <- phylum_hosts_p_value[,c(1,3,4,8,11,15,16)]
write.csv(phylum_hosts_p_value, file = "phylum_hosts_p_value.csv")

df_p_val2 <- phylum_final_new %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Group1) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val2$y.position <- 105
phylum_Group1_p_value <- as.data.frame(df_p_val2)
phylum_Group1_p_value <- phylum_Group1_p_value[,c(1,3,4,8,9,13,14)]
write.csv(phylum_Group1_p_value, file = "phylum_Group1_p_value.csv")

phylum_final3 <- phylum_final_new[phylum_final_new$Group2!="TA",]
phylum_final3 <- phylum_final3[phylum_final3$Group2!="TH",]
df_p_val3 <- phylum_final3 %>% 
  group_by(phylum) %>% 
  wilcox_test(formula = profile ~ Group2) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='phylum')
df_p_val3$y.position <- 100
phylum_Group2_p_value <- as.data.frame(df_p_val3)
phylum_Group2_p_value <- phylum_Group2_p_value[,c(1,3,4,8,9,13,14)]
write.csv(phylum_Group2_p_value, file = "phylum_Group2_p_value.csv")

plot_phylum <- ggplot(phylum_final_new, aes(phylum, profile)) + 
  stat_boxplot(aes(color = Host), geom = "errorbar",width=0.6,
               position=position_dodge(0.8)) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  labs(x = '', y = 'Relativate abudance (%)')+
  geom_boxplot(aes(fill = Host),position=position_dodge(0.8),
               outlier.size = 0.1,width=0.6)+
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                              "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  stat_pvalue_manual(df_p_val1_new,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val2,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val3,label = '{p.signif}',tip.length = 0.01) +
  coord_flip() +
  theme_classic() 

ggsave(plot_phylum, filename = "1.result/phylum_7_plot.pdf", height = 8, width = 8)

### genus
genus_profile <- all_taxo_profile[, -c(1, 2, 3, 4, 5, 6, 8)]
genus_profile_sum <- aggregate(genus_profile[, -1], by=list(genus_profile$Genus), sum)
### plot
genus_profile_sum_str <- as.data.frame(t(genus_profile_sum))
colnames(genus_profile_sum_str) <- genus_profile_sum_str[1,]
genus_profile_sum_str <- genus_profile_sum_str[-1,]
genus_profile_sum_str$Sample <- rownames(genus_profile_sum_str)
genus_profile_sample <- merge(genus_profile_sum_str, Sample, by.x = "Sample", by.y = "Sample_ID")
#write.csv(genus_profile_sample, file = "genus_profile_sample.csv")
genus_profile_sample <- genus_profile_sample[, -1]
genus_profile_sample_mean <- aggregate(as.data.frame(lapply(genus_profile_sample[, -960],as.numeric)), by=list(type=genus_profile_sample$Host), mean)
#write.csv(genus_profile_sample_mean, file = "genus_profile_sample_mean.csv")
genus_final <- data.frame(matrix(ncol = 0, nrow = 1412))
genus_final$Host <- genus_profile_sample$Host
genus_final$g__Cryptobacteroides <- genus_profile_sample$g__Cryptobacteroides
genus_final$g__Alistipes <- genus_profile_sample$g__Alistipes
genus_final$g__RF16 <- genus_profile_sample$g__RF16
genus_final$g__HGM04593 <- genus_profile_sample$g__HGM04593
genus_final$g__UBA1189 <- genus_profile_sample$g__UBA1189
genus_final$g__Faecousia <- genus_profile_sample$g__Faecousia
genus_final$g__Phocaeicola <- genus_profile_sample$g__Phocaeicola
genus_final$g__Treponema_D <- genus_profile_sample$g__Treponema_D
genus_final$g__UBA1258 <- genus_profile_sample$g__UBA1258
genus_final$g__NK4A136 <- genus_profile_sample$g__NK4A136
genus_final$g__UBA1740 <- genus_profile_sample$g__UBA1740
genus_final$g__Prevotella <- genus_profile_sample$g__Prevotella
genus_final$g__UBA3663 <- genus_profile_sample$g__UBA3663
#write.csv(genus_final, file = "genus_final.csv")
#genus_final <- read.csv("genus_final.csv")
#genus_final <- genus_final[,-1]
genus_final_new <- gather(genus_final, genus, profile, g__Cryptobacteroides:g__UBA3663)
genus_final_new$Group1 <- genus_final_new$Host
genus_final_new$Group1[genus_final_new$Group1=="TA"] <- "Perissodactyla"
genus_final_new$Group1[genus_final_new$Group1=="TH"] <- "Perissodactyla"
genus_final_new$Group1[genus_final_new$Group1=="TAN"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="TS"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="Yak"] <- "Artiodactyla"
genus_final_new$Group1[genus_final_new$Group1=="TC"] <- "Artiodactyla"
genus_final_new$Group2 <- genus_final_new$Host
genus_final_new$Group2[genus_final_new$Group2=="TAN"] <- "Caprinae"
genus_final_new$Group2[genus_final_new$Group2=="TS"] <- "Caprinae"
genus_final_new$Group2[genus_final_new$Group2=="Yak"] <- "Bovinae"
genus_final_new$Group2[genus_final_new$Group2=="TC"] <- "Bovinae"
genus_final_new$Host <- factor(genus_final_new$Host, levels = c("TA", "TH",
                                                                "TAN", "TS",
                                                                "Yak", "TC"))
genus_final_new$profile <- as.numeric(genus_final_new$profile)
gene_sort <- aggregate(genus_final_new$profile, by=list(type=genus_final_new$genus),mean)
genus_final_new$genus <- factor(genus_final_new$genus, levels = c("g__Cryptobacteroides", "g__Alistipes",
                                                                  "g__Faecousia", "g__HGM04593",
                                                                  "g__RF16", "g__UBA1189",
                                                                  "g__Phocaeicola", "g__UBA1740",
                                                                  "g__Treponema_D", "g__Prevotella",
                                                                  "g__NK4A136", "g__UBA3663", "g__UBA1258"))

genus_final_new <- genus_final_new[genus_final_new$profile>0,]
df_p_val1 <- genus_final_new %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Host) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')

df_p_val1_new1 <- df_p_val1[(df_p_val1$group1=="TA"&df_p_val1$group2=="TH"),]
df_p_val1_new2 <- df_p_val1[(df_p_val1$group1=="TAN"&df_p_val1$group2=="TS"),]
df_p_val1_new3 <- df_p_val1[(df_p_val1$group1=="Yak"&df_p_val1$group2=="TC"),]
df_p_val1_new <- rbind(df_p_val1_new1,df_p_val1_new2,df_p_val1_new3)
df_p_val1_new$y.position <- 40
genus_host_p_value <- as.data.frame(df_p_val1_new)
genus_host_p_value <- genus_host_p_value[,c(1,3,4,8,11,15,16)]
write.csv(genus_host_p_value, file = "genus_host_p_value.csv")

genus_final2 <- genus_final_new[genus_final_new$genus!="g__Alistipes",]
genus_final2 <- genus_final2[genus_final2$genus!="g__HGM04593",]
genus_final2 <- genus_final2[genus_final2$genus!="g__Phocaeicola",]
genus_final2 <- genus_final2[genus_final2$genus!="g__UBA1258",]
df_p_val2 <- genus_final2 %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Group1) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')
df_p_val2$y.position <- 46
genus_Group1_p_value <- as.data.frame(df_p_val2)
genus_Group1_p_value <- genus_Group1_p_value[,c(1,3,4,8,9,13,14)]
write.csv(genus_Group1_p_value, file = "genus_Group1_p_value.csv")

genus_final3 <- genus_final_new[genus_final_new$Group2!="TA",]
genus_final3 <- genus_final3[genus_final3$Group2!="TH",]
genus_final3 <- genus_final3[genus_final3$genus!="g__UBA1258",]
df_p_val3 <- genus_final3 %>% 
  group_by(genus) %>% 
  wilcox_test(formula = profile ~ Group2) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='genus')
df_p_val3$y.position <- 43
genus_Group2_p_value <- as.data.frame(df_p_val3)
genus_Group2_p_value <- genus_Group2_p_value[,c(1,3,4,8,9,13,14)]
write.csv(genus_Group2_p_value, file = "genus_Group2_p_value.csv")

plot_genus <- ggplot(genus_final_new, aes(genus, profile)) + 
  stat_boxplot(aes(color = Host), geom = "errorbar",width=0.6,
               position=position_dodge(0.8)) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  geom_boxplot(aes(fill = Host), width=0.6,
               position=position_dodge(0.8), outlier.size = 0.1)+
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  labs(x = '', y = 'Relativate abudance (%)')+
  stat_pvalue_manual(df_p_val1_new,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val2,label = '{p.signif}',tip.length = 0.01) +
  stat_pvalue_manual(df_p_val3,label = '{p.signif}',tip.length = 0.01) +
  coord_flip() +
  theme_classic() 

ggsave(plot_genus, filename = "1.result/genus_13_plot.pdf", height = 10, width = 8)

### 0,1
genus_profile_01 <- genus_profile_sum
rownames(genus_profile_01) <- genus_profile_01$Group.1
genus_profile_01 <- genus_profile_01[,-1]
genus_profile_01[genus_profile_01>0] <- 1
genus_profile_01_str <- as.data.frame(t(genus_profile_01))
genus_profile_01_str$Sample <- rownames(genus_profile_01_str)
genus_profile_sample <- merge(genus_profile_01_str, Sample, by.x = "Sample", by.y = "Sample_ID")
genus_profile_sample <- genus_profile_sample[,-1]
genus_profile_sample_sum <- aggregate(genus_profile_sample[, -960], by=list(genus_profile_sample$Host), sum)
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
