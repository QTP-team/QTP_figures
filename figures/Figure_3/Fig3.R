#Fig.3 Gut microbial diversity features of six animals host species
library(tidyverse)
library(ggpubr)
library(vegan)
library(ade4)
library(ggtree)
library(readr)
library(cowplot)
library(aplot)

## a. Alpha diversity
cal_shannon_richness <- function(infile){
  #infile <- "0.data/02.profile/species/Ass.s48.SGBs.Rel_Ab.profile.csv"
  profile.dset <- read.csv(infile, header = T, row.names = 1, stringsAsFactors=FALSE, check.names = F)
  shannon.res <- diversity(profile.dset, index = "shannon")
  simpson.res <- diversity(profile.dset, index = "simpson")
  
  res.dset <- tibble(Sample_ID = rownames(profile.dset), Richness = Reduce(`+`, as.data.frame(profile.dset > 0))) %>%
    left_join(tibble(Sample_ID = names(shannon.res), Shannon_index = shannon.res), by = "Sample_ID") %>%
    left_join(tibble(Sample_ID = names(simpson.res), Simpson_index = simpson.res), by = "Sample_ID")
  return(res.dset) 
}

d.res <- cal_shannon_richness("0.data/02.profile/SGBs_14062_profile.csv")
sample.dset <- read_csv("0.data/01.summary/TableS0_samples_info.csv")
sample.dset <- sample.dset[,c(2,3)]
meta.alpha <- left_join(sample.dset, d.res, by = "Sample_ID") %>%
  pivot_longer(Richness : Simpson_index,
               names_to = "Diversity",
               values_to = "value")

meta.alpha$Host[meta.alpha$Host=="Tibetan_Ass"] <- "TA"
meta.alpha$Host[meta.alpha$Host=="Tibetan_Horse"] <- "TH"
meta.alpha$Host[meta.alpha$Host=="Tibetan_Cattle"] <- "TC"
meta.alpha$Host[meta.alpha$Host=="Tibetan_Yak"] <- "Yak"
meta.alpha$Host[meta.alpha$Host=="Tibetan_Sheep"] <- "TS"
meta.alpha$Host[meta.alpha$Host=="Tibetan_Antelope"] <- "TAN"
meta.alpha$Host <- factor(meta.alpha$Host, levels = c("TA", "TH", "TC", "Yak", "TS", "TAN"))

### host alpha
host_comparisons <- list(c("TA", "TH"), c("TC", "Yak"), c("TS", "TAN"))

richness.plot <- ggviolin(filter(meta.alpha, Diversity == "Richness"), x = "Host", y = "value", fill = "Host") +
  stat_compare_means(comparisons = host_comparisons,label.y = 2000) +
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("Richness") +
  xlab("")

ggsave("1.result/fig3A.richness.pdf", richness.plot, width = 5, height = 3)

Shannon.plot <- ggviolin(filter(meta.alpha, Diversity == "Shannon_index"), x = "Host", y = "value", fill = "Host") +
  stat_compare_means(comparisons = host_comparisons, label.y = 7) +
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("Shannon index") +
  xlab("")
ggsave("1.result/fig3A.Shannon.pdf", Shannon.plot, width = 5, height = 3)

## c. Beta diversity
### all hosts
merge_profile <- function(profile.dir, file.suffix){
  RA.file.lst = list()
  n <-  1
  for (RA.file in dir(profile.dir)) {
    if (endsWith(RA.file, file.suffix)) {
      RA.file.lst[n] <- RA.file
      n = n + 1
    }
  }
  
  RA.file.lst = unlist(RA.file.lst)
  
  all.merge.RA <- read_csv(paste0(profile.dir, RA.file.lst[1]))
  
  for (RA.file in RA.file.lst[-1]){
    tmp.RA <- read_csv(paste0(profile.dir, RA.file))
    all.merge.RA <- full_join(all.merge.RA, tmp.RA)
  }
  
  all.merge.RA[is.na(all.merge.RA)] <- 0 
  
  return(all.merge.RA)
}

p.merge.dset <- merge_profile("0.data/02.profile/", "profile.csv")
all.profile.mat <- as.matrix(p.merge.dset[-1])
rownames(all.profile.mat) <- p.merge.dset$Sample
all.bc <- vegdist(all.profile.mat, method="bray")
all.bc.pcoa <- dudi.pco(all.bc, scann = FALSE, nf = 2) 
pcoa_eig <- (all.bc.pcoa$eig)[1:2] / sum(all.bc.pcoa$eig) #解释度
data <- all.bc.pcoa$li
data[3] <- p.merge.dset[[1]]
colnames(data) <- c("PCoA1", "PCoA2", "Sample_ID")
host <- read_csv("0.data/01.summary/TableS0_samples_info.csv")
merge_data <- merge(data,host,by = "Sample_ID")
merge_data$Host[merge_data$Host=="Tibetan_Ass"] <- "TA"
merge_data$Host[merge_data$Host=="Tibetan_Horse"] <- "TH"
merge_data$Host[merge_data$Host=="Tibetan_Cattle"] <- "TC"
merge_data$Host[merge_data$Host=="Tibetan_Yak"] <- "Yak"
merge_data$Host[merge_data$Host=="Tibetan_Sheep"] <- "TS"
merge_data$Host[merge_data$Host=="Tibetan_Antelope"] <- "TAN"
#write.csv(merge_data, file = "0.data/03.pcoa/all_sample_pcoa.csv")
host.pcoa.p <- merge_data %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Host, fill = Host)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                              "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) + 
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 

ggsave("1.result/fig3C.host_pcoa.pdf", host.pcoa.p, width = 5, height = 5)

# Ass.Horse
Ass.Horse.sample <- merge_data[which(merge_data$Host=="TA"|merge_data$Host=="TH"),1]
Ass.Horse.profile.dset <- as.data.frame(p.merge.dset[which(p.merge.dset$Sample%in%Ass.Horse.sample),])
Ass.Horse.profile.mat <- as.matrix(Ass.Horse.profile.dset[-1])
rownames(Ass.Horse.profile.mat) <- Ass.Horse.profile.dset$Sample
Ass.Horse.bc <- vegdist(Ass.Horse.profile.mat, method="bray")
Ass.Horse.bc.pcoa <- dudi.pco(Ass.Horse.bc, scann = FALSE, nf = 2) 
Ass.Horse_pcoa_eig <- (Ass.Horse.bc.pcoa$eig)[1:2] / sum(Ass.Horse.bc.pcoa$eig)
Ass.Horse.data <- Ass.Horse.bc.pcoa$li
Ass.Horse.data[3] <- Ass.Horse.profile.dset[[1]]
colnames(Ass.Horse.data) <- c("PCoA1", "PCoA2", "Sample_ID")
Ass.Horse.merge_data <- merge(Ass.Horse.data,host,by = "Sample_ID")
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Ass"] <- "TA"
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Horse"] <- "TH"
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Cattle"] <- "TC"
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Yak"] <- "Yak"
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Sheep"] <- "TS"
Ass.Horse.merge_data$Host[Ass.Horse.merge_data$Host=="Tibetan_Antelope"] <- "TAN"
Ass.Horse.host.pcoa.p <- Ass.Horse.merge_data %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Host, fill = Host)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                              "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) + 
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 
ggsave("1.result/fig3C.Ass.Horse.host_pcoa.pdf", Ass.Horse.host.pcoa.p, width = 5, height = 5)

# Antelope.Sheep
Antelope.Sheep.sample <- merge_data[which(merge_data$Host=="TS"|merge_data$Host=="TAN"),1]
Antelope.Sheep.profile.dset <- as.data.frame(p.merge.dset[which(p.merge.dset$Sample%in%Antelope.Sheep.sample),])
Antelope.Sheep.profile.mat <- as.matrix(Antelope.Sheep.profile.dset[-1])
rownames(Antelope.Sheep.profile.mat) <- Antelope.Sheep.profile.dset$Sample
Antelope.Sheep.bc <- vegdist(Antelope.Sheep.profile.mat, method="bray")
Antelope.Sheep.bc.pcoa <- dudi.pco(Antelope.Sheep.bc, scann = FALSE, nf = 2) 
Antelope.Sheep_pcoa_eig <- (Antelope.Sheep.bc.pcoa$eig)[1:2] / sum(Antelope.Sheep.bc.pcoa$eig)
Antelope.Sheep.data <- Antelope.Sheep.bc.pcoa$li
Antelope.Sheep.data[3] <- Antelope.Sheep.profile.dset[[1]]
colnames(Antelope.Sheep.data) <- c("PCoA1", "PCoA2", "Sample_ID")
Antelope.Sheep.merge_data <- merge(Antelope.Sheep.data,host,by = "Sample_ID")
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Ass"] <- "TA"
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Horse"] <- "TH"
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Cattle"] <- "TC"
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Yak"] <- "Yak"
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Sheep"] <- "TS"
Antelope.Sheep.merge_data$Host[Antelope.Sheep.merge_data$Host=="Tibetan_Antelope"] <- "TAN"
Antelope.Sheep.host.pcoa.p <- Antelope.Sheep.merge_data %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Host, fill = Host)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                              "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) + 
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 
ggsave("1.result/fig3C.Antelope.Sheep.host_pcoa.pdf", Antelope.Sheep.host.pcoa.p, width = 5, height = 5)

# Cattle.Yak
Cattle.Yak.sample <- merge_data[which(merge_data$Host=="TC"|merge_data$Host=="Yak"),1]
Cattle.Yak.profile.dset <- as.data.frame(p.merge.dset[which(p.merge.dset$Sample%in%Cattle.Yak.sample),])
Cattle.Yak.profile.mat <- as.matrix(Cattle.Yak.profile.dset[-1])
rownames(Cattle.Yak.profile.mat) <- Cattle.Yak.profile.dset$Sample
Cattle.Yak.bc <- vegdist(Cattle.Yak.profile.mat, method="bray")
Cattle.Yak.bc.pcoa <- dudi.pco(Cattle.Yak.bc, scann = FALSE, nf = 2) 
Cattle.Yak_pcoa_eig <- (Cattle.Yak.bc.pcoa$eig)[1:2] / sum(Cattle.Yak.bc.pcoa$eig)
Cattle.Yak.data <- Cattle.Yak.bc.pcoa$li
Cattle.Yak.data[3] <- Cattle.Yak.profile.dset[[1]]
colnames(Cattle.Yak.data) <- c("PCoA1", "PCoA2", "Sample_ID")
Cattle.Yak.merge_data <- merge(Cattle.Yak.data,host,by = "Sample_ID")
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Ass"] <- "TA"
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Horse"] <- "TH"
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Cattle"] <- "TC"
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Yak"] <- "Yak"
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Sheep"] <- "TS"
Cattle.Yak.merge_data$Host[Cattle.Yak.merge_data$Host=="Tibetan_Antelope"] <- "TAN"
Cattle.Yak.host.pcoa.p <- Cattle.Yak.merge_data %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Host, fill = Host)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                              "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) + 
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 
ggsave("1.result/fig3C.Cattle.Yak.host_pcoa.pdf", Cattle.Yak.host.pcoa.p, width = 5, height = 5)
#pcoa.Ass.Horse.species <- merge_data[which(merge_data$Host=="TA"|merge_data$Host=="TH"),c(1,2,3,5)]
#write.csv(pcoa.Ass.Horse.species, file = "../Figure S4.diversity_3/0.data/pcoa.Ass.Horse.speciese_data.csv")
#pcoa.Sheep.Antelope.species <- merge_data[which(merge_data$Host=="TS"|merge_data$Host=="TAN"),c(1,2,3,5)]
#write.csv(pcoa.Sheep.Antelope.species, file = "../Figure S4.diversity_3/0.data/pcoa.Sheep.Antelope.speciese_data.csv")
#pcoa.Cattle.Yak.species <- merge_data[which(merge_data$Host=="TC"|merge_data$Host=="Yak"),c(1,2,3,5)]
#write.csv(pcoa.Cattle.Yak.species, file = "../Figure S4.diversity_3/0.data/pcoa.Cattle.Yak.speciese_data.csv")

## adonis2
# all
set.seed(1)
dune.div <- adonis2(all.profile.mat ~ Host, data = host, permutations = 999, method="bray")
# TA TH
set.seed(1)
Host_TA_TH <- host[which(host$Host=="Tibetan_Ass" | host$Host=="Tibetan_Horse"),]
p.merge.dset_TA_TH <- p.merge.dset[which(p.merge.dset$Sample%in%Host_TA_TH$Sample_ID),]
all.profile.mat_TA_TH <- as.matrix(p.merge.dset_TA_TH[-1])
rownames(all.profile.mat_TA_TH) <- p.merge.dset_TA_TH$Sample
dune.div.TA_TH <- adonis2(all.profile.mat_TA_TH ~ Host, data = Host_TA_TH, permutations = 999, method="bray")
# TS TAN
set.seed(1)
Host_TS_TAN <- host[which(host$Host=="Tibetan_Sheep" | host$Host=="Tibetan_Antelope"),]
p.merge.dset_TS_TAN <- p.merge.dset[which(p.merge.dset$Sample%in%Host_TS_TAN$Sample_ID),]
all.profile.mat_TS_TAN <- as.matrix(p.merge.dset_TS_TAN[-1])
rownames(all.profile.mat_TS_TAN) <- p.merge.dset_TS_TAN$Sample
dune.div.TS.TAN <- adonis2(all.profile.mat_TS_TAN ~ Host, data = Host_TS_TAN, permutations = 999, method="bray")
# TC Yak
set.seed(1)
Host_TC_Yak <- host[which(host$Host=="Tibetan_Cattle" | host$Host=="Tibetan_Yak"),]
p.merge.dset_TC_Yak <- p.merge.dset[which(p.merge.dset$Sample%in%Host_TC_Yak$Sample_ID),]
all.profile.mat_TC_Yak <- as.matrix(p.merge.dset_TC_Yak[-1])
rownames(all.profile.mat_TC_Yak) <- p.merge.dset_TC_Yak$Sample
dune.div.TC.Yak <- adonis2(all.profile.mat_TC_Yak ~ Host, data = Host_TC_Yak, permutations = 999, method="bray")

## b. Pairwise Bray-Curtis dissimilarity between pairs of samples within the host taxonomic rank
all.bc.mat <- as.matrix(all.bc)

pair.bc.dist <- function(bc.mat, taxon_group, type_name){
  sample.set <- tibble(Host = taxon_group) %>% left_join(host) %>% .$Sample_ID
  bc.mat.sel <- bc.mat[sample.set, sample.set]
  bc.mat.sel[upper.tri(bc.mat.sel, diag = T)] = 0
  
  res.dset <- as_tibble(bc.mat.sel)
  res.dset$sample1 <- names(res.dset)
  res.dset <- res.dset %>% pivot_longer(-sample1, names_to = "sample2", values_to = "pair_bc") %>%
    filter(pair_bc != 0) %>%
    mutate(type = type_name)
  return(res.dset)
}
### inner_species
horse <- pair.bc.dist(all.bc.mat, "Tibetan_Horse", "inner_species")
ass <- pair.bc.dist(all.bc.mat, "Tibetan_Ass", "inner_species")
cattle <- pair.bc.dist(all.bc.mat, "Tibetan_Cattle", "inner_species")
yak <- pair.bc.dist(all.bc.mat, "Tibetan_Yak", "inner_species")
antelope <- pair.bc.dist(all.bc.mat, "Tibetan_Antelope", "inner_species")
sheep <- pair.bc.dist(all.bc.mat, "Tibetan_Sheep", "inner_species")

inner_species <- bind_rows(horse, ass, cattle, yak, antelope, sheep)

### plot species
horse <- pair.bc.dist(all.bc.mat, "Tibetan_Horse", "TH")
ass <- pair.bc.dist(all.bc.mat, "Tibetan_Ass", "TA")
cattle <- pair.bc.dist(all.bc.mat, "Tibetan_Cattle", "TC")
yak <- pair.bc.dist(all.bc.mat, "Tibetan_Yak", "Yak")
antelope <- pair.bc.dist(all.bc.mat, "Tibetan_Antelope", "TAN")
sheep <- pair.bc.dist(all.bc.mat, "Tibetan_Sheep", "TS")

inner_species <- bind_rows(horse, ass, cattle, yak, antelope, sheep)
inner_species$type <- factor(inner_species$type, levels = c("TA", "TH", "TC", "Yak", "TS", "TAN"))
p_species <- ggviolin(inner_species, x = "type", y = "pair_bc", fill = "type") +
  stat_compare_means(comparisons = list(c("TA", "TH"), c("TC", "Yak"), c("TS", "TAN")), label.y = 1.2) +
  scale_fill_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                             "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("Pairwise BC dissimilarity") +
  xlab("")
ggsave(p_species, filename = "1.result/fig3B.species_bc.pdf", width = 5, height = 3)
### inner_genues
Bos <- pair.bc.dist(all.bc.mat, c("Tibetan_Cattle", "Tibetan_Yak"), "inner_genues")
Equus <- pair.bc.dist(all.bc.mat, c("Tibetan_Ass", "Tibetan_Horse"), "inner_genues")
Ovis <- pair.bc.dist(all.bc.mat, "Tibetan_Sheep", "inner_genues")
Pantholops <- pair.bc.dist(all.bc.mat, "Tibetan_Antelope", "inner_genues")

inner_genues <- bind_rows(Bos, Equus, Ovis, Pantholops)

### inner_family
Bovidae <- pair.bc.dist(all.bc.mat, c("Tibetan_Cattle", "Tibetan_Yak", "Tibetan_Sheep", "Tibetan_Antelope"), "inner_family")
Equidae <- pair.bc.dist(all.bc.mat, c("Tibetan_Ass", "Tibetan_Horse"), "inner_family")
inner_family <- bind_rows(Bovidae, Equidae)

### inner_class
Mammalia <- pair.bc.dist(all.bc.mat, c("Tibetan_Cattle", "Tibetan_Yak", "Tibetan_Sheep", "Tibetan_Antelope", "Tibetan_Ass", "Tibetan_Horse"), "inner_class")

### plot inner
inner.dset <- rbind(inner_species, inner_genues, inner_family, Mammalia)
inner.dset$type <- factor(inner.dset$type, levels = c("inner_species", "inner_genues", "inner_family", "inner_class"))

my_comparisons <- list(c("inner_species", "inner_genues"), c("inner_genues", "inner_family"), c("inner_family", "inner_class"))
#ggviolin(inner.dset, x = "type", y = "pair_bc", fill = "type") + stat_compare_means(comparisons = my_comparisons)

### plot inter
inner_genues_inter <- inner_genues %>% anti_join(inner_species, by = c("sample1", "sample2"))
inner_family_inter <- inner_family %>% anti_join(inner_genues, by = c("sample1", "sample2"))
inner_class_inter <- Mammalia %>% anti_join(inner_family, by = c("sample1", "sample2"))

inter.dset <- rbind(inner_species, inner_genues_inter, inner_family_inter, inner_class_inter)
inner.p <- ggboxplot(inter.dset, x = "type", y = "pair_bc", fill = "type", outlier.size = 0.5) + 
  stat_compare_means(comparisons = my_comparisons) + 
  theme(legend.position="none")+
  xlab("")
ggsave(plot = inner.p, "1.result/fig3B.inner_inter_bc.pdf", width = 5, height = 3)

## d.all ggtree
### all
all.hc <- hclust(all.bc, method="complete")

### anno
meta.dset <- read_csv("0.data/01.summary/TableS0_samples_info.csv")

hc.anno <- tibble(Sample_ID = all.hc$labels, order = all.hc$order) %>%
  arrange(order) %>%
  left_join(meta.dset, by = "Sample_ID")

hc.anno.otu <- list(TAN = hc.anno %>% filter(Host == "Tibetan_Antelope") %>% .$Sample_ID,
                    TA = hc.anno %>% filter(Host == "Tibetan_Ass") %>% .$Sample_ID,
                    TC = hc.anno %>% filter(Host == "Tibetan_Cattle") %>% .$Sample_ID,
                    TH = hc.anno %>% filter(Host == "Tibetan_Horse") %>% .$Sample_ID,
                    TS = hc.anno %>% filter(Host == "Tibetan_Sheep") %>% .$Sample_ID,
                    Yak = hc.anno %>% filter(Host == "Tibetan_Yak") %>% .$Sample_ID)

p.hc <- ggtree(all.hc, branch.length='none') + geom_rootedge(rootedge = 0.5)
p.hc.anno <- groupOTU(p.hc, hc.anno.otu, 'Host') + aes(color = Host) + scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B", "TC" = "#00FF00",
                                                                                                   "Yak" = "#336600", "TS" = "#33FFCC", "TAN" = "#0000CD")) #+ geom_label(aes(label=node))
p.hc.anno.done <- rotate(p.hc.anno, 1417) %>% rotate(1488) %>% rotate(1550) %>% rotate(1431) %>% rotate(1445) %>% rotate(1416) %>% rotate(1417) %>% rotate(1420) %>% rotate(1440) %>% rotate(1446)
#ggsave("1.result/fig3D.tree.tmp.pdf", p.hc.anno.done, width = 20, height = 50,limitsize = FALSE)
ggsave("1.result/fig3D.tree.pdf", p.hc.anno.done, width = 6, height = 10)

