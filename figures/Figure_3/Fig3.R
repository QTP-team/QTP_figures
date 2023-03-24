# Figure 3. Gut microbial diversity features of the six animals host species.
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
  profile.dset <- t(read.csv(infile, header = T, row.names = 1, stringsAsFactors=FALSE, check.names = F))
  shannon.res <- diversity(profile.dset, index = "shannon")
  simpson.res <- diversity(profile.dset, index = "simpson")
  
  res.dset <- tibble(Sample_ID = rownames(profile.dset), Richness = Reduce(`+`, as.data.frame(profile.dset > 0))) %>%
    left_join(tibble(Sample_ID = names(shannon.res), Shannon_index = shannon.res), by = "Sample_ID") %>%
    left_join(tibble(Sample_ID = names(simpson.res), Simpson_index = simpson.res), by = "Sample_ID")
  return(res.dset) 
}

alpha_diversity_lst = list()
n = 1
for (profile in dir("0.data/02.profile/species/")) {
  if (endsWith(profile, "Rel_Ab.profile.csv")) {
    d.res <- cal_shannon_richness(paste0("0.data/02.profile/species/", profile))
    alpha_diversity_lst[[n]] <-  d.res
    n = n +1
  }
}

Tibetan.alpha_diversity.dset <- bind_rows(alpha_diversity_lst)
write_csv(Tibetan.alpha_diversity.dset, "Tibetan.alpha_diversity.csv")

alpha_diversity.dset <- read_csv("Tibetan.alpha_diversity.csv")

meta.dset <- read_csv("0.data/01.summary/TableS0_samples_info.csv")

meta.alpha <- left_join(meta.dset, alpha_diversity.dset, by = "Sample_ID") %>%
  pivot_longer(Richness : Simpson_index,
               names_to = "Diversity",
               values_to = "value")

meta.alpha$Host <- factor(meta.alpha$Host, levels = c("Tibetan_Ass", "Tibetan_Horse", "Tibetan_Cattle", "Tibetan_Yak", "Tibetan_Sheep", "Tibetan_Antelope"))

### host alpha
host_comparisons <- list(c("Tibetan_Ass", "Tibetan_Horse"), c("Tibetan_Cattle", "Tibetan_Yak"), c("Tibetan_Sheep", "Tibetan_Antelope"))

ggboxplot(inter.dset, x = "type", y = "pair_bc", fill = "type") + stat_compare_means(comparisons = host_comparisons)

richness.plot <- ggviolin(filter(meta.alpha, Diversity == "Richness"), x = "Host", y = "value", fill = "Host") +
  stat_compare_means(comparisons = host_comparisons) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  ylab("Richness")

ggsave("1.result/Fig3A.richness.pdf", richness.plot, width = 5.5, height = 3)

Shannon.plot <- ggviolin(filter(meta.alpha, Diversity == "Shannon_index"), x = "Host", y = "value", fill = "Host") +
  stat_compare_means(comparisons = host_comparisons) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position="none") +
  ylab("Shannon index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("1.result/Fig3A.Shannon.pdf", Shannon.plot, width = 5.5, height = 4)

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
ggviolin(inner.dset, x = "type", y = "pair_bc", fill = "type") + stat_compare_means(comparisons = my_comparisons)

### plot inter
inner_genues_inter <- inner_genues %>% anti_join(inner_species, by = c("sample1", "sample2"))
inner_family_inter <- inner_family %>% anti_join(inner_genues, by = c("sample1", "sample2"))
inner_class_inter <- Mammalia %>% anti_join(inner_family, by = c("sample1", "sample2"))

inter.dset <- rbind(inner_species, inner_genues_inter, inner_family_inter, inner_class_inter)
inner.p <- ggboxplot(inter.dset, x = "type", y = "pair_bc", fill = "type") + stat_compare_means(comparisons = my_comparisons) + theme(legend.position="none")
ggsave(plot = inner.p, "1.result/Fig3B.inner_inter_bc.pdf", width = 5.5, height = 4)

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

p.merge.dset <- merge_profile("0.data/02.profile/species/", "profile.csv")
all.profile.mat <- t(as.matrix(p.merge.dset[-1]))
colnames(all.profile.mat) <- p.merge.dset[[1]]
all.bc <- vegdist(all.profile.mat, method="bray")
all.bc.pcoa <- dudi.pco(all.bc, scann = FALSE, nf = 2) 

all.bc.pcoa <- read_csv("0.data/03.pcoa/pcoa.all.hosts.species.csv")
host <- read_csv("0.data/01.summary/TableS0_samples_info.csv")

host.pcoa.p <- all.bc.pcoa %>%
  left_join(host) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(alpha = 0.75, size = 1) +
  scale_color_brewer(palette = "Set2") + 
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("1.result/Fig3C.host_pcoa.pdf", host.pcoa.p, width = 5.5, height = 4)

## d.all ggtree
### all
all.hc <- hclust(all.bc, method="complete")

### anno
meta.dset <- read_csv("0.data/01.summary/TableS0_samples_info.csv")

hc.anno <- tibble(Sample_ID = all.hc$labels, order = all.hc$order) %>%
  arrange(order) %>%
  left_join(meta.dset, by = "Sample_ID")

hc.anno.otu <- list(Antelope = hc.anno %>% filter(Host == "Tibetan_Antelope") %>% .$Sample_ID,
                    Ass = hc.anno %>% filter(Host == "Tibetan_Ass") %>% .$Sample_ID,
                    Cattle = hc.anno %>% filter(Host == "Tibetan_Cattle") %>% .$Sample_ID,
                    Horse = hc.anno %>% filter(Host == "Tibetan_Horse") %>% .$Sample_ID,
                    Sheep = hc.anno %>% filter(Host == "Tibetan_Sheep") %>% .$Sample_ID,
                    Yak = hc.anno %>% filter(Host == "Tibetan_Yak") %>% .$Sample_ID)

p.hc <- ggtree(all.hc, branch.length='none') + geom_rootedge(rootedge = 0.5)
p.hc.anno <- groupOTU(p.hc, hc.anno.otu, 'Host') + aes(color = Host) + scale_color_manual("Host", values = c("#0000CD", "#DC143B", "#00FF00", "#FF99FF", "#33FFCC", "#336600")) #+ geom_label(aes(label=node))
p.hc.anno.done <- rotate(p.hc.anno, 1440) %>% rotate(1417) %>% rotate(1414) %>% rotate(1434) %>% rotate(1436)
ggsave("1.result/Fig3D.tree.pdf", p.hc.anno.done, width = 6, height = 10)
