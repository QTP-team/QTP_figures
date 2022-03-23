# Fig.1 The characteristics of Plateau Animal Microbiome Project (PAMP) database
library(tidyverse)
library(ggplot2)

## a.A rarefaction curve for the assessment of SGB obtained.
dset.a <- read_csv("0.data/dset.fig1A.csv")

plot.a <- ggplot(dset.a, aes(x = Number_of_samples, y = Number_of_SGBs, color = taxo)) +
  geom_point() +
  geom_smooth() +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1500)) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = c(0,1), legend.justification = c(0, 1))

ggsave("1.result/fig1A.curve.pdf", plot.a, width = 3.5, height = 3)

## b.The mapping rates of our sample reads to different metagenomic datasets.
dset.b <- read_csv("0.data/dset.fig1B.csv")
dset.b$Database <- factor(dset.b$Database, levels = c("RUG","GTDB","Earth","All_SGBs"))

plot.b <- ggplot(dset.b, aes(x = Host, y = Mapping_rate, fill = Database)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = 'Host', y = 'Mapping rate(%)') + 
  scale_x_discrete(limits=c("Tibetan horse","Tibetan cattle", "Tibetan ass","Yak","Tibetan sheep","Tibetan antelope")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.key.size = unit(10, "pt"))
  

ggsave("1.result/fig1B.mapping_rate.pdf", plot.b, width = 5.5, height = 3)

## c.The Mash distance among our identified SGBs and other known SGBs released by Genome Taxonomy database (GTDB 05-RS95) and Earth Microbiome Project
dset.c <- read_csv("0.data/dset.fig1C.csv")

dist.p <- ggplot(dset.c) + 
  geom_line(aes(x = mash_dist, y = count, color = database)) +
  geom_point(aes(x = 0.05, 153), color = "#13b6bb") +
  geom_point(aes(x = 0.05, 334), color = "#f56b64") +
  geom_point(aes(x = 0.15, 431), color = "#13b6bb") +
  geom_point(aes(x = 0.15, 1331), color = "#f56b64") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3600), breaks = seq(0,3600,600)) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = c(0,1), legend.justification = c(0, 1))

ggsave("1.result//fig1C.mash.pdf", dist.p, width = 3.5, height = 3)

## d. Compared with GTDB database, assessing the  percentage rate of expanding the microbial genomic diversity at the phylum level.
dset.d <- read_csv("0.data/dset.fig1D.csv") %>% 
  filter(improve_rate > 20 & TP_freq > 10)

dset.d$GTDB_phylum <- factor(dset.d$GTDB_phylum, levels = dset.d$GTDB_phylum)
uc.p <- ggplot(dset.d, aes(x = improve_rate, y = GTDB_phylum, fill = GTDB_phylum)) + 
  geom_col() + 
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 550), breaks = seq(0, 550, 100)) + 
  theme_classic(base_size = 7, base_family = "Helvetica") + 
  xlab("Improve rate(%)") + 
  guides(fill = FALSE)

ggsave("1.result/fig1D.Imporve_GTDB.pdf", uc.p, width = 3.5, height = 3)
