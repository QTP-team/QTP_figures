#Fig.S2 Gut microbial diversity features of six animals host species
library(tidyverse)
library(ggpubr)
library(vegan)
library(ade4)

### Ass vs Horse
AH <- read_csv("0.data/pcoa.Ass.Horse.species.csv")

AH.pcoa.p <- ggplot(AH, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(alpha = 0.75, size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("1.result/figS2A.AH_pcoa.pdf", AH.pcoa.p, width = 5.5, height = 4)

### Sheep vs Antelope
SA <- read_csv("0.data/pcoa.Sheep.Antelope.species.csv")
SA.pcoa.p <- ggplot(SA, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(alpha = 0.75, size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("1.result/figS2B.SA_pcoa.pdf", SA.pcoa.p, width = 5.5, height = 4)

### Cattle vs Yak
CY <- read_csv("0.data/pcoa.Cattle.Yak.species.csv")
CY.pcoa.p <- ggplot(CY, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(alpha = 0.75, size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave("1.result/figS2C.CY_pcoa.pdf", CY.pcoa.p, width = 5.5, height = 4)