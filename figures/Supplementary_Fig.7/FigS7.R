#Fig.S7 Gut microbial diversity features of six animals host species
library(tidyverse)
library(ggpubr)
library(vegan)
library(ade4)

### Ass vs Horse
AH <- read.csv("0.data/pcoa.Ass.Horse.speciese_data.csv")

AH.pcoa.p <- ggplot(AH, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
  scale_color_manual(values=c("TH" = "#FF99FF", "TA" = "#DC143B"))
ggsave("1.result/figS4A.AH_pcoa.pdf", AH.pcoa.p, width = 5, height = 5)

### Sheep vs Antelope
SA <- read_csv("0.data/pcoa.Sheep.Antelope.speciese_data.csv")
SA.pcoa.p <- ggplot(SA, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point(size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))+
  scale_color_manual(values=c("TS" = "#33FFCC", "TAN" = "#0000CD"))
ggsave("1.result/fiS4B.SA_pcoa.pdf", SA.pcoa.p, width = 5, height = 5)

### Cattle vs Yak
CY <- read_csv("0.data/pcoa.Cattle.Yak.speciese_data.csv")
CY.pcoa.p <- ggplot(CY, aes(x = PCoA1, y = PCoA2, color = Host)) +
  geom_point( size = 1) +
  theme_classic() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
  scale_color_manual(values=c("TC" = "#00FF00", "Yak" = "#336600"))
ggsave("1.result/fiS4C.CY_pcoa.pdf", CY.pcoa.p, width = 5, height = 5)
