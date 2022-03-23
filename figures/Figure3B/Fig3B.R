# Fig.3B Evolutionary dynamics of gut core microbiome  
library(tidyverse)
library(ggpubr)
library(ggplot2)

## b.The difference of samxple occurrence frequency between common and gained SGBs.(statistical method and p value)
dset.b <- read_csv("0.data/dset.fig3B.csv")
dset.b$Host <- factor(dset.b$Host, levels = c("Tibetan ass","Tibetan horse", "Tibetan antelope","Tibetan sheep","Yak", "Tibetan cattle"))

plot.b <- ggplot(dset.b, aes(Type, Freq)) +
  geom_boxplot(aes(color = Type), outlier.size = 0.5) +
  facet_wrap(~Host, nrow = 3) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  stat_compare_means(label = "p.signif") + 
  stat_boxplot(geom = "errorbar",width = 0.4,aes(color = Type)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.4, 1.05)) +
  theme_classic(base_family = "Helvetica") 
  #+ theme(legend.position = "none")

ggsave("1.result/fig3B.pdf", plot.b, width = 4, height = 5)
