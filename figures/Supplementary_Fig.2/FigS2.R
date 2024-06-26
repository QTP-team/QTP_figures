#Fig.S2 Quality assessing of 14062 SGBs
library(tidyverse)
library(reshape2)

## A. The completeness and contamination of SGBs
SGB.stat <- read_csv("0.data/SGBs_14062_quality.csv")

SGB.stat$Quality <- factor(SGB.stat$Quality, levels = c("High", "Medium"))
cont.comp.plot <- ggplot(SGB.stat, aes(x = completeness, y = contamination, color = Quality)) +
  geom_point(alpha = 0.5) +
  ylim(0,10) +
  theme_classic() +
  guides(color = F) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
  scale_x_continuous(expand = c(0, 0), limits = c(50, 100))

ggsave("1.result/figS2A.cont_comp_plot.pdf", width = 5.5, height = 5)

cont.bar.plot <- ggplot(SGB.stat, aes(y = contamination, fill = Quality)) +
  geom_histogram(bins = 40) +
  theme_classic() +
  guides(fill = F)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2000))
ggsave("1.result/figS2A.cont_plot.pdf", cont.bar.plot, width = 2, height = 5)

comp.bar.plot <- ggplot(SGB.stat, aes(x = completeness, fill = Quality)) +
  geom_histogram(bins = 40) +
  theme_classic() +
  scale_y_continuous(position = "right", expand = c(0, 0), limits = c(0, 1500)) +
  guides(fill = F)
ggsave("1.result/fig2A.comp_plot.pdf", comp.bar.plot, width = 5.5, height = 1.5)

## B. Genome sizes
Gszie.plot <- ggplot(SGB.stat, aes(x = Quality, y = length/1000000, fill = Quality)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_y_log10() +
  ylab("Genome size (Mbp)") +
  guides(fill = F) +
  coord_flip() 
ggsave("1.result/figS2B.Gsize_plot.pdf", Gszie.plot, width = 5.5, height = 1)

## C. N50 statistics
N50.plot <- ggplot(SGB.stat, aes(x = Quality, y = N50/1000, fill = Quality)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_y_log10() +
  ylab("Contig N50 (Kbp)") +
  guides(fill = F) +
  coord_flip()
ggsave("1.result/figS2C.N50_plot.pdf", N50.plot, width = 5.5, height = 1)

## D. Distribution of MAG numbers supporting SGBs.

plot.d.dset <- SGB.stat %>% 
  group_by(MAGs_in_SGB, Quality) %>%
  summarise(count = n()) %>%
  ungroup()

plot.d.5 <- filter(plot.d.dset, MAGs_in_SGB <= 5) %>%
  mutate(MAGs_in_SGB = as.character(MAGs_in_SGB))

plot.d.6_10 <- filter(plot.d.dset, MAGs_in_SGB <= 10 & MAGs_in_SGB >= 6) %>%
  group_by(Quality) %>%
  summarise(count = sum(count),
            MAGs_in_SGB = "6-10")

plot.d.11_20 <- filter(plot.d.dset, MAGs_in_SGB <= 20 & MAGs_in_SGB >= 11) %>%
  group_by(Quality) %>%
  summarise(count = sum(count),
            MAGs_in_SGB = "11-20")

plot.d.21_50 <- filter(plot.d.dset, MAGs_in_SGB <= 50 & MAGs_in_SGB >= 21) %>%
  group_by(Quality) %>%
  summarise(count = sum(count),
            MAGs_in_SGB = "21-50")

plot.d.51 <- filter(plot.d.dset, MAGs_in_SGB >= 51) %>%
  group_by(Quality) %>%
  summarise(count = sum(count),
            MAGs_in_SGB = ">51")

plot.d <- bind_rows(plot.d.5, plot.d.6_10, plot.d.11_20, plot.d.21_50, plot.d.51) %>%
  mutate(MAGs_in_SGB = as_factor(MAGs_in_SGB)) %>%
  ggplot(aes(x = MAGs_in_SGB, y = count, fill = Quality)) +
  geom_bar(stat="identity") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000))

ggsave("1.result/figS2D.SGBs_size.pdf", plot.d,  width = 5.5, height = 4)
