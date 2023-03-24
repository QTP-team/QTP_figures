# Supplementary Fig. 1. The effective estimation of co-binned MAGs and SGBs using all 79 Tibetan horse samples.
library(tidyverse)
library(ggpubr)
library(reshape2)
library(UpSetR)

## a. The effect of co-binned sample size on the total MAG number
dset.a <- read_tsv("0.data/dset.figA.tsv")

dset.a.plot <- melt(dset.a, id.vars = "sample_num", variable.name = "Type", value.name = "num")
dset.a.plot$sample_num <- factor(dset.a.plot$sample_num, levels = c("s1","s4","s6","s8","s10", "s15", "s20"))
dset.a.plot$Type <- factor(dset.a.plot$Type, levels = c("HQ_MAGs", "HMQ_MAGs", "Bins"))

plot.a <- ggplot(dset.a.plot, aes(x = sample_num, y = num, fill = Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=num), vjust=-0.3, size=2, position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30000)) +
  scale_fill_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A")) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = c(0,1), legend.justification = c(0, 1))

ggsave("1.result/FigS1A.MAGs.pdf", plot.a, width = 6, height = 4)

## b. The effect of co-binned sample size on the QS50 MAGs per sample
dset.b <- read_csv("0.data/dset.figB.csv")
plot.b <- ggpaired(dset.b, x = "Sample_n", y = "HMQ_MAGs", id = 'Sample_ID', color = "Sample_n", line.color = "gray90") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350), breaks = seq(0,350,50)) + 
  theme_classic(base_family = "Helvetica") + theme(legend.position = "none")

ggsave("1.result/FigS1B.QS50_per_Sample.pdf", plot.b, width = 6, height = 4)

## c. Estimation of effective SGB improvement in different sample size
dset.c <- read_csv("0.data/dset.figC.csv")
dset.c$Sample_n <- factor(dset.c$Sample_n, levels = c("s1","s4","s6","s8","s10","s15","s20"))
dset.c$Type <- factor(dset.c$Type, levels = c("MQ_SGBs_num", "HQ_SGBs_num"))

plot.c <- ggplot(dset.c, aes(x = Sample_n, y = num, fill = Type)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=label_ypos, label=num), vjust=1.6, color="white", size=3.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1800), breaks = seq(0,1800,600)) + 
  scale_fill_manual(values = c("#FDBF6F", "#FF7F00")) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = c(0,1), legend.justification = c(0, 1))

ggsave("1.result/FigS1C.SGBs.pdf", plot.c, width = 6, height = 4)

## d. Overlapping of effective SGBs obtained from different sample size settings
dset.d <- read.csv("0.data/dset.figD.csv",header=TRUE,row.names=1)
sample_n_lst = c("s1","s4","s6","s8","s10","s15","s20")

plot.d <- upset(dset.d, sets = sample_n_lst, mb.ratio = c(0.55, 0.45), 
      order.by = c("freq"),
      queries = list(list(query = intersects, params = list("s1","s4","s6","s8","s10","s15","s20"), color = "#1F78B4", active = T)
      ))

