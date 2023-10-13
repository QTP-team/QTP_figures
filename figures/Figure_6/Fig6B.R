library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(cowplot)
data <- read.csv("0.data/genome_coverage.csv")
ggplot(data, aes(x = Coverage, fill = Group)) +
  geom_histogram(position = "identity", alpha = 0.4)

ggplot(data, aes(x = Coverage, fill = Group)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  scale_x_continuous(limits = c(0, 5))

data_filter <- data[which(data$Coverage>=3),]
sample_data <- read.csv("0.data/TableS0_samples_info.csv")
sample_data <- sample_data[,c(3,2)]
sample_data$Host[sample_data$Host=="Tibetan_Ass"] <- "TA"
sample_data$Host[sample_data$Host=="Tibetan_Horse"] <- "TH"
sample_data$Host[sample_data$Host=="Tibetan_Cattle"] <- "TC"
sample_data$Host[sample_data$Host=="Tibetan_Yak"] <- "Yak"
sample_data$Host[sample_data$Host=="Tibetan_Sheep"] <- "TS"
sample_data$Host[sample_data$Host=="Tibetan_Antelope"] <- "TAN"
sample_data$Group1 <- sample_data$Host
sample_data$Group1[sample_data$Group1=="TA"] <- "Perissodactyla"
sample_data$Group1[sample_data$Group1=="TH"] <- "Perissodactyla"
sample_data$Group1[sample_data$Group1=="TAN"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="TS"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="Yak"] <- "Artiodactyla"
sample_data$Group1[sample_data$Group1=="TC"] <- "Artiodactyla"
sample_data$Group2 <- sample_data$Host
sample_data$Group2[sample_data$Group2=="TAN"] <- "Caprinae"
sample_data$Group2[sample_data$Group2=="TS"] <- "Caprinae"
sample_data$Group2[sample_data$Group2=="Yak"] <- "Bovinae"
sample_data$Group2[sample_data$Group2=="TC"] <- "Bovinae"
data_filter_sample <- merge(data_filter, sample_data, by.x = "Sample", by.y = "Sample_ID")
#write.csv(data_filter_sample, file = "data_filter_sample.csv")
all_gene_ko <- read.table("0.data/all.gene2KO.lst.Gene2KO.xls", sep = "\t", header = T)
all_gene_ko$SGB <- str_split_fixed(all_gene_ko$qaccver, "_k141", 2)[,1]
ko_pathway <- read.csv("ko2pathway.csv")

## data
SIG77_data_Yak <- data_filter_sample[which(data_filter_sample$Host=="Yak"&data_filter_sample$Group=="SGB_6458"),]
SIG77_data_all <- read.csv("0.data/all_sample_g_SIG77_SGB_6458_pN_pS_all.csv")
SIG77_data_all_Yak <- SIG77_data_all[which(SIG77_data_all$sample%in%SIG77_data_Yak$Sample),]
SIG77_data_filter <- SIG77_data_all[which(SIG77_data_all$coverage>=3),]
SIG77_data_filter_Yak <- SIG77_data_filter[which(SIG77_data_filter$sample%in%SIG77_data_Yak$Sample),]

Treponema_F_data_Yak <- data_filter_sample[which(data_filter_sample$Host=="Yak"&data_filter_sample$Group=="SGB_8084"),]
Treponema_F_data_all <- read.csv("0.data/all_sample_g_Treponema_F_SGB_8084_pN_pS_all.csv")
Treponema_F_data_all_Yak <- Treponema_F_data_all[which(Treponema_F_data_all$sample%in%Treponema_F_data_Yak$Sample),]
Treponema_F_data_filter <- Treponema_F_data_all[which(Treponema_F_data_all$coverage>=3),]
Treponema_F_data_filter_Yak <- Treponema_F_data_filter[which(Treponema_F_data_filter$sample%in%Treponema_F_data_Yak$Sample),]

## pN/pS medium
SIG77_data_median <- aggregate(SIG77_data_filter_Yak$pNpS_variants, by=list(SIG77_data_filter_Yak$sample), median)
colnames(SIG77_data_median) <- c("Sample", "pNpS_variants_median")
SIG77_data_median$Group <- "SGB_6458"
Treponema_F_data_median <- aggregate(Treponema_F_data_filter_Yak$pNpS_variants, by=list(Treponema_F_data_filter_Yak$sample), median)
colnames(Treponema_F_data_median) <- c("Sample", "pNpS_variants_median")
Treponema_F_data_median$Group <- "SGB_8084"

Two_SGBs_Yak_median <- rbind(SIG77_data_median,Treponema_F_data_median)
write.csv(Two_SGBs_Yak_median, file = "Yak/Two_SGBs_Yak_median.csv")
p_Two_SGBs_Yak_median <- ggplot(Two_SGBs_Yak_median,aes(x=Group,y=pNpS_variants_median,fill=Group))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  geom_signif(comparisons = list(c("SGB_6458", "SGB_8084")), y_position = 0.5,
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442"))+  #设置填充的颜色
  scale_color_manual(values=c("black","black","black"))+ #设置散点图的圆圈的颜色为黑色
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("Median value of pN/pS per sample")+xlab("") 


ggsave(p_Two_SGBs_Yak_median, filename = "Yak/Two_SGBs_Yak_median.pdf", width = 5, height = 5)

## pN/pS mean
SIG77_data_average <- aggregate(SIG77_data_filter_Yak$pNpS_variants, by=list(SIG77_data_filter_Yak$sample), mean)
colnames(SIG77_data_average) <- c("Sample", "pNpS_variants_average")
SIG77_data_average$Group <- "SGB_6458"
Treponema_F_data_average <- aggregate(Treponema_F_data_filter_Yak$pNpS_variants, by=list(Treponema_F_data_filter_Yak$sample), mean)
colnames(Treponema_F_data_average) <- c("Sample", "pNpS_variants_average")
Treponema_F_data_average$Group <- "SGB_8084"

Two_SGBs_Yak_mean <- rbind(SIG77_data_average,Treponema_F_data_average)
write.csv(Two_SGBs_Yak_mean, file = "Yak/Two_SGBs_Yak_mean.csv")
p_Two_SGBs_Yak_mean <- ggviolin(Two_SGBs_Yak_mean, x = "Group", y = "pNpS_variants_average", fill = "Group") +
  geom_boxplot(aes(fill = Group),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("SGB_6458", "SGB_8084")), y_position = 1,
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none") +
  xlab("") +
  ylab("pN/pS average value")
ggsave(p_Two_SGBs_Yak_mean, filename = "Yak/Two_SGBs_Yak_mean.pdf", width = 5, height = 5)

## pN/pS >1
SIG77_data_Yak_pN <- as.data.frame(table(SIG77_data_filter_Yak[which(SIG77_data_filter_Yak$pNpS_variants>1),1]))
colnames(SIG77_data_Yak_pN) <- c("Sample", "Gene_num")
SIG77_data_Yak_pN$Group <- "SGB_6458"
Treponema_F_data_Yak_pN <- as.data.frame(table(Treponema_F_data_filter_Yak[which(Treponema_F_data_filter_Yak$pNpS_variants>1),1]))
colnames(Treponema_F_data_Yak_pN) <- c("Sample", "Gene_num")
Treponema_F_data_Yak_pN$Group <- "SGB_8084"
Two_SGBs_Yak_pN <- rbind(SIG77_data_Yak_pN,Treponema_F_data_Yak_pN)
p_Two_SGBs_Yak_pN <- ggviolin(Two_SGBs_Yak_pN, x = "Group", y = "Gene_num", fill = "Group") +
  geom_boxplot(aes(fill = Group),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("SGB_6458", "SGB_8084")), y_position = 40,
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none") +
  xlab("") +
  ylab("pN/pS>1 gene num")
ggsave(p_Two_SGBs_Yak_pN, filename = "Yak/Two_SGBs_Yak_pN.pdf", width = 5, height = 5)

## pN/pS <1
SIG77_data_Yak_pS <- as.data.frame(table(SIG77_data_filter_Yak[which(SIG77_data_filter_Yak$pNpS_variants<1),1]))
colnames(SIG77_data_Yak_pS) <- c("Sample", "Gene_num")
SIG77_data_Yak_pS$Group <- "SGB_6458"
Treponema_F_data_Yak_pS <- as.data.frame(table(Treponema_F_data_filter_Yak[which(Treponema_F_data_filter_Yak$pNpS_variants<1),1]))
colnames(Treponema_F_data_Yak_pS) <- c("Sample", "Gene_num")
Treponema_F_data_Yak_pS$Group <- "SGB_8084"
Two_SGBs_Yak_pS <- rbind(SIG77_data_Yak_pS,Treponema_F_data_Yak_pS)
p_Two_SGBs_Yak_pS <- ggviolin(Two_SGBs_Yak_pS, x = "Group", y = "Gene_num", fill = "Group") +
  geom_boxplot(aes(fill = Group),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("SGB_6458", "SGB_8084")), y_position = 1200,
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none") +
  xlab("") +
  ylab("pN/pS<1 gene num")
ggsave(p_Two_SGBs_Yak_pS, filename = "Yak/Two_SGBs_Yak_pS.pdf", width = 5, height = 5)

## gene pN/pS medium
SIG77_data_median <- aggregate(SIG77_data_filter_Yak$pNpS_variants, by=list(SIG77_data_filter_Yak$gene), median)
colnames(SIG77_data_median) <- c("Gene", "pNpS_variants_median")
SIG77_data_median$Group <- "SGB_6458"
Treponema_F_data_median <- aggregate(Treponema_F_data_filter_Yak$pNpS_variants, by=list(Treponema_F_data_filter_Yak$gene), median)
colnames(Treponema_F_data_median) <- c("Gene", "pNpS_variants_median")
Treponema_F_data_median$Group <- "SGB_8084"

Two_SGBs_Yak_medium <- rbind(SIG77_data_median,Treponema_F_data_median)
p_Two_SGBs_Yak_medium <- ggviolin(Two_SGBs_Yak_medium, x = "Group", y = "pNpS_variants_median", fill = "Group") +
  geom_boxplot(aes(fill = Group),width=0.1, outlier.size = 0.5) +
  geom_signif(comparisons = list(c("SGB_6458", "SGB_8084")), y_position = 3,
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none") +
  xlab("") +
  ylab("pN/pS median value")
ggsave(p_Two_SGBs_Yak_medium, filename = "Yak/Two_SGBs_Yak_medium.pdf", width = 5, height = 5)

## Fisher test
# SGB_6458
# all
SIG77_data_all_ko <- all_gene_ko[which(all_gene_ko$SGB=="SGB_6458"),]
SIG77_data_all_ko <- unique(SIG77_data_all_ko[,c(1,2)])
SIG77_data_all_ko_num <- as.data.frame(table(SIG77_data_all_ko$KO))
colnames(SIG77_data_all_ko_num) <- c("KO", "SGB_6458_all")

# pN/pS >1
SIG77_data_Yak_pN_all <- SIG77_data_filter_Yak[which(SIG77_data_filter_Yak$pNpS_variants>1),]
SIG77_data_filter_pN_ko_all <- merge(SIG77_data_Yak_pN_all, all_gene_ko, by.x = "gene", by.y = "qaccver", all.x = T)
SIG77_data_filter_pN_function_all <- merge(SIG77_data_filter_pN_ko_all, ko_pathway, by.x = "KO", by.y = "ModuleID", all.x = T)
write.csv(SIG77_data_filter_pN_function_all, file = "Yak/SIG77_data_filter_pN_function_all.csv")
SIG77_data_filter_pN_ko <- merge(SIG77_data_Yak_pN_all, all_gene_ko, by.x = "gene", by.y = "qaccver")
SIG77_data_filter_pN_ko <- unique(SIG77_data_filter_pN_ko[,c(1,7)])
SIG77_data_filter_pN_ko_num <- as.data.frame(table(SIG77_data_filter_pN_ko$KO))
colnames(SIG77_data_filter_pN_ko_num) <- c("KO", "SGB_6458_pN")

KEGG_ko_merge_pN <- merge(SIG77_data_all_ko_num, SIG77_data_filter_pN_ko_num, by = "KO", all = T)
KEGG_ko_merge_pN[is.na(KEGG_ko_merge_pN)] <- 0
KEGG_ko_merge_pN$SGB_6458_all_no <- 2175-KEGG_ko_merge_pN$SGB_6458_all
KEGG_ko_merge_pN$SGB_6458_pN_no <- 181-KEGG_ko_merge_pN$SGB_6458_pN
write.csv(KEGG_ko_merge_pN, file = "Yak/KEGG_ko_merge_pN_SGB_6458.csv")
mydat<-read.xlsx("Yak/KEGG_ko_merge_pN_SGB_6458.xlsx",sheet=1,colNames = T)
result<-c()
for (i in (2:ncol(mydat))){
  newdat<-mydat[,i]
  test<-matrix(newdat,nrow=2,ncol=2,byrow=TRUE)
  qq<-fisher.test(test)   #fisher.test()   chisq.test
  Odds_ratio<-qq[["estimate"]][["odds ratio"]]
  p.result<-qq$p.value
  result.linshi<-cbind(i,Odds_ratio,p.result)
  result<-rbind(result,result.linshi)
}

KO<-names(mydat[,-1])  

result<-cbind(KO,result)
result <- as.data.frame(result)
result$FDR <- p.adjust(result$p.result,"fdr")

result$sig <- "NS"
result$sig[result$FDR<0.05]="*"
result$sig[result$FDR<0.01]="**"
result$sig[result$FDR<0.001]="***"
result$sig[result$FDR<0.0001]="****"
write.table(result, file="fisher_result_pN_SGB_6458.xls",col.names = T,row.names = F,sep="\t",quote=FALSE)

# pN/pS <1
SIG77_data_Yak_pS_all <- SIG77_data_filter_Yak[which(SIG77_data_filter_Yak$pNpS_variants<1),]
SIG77_data_filter_pS_ko_all <- merge(SIG77_data_Yak_pS_all, all_gene_ko, by.x = "gene", by.y = "qaccver", all.x = T)
SIG77_data_filter_pS_function_all <- merge(SIG77_data_filter_pS_ko_all, ko_pathway, by.x = "KO", by.y = "ModuleID", all.x = T)
write.csv(SIG77_data_filter_pS_function_all, file = "Yak/SIG77_data_filter_pS_function_all.csv")
SIG77_data_filter_pS_ko <- merge(SIG77_data_Yak_pS_all, all_gene_ko, by.x = "gene", by.y = "qaccver")
SIG77_data_filter_pS_ko <- unique(SIG77_data_filter_pS_ko[,c(1,7)])
SIG77_data_filter_pS_ko_num <- as.data.frame(table(SIG77_data_filter_pS_ko$KO))
colnames(SIG77_data_filter_pS_ko_num) <- c("KO", "SGB_6458_pS")

KEGG_ko_merge_pS <- merge(SIG77_data_all_ko_num, SIG77_data_filter_pS_ko_num, by = "KO", all = T)
KEGG_ko_merge_pS[is.na(KEGG_ko_merge_pS)] <- 0
KEGG_ko_merge_pS$SGB_6458_all_no <- 2175-KEGG_ko_merge_pS$SGB_6458_all
KEGG_ko_merge_pS$SGB_6458_pS_no <- 1754-KEGG_ko_merge_pS$SGB_6458_pS
write.csv(KEGG_ko_merge_pS, file = "Yak/KEGG_ko_merge_pS_SGB_6458.csv")
mydat<-read.xlsx("Yak/KEGG_ko_merge_pS_SGB_6458.xlsx",sheet=1,colNames = T)
result<-c()
for (i in (2:ncol(mydat))){
  newdat<-mydat[,i]
  test<-matrix(newdat,nrow=2,ncol=2,byrow=TRUE)
  qq<-fisher.test(test)   #fisher.test()   chisq.test
  Odds_ratio<-qq[["estimate"]][["odds ratio"]]
  p.result<-qq$p.value
  result.linshi<-cbind(i,Odds_ratio,p.result)
  result<-rbind(result,result.linshi)
}

KO<-names(mydat[,-1])  

result<-cbind(KO,result)
result <- as.data.frame(result)
result$FDR <- p.adjust(result$p.result,"fdr")

result$sig <- "NS"
result$sig[result$FDR<0.05]="*"
result$sig[result$FDR<0.01]="**"
result$sig[result$FDR<0.001]="***"
result$sig[result$FDR<0.0001]="****"
write.table(result, file="fisher_result_pS_SGB_6458.xls",col.names = T,row.names = F,sep="\t",quote=FALSE)

# SGB_8084
# all
Treponema_F_data_all_ko <- all_gene_ko[which(all_gene_ko$SGB=="SGB_8084"),]
Treponema_F_data_all_ko <- unique(Treponema_F_data_all_ko[,c(1,2)])
Treponema_F_data_all_ko_num <- as.data.frame(table(Treponema_F_data_all_ko$KO))
colnames(Treponema_F_data_all_ko_num) <- c("KO", "SGB_8084_all")

# pN/pS >1
Treponema_F_data_Yak_pN_all <- Treponema_F_data_filter_Yak[which(Treponema_F_data_filter_Yak$pNpS_variants>1),]
Treponema_F_data_filter_pN_ko_all <- merge(Treponema_F_data_Yak_pN_all, all_gene_ko, by.x = "gene", by.y = "qaccver", all.x = T)
Treponema_F_data_filter_pN_function_all <- merge(Treponema_F_data_filter_pN_ko_all, ko_pathway, by.x = "KO", by.y = "ModuleID", all.x = T)
write.csv(Treponema_F_data_filter_pN_function_all, file = "Yak/Treponema_F_data_filter_pN_function_all.csv")
Treponema_F_data_filter_pN_ko <- merge(Treponema_F_data_Yak_pN_all, all_gene_ko, by.x = "gene", by.y = "qaccver")
Treponema_F_data_filter_pN_ko <- unique(Treponema_F_data_filter_pN_ko[,c(1,7)])
Treponema_F_data_filter_pN_ko_num <- as.data.frame(table(Treponema_F_data_filter_pN_ko$KO))
colnames(Treponema_F_data_filter_pN_ko_num) <- c("KO", "SGB_8084_pN")

KEGG_ko_merge_pN <- merge(Treponema_F_data_all_ko_num, Treponema_F_data_filter_pN_ko_num, by = "KO", all = T)
KEGG_ko_merge_pN[is.na(KEGG_ko_merge_pN)] <- 0
KEGG_ko_merge_pN$SGB_8084_all_no <- 2675-KEGG_ko_merge_pN$SGB_8084_all
KEGG_ko_merge_pN$SGB_8084_pN_no <- 179-KEGG_ko_merge_pN$SGB_8084_pN
write.csv(KEGG_ko_merge_pN, file = "Yak/KEGG_ko_merge_pN_SGB_8084.csv")
mydat<-read.xlsx("Yak/KEGG_ko_merge_pN_SGB_8084.xlsx",sheet=1,colNames = T)
result<-c()
for (i in (2:ncol(mydat))){
  newdat<-mydat[,i]
  test<-matrix(newdat,nrow=2,ncol=2,byrow=TRUE)
  qq<-fisher.test(test)   #fisher.test()   chisq.test
  Odds_ratio<-qq[["estimate"]][["odds ratio"]]
  p.result<-qq$p.value
  result.linshi<-cbind(i,Odds_ratio,p.result)
  result<-rbind(result,result.linshi)
}

KO<-names(mydat[,-1])  

result<-cbind(KO,result)
result <- as.data.frame(result)
result$FDR <- p.adjust(result$p.result,"fdr")

result$sig <- "NS"
result$sig[result$FDR<0.05]="*"
result$sig[result$FDR<0.01]="**"
result$sig[result$FDR<0.001]="***"
result$sig[result$FDR<0.0001]="****"
write.table(result, file="fisher_result_pN_SGB_8084.xls",col.names = T,row.names = F,sep="\t",quote=FALSE)

# pN/pS <1
Treponema_F_data_Yak_pS_all <- Treponema_F_data_filter_Yak[which(Treponema_F_data_filter_Yak$pNpS_variants<1),]
Treponema_F_data_filter_pS_ko_all <- merge(Treponema_F_data_Yak_pS_all, all_gene_ko, by.x = "gene", by.y = "qaccver", all.x = T)
Treponema_F_data_filter_pS_function_all <- merge(Treponema_F_data_filter_pS_ko_all, ko_pathway, by.x = "KO", by.y = "ModuleID", all.x = T)
write.csv(Treponema_F_data_filter_pS_function_all, file = "Yak/Treponema_F_data_filter_pS_function_all.csv")
Treponema_F_data_filter_pS_ko <- merge(Treponema_F_data_Yak_pS_all, all_gene_ko, by.x = "gene", by.y = "qaccver")
Treponema_F_data_filter_pS_ko <- unique(Treponema_F_data_filter_pS_ko[,c(1,7)])
Treponema_F_data_filter_pS_ko_num <- as.data.frame(table(Treponema_F_data_filter_pS_ko$KO))
colnames(Treponema_F_data_filter_pS_ko_num) <- c("KO", "SGB_8084_pS")

KEGG_ko_merge_pS <- merge(Treponema_F_data_all_ko_num, Treponema_F_data_filter_pS_ko_num, by = "KO", all = T)
KEGG_ko_merge_pS[is.na(KEGG_ko_merge_pS)] <- 0
KEGG_ko_merge_pS$SGB_8084_all_no <- 2675-KEGG_ko_merge_pS$SGB_8084_all
KEGG_ko_merge_pS$SGB_8084_pS_no <- 1984-KEGG_ko_merge_pS$SGB_8084_pS
write.csv(KEGG_ko_merge_pS, file = "Yak/KEGG_ko_merge_pS_SGB_8084.csv")
mydat<-read.xlsx("Yak/KEGG_ko_merge_pS_SGB_8084.xlsx",sheet=1,colNames = T)
result<-c()
for (i in (2:ncol(mydat))){
  newdat<-mydat[,i]
  test<-matrix(newdat,nrow=2,ncol=2,byrow=TRUE)
  qq<-fisher.test(test)   #fisher.test()   chisq.test
  Odds_ratio<-qq[["estimate"]][["odds ratio"]]
  p.result<-qq$p.value
  result.linshi<-cbind(i,Odds_ratio,p.result)
  result<-rbind(result,result.linshi)
}

KO<-names(mydat[,-1])  

result<-cbind(KO,result)
result <- as.data.frame(result)
result$FDR <- p.adjust(result$p.result,"fdr")

result$sig <- "NS"
result$sig[result$FDR<0.05]="*"
result$sig[result$FDR<0.01]="**"
result$sig[result$FDR<0.001]="***"
result$sig[result$FDR<0.0001]="****"
write.table(result, file="fisher_result_pS_SGB_8084.xls",col.names = T,row.names = F,sep="\t",quote=FALSE)






