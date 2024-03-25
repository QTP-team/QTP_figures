setwd("~/Desktop/BGI/高原项目/QTP_redo/Figure 6.pN_pS/Genus_38/")
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(ggpubr)

sample_data <- read.csv("0.data/TableS0_samples_info.csv")
SGB_data <- read.csv("0.data/SGBs_data.csv")

process <- function(SGBID) {
  pwd <- "1.pN_pS_result/"
  pN_pS_data <- read.csv(paste0(pwd,SGBID,"_pNpS_value.csv"), header = F)
  colnames(pN_pS_data) <- c("Sample", "Contig", "Gene", "Coverage", "Breadth", "pNpS_variant")
  pN_pS_data$SGB <- SGBID
  pN_pS_data_merge <- merge(pN_pS_data, sample_data, by.x = "Sample", by.y = "Sample_ID")
  pN_pS_data_merge1 <- merge(pN_pS_data_merge, SGB_data, by.x = "SGB", by.y = "SGB_ID")
  pN_pS_result <- pN_pS_data_merge1[,c(1,2,3,4,5,6,7,9,17)]
  return(pN_pS_result)
}

### data

## method2
setwd("~/Desktop/BGI/高原项目/QTP_redo/Figure 6.pN_pS/Genus_38/1.pN_pS_result/")
all_file <- list.files()
n <- length(all_file)
method2.empty <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(method2.empty) <- c("SGB_ID", "Tibetan_Ass", "Tibetan_Horse", "Tibetan_Antelope", "Tibetan_Sheep", "Tibetan_Yak", "Tibetan_Cattle")
host.process <- function(host) {
  pN_pS_data <- pN_pS_result[which(pN_pS_result$Host==host),]
  pN_pS_median <- aggregate(pN_pS_data$pNpS_variant, by=list(pN_pS_data$Gene), median)
  colnames(pN_pS_median) <- c("Gene", "pNpS_variants_median")
  return(pN_pS_median)
}

for (i in 1:376){
  SGBID <- strsplit(all_file[i],split = "_pNpS")[[1]][1]
  setwd("~/Desktop/BGI/高原项目/QTP_redo/Figure 6.pN_pS/Genus_38")
  pN_pS_result <- process(SGBID)
  sample_result <- as.data.frame(table(unique(pN_pS_result[,c(2,8)])$Host))
  colnames(sample_result) <- c("Host", "Freq")
  host <- as.character(sample_result[,1])
  method2.empty[i,1] <- SGBID
  for(j in 1:length(host)){
    assign(paste0(host[j],"_data"), host.process(host[j]))
    method2.empty[i,which(names(method2.empty)%in%host[j])] <- median(get(paste0(host[j],"_data"))$pNpS_variants_median)
  }
}
colnames(method2.empty) <- c("SGB_ID", "TA", "TH", "TAN", "TS", "Yak", "TC")

write.csv(method2.empty, file = "method2_pNpS.csv")


### r2
data <- read.csv("0.data/SGBs_data_final.csv")

# Co.phylogeny all
linear_model <- lm(pNpS_ratio ~ Co.phylogeny_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_pnps_all <- ggplot(data, aes(x = Co.phylogeny_rate, y = pNpS_ratio)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Co-phylogeny rate",
       y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_co_phylogeny_pnps_all, filename = "2.plot_result/co_phylogeny_pnps_r2.pdf", width = 5, height = 4)

linear_model <- lm(More ~ Co.phylogeny_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_pnps_all <- ggplot(data, aes(x = Co.phylogeny_rate, y = More)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Co-phylogeny rate",
       y = "Ka/Ks > 1") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_co_phylogeny_pnps_all, filename = "2.plot_result/co_phylogeny_kaks_more_r2.pdf", width = 5, height = 4)

linear_model <- lm(Less ~ Co.phylogeny_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_pnps_all <- ggplot(data, aes(x = Co.phylogeny_rate, y = Less)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Co-phylogeny rate",
       y = "Ka/Ks < 1") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_co_phylogeny_pnps_all, filename = "2.plot_result/co_phylogeny_kaks_less_r2.pdf", width = 5, height = 4)
# Host.swap all
linear_model <- lm(pNpS_ratio ~ Host.swap_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_pnps_all <- ggplot(data, aes(x = Host.swap_rate, y = pNpS_ratio)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Host-swap rate",
       y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_host_swap_pnps_all, filename = "2.plot_result/host_swap_pnps_r2.pdf", width = 5, height = 4)

linear_model <- lm(More ~ Host.swap_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_pnps_all <- ggplot(data, aes(x = Host.swap_rate, y = More)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Host-swap rate",
       y = "Ka/Ks > 1") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_host_swap_pnps_all, filename = "2.plot_result/host_swap_kaks_more_r2.pdf", width = 5, height = 4)

linear_model <- lm(Less ~ Host.swap_rate, data = data)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 10))
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_pnps_all <- ggplot(data, aes(x = Host.swap_rate, y = Less)) +
  geom_point(shape = 1, color = "#000000") +
  geom_smooth(method = "lm", se = TRUE, color = "#808080") +
  labs(title = title_name,
       x = "Host-swap rate",
       y = "Ka/Ks < 1") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')
ggsave(p_host_swap_pnps_all, filename = "2.plot_result/host_swap_kaks_less_r2.pdf", width = 5, height = 4)

