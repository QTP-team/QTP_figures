library(readxl)
library(ggplot2)
library(cowplot)
library(ggpubr)

### Figure 6
## pN/pS
# 读取Excel文件
file_path <- "0.data/SGBs_90_pN_pS_values.xlsx"
sheet_name <- "SGBs_90_pN_pS_values"
data <- read_excel(file_path, sheet = sheet_name)

## "Host-swap" 
data_HS <- data[which(data$Type=="Host-swap"),]
data_HS_filter <- data_HS
# 执行线性回归分析
linear_model <- lm(Rate ~ Median, data = data_HS_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_pnps_all <- ggplot(data_HS_filter, aes(x = Rate, y = Median, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Host-swap rate",
    y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

### Group
cal <- function(data_HS_group){linear_model <- lm(Rate ~ Median, data = data_HS_group)
data_HS_group_name <- as.character(substitute(data_HS_group))
group <- strsplit(data_HS_group_name, 'filter_')[[1]][2]
summary(linear_model)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))
title_name <- paste0(group, ":\n", Linear_Equation, ", ", R_squared, ", ", p_value)
return(title_name)}
## phylum
data_HS_filter_Bacteroidota <- data_HS_filter[which(data_HS_filter$Phylum=="p__Bacteroidota"),]
title_Bacteroidota <- cal(data_HS_filter_Bacteroidota)
Phylum.labs <- c("p__Bacteroidota") # 这个是我们希望展示出来的标签名
names(Phylum.labs) <- c("p__Bacteroidota") # 这个是我们希望隐藏的标签名
p_host_swap_pnps_Bacteroidota <- ggplot(data_HS_filter_Bacteroidota, aes(x = Rate, y = Median, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(x = "Host-swap rate",
       y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_grid(. ~ Phylum,labeller = labeller(Phylum = Phylum.labs))+
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) 
## host
data_HS_filter_TA <- data_HS_filter[which(data_HS_filter$Host=="Tibetan_Ass"),]
title_TA <- cal(data_HS_filter_TA)
data_HS_filter_TC <- data_HS_filter[which(data_HS_filter$Host=="Tibetan_Cattle"),]
title_TC <- cal(data_HS_filter_TC)
data_HS_filter_TH <- data_HS_filter[which(data_HS_filter$Host=="Tibetan_Horse"),]
title_TH <- cal(data_HS_filter_TH)
data_HS_filter_TS <- data_HS_filter[which(data_HS_filter$Host=="Tibetan_Sheep"),]
title_TS <- cal(data_HS_filter_TS)
data_HS_filter_Yak <- data_HS_filter[which(data_HS_filter$Host=="Tibetan_Yak"),]
title_Yak <- cal(data_HS_filter_Yak)
Host.labs <- c("TA", "TC", "TH", "TS", "Yak") # 这个是我们希望展示出来的标签名
names(Host.labs) <- c("Tibetan_Ass", "Tibetan_Cattle", "Tibetan_Horse", "Tibetan_Sheep", "Tibetan_Yak") # 这个是我们希望隐藏的标签名
p_host_swap_pnps_host <- ggplot(data_HS_filter, aes(x = Rate, y = Median)) +
  geom_point(aes(color = Host)) +
  geom_smooth(method = "lm", se = TRUE, aes(color = Host)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Host-swap rate",
       y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_grid(. ~ Host,labeller = labeller(Host = Host.labs))+
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8))

## "Co-phylogeny" 
data_CP <- data[which(data$Type=="Co-phylogeny"),]
data_CP_filter <- data_CP
#data_CP_filter <- data_CP_filter[which(data_CP_filter$Phylum=="p__Bacteroidota"),]
#data_CP_filter <- data_CP_filter[which(data_CP_filter$Host=="Tibetan_Sheep"),]
# 执行线性回归分析
linear_model <- lm(Rate ~ Median, data = data_CP_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_pnps_all <- ggplot(data_CP_filter, aes(x = Rate, y = Median, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Co-phylogeny rate",
    y = "pN/pS ratio") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

## ka/ks
# 读取Excel文件
file_path <- "0.data/kaks.xlsx"
sheet_name <- "Sheet1"
data <- read_excel(file_path, sheet = sheet_name)

## "Host-swap" ka.ks>1
data_HS <- data[which(data$Event=="Host-swap"),]
data_HS_filter <- data_HS
# 执行线性回归分析
linear_model <- lm(Rate ~ CountGreaterThan1, data = data_HS_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_kaks_all_CountGreaterThan1 <- ggplot(data_HS_filter, aes(x = Rate, y = CountGreaterThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Host-swap rate",
    y = "Frequency of single copy gene (Ka/Ks>1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

### Group
cal <- function(data_HS_group){linear_model <- lm(Rate ~ CountGreaterThan1, data = data_HS_group)
data_HS_group_name <- as.character(substitute(data_HS_group))
group <- strsplit(data_HS_group_name, 'filter_')[[1]][2]
summary(linear_model)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))
title_name <- paste0(group, ":\n", Linear_Equation, ", ", R_squared, ", ", p_value)
return(title_name)}
## phylum
data_HS_filter_Bacteroidota <- data_HS_filter[which(data_HS_filter$Phyla=="p__Bacteroidota"),]
title_Bacteroidota <- cal(data_HS_filter_Bacteroidota)
Phyla.labs <- c("p__Bacteroidota") # 这个是我们希望展示出来的标签名
names(Phyla.labs) <- c("p__Bacteroidota") # 这个是我们希望隐藏的标签名
p_host_swap_kaks_CountGreaterThan1_Bacteroidota <- ggplot(data_HS_filter_Bacteroidota, aes(x = Rate, y = CountGreaterThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(x = "Host-swap rate",
       y = "Frequency of single copy gene (Ka/Ks>1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_grid(. ~ Phyla,labeller = labeller(Phyla = Phyla.labs))+
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) 

## "Host-swap" ka.ks<1
data_HS <- data[which(data$Event=="Host-swap"),]
data_HS_filter <- data_HS
# 执行线性回归分析
linear_model <- lm(Rate ~ CountLessThan1, data = data_HS_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_host_swap_kaks_all_CountLessThan1 <- ggplot(data_HS_filter, aes(x = Rate, y = CountLessThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Host-swap rate",
    y = "Frequency of single copy gene (Ka/Ks<1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

### Group
cal <- function(data_HS_group){linear_model <- lm(Rate ~ CountLessThan1, data = data_HS_group)
data_HS_group_name <- as.character(substitute(data_HS_group))
group <- strsplit(data_HS_group_name, 'filter_')[[1]][2]
summary(linear_model)
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))
title_name <- paste0(group, ":\n", Linear_Equation, ", ", R_squared, ", ", p_value)
return(title_name)}
## phylum
data_HS_filter_Bacteroidota <- data_HS_filter[which(data_HS_filter$Phyla=="p__Bacteroidota"),]
title_Bacteroidota <- cal(data_HS_filter_Bacteroidota)
Phyla.labs <- c("p__Bacteroidota") # 这个是我们希望展示出来的标签名
names(Phyla.labs) <- c("p__Bacteroidota") # 这个是我们希望隐藏的标签名
p_host_swap_kaks_CountLessThan1_Bacteroidota <- ggplot(data_HS_filter_Bacteroidota, aes(x = Rate, y = CountLessThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(x = "Host-swap rate",
       y = "Frequency of single copy gene (Ka/Ks<1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_grid(. ~ Phyla,labeller = labeller(Phyla = Phyla.labs))+
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8))

## "Co-phylogeny" ka.ks>1
data_CP <- data[which(data$Event=="Co-phylogeny"),]
data_CP_filter <- data_CP
# 执行线性回归分析
linear_model <- lm(Rate ~ CountGreaterThan1, data = data_CP_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_kaks_all_CountGreaterThan1 <- ggplot(data_CP_filter, aes(x = Rate, y = CountGreaterThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Co-phylogeny rate",
    y = "Frequency of single copy gene (Ka/Ks>1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

## "Co-phylogeny" ka.ks<1
data_CP <- data[which(data$Event=="Co-phylogeny"),]
data_CP_filter <- data_CP
# 执行线性回归分析
linear_model <- lm(Rate ~ CountLessThan1, data = data_CP_filter)

# 打印回归分析摘要
summary(linear_model)

# 提取回归方程的系数和R平方
coefficients <- coef(linear_model)
intercept <- coefficients[1]
slope <- coefficients[2]
r_squared <- summary(linear_model)$r.squared

# 打印回归方程和R平方
Linear_Equation <- paste0("y=", round(intercept, 3), "+", round(slope, 3), "*x")
R_squared <- paste0("R-squared=", round(r_squared, 3))

# 进行F检验
anova_result <- anova(linear_model)
p_value <- paste0("p-value=", round(anova_result$`Pr(>F)`[1], 5))

# 绘制散点图和回归线
title_name <- paste0(Linear_Equation, ", ", R_squared, ", ", p_value)
p_co_phylogeny_kaks_all_CountLessThan1 <- ggplot(data_CP_filter, aes(x = Rate, y = CountLessThan1, color = "#FF8080")) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#FF8080") +
  labs(#title = title_name,
    x = "Co-phylogeny rate",
    y = "Frequency of single copy gene (Ka/Ks<1)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_text(size = 10)) + theme(axis.title.y=element_text(size = 10)) +
  theme(axis.text.x=element_text(size = 8)) + theme(axis.text.y=element_text(size = 8)) +
  theme(legend.position = 'none')

p_merge1 <- plot_grid(p_host_swap_kaks_CountGreaterThan1_Bacteroidota, p_host_swap_kaks_CountLessThan1_Bacteroidota, 
                      p_host_swap_pnps_Bacteroidota, 
                      align = "h", labels = c("A", "B", "C"), nrow = 1)
p_merge2 <- plot_grid(p_merge1, p_host_swap_pnps_host, 
                      labels = c("", "D"), nrow = 2)
ggsave(p_merge2, filename = "FigureS.pdf", height = 8, width = 10)