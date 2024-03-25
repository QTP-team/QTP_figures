library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(ggpubr)

sample_data <- read.csv("0.data/TableS0_samples_info.csv")
SGB_data <- read.csv("0.data/SGBs_data.csv")
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

