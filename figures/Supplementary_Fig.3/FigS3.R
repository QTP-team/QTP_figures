library(ggplot2)
library(cowplot)

#Horse
Horse_data <- read.table("0.data/Horse_result.txt")
Horse_data[,3] <- rep(c(seq(from = 0, to = 75, by =5), 79), each = 10)
colnames(Horse_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Horse_p <- ggplot(Horse_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(300, 1500)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 80)) +
  labs(x = 'Numer of Horse samples', y = 'Number of SGBs')

#Ass
Ass_data <- read.table("0.data/Ass_result.txt")
Ass_data[,3] <- rep(c(seq(from = 0, to = 45, by =5), 48), each = 10)
colnames(Ass_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Ass_p <- ggplot(Ass_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(400, 1500)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 50)) +
  labs(x = 'Numer of Ass samples', y = 'Number of SGBs')

#Cattle
Cattle_data <- read.table("0.data/Cattle_result.txt")
Cattle_data[,3] <- rep(c(seq(from = 0, to = 180, by =20), 196), each = 10)
colnames(Cattle_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Cattle_p <- ggplot(Cattle_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(400, 2500)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 200)) +
  labs(x = 'Numer of Cattle samples', y = 'Number of SGBs')

#Yak
Yak_data <- read.table("0.data/Yak_result.txt")
Yak_data[,3] <- rep(c(seq(from = 0, to = 350, by =50), 388), each = 10)
colnames(Yak_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Yak_p <- ggplot(Yak_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(1000, 5000)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 400)) +
  labs(x = 'Numer of Yak samples', y = 'Number of SGBs')

#Sheep
Sheep_data <- read.table("0.data/Sheep_result.txt")
Sheep_data[,3] <- rep(c(seq(from = 0, to = 420, by =30), 446), each = 10)
colnames(Sheep_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Sheep_p <- ggplot(Sheep_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(1000, 7000)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 450)) +
  labs(x = 'Numer of Sheep samples', y = 'Number of SGBs')

#Antelope
Antelope_data <- read.table("0.data/Antelope_result.txt")
Antelope_data[,3] <- rep(c(seq(from = 0, to = 240, by =20), 255), each = 10)
colnames(Antelope_data) <- c("taxo", "Number_of_SGBs", "Number_of_samples")
Antelope_p <- ggplot(Antelope_data, aes(x = Number_of_samples, y = Number_of_SGBs, colour = taxo)) +
  geom_point(size = 1) +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth()+ 
  theme_classic() +
  theme(legend.position = "None") +
  #scale_y_continuous(expand = c(0, 0), limits = c(1000, 6000)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 260)) +
  labs(x = 'Numer of Antelope samples', y = 'Number of SGBs')


merge_plot <- plot_grid(Ass_p, Horse_p, Antelope_p, Sheep_p, Yak_p, Cattle_p, 
                        nrow = 3, ncol = 2, 
                        label_x = 0,
                        label_y = 1,
                        hjust = -0.5,
                        vjust = 1.5,
                        labels = c("A", "B", "C", "D", "E", "F"))
ggsave(merge_plot, filename = "Merge.pdf", width = 10, height = 10)
