library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(devtools)
library(vegan)
library(patchwork)

#citation("vegan")

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

#### Read in BestConfig files #### (manually cleaned them up in notepad++)
config_2019 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2020/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2020 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2021/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2021 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2022/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

#### Function from Bobby ####
Ns_calc <- function(family){
  require(vegan)
  family$FatherID <- paste0("Dad", family$FatherID)
  family$MotherID <- paste0("Mom", family$MotherID)
  #making matrix for dads
  dads <- data.frame(matrix(0,nrow = length(family$OffspringID), ncol = length(unique(family$FatherID))))
  colnames(dads) <- unique(family$FatherID)
  rownames(dads) <- family$OffspringID
  #making matrix for moms
  moms <- data.frame(matrix(0, nrow = length(family$OffspringID), ncol = length(unique(family$MotherID))))
  colnames(moms) <- unique(family$MotherID)
  rownames(moms) <- family$OffspringID
  
  #loop to fill in matrix with parents
  for (i in 1:length(family$OffspringID)) {
    off <- family[i,]
    dadn <- which(off$FatherID == colnames(dads))
    dads[i, dadn] <- 1
  }
  for (i in 1:length(family$OffspringID)) {
    off <- family[i,]
    momn <- which(off$MotherID == colnames(moms))
    moms[i, momn] <- 1
  }
  parents <- cbind(moms, dads)
  Ns <- as.matrix(ncol(parents))
  colnames(Ns) <- "Ns"
  ns_points <- specaccum(parents, method = "random")
  asymp <- specpool(parents)
  output <- list(ns_points, asymp, Ns)
  output
}

#### Run config files through the function ####
# 2020 Offspring
Ns_2019 <- Ns_calc(family = config_2019)

Ns_2019_data <- tibble(y = Ns_2019[[1]]$richness,
                       y_sd = Ns_2019[[1]]$sd,
                       x = Ns_2019[[1]]$sites)

# 2021 Offspring
Ns_2020 <- Ns_calc(family = config_2020)

Ns_2020_data <- tibble(y = Ns_2020[[1]]$richness,
                       y_sd = Ns_2020[[1]]$sd,
                       x = Ns_2020[[1]]$sites)
# 2022 Offspring
Ns_2021 <- Ns_calc(family = config_2021)

Ns_2021_data <- tibble(y = Ns_2021[[1]]$richness,
                       y_sd = Ns_2021[[1]]$sd,
                       x = Ns_2021[[1]]$sites)

#### Plot #### (jackknife commented out because untrustworthy)
# 2020 Offspring
Plot_2019 <- Ns_2019_data %>% 
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y-y_sd, ymax = y + y_sd), 
              alpha = 0.95, 
              fill = "lightgrey") +
  geom_point(size = 0.5) +
  geom_hline(yintercept = Ns_2019[[2]]$chao, 
             linetype = "dashed") +
  geom_text(data = data.frame(x = 0, y = Ns_2019[[2]]$chao), 
            aes(x, y), 
            label = paste("Chao estimate =", round(Ns_2019[[2]]$chao, digits = 2)),
            hjust = -0.25, 
            vjust = 1.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = NULL,
       y = bquote(N[S]),
       title = bquote("(A)  " ~ 2019 ~ (N[offspring] ~ "= 145"))) +
  theme_classic()

# 2021 Offspring
Plot_2020 <- Ns_2020_data %>% 
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y-y_sd, ymax = y + y_sd), 
              alpha = 0.95, 
              fill = "lightgrey") +
  geom_point(size = 0.5) +
  geom_hline(yintercept = Ns_2020[[2]]$chao, 
             linetype = "dashed") +
  geom_text(data = data.frame(x = 0, y = Ns_2020[[2]]$chao), 
            aes(x, y), 
            label = paste("Chao estimate =", round(Ns_2020[[2]]$chao, digits = 2)),
            hjust = -0.25, 
            vjust = 1.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = bquote(N["offspring sampled"]),
       y = NULL,
       title = bquote("(B)  " ~ 2020 ~ (N[offspring] ~ "= 72"))) +
  theme_classic() +
  theme(axis.text.y = element_blank())

# 2022 Offspring
Plot_2021 <- Ns_2021_data %>% 
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y-y_sd, ymax = y + y_sd), 
              alpha = 0.95, 
              fill = "lightgrey") +
  geom_point(size = 0.5) +
  geom_hline(yintercept = Ns_2021[[2]]$chao, 
             linetype = "dashed") +
  geom_text(data = data.frame(x = 0, y = Ns_2021[[2]]$chao), 
            aes(x, y), 
            label = paste("Chao estimate =", round(Ns_2021[[2]]$chao, digits = 2)),
            hjust = -0.25, 
            vjust = 1.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  labs(x = NULL,
       y = NULL,
       title = bquote("(C)  " ~ 2021 ~ (N[offspring] ~ "= 41"))) +
  theme_classic() +
  theme(axis.text.y = element_blank())
  

# Combine the plots
all_plots <- Plot_2019 + Plot_2020 + Plot_2021

ggsave(filename = "PAA_plots_all.pdf",
       plot = all_plots,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Pedigree_accum_3runs",
       height = 4,
       width = 9,
       units = "in")

ggsave(filename = "PAA_plots_all.png",
       plot = all_plots,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Pedigree_accum_3runs",
       height = 4,
       width = 9,
       units = "in")
