library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(devtools)
#library(cowplot)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(forcats)

Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

#### Read in BestConfig files #### (manually cleaned them up in notepad++)
config_2020 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2020/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2021 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2021/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2022 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2022/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

# Make df to spot check bestconfig parentage assignments
temp_wild <- Samples_2111 %>%
  select(SampleID, Cohort_year, Cohort) %>%
  filter(Cohort == "Wild") %>% 
  rename(OffspringID = SampleID,
         Offspring_cohort = Cohort_year) %>% 
  select(-Cohort)

temp_dads <- Samples_2111 %>%
  select(SampleID, Cohort_year, Cohort) %>%
  filter(Cohort != "Wild") %>% 
  rename(FatherID = SampleID,
         Father_cohort = Cohort_year) %>% 
  select(-Cohort)

temp_moms <- Samples_2111 %>%
  select(SampleID, Cohort_year, Cohort) %>%
  filter(Cohort != "Wild") %>% 
  rename(MotherID = SampleID,
         Mother_cohort = Cohort_year) %>% 
  select(-Cohort)

spot_check <- config_2020 %>% 
  bind_rows(config_2021) %>% 
  bind_rows(config_2022) %>% 
  left_join(temp_wild) %>% 
  left_join(temp_dads) %>% 
  left_join(temp_moms) %>% 
  select(Offspring_cohort, Father_cohort, Mother_cohort, OffspringID, FatherID, MotherID) %>% 
  arrange(Offspring_cohort)