library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(remotes)
#install_github("thierrygosselin/radiator")
#library(radiator)

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) #%>% 
  #mutate(GD_groups = case_when(Cohort == "Wild" ~ "Wild",
  #                             .default = Cohort_year))

Samples_2111 %>% 
  count(Cohort_year)

# Fill the pop slots
Data_2111@pop <- as_factor(Samples_2111$Cohort_year)

Data_2111 <- popsub(Data_2111,
                    exclude = c("F2_2019",
                                "W_2019",
                                "W_2020",
                                "W_2021"))

#### Calculate genetic diversity measures ####
library(hierfstat)

gd <- basic.stats(Data_2111)

H_expected <- as.data.frame.matrix(gd$Hs) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:9,
               names_to = "pop",
               values_to = "He") %>% 
  drop_na(He) %>% 
  group_by(pop) %>% 
  summarize("He" = round(mean(He), 3))

H_observed <- as.data.frame.matrix(gd$Ho) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:9,
               names_to = "group",
               values_to = "Ho") %>% 
  drop_na(Ho) %>% 
  group_by(group) %>% 
  summarize("Ho" = round(mean(Ho), 3))

F_is <- as.data.frame.matrix(gd$Fis) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:9,
               names_to = "group",
               values_to = "Fis") %>% 
  drop_na(Fis) %>% 
  group_by(group) %>% 
  summarize("Fis" = round(mean(Fis), 3))

Allelic_richness <- as.data.frame(allelic.richness(Data_2111)) %>% 
  as_tibble(rownames = "locus") %>% 
  select(-min.all) %>% 
  pivot_longer(2:9,
               names_to = "group",
               values_to = "ar") %>% 
  drop_na(ar) %>% 
  group_by(group) %>% 
  summarize("ar" = round(mean(ar), 3))


cohorts_years <- Samples_2111 %>% 
  select(Cohort, Year, Cohort_year) %>% 
  distinct()

GD_tibble <- bind_cols(H_expected$pop,
                       Allelic_richness$ar, 
                       H_observed$Ho, 
                       H_expected$He,
                       F_is$Fis) %>% 
  rename("Cohort_year" = 1,
         "Ar" = 2,
         "Ho" = 3, 
         "He" = 4,
         "Fis" = 5) %>% 
  left_join(cohorts_years)

write_csv(GD_tibble, "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2111.csv")

#############################################
#### Table of offspring per group and Ho ####

#### Read in BestConfig files #### (manually cleaned them up in notepad++)
config_2019 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2020/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2020 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2021/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2021 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2022/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

# Bind them together to make tidy config df
temp_wild <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort == "Wild") %>% 
  rename(OffspringID = SampleID,
         Offspring_cohort = Cohort,
         Offspring_year = Year)

temp_parents <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort != "Wild") %>% 
  rename(ParentID = SampleID,
         Parent_cohort = Cohort,
         Parent_year = Year)

config_all <- config_2019 %>% 
  bind_rows(config_2020) %>% 
  bind_rows(config_2021) %>% 
  pivot_longer(cols = c(FatherID, MotherID), 
               names_to = "Parent_type", 
               values_to = "ParentID") %>% 
  left_join(temp_parents) %>% 
  left_join(temp_wild) %>% 
  select(-Parent_type) %>% 
  drop_na()
#mutate(Parent_cohort = case_when(is.na(Parent_cohort) ~ "Unknown", .default = Parent_cohort))

N_offspring <- config_all %>% 
  group_by(Parent_cohort, Parent_year) %>% 
  summarize(n_offspring = n_distinct(OffspringID)) %>% 
  rename(Cohort = Parent_cohort,
         Year = Parent_year)

GD_offspring <- GD_tibble %>% 
  select(-Cohort_year) %>% 
  left_join(N_offspring) %>% 
  select(Year, Cohort, n_offspring, Ar, Ho, He, Fis)

write_csv(GD_offspring, "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Genetic_diversity/GD_n_offspring.csv")

