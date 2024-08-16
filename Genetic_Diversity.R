# Load packages
library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(remotes)
library(hierfstat)

#install_github("thierrygosselin/radiator")
#library(radiator)

######################################################################
#### Estimate measures of genetic diversity and write as csv file ####
######################################################################

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

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
gd <- basic.stats(Data_2111)

# Extract each measure of genetic diversity, clean it up, and bind them back together at the end
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

# Create tibble to bring metadata into genetic diversity tibble
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

# Write as csv
write_csv(GD_tibble, "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2111.csv")
