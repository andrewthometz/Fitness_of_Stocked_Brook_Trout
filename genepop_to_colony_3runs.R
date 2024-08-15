library(tidyverse)
library(readxl)
library(ggrepel)
library(poppr)
library(adegenet)

#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("thierrygosselin/radiator")
library(radiator)

#install.packages("remotes")
#remotes::install_github("jaredhomola/MCGLfxns")
library(MCGLfxns)

#### Read in genetic data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

#### Read in 2111 metadata ####
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  arrange(Cohort) %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

cohort_view <- Samples_2111 %>% 
  select(SampleID, Cohort)

#### subset to genetic data by Cohort ####
Data_2111@pop <- as_factor(Samples_2111$Cohort)

# Write COLONY file for 2020 wild offspring, only including possible parents groups
W_2020 <- popsub(Data_2111, sublist = c("W_2020", 
                                        "D_2018",
                                        "F1_2018",
                                        "F2_2018")) %>% 
  tidy_genind()

#write_colony(W_2020, filename = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Input_files/Wild_2020_input.colony")

# Write COLONY file for 2021 wild offspring, only including possible parents groups
W_2021 <- popsub(Data_2111, sublist = c("W_2021", 
                                        "D_2018",
                                        "D_2019",
                                        "F1_2018",
                                        "F1_2019",
                                        "F2_2018"
                                        #"F2_2019"
                                        )) %>% 
  tidy_genind()

#write_colony(W_2021, filename = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Input_files/Wild_2021_input.colony")

# Write COLONY file for 2022 wild offspring, only including possible parents groups
W_2022 <- popsub(Data_2111, sublist = c("W_2022",
                                        "D_2018",
                                        "D_2019",
                                        "D_2020",
                                        "F1_2018",
                                        "F1_2019",
                                        "F1_2020",
                                        "F2_2018",
                                        #"F2_2019",
                                        "F2_2020")) %>% 
  tidy_genind()

#write_colony(W_2022, filename = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Input_files/Wild_2022_input.colony")

#### Likelihood of parentage for W_2020 cohort ####
Total_n_stocked <- 3 * 500 * 1  # 3 groups, 500 fish each, 1 year

potential_parents <- 3 * 200

potential_parents / Total_n_stocked # 0.40

#### Likelihood of parentage for W_2021 cohort ####
Total_n_stocked <- 3 * 500 * 2  # 3 groups, 500 fish each, 2 years

potential_parents <- 3 * 200 + 3 * 50 - 50 # minus 50 for 2019 F2's that are aunts/uncles

potential_parents / Total_n_stocked # 0.2333

#### Likelihood of parentage for W_2022 cohort ####
Total_n_stocked <- 3 * 500 * 3  # 3 groups, 500 fish each, 3 years

potential_parents <- 3 * 2 * 200 + 3 * 50 - 50 # minus 50 for 2019 F2's that are aunts/uncles

potential_parents / Total_n_stocked # 0.2888
