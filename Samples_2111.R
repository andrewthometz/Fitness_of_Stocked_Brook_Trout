library(tidyverse)
library(readxl)

### Reading in MCGL Data ###
MCGL1 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Archive (04-00001 to 10-18506).xlsx", 
                    guess_max = 1000000) %>%
  select(SampleID, TotalLength, ProjectID, PlateID, CommonName, YearCollected, MonthCollected, DayCollected, WaterbodyName, State, Comments) %>% 
  mutate(DayCollected = as.numeric(DayCollected))
         
MCGL2 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (11-00001 to 17-15647).xlsx",
                    guess_max = 1000000) %>%
  select(SampleID, TotalLength, ProjectID, PlateID, CommonName, YearCollected, MonthCollected, DayCollected, WaterbodyName, State, Comments) %>% 
  mutate(DayCollected = as.numeric(DayCollected),
         TotalLength = as.numeric(TotalLength))

MCGL3 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (18-00001 to 19-28271).xlsx",
                    guess_max = 1000000) %>%
  select(SampleID, TotalLength, ProjectID, PlateID, CommonName, YearCollected, MonthCollected, DayCollected, WaterbodyName, State, Comments) %>% 
  mutate(DayCollected = as.numeric(DayCollected),
         TotalLength = as.numeric(TotalLength))

MCGL4 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (21-00001 to 22-19350).xlsx",
                    guess_max = 1000000) %>% 
  select(SampleID, TotalLength, ProjectID, PlateID, CommonName, YearCollected, MonthCollected, DayCollected, WaterbodyName, State, Comments) %>% 
  mutate(YearCollected = as.numeric(YearCollected),
         DayCollected = as.numeric(DayCollected),
         TotalLength = as.numeric(TotalLength))

MCGL_all <- bind_rows(MCGL1, MCGL2, MCGL3, MCGL4) %>% 
  mutate(TotalLength = round(TotalLength, digits = 0))

### Filtering to all 2111 samples ### (1608 total fish)
Samples_2111 <- MCGL_all %>% 
  filter(str_detect(ProjectID, "2111"),
         SampleID != "18-10400",
         SampleID != "18-10398",
         SampleID != "18-10399",
         SampleID != "21-07852",
         SampleID != "21-07851",
         SampleID != "21-07921",
         WaterbodyName != "Cutler Creek") %>% 
  mutate(Cohort_year = case_when(YearCollected == 2018 & WaterbodyName == "Hay River" ~ "F1_2018",
                                 YearCollected == 2018 & WaterbodyName == "Melancthon Creek" ~ "F2_2018",
                                 YearCollected == 2018 & str_detect(WaterbodyName, "Croix") ~ "D_2018",
                                 YearCollected == 2019 & WaterbodyName == "Cady Creek" ~ "F1_2019",
                                 YearCollected == 2019 & WaterbodyName == "West Fork Kickapoo" ~ "F2_2019",
                                 YearCollected == 2019 & str_detect(WaterbodyName, "Croix") ~ "D_2019",
                                 YearCollected == 2020 & WaterbodyName == "Lowery Creek" ~ "F1_2020",
                                 YearCollected == 2020 & WaterbodyName == "SFK Hay River x (Lowery or Melancthon)" ~ "F2_2020",
                                 YearCollected == 2020 & str_detect(WaterbodyName, "Croix") ~ "D_2020",
                                 YearCollected == 2020 & MonthCollected == "Sep" & DayCollected == 8 ~ "W_2019",
                                 YearCollected == 2021 & MonthCollected == "Apr" & DayCollected == 5 ~ "W_2019",
                                 YearCollected == 2021 & MonthCollected == "Sep" & DayCollected == 8 ~ "W_2020",
                                 YearCollected == 2022 & MonthCollected == "Apr" & DayCollected == 4 ~ "W_2020",
                                 YearCollected == 2022 & MonthCollected == "Oct" & DayCollected == 6 ~ "W_2021"),
         Cohort = case_when(str_detect(Cohort_year, "F1") ~ "F1",
                            str_detect(Cohort_year, "F2") ~ "F2",
                            str_detect(Cohort_year, "D") ~ "Domestic",
                            str_detect(Cohort_year, "W") ~ "Wild"),
         Year = case_when(str_detect(Cohort_year, "2018") ~ "2018",
                          str_detect(Cohort_year, "2019") ~ "2019",
                          str_detect(Cohort_year, "2020") ~ "2020",
                          str_detect(Cohort_year, "2021") ~ "2021")) %>%
                          #str_detect(Cohort_year, "2022") ~ "2022")
  write_csv("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv")
         
Samples_2111 %>% 
  group_by(Cohort) %>% 
  count(Year)
