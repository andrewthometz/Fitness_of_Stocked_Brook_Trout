library(readxl)
library(tidyverse)

# My samples
Samples_2205 <- read_delim("C:/Users/27tho/UWSP Backup/2205/Thometz_scripts/Samples_2205.csv") %>% 
  select(WaterbodyName, WBIC, Latitude, Longitude, HUC_4, HUC_8) %>% 
  distinct(WaterbodyName, .keep_all = TRUE)

# Read in Erdman pops
erdmanPops <- read_excel("C:/Users/27tho/UWSP Backup/2205/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx") %>% 
  select(WaterbodyName, WBIC, Latitude, Longitude, HUC_4, HUC_8) %>% 
  distinct(WaterbodyName, .keep_all = TRUE) #%>% 
  #arrange(WaterbodyName)

pops_in_both <- Samples_2205 %>% 
  filter(WaterbodyName %in% erdmanPops$WaterbodyName)

write.csv(pops_in_both, 
          "C:/Users/27tho/UWSP Backup/2205/Streams_moving_forward/Streams_reanalyzed.csv",
          row.names = FALSE)

erdman_exclusives <- erdmanPops %>% 
  filter(!(WaterbodyName %in% Samples_2205$WaterbodyName))

write.csv(erdman_exclusives, 
          "C:/Users/27tho/UWSP Backup/2205/Streams_moving_forward/Erdman_exclusives.csv",
          row.names = FALSE)

# Manually combined the two csv file into one excel spreadsheet






