library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(devtools)
#library(cowplot)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(forcats)

#### Read in BestConfig file #### (manually cleaned it up in notepad++)
Best_config <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_pedigree_simulation/PA_postSim_2111/PA_postSim_2111.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

#### Identify incorrect assignments ####
Best_config %>% 
  mutate(OffspringID = str_remove(OffspringID, "C1"),
         OffspringID = str_remove(OffspringID, "C2"),
         OffspringID = str_remove(OffspringID, "C3"),
         OffspringID = str_remove(OffspringID, "C4"),
         OffspringID = str_remove(OffspringID, "C5"),
         OffspringID = str_remove(OffspringID, "C6")) %>% 
  unite(col = "Assignment", 2:3, sep = "") %>% 
  filter(OffspringID != Assignment) # No incorrect assignments found...
         
