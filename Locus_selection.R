
################ This has now been implemented into the 2111 and 2205 megasat_to_genepop scripts ###############

library(tidyverse)
library(readxl)
library(adegenet)

Locus_data <- read_excel("X:/2111_F1F2D_BKT/BKT_Locus_Evaluation.xlsx") %>% 
  select(1:9)

Locus_data %>% count(Legacy_or_New)

hwe_cutoff <- 0.15
naf_cutoff <- 0.2

Locus_data %>% 
  filter(naf_2205 < naf_cutoff &
         HWE_prop_2205 < hwe_cutoff |
         Locus == "SfoC38" | 
         Locus == "SFOC86",
         Locus != "Salv_1_Di_30_1186") %>% 
  count(Legacy_or_New)

Locus_data %>% 
  ggplot(aes(x = HWE_prop_2205)) +
  geom_vline(xintercept = hwe_cutoff) +
  geom_histogram()

Final_loci <- Locus_data %>% 
  filter(naf_2205 < naf_cutoff &
         HWE_prop_2205 < hwe_cutoff |
         Locus == "SfoC38" | 
         Locus == "SFOC86",
         Locus != "Salv_1_Di_30_1186")
  mutate(Locus = str_replace_all(Locus, c("\\." = "_", "-" = "_")), .keep = "unused") %>% 
  select(Locus)

write.csv(Final_loci, file = "X:/2111_F1F2D_BKT/Final_loci_2111_2205.csv")

dat.genind[loc = Locus_vector]
