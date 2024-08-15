library(tidyverse)
library(ggforce)

### Function to make length_dist files tidy ###

length_dist_tidy <- function(x){
  
  BKT_ID <- deparse(substitute(x))
  
  output <- x %>% 
    mutate(uSat_locus = str_replace_all(Microsatellite, "-", "_"), .keep = "unused") %>% 
    select(-sum) %>%
    pivot_longer(-c(uSat_locus, scores), names_to = "Length", values_to = "Read_count") %>% 
    add_column(SampleID = BKT_ID) %>% 
    mutate(SampleID = str_remove_all(SampleID, "BKT_")) %>% 
    mutate(Length = as.numeric(Length))
  
  print(output)
}

### Read in 12 BKT ### Grabbed largest few files from each cat folder. Used some from each folder to get "even" representation

BKT_21_08632 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_1of5_output/length_distribution/Genotype_21-08632.txt")
BKT_21_08633 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_1of5_output/length_distribution/Genotype_21-08633.txt")
BKT_21_08624 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_1of5_output/length_distribution/Genotype_21-08624.txt")

BKT_21_09151 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_2of5_output/length_distribution/Genotype_21-09151.txt")
BKT_21_09203 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_2of5_output/length_distribution/Genotype_21-09203.txt")
BKT_21_09148 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_2of5_output/length_distribution/Genotype_21-09148.txt")

BKT_21_09303 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_3of5_output/length_distribution/Genotype_21-09303.txt")
BKT_21_09425 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_3of5_output/length_distribution/Genotype_21-09425.txt")

BKT_22_11471 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_4of5_output/length_distribution/Genotype_22-11471.txt")
BKT_22_11380 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_4of5_output/length_distribution/Genotype_22-11380.txt")

BKT_22_11777 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_5of5_output/length_distribution/Genotype_22-11777.txt")
BKT_22_11658 <- read_delim("X:/2111_F1F2D_BKT/MEGAsat_outputs/MEGAsat_2111_cat_5of5_output/length_distribution/Genotype_22-11658.txt")

### Run each fish's length_dist file through my tidying function while joining them ###

All_LD_data <- bind_rows(length_dist_tidy(BKT_21_08632),
                         length_dist_tidy(BKT_21_08633),
                         length_dist_tidy(BKT_21_08624),
                         length_dist_tidy(BKT_21_09151),
                         length_dist_tidy(BKT_21_09203),
                         length_dist_tidy(BKT_21_09148),
                         length_dist_tidy(BKT_21_09303),
                         length_dist_tidy(BKT_21_09425),
                         length_dist_tidy(BKT_22_11471),
                         length_dist_tidy(BKT_22_11380),
                         length_dist_tidy(BKT_22_11777),
                         length_dist_tidy(BKT_22_11658))

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.25) %>% 
  ungroup() %>% 
  summarize(x = n_distinct(uSat_locus))

### Make peak morphology plots, one locus and 12 BKT per page ### Saves directly as pdf

# Beware, this takes about 40min to run #

pdf("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_PeakMorph_Thometz.pdf", paper = "a4r", width = 11, height = 9)

ProgressBar <- txtProgressBar(min = 0, max = 91, style = 3)

for(i in 1:91){
  
print(All_LD_data %>% 
        group_by(uSat_locus, SampleID) %>% 
        #slice_max(Read_count, n = 8) %>% ### This seems to determine run time (4 ~ 40min, 10 ~ 1.5hrs)
        filter(Read_count >= max(Read_count)*0.05) %>% 
        ggplot(aes(x = Length, y = Read_count)) +
        geom_col(fill = "seagreen") +
        geom_text(aes(label = Length), 
                  hjust = 0.5, 
                  vjust = 1.5, 
                  colour = "black") +
        ylab("Read count") +
        xlab("Allele length") +
        scale_y_continuous(expand = c(0,0)) +
        #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
        theme_classic() +
        facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                            nrow = 3,
                            ncol = 4,
                            page = i, 
                            scales = "free"))
  Sys.sleep(0.1)
  setTxtProgressBar(ProgressBar, i)
}
close(ProgressBar)
dev.off()

### Single page of peak morphs for testing plot changes ###

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.05) %>% 
  ggplot(aes(x = Length, y = Read_count)) +
  geom_col(fill = "seagreen") +
  geom_text(aes(label = Length), hjust = 0.5, vjust = 1.5, colour = "black") +
  ylab("Read count") +
  xlab("Allele length") +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
  theme_classic() +
  facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                      nrow = 3,
                      ncol = 4,
                      page = 49, 
                      scales = "free")
