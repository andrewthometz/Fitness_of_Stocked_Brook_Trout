library(tidyverse)
library(readxl)
library(miscTools)

#### Reading in 2111 metadata ####
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  arrange(Cohort)

Samples_2111 %>% count(Cohort)

#### Final locus selections ####
Locus_data <- read_excel("X:/2111_F1F2D_BKT/BKT_Locus_Evaluation.xlsx") %>% 
  select(1:9)

hwe_cutoff <- 0.15
naf_cutoff <- 0.2

Final_loci <- Locus_data %>% 
  filter(naf_2205 < naf_cutoff &
         HWE_prop_2205 < hwe_cutoff |
         Locus == "SfoC38" | 
         Locus == "SFOC86",
         Locus != "Salv_1_Di_30_1186") %>% 
  select(Locus)

#### Read in genotypes, merge with sample data, organize by population, remove negative controls ####
genotype_2111 <- read_excel("X:/2111_F1F2D_BKT/MEGAsat_outputs/Genotype_2111_SP.xlsx", 
                           col_types = "text") %>% 
  rename(SampleID = Sample_idx1_idx2)

# Need genotype data from the 2111 fish in my 2205 sequence run
genotype_2205 <- read_excel("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/Genotype_2205.xlsx",
                             col_types = "text") %>% 
  rename(SampleID = Sample_idx1_idx2)

genotypes_all <- genotype_2111 %>% 
  bind_rows(genotype_2205)

genotype_data_1 <- Samples_2111 %>% 
  select(SampleID, Cohort) %>% 
  left_join(genotypes_all, by = "SampleID") #%>% 

colnames(genotype_data_1) <- genotype_data_1 %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

genotype_data_1 <- genotype_data_1 %>% 
  select(SampleID, Cohort, matches(Final_loci$Locus))

#### Correct names for each locus #### and filter by final loci
locus_names <- colnames(genotype_data_1) %>%
  as_tibble() %>%
  #mutate(Locus = str_replace_all(value, c("\\." = "_", "-" = "_")), .keep = "unused") %>% 
  slice(-c(1:2)) %>% 
  rename("Locus" = value) %>% 
  filter(!str_detect(Locus, "_b"))

#### Change MegaSat notation to be missing genepop calls (000) ####
genotype_data_2 <- genotype_data_1 %>%
  select(-c(1:2)) %>%
  mutate(across(everything(), ~replace(., . ==  0, "000")),
         across(everything(), ~replace(., . ==  "X" , "000")),
         across(everything(), ~replace(., . ==  "Unscored" , "000")),
         across(everything(), ~str_pad(., 3, pad = "0")))

#### Unite the alleles for each locus ####
united_alleles <- genotype_data_1 %>% 
  select(SampleID) %>%   
  mutate(SampleID = paste(genotype_data_1$SampleID, ","))

# 135 because 68 final selected loci
odds <- seq(1, 135, by = 2)

for(i in odds){
  z <- i + 1 
  locus <- colnames(genotype_data_2)[i]
  united_locus <- genotype_data_2 %>% select(i, z) %>% unite({{locus}}, sep = "")
  united_alleles[, ncol(united_alleles) + 1] <- united_locus
}

united_alleles

################### Construct the Genepop file ################################

genotype_matrix <- as.matrix(united_alleles)

#### Make Genepop header ####
file_date <- format(Sys.time(), "%Y%m%d@%H%M") # date and time
header <- paste("Genepop file format", "MCGL2111 ", file_date)

#### List of locus names separated by commas ####
locus_names_2 <- paste(locus_names$Locus, collapse = ",")

#### Generate appropriate "Pop" rows ####
# Pop label that will separate each population
pop_line <- c("Pop", rep("", ncol(genotype_matrix)-1))

#### Count the number of individuals in each population ####
pop_counts <- data.frame(Counts = count(genotype_data_1, Cohort))

#### Add a column totalling the cumulative sum ####
pop_counts <- pop_counts %>% mutate(Sum = cumsum(pop_counts$Counts.n))

#### Insert a Pop row between each population ####
for (i in 1:nrow(pop_counts)){
  # i is the row number and increases by 1 after each iteration to compensate
  # for the extra row being inserted each run through the loop
  pop.row <- rep(NA, nrow(pop_counts))
  pop.row[i] <- pop_counts$Sum[i] + i
  genotype_matrix <- insertRow(genotype_matrix, pop.row[i], pop_line)
}

genotype_matrix <- genotype_matrix[-nrow(genotype_matrix), ]

# Insert title, locus and pop rows at the beginning
genotype_matrix <- insertRow(genotype_matrix, 1, c(header, rep("", ncol(genotype_matrix)-1 )))
genotype_matrix <- insertRow(genotype_matrix, 2, c(locus_names_2, rep("", ncol(genotype_matrix)-1 )))
genotype_matrix <- insertRow(genotype_matrix, 3, pop_line)

# Export file
write.table(genotype_matrix, file = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
