library(tidyverse)
library(readxl)
library(devtools)
#install_github(repo = "nicksard/mater")
library(mater)
library(adegenet)

#### Prepare mating matrix input file ####
# Create breeding matrix
matrix_1 <- brd.mat(moms = 50, dads = 50, lambda.low = 4, lambda.high = 4) # lambda of 4 = Poisson distribution?

# Define the number of offspring produced by each mate pair
matrix_1_offspring <- brd.mat.fitness(mat = matrix_1, min.fert = 250, max.fert = 900, type = "uniform")
                                                  # This fecundity info is from McFadden (1961)

# To get basic stats of the matrix
mat.stats(mat = matrix_1_offspring, id.col = FALSE)

# Take subsample of the matrix
subsample_df <- mat.sub.sample(mat = matrix_1_offspring, noff = 100)
head(subsample_df)

# Convert the subsample df to a "typical three column pedigree"
pedigree_1 <- convert2ped(df = subsample_df)
head(pedigree_1)

# Convert the pedigree back into a matrix for COLONY and write as .txt
subsample_matrix <- ped2mat(ped = pedigree_1)
mat.stats(mat = subsample_matrix)

subsample_matrix %>% 
write.table("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_pedigree_simulation/Mating_matrix.txt", 
            #sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)

#### Create marker input file ####
# Read in genind file
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

library(OncoBayes2) # This packages has the bind_rows_0() function

Marker_types <- tibble(x = 0, y = locNames(Data_2111)) %>% pivot_wider(names_from = y, values_from = x)
Allelic_dropout <- tibble(x = 0, y = locNames(Data_2111)) %>% pivot_wider(names_from = y, values_from = x)
False_allele <- tibble(x = 0.02, y = locNames(Data_2111)) %>% pivot_wider(names_from = y, values_from = x)

Allele_counts <- nAll(Data_2111) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Marker") %>% 
  as_tibble() %>% 
  rename("n_alleles" = ".") %>% 
  pivot_wider(names_from = Marker, values_from = n_alleles) 

# Average number of alleles per locus
# nAll(Data_2111) %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "Marker") %>% 
#   as_tibble() %>% 
#   rename("n_alleles" = ".") %>% 
#   count(n_alleles) %>% 
#   summarize(mean_alleles = mean(n_alleles))

Marker_types %>%
  bind_rows(Allele_counts) %>% 
  bind_rows_0(Allelic_dropout) %>% 
  bind_rows_0(False_allele) %>% 
  write.table("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_pedigree_simulation/Marker_input.txt", 
              #sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE,
              quote = FALSE,
              eol = "\n")

#### Gather allele frequency data ####
#library(PopGenReport)
#library(poppr)

ProgressBar <- txtProgressBar(min = 0, max = 68, style = 3)
for (i in 1:68){

df <- allele.dist(Data_2111, mk.figures = FALSE)$frequency[[i]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Allele") %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  pivot_longer(cols = c(2:ncol(.)), names_to = "Pop", values_to = "freq") %>% 
  group_by(Allele) %>% 
  dplyr::summarize(mean_freq = round(mean(freq), digits = 3)) %>%  # call dplyr to make sure summarize is recognizing groupings
  pivot_wider(names_from = Allele, values_from = mean_freq)

  vector <- as.numeric(df[1,])
  
  vector %>% write(file = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_pedigree_simulation/Allele_frequencies_new.txt",
                   append = TRUE,
                   ncolumns = length(vector))
  
  setTxtProgressBar(ProgressBar, i)
}
close(ProgressBar)
