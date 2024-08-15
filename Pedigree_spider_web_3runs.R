library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(devtools)
#library(cowplot)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(forcats)

#### Create spiderweb ("pedigree plot") ####

Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

#### Read in BestConfig files #### (manually cleaned them up in notepad++)
config_2020 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony/Output_files/2020/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2021 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony/Output_files/2021/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2022 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony/Output_files/2022/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

# Bind them together to make tidy config df
# temp_wild <- Samples_2111 %>%
#   select(SampleID, Cohort, Year) %>%
#   filter(Cohort == "Wild") %>% 
#   rename(OffspringID = SampleID,
#          Offspring_cohort = Cohort,
#          Offspring_year = Year)

# Create MotherID and FatherID, (need to use this language for plotting script)
temp_parents_1 <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort != "Wild") %>% 
  rename(FatherID = SampleID,
         Parent_cohort_1 = Cohort,
         Parent_year_1 = Year)

temp_parents_2 <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort != "Wild") %>% 
  rename(MotherID = SampleID,
         Parent_cohort_2 = Cohort,
         Parent_year_2 = Year)

config_all <- config_2020 %>% 
  bind_rows(config_2021) %>% 
  bind_rows(config_2022) %>% 
  #rename(Parent_1 = FatherID,
  #       Parent_2 = MotherID) %>% 
  left_join(temp_parents_1) %>% 
  left_join(temp_parents_2) %>% 
  mutate(Parent_cohort_1 = case_when(is.na(Parent_cohort_1) ~ "Unknown",
                                     .default = Parent_cohort_1),
         Parent_cohort_2 = case_when(is.na(Parent_cohort_2) ~ "Unknown",
                                     .default = Parent_cohort_2)) %>% 
  select(-c(Parent_year_1, Parent_year_2))

test <- config_all %>% 
  filter(MotherID %in% FatherID)

#############################################
#### Bring in code for pedigree plotting ####

pedigree.plot <- function(family,title = "Pedigree Plot", cohortbox=T){
  if(cohortbox == T) {
    family <- family[order(family$cohort,family$FatherID,family$MotherID),]
  } else {
    family <- family[order(family$FatherID,family$MotherID),]
  }
  moms <- unique(family$MotherID)
  dads <- unique(family$FatherID)
  
  momdots <- data.frame(ID = moms, y = seq(from=1,to=length(family$OffspringID),length.out=length(moms)))
  daddots <- data.frame(ID = dads, y = seq(from=1,to=length(family$OffspringID),length.out=length(dads)))
  
  #make the plot points
  plot(x = c(1,2,3), 
       y = c(0, length(family$OffspringID), length(family$OffspringID)), 
       pch = "", axes = FALSE, 
       xlab = "", ylab = "",
       main = title)
  points(x = rep(2,length(family$OffspringID)), y = c(1:length(family$OffspringID)), pch = "-", cex = 0.25)
  points(x = rep(1,length(moms)), y = momdots$y, pch = 19, cex = 0.5)
  points(x = rep(3,length(dads)), y = daddots$y, pch = 19, cex = 0.5)
  
  #make the lines in the plot
  for(r in 1:length(family$OffspringID)) {
    #mom line
    lines(x=c(2,1), y = c(r, momdots$y[which(family$MotherID[r] == momdots$ID)]), lwd = 0.5,col = "grey")
    
    #dad line
    lines(x =c(2,3), y = c(r, daddots$y[which(family$FatherID[r] == daddots$ID)]), lwd = 0.5,col = "grey")
  }
  
  #adding labels
  mtext("Parent 1", side = 3, line = 0, at = 1, cex = 0.9)
  mtext("Parent 2", side = 3, line = 0, at = 3, cex = 0.9)
  mtext("Offspring (n = 258)", side = 3, line = 0, at = 2, cex = 0.9)
  #adding rectangles (optional)
  if(cohortbox == T){
    if(length(unique(family$loc > 1))){
      cohorts <- unique(family$loc)
      for (r in 1:length(cohorts)) {
        c <- cohorts[r]
        rect(xleft=1.75,
             xright = 2.25,
             ybottom = max(which(family$clust == c)), 
             ytop = min(which(family$clust == c)),
             border = "black", lwd = 2)
        text(x = 1.7, y = mean(which(family$clust == c)), c, pos =2, col = "black")
      }
    }
    cohorts <- unique(family$cohort)
    for (r in 1:length(cohorts)) {
      c <- cohorts[r]
      rect(xleft=1.75,
           xright = 2.25,
           ybottom = max(which(family$cohort == c)), 
           ytop = min(which(family$cohort == c)),
           border = "black", lwd = 2)
      text(x = 1.7, y = mean(which(family$cohort == c)), paste("inferred\n",c), pos =2, col = "black")
    }
  }
}

pedigree.plot(config_all, title = "test_title", cohortbox = F)
