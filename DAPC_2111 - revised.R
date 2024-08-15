library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(ggrepel)

# DAPC tests a hypothesis, PCA does not

# DAPC guidelines from Thia 2022:
# n.da = k groups (should be determined a priori, # of sample pops)
# n.pca must be =< k-1 (only k-1 PCs are biologically informative)

#### Read in genetic data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

#### Read in 2111 metadata ####
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "St. Croix",
                                   .default = WaterbodyName))

#### subset to genetic data by Cohort ####
Data_2111@pop <- as_factor(Samples_2111$Cohort)

Stocked_fish <- popsub(Data_2111, 
                       sublist = c("F1", "F2", "Domestic"),
                       drop = FALSE)

Wild_offspring <- popsub(Data_2111, 
                         sublist = c("Wild"),
                         drop = FALSE)

#### Run DAPC using just stocked fish ####
set.seed(27)
clusters_stocked <- find.clusters.genind(Stocked_fish, 
                                         max.n.clust = 5,
                                         n.pca = 250,
                                         n.clust = 3)
                                         
# Select 250 to retain all PCs
# Selecting 3 clusters, lowest BIC

table(Stocked_fish$pop, clusters_stocked$grp)

DAPC_stocked <- dapc.genind(Stocked_fish, 
                            #pop = clusters_hybrids$grp,
                            #n.clust = 4,
                            n.pca = nPop(Stocked_fish) - 1,
                            n.da = nPop(Stocked_fish))

summary(DAPC_stocked)

# Create figure for reassignment percentages
reassignment_plot <- summary(DAPC_stocked)$assign.per.pop %>% 
  as_tibble(rownames = "Cohort") %>% 
  mutate(reassignment = round(value*100, 0),
         Cohort = fct_relevel(Cohort, "F1", "F2", "Domestic")) %>% 
  ggplot(aes(x = Cohort, y = reassignment)) +
  geom_col() +
  coord_flip() +
  labs(x = "Stocked group",
       y = "% reassignment to correct group") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(plot.margin = margin(t = 5,  # Top margin
                             r = 10,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))

ggsave(filename = "reassignment_plot.pdf",
       plot = reassignment_plot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/",
       height = 4,
       width = 6,
       units = "in")

ggsave(filename = "reassignment_plot.png",
       plot = reassignment_plot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/",
       height = 4,
       width = 6,
       units = "in")
  
A_score_stocked <- optim.a.score(DAPC_stocked) # Optimal n.pca = 2

# Bring wild offspring in as supplementary individuals to assign back to stocked cohorts
prediction <- predict.dapc(DAPC_stocked, newdata = Wild_offspring)

prediction$assign
prediction$posterior
prediction$ind.scores

Assignment_probs <- prediction$posterior %>%
  as_tibble(rownames = "SampleID") %>% 
  pivot_longer(cols = 2:4, 
               names_to = "Cluster", 
               values_to = "Probability")

# Plot all assignment probabilities
Ass_prob_plot <- Assignment_probs %>% 
  mutate(Cluster = fct_relevel(Cluster, "F1", "F2", "Domestic")) %>% 
  ggplot(aes(x = SampleID, y = Probability, fill = Cluster)) +
  geom_col(color = "darkgrey", linewidth = 0.1) +
    #facet_grid(#~fct_inorder(Offspring_year), 
    #           switch = "x", 
    #           scales = "free", 
    #           space = "free") +
  labs(x = "Wild-caught offspring (n = 258)", 
       y = "Assignment probability",
       fill = "Stocked group") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("#1E88E5","#FFC107", "#D81B60")) +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top")

ggsave(filename = "DAPC_ass_prob_plot.pdf",
       plot = Ass_prob_plot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 8,
       units = "in")

ggsave(filename = "DAPC_ass_prob_plot.png",
       plot = Ass_prob_plot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 8,
       units = "in")

#### Try plotting like PCA ####
# Make tibble from original DAPC results
row_names <- row.names(DAPC_stocked$grp.coord)

dapc_centroids <- DAPC_stocked$grp.coord %>% 
  as_tibble()

dapc_centroids$Cohort <- row_names

stocked_centroids <- dapc_centroids %>% 
  select(Cohort) %>% 
  right_join(aggregate(cbind(LD1, LD2) ~ Cohort, dapc_centroids, mean)) %>% 
  distinct(.keep_all = TRUE) #%>% 
  #mutate(Cohort = str_replace_all(Cohort, "_", " "), .keep = "unused")

# Make tibble from prediction results
row_names <- row.names(prediction$ind.score)

wild_data_points <- prediction$ind.scores %>% 
  as_tibble()

wild_data_points$SampleID <- row_names

wild_data_points <- wild_data_points %>% 
  left_join(Samples_2111)

# Make tibble for stocked cohort ellipses 
row_names <- row.names(DAPC_stocked$ind.coord)

stocked_data_points <- DAPC_stocked$ind.coord %>% 
  as_tibble()

stocked_data_points$SampleID <- row_names

stocked_data_points <- stocked_data_points %>% 
  left_join(Samples_2111) %>% 
  mutate(Cohort = fct_relevel(Cohort, "F1", "F2", "Domestic"))

# Plot
DAPC_ord_plot <- stocked_data_points %>% 
  ggplot(aes(x = LD1, y = LD2, color = Cohort)) + 
  geom_point(alpha = 0.25,
             show.legend = FALSE) +
  #geom_point(data = stocked_centroids,
  #           size = 6) +
             #color = Cohort) +
  stat_ellipse(data = stocked_data_points,
               size = 0.75,
               alpha = 0.9,
               show.legend = FALSE) +
  geom_point(data = wild_data_points,
             color = "black") +
  annotate("text",
           label = "Wild-caught offspring",
           x = -1.5,
           y = 2.1) +
  annotate("text",
           label = "F1",
           x = 4,
           y = -2.5,
           color = "#1E88E5") +
  annotate("text",
           label = "F2",
           x = 4,
           y = 3,
           color = "#FFC107") +
  annotate("text",
           label = "Domestic",
           x = -5,
           y = -1.5,
           color = "#D81B60") +
  labs(x = "Axis 1",
       y = "Axis 2",
       color = "Stocked group") +
  scale_color_manual(values = c("#1E88E5","#FFC107", "#D81B60")) +
  theme_classic() + 
  theme(legend.position = "top")

ggsave(filename = "DAPC_ord_by_cohort.pdf",
       plot = DAPC_ord_plot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 8,
       units = "in")

ggsave(filename = "DAPC_ord_by_cohort.png",
       plot = DAPC_ord_plot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 8,
       units = "in")

##################################################### Not feasible ##############
#### Re-run DAPC by source stream ####
Data_2111@pop <- as_factor(Samples_2111$WaterbodyName)

Stocked_fish <- popsub(Data_2111, 
                       exclude = "Strutt Creek",
                       drop = FALSE)

Wild_offspring <- popsub(Data_2111, 
                         sublist = "Strutt Creek",
                         drop = FALSE)

#### Run DAPC using just stocked fish ####
set.seed(27)
clusters_stocked <- find.clusters.genind(Stocked_fish, 
                                         max.n.clust = 12,
                                         n.pca = 250,
                                         n.clust = 4)

# Select 250 to retain all PCs
# Selecting 4 clusters, lowest BIC

table(Stocked_fish$pop, clusters_stocked$grp)

DAPC_stocked <- dapc.genind(Stocked_fish, 
                            #pop = clusters_hybrids$grp,
                            #n.clust = 4,
                            n.pca = nPop(Stocked_fish) - 1,
                            n.da = nPop(Stocked_fish))

summary(DAPC_stocked)

A_score_stocked <- optim.a.score(DAPC_stocked) # Optimal n.pca = 6

# Bring wild offspring in as supplementary individuals to assign back to stocked cohorts
prediction <- predict.dapc(DAPC_stocked, newdata = Wild_offspring)

prediction$assign
prediction$posterior
prediction$ind.scores

Assignment_probs <- round(prediction$posterior, 2) %>%
  as_tibble(rownames = "Offspring_ID") %>% 
  #mutate(Offspring_year = Wild_offspring@pop) %>% 
  pivot_longer(cols = 2:8, 
               names_to = "Cluster", 
               values_to = "Probability")

# Plot all assignment probabilities
Assignment_probs %>% 
  ggplot(aes(Offspring_ID, Probability, fill = factor(Cluster))) +
  geom_col(color = "gray", linewidth = 0.1) +
  #facet_grid(~fct_inorder(Offspring_year), 
  #           switch = "x", 
  #           scales = "free", 
  #           space = "free") +
  labs(x = "Offspring", y = "Assignment probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_discrete("Dark2") +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  guides(fill = guide_legend(title = "Genetic lineage")) +
  theme(strip.text.x = element_text(angle = -90))

#### Plotting like PCA ####
# Make tibble from prediction results
row_names <- row.names(prediction$ind.score)

prediction_data <- prediction$ind.scores %>% 
  as_tibble()

prediction_data$SampleID <- row_names

prediction_data <- prediction_data %>% 
  left_join(Samples_2111)

# Make tibble from original DAPC results
row_names <- row.names(DAPC_stocked$grp.coord)

dapc_data <- DAPC_stocked$grp.coord %>% 
  as_tibble()

dapc_data$WaterbodyName <- row_names

# Plot
centroids <- dapc_data %>% 
  select(WaterbodyName) %>% 
  right_join(aggregate(cbind(LD1, LD2) ~ WaterbodyName, dapc_data, mean)) %>% 
  distinct(.keep_all = TRUE)

prediction_data %>% 
  ggplot(aes(x = LD1, y = LD2, color = WaterbodyName)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids, 
                  aes(label = WaterbodyName), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")
