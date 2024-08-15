library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(devtools)
#library(cowplot)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(forcats)

Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

nAll(Data_2111) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Marker") %>% 
  as_tibble() %>% 
  rename("n_alleles" = ".") %>% 
  count(n_alleles) %>% 
  summarize(mean_alleles = mean(n_alleles),
            max_alleles = max(n_alleles),
            min_alleles = min(n_alleles))

#### Read in BestConfig files #### (manually cleaned them up in notepad++)
config_2019 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2020/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2020 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2021/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

config_2021 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony_3_runs/Output_files/2022/cleaned_up.BestConfig_Ordered") %>% 
  select(-ClusterIndex)

# Bind them together to make tidy config df
temp_wild <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort == "Wild") %>% 
  rename(OffspringID = SampleID,
         Offspring_cohort = Cohort,
         Offspring_year = Year)

temp_parents <- Samples_2111 %>%
  select(SampleID, Cohort, Year) %>%
  filter(Cohort != "Wild") %>% 
  rename(ParentID = SampleID,
         Parent_cohort = Cohort,
         Parent_year = Year)

config_all <- config_2019 %>% 
  bind_rows(config_2020) %>% 
  bind_rows(config_2021) %>% 
  pivot_longer(cols = c(FatherID, MotherID), 
               names_to = "Parent_type", 
               values_to = "ParentID") %>% 
  left_join(temp_parents) %>% 
  left_join(temp_wild) %>% 
  select(-Parent_type) %>% 
  drop_na()

N_offspring <- config_all %>% 
  group_by(Parent_cohort, Parent_year) %>% 
  summarize(n_offspring = n_distinct(OffspringID)) %>% 
  rename(Cohort = Parent_cohort,
         Year = Parent_year)

test <- config_all %>% 
  distinct(OffspringID)

#write.csv(config_all, "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Analyses/Colony/Best_Config_All.csv")

#### Create info describing... ####
# Whether a fish spawned or not (TRUE, FALSE)
successful_spawners <- config_all %>% 
  select(-c(OffspringID, Offspring_cohort)) %>% 
  group_by(ParentID) %>% 
  summarize(n_spawn_yrs = n_distinct(Offspring_year)) %>% 
  filter(str_detect(ParentID, "-"),
         n_spawn_yrs >= 1)

# How many years each fish spawned
spawned_once <- successful_spawners %>% 
  filter(n_spawn_yrs == 1)

spawned_twice <- successful_spawners %>% 
  filter(n_spawn_yrs == 2)

spawned_thrice <- successful_spawners %>% 
  filter(n_spawn_yrs == 3)

# Whether a fish produced more than 1 offspring (TRUE, FALSE)
polyparents <- config_all %>% 
  group_by(ParentID) %>% 
  count() %>%
  filter(n > 1)

# Bring them into the sample dataframe
Samples_revised <- Samples_2111 %>% 
  mutate(Spawned = case_when(Cohort == "Wild" ~ NA,
                             SampleID %in% successful_spawners$ParentID ~ "TRUE",
                             .default = "FALSE"),
         Polyparent = case_when(SampleID %in% polyparents$ParentID ~ "TRUE",
                                .default = "FALSE"),
         n_spawn_yrs = case_when(Cohort == "Wild" ~ NA,
                                 !(SampleID %in% successful_spawners$ParentID) ~ 0,
                                 SampleID %in% spawned_once$ParentID ~ 1,
                                 SampleID %in% spawned_twice$ParentID ~ 2,
                                 SampleID %in% spawned_thrice$ParentID ~ 3))

######################################################################
#### Begin calculating relative survival and reproductive success ####

# Bring in Mitro's catch data
Mitro_catch <- read_excel("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Offspring_catch_data.xlsx") 

# RELATIVE SURVIVAL: relative survival proportions
Survival <- Mitro_catch %>% 
  select(n_caught, Cohort) %>%
  group_by(Cohort) %>% 
  summarize(fall_recaps = sum(n_caught)) %>% 
  ungroup() %>% 
  mutate(RS = round(fall_recaps / 69, 3)) # Scale to F1s (69)

# Relative reproductive success (likelihood of producing at least one offspring relative to max value) (F2s corrected for sample size)
RRS <- Samples_revised %>% 
  select(Cohort, Spawned) %>% 
  filter(Spawned == "TRUE") %>% 
  group_by(Cohort) %>% 
  count() %>% 
  rename(n_parents = n) %>% 
  ungroup() %>% 
  mutate(subsample = case_when(Cohort == "F1" ~ 450,
                               Cohort == "F2" ~ 400, # F2s have 50 less because 2019 F2s cannot be parents
                               Cohort == "Domestic" ~ 450),
         RRS = round(n_parents / (subsample/max(subsample)) / max(n_parents), 3)) # Scale to max value

# Identify quantities of interannual (IA) parents per cohort #### Not using because can't confirm fish ages
#IA_parent_counts <- Samples_revised %>% 
#  select(Cohort, n_spawn_yrs) %>% 
#  filter(n_spawn_yrs > 1) %>% 
#  group_by(Cohort) %>% 
#  count() %>% 
#  rename(IA_parents = n)

# Join IA parents to RRS df and perform calculation
#RRS <- RRS %>% 
#  left_join(IA_parent_counts) %>% 
#  mutate(IA_parents = case_when(is.na(IA_parents) ~ 0,
#                                .default = IA_parents),
#         IA_RRS = round(IA_parents/ (subsample/max(subsample)) / 8, 3)) # Scale to F1s (8)

# Identify quantities of polyparents per cohort
polyparent_counts <- Samples_revised %>% 
  select(Cohort, Polyparent) %>% 
  filter(Polyparent == "TRUE") %>% 
  group_by(Cohort) %>% 
  count() %>% 
  rename(n_polyparents = n)

# Join polyparents (produced 2 or more offspring) to RRS df and perform calculation
RRS <- RRS %>% 
  left_join(polyparent_counts) %>% 
  mutate(n_polyparents = case_when(is.na(n_polyparents) ~ 0,
                                   .default = n_polyparents),
         Poly_RRS = round(n_polyparents/ (subsample/max(subsample)) / max(n_polyparents), 3)) # Scale to max value

# Implement relative survival to into RRS df
RRS <- RRS %>% 
  left_join(Survival) %>% 
  mutate(SI_RRS = round((RRS / RS), 3),
         SI_RRS_scaled = round((SI_RRS / max(SI_RRS)), 3),
         Cohort = fct_relevel(Cohort, "F1", "Domestic", "F2"))

######################################################
#### Create plots for each value calculated above ####

# Plot relative survival
RS_plot <- RRS %>% 
  filter(Cohort != "F1") %>% 
  ggplot(aes(x = Cohort, y = RS)) + 
  geom_col(fill = "grey50") +
  geom_hline(yintercept = 1, 
             #color = "red",
             linetype = "dashed") +
  geom_text(aes(label = "F1"),
            color = "black",
            #hjust = -0.1,
            x = 2.4, 
            y = 1.04,
            size = 3) +
  geom_text(aes(label = round(RS, 3)),
            color = "black",
            hjust = -0.2,
            size = 3) +
  labs(x = NULL,
       y = "Relative survival",
       title = "(A)") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     breaks = seq(0, 1.2, by = 0.2)) +
  coord_flip() +
  theme_classic()

# Plot RRS
RRS_plot <- RRS %>% 
  filter(Cohort != "F1") %>% 
  ggplot(aes(x = Cohort, y = RRS)) +
  geom_col(fill = "grey50") +
  geom_hline(yintercept = 1, 
             #color = "red",
             linetype = "dashed") +
  geom_text(aes(label = "F1"),
            color = "black",
            #hjust = -0.1,
            x = 2.4, 
            y = 1.04,
            size = 3) +
  geom_text(aes(label = RRS),
            color = "black",
            hjust = -0.2,
            size = 3) +
  labs(x = NULL,
       y = "RRS (\u22651 offspring)",
       title = "(B)") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     breaks = seq(0, 1.2, by = 0.2)) +
  coord_flip() +
  theme_classic()

# Plot interannual RRS and poly RRS (Not plotting because it's just zero for F2 and D)
# IA_RRS_plot <- RRS %>% 
#   ggplot(aes(x = Cohort, y = IA_RRS)) +
#   geom_col() +
#   #geom_text(aes(label = IA_RRS),
#   #          color = "black",
#   #          hjust = -0.1,
#   #          size = 3) +
#   geom_hline(yintercept = 1, 
#              #color = "red",
#              linetype = "dashed") +
#   geom_text(aes(label = "F1"),
#             color = "black",
#             #hjust = -0.1,
#             x = 2.5, 
#             y = 1.02,
#             size = 4,
#             inherit.aes = FALSE) +
#   labs(x = NULL,
#        y = "Interannual RRS",
#        title = "(D)") +
#   scale_y_continuous(expand = c(0, 0),
#                      limits = c(0, 1.2),
#                      breaks = seq(0, 1.2, by = 0.2)) +
#   coord_flip() +
#   theme_classic()

poly_RRS_plot <- RRS %>% 
  filter(Cohort != "F1") %>% 
  ggplot(aes(x = Cohort, y = Poly_RRS)) +
  geom_col(fill = "grey50") +
  geom_hline(yintercept = 1, 
             #color = "red",
             linetype = "dashed") +
  geom_text(aes(label = "F1"),
            color = "black",
            #hjust = -0.1,
            x = 2.4, 
            y = 1.04,
            size = 3) +
  geom_text(aes(label = Poly_RRS),
            color = "black",
            hjust = -0.2,
            size = 3) +
  labs(x = NULL,
       y = "Multiple offspring RRS (\u22652 offspring)",
       title = "(C)") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     breaks = seq(0, 1.2, by = 0.2)) +
  coord_flip() +
  theme_classic()

# Plot survival-independent RRS
SI_RRS_plot <- RRS %>% 
  filter(Cohort != "Domestic") %>% 
  ggplot(aes(x = Cohort, y = SI_RRS_scaled)) +
  geom_col(fill = "grey50") +
  geom_hline(yintercept = 1, 
             #color = "red",
             linetype = "dashed") +
  geom_text(aes(label = "Domestic"),
            color = "black",
            #hjust = -0.1,
            x = 2.4, 
            y = 1.09,
            size = 3) +
  geom_text(aes(label = SI_RRS_scaled),
            color = "black",
            hjust = -0.2,
            size = 3) +
  labs(x = NULL,
       y = "Survival independent RRS",
       title = "(D)") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.2),
                     breaks = seq(0, 1.2, by = 0.2)) +
  #scale_y_continuous(expand = c(0, 0),
  #                   limits = c(0, 14),
  #                   breaks = seq(0, 14, by = 2)) +
  coord_flip() +
  theme_classic()

all_plots <- RS_plot / RRS_plot / poly_RRS_plot / SI_RRS_plot

ggsave(filename = "allplots_horizontal.pdf",
       plot = all_plots,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 6,
       width = 5,
       units = "in")

ggsave(filename = "allplots_horizontal.png",
       plot = all_plots,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 6,
       width = 5,
       units = "in")

######################################################
#### Plot the number of offspring per parent fish ####

# Create F1 df to order the bars in desc order then plot
F1_bar_order <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  filter(Parent_cohort == "F1") %>% 
  arrange(desc(n))

wild_per_F1_plot <- config_all %>% 
  filter(Parent_cohort == "F1") %>% 
  ggplot() +
  geom_bar(aes(x = ParentID 
               #fill = as_factor(Offspring_year)
               ),
          show.legend = FALSE) +
  labs(y = bquote(N[offspring]),
       x = NULL,
       fill = "Year",
       title = "(A)   F1") +
  scale_y_continuous(limits = c(0, 55), 
                     expand = c(0, 0),
                     breaks = seq(0, 55, by = 5)) +
  scale_x_discrete(limits = c(F1_bar_order$ParentID)) +
  #scale_fill_manual(values = c("grey65","grey40", "black")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Create F1 df to order the bars in desc order the plot
F2_bar_order <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  filter(Parent_cohort == "F2") %>% 
  arrange(desc(n))

wild_per_F2_plot <- config_all %>% 
  filter(Parent_cohort == "F2") %>% 
  ggplot() +
  geom_bar(aes(x = ParentID 
               #fill = as_factor(Offspring_year)
               ),
           show.legend = FALSE) +
  labs(y = NULL,
       x = "Inferred parent fish",
       fill = "Year",
       title = "(B)   F2") +
  scale_y_continuous(limits = c(0, 55), 
                     expand = c(0, 0),
                     breaks = seq(0, 55, by = 5)) +
  scale_x_discrete(limits = c(F2_bar_order$ParentID)) +
  #scale_fill_manual(values = c("grey65","grey40", "black")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks.x = element_blank())

# Create domestic df to order the bars in desc order then plot
D_bar_order <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  filter(Parent_cohort == "Domestic") %>% 
  arrange(desc(n))

wild_per_D_plot <- config_all %>% 
  filter(Parent_cohort == "Domestic") %>% 
  ggplot() +
  geom_bar(aes(x = ParentID 
               #fill = as_factor(Offspring_year)
               ),
           show.legend = TRUE) +
  labs(y = NULL,
       x = NULL,
       fill = "Year",
       title = "(C)   Domestic") +
  scale_y_continuous(limits = c(0, 55), 
                     expand = c(0, 0),
                     breaks = seq(0, 55, by = 5)) +
  scale_x_discrete(limits = c(D_bar_order$ParentID)) +
  #scale_fill_manual(values = c("grey65","grey40", "black")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks.x = element_blank())

# Combine the 3 plots and save
wild_per_cohort_plot <- wild_per_F1_plot + wild_per_F2_plot + wild_per_D_plot #+ plot_layout(guides = 'collect'))

ggsave(filename = "Wild_per_cohort3panel.pdf",
       plot = wild_per_cohort_plot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 3,
       width = 7,
       units = "in")

ggsave(filename = "Wild_per_cohort3panel.png",
       plot = wild_per_cohort_plot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 3,
       width = 7,
       units = "in")

#################################################################################
#### Run calculation to find median number of offspring per parent by cohort ####

# Isolate number of offspring per parent for each cohort
Wild_per_parent <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  mutate(n_offspring = n,
         Cohort = fct_relevel(Parent_cohort, "F1", "F2", "Domestic"), .keep = "unused")

Wild_per_parent %>% group_by(Cohort) %>% 
  summarize(n_parents = n(),
            total_offspring_n = sum(n_offspring),
            mean = mean(n_offspring),
            sd = sd(n_offspring),
            median = median(n_offspring))

# Run Kruskal-wallis test (non-parametric anova)
model_1 <- kruskal.test(n_offspring ~ Cohort, Wild_per_parent)

# Plot boxplot
boxplot <- Wild_per_parent %>% 
  ggplot(aes(x = Cohort, y = n_offspring)) +
  geom_boxplot() +
  annotate("text",
           x = "Domestic",
           y = 35,
           label = "p-value = 0.071") +
  labs(x = "Stocked group",
       y = bquote(N[offspring])) +
  scale_y_continuous(limits = c(0,60),
                     expand = c(0,0)) +
  theme_classic()

ggsave(filename = "Wild_per_cohort_boxplot.pdf",
       plot = boxplot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "Wild_per_cohort_boxplot.png",
       plot = boxplot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 5,
       width = 5,
       units = "in")

# Barplot ################## Unnecessary without unknown fish # Not using this single colony run
barplot_df <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  mutate(n_offspring = n,
         Cohort = fct_relevel(Parent_cohort, c("F1", "F2", "Domestic")), .keep = "unused") %>% 
  arrange(desc(n_offspring))

barplot <- ggplot(barplot_df) +
  geom_col(aes(x = ParentID, y = n_offspring, fill = Cohort)) +
  labs(x = "Inferred parent fish",
       y = bquote(N[offspring]),
       fill = "Stocked\ngroup") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,55),
                     breaks = seq(0,55, by = 5)) +
  scale_x_discrete(limits = c(barplot_df$ParentID)) +
  scale_fill_manual(values = c("#1E88E5","#FFC107", "#D81B60")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "Wild_per_cohort_barplot.pdf",
       plot = barplot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_1_run",
       height = 4,
       width = 5,
       units = "in")

ggsave(filename = "Wild_per_cohort_barplot.png",
       plot = barplot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_1_run",
       height = 4,
       width = 5,
       units = "in")

#######################################################################
#### Alternative plotting ideas ############### Not so many barplots...

scatter_df <- config_all %>%
  group_by(ParentID) %>%
  count() %>% 
  left_join(temp_parents) %>% 
  ungroup() %>% 
  mutate(n_offspring = n,
         Cohort = fct_relevel(Parent_cohort, c("F1", "F2", "Domestic")), .keep = "unused") %>% 
  arrange(desc(n_offspring))

scatter_plot <- scatter_df %>% 
  ggplot(aes(x = Parent_year, y = n_offspring)) +
  geom_point(aes(color = Cohort),
             position = position_jitter(h = 0, w = 0.15),
             alpha = 0.5,
             size = 1.5) +
  scale_colour_manual(name = "Stocked group", # Can use this or scale_color_discrete()
                      values = c("#1E88E5","#FFC107", "#D81B60")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,55),
                     breaks = seq(0,55, by = 5)) +
  scale_x_continuous(breaks = seq(2018, 2020, by = 1)) +
  labs(x = "Stocking year",
       y = bquote(N[offspring])) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top")

ggsave(filename = "parental_scatterplot.pdf",
       plot = scatter_plot,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 3,
       width = 5.5,
       units = "in")

ggsave(filename = "parental_scatterplot.png",
       plot = scatter_plot,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures/Colony_3_runs",
       height = 3,
       width = 5,
       units = "in")
