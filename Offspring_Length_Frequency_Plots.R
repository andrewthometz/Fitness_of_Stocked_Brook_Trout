# Load packages
library(tidyverse)
library(readxl)
library(devtools)
devtools::install_github("thomasp85/patchwork")
library(patchwork)

###############################################################################
#### Create length-frequency distributions for each year of wild offspring ####
###############################################################################

# Wild offspring cohorts indicate the year a fish was believed to have been naturally produced
# For example, wild offspring captured in 2020 are in the W_2019 cohort because the natural reproduction would have occured in Fall of 2019

# Read in metadata
Offspring_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort_year == "W_2019" |
         Cohort_year == "W_2020" |
         Cohort_year == "W_2021") %>% 
  unite(col = "Capture_date", c(MonthCollected, DayCollected, YearCollected), sep = "-") %>% 
  mutate(Capture = case_when(Capture_date == "Sep-8-2020" ~ "Fall 2020",
                             Capture_date == "Apr-5-2021" ~ "Spring 2021",
                             Capture_date == "Sep-8-2021" ~ "Fall 2021",
                             Capture_date == "Apr-4-2022" ~ "Spring 2022",
                             Capture_date == "Oct-6-2022" ~ "Fall 2022")) # Simplify capture date to season + year

Offspring_2111 %>% count(Capture_date)
Offspring_2111 %>% count(Cohort)

# Plot for offspring from W_2019 cohort (We don't have length data for fish captured on 4/5/2020 (n = 28)... not a big deal though)
plot_2019 <- Offspring_2111 %>% 
  filter(Capture_date == "Sep-8-2020") %>%
  ggplot(aes(x = TotalLength, fill = Capture)) +
  geom_histogram(binwidth = 20,
                 center = 10,
                 color = "black",
                 linewidth = 0.25) +
  labs(x = NULL,
       y = "Count",
       title = "(A)   Colony Run 1 (n = 145)") +
  scale_y_continuous(limits = c(0,70),
                     expand = c(0,0),
                     breaks = seq(0, 70, by = 10)) +
  scale_x_continuous(limits = c(40,200),
                     expand = c(0,0),
                     breaks = seq(40, 200, by = 20)) +
  scale_fill_manual(values = c("grey65")) +
  theme_classic() +
  theme(legend.position = "top",
        plot.margin = margin(0, 10, 0, 0))

# Plot for offspring from W_2020 cohort
plot_2020 <- Offspring_2111 %>% 
  filter(Cohort_year == "W_2020") %>% 
  ggplot(aes(x = TotalLength, fill = Capture)) +
  geom_histogram(binwidth = 20,
                 center = 10, 
                 color = "black",
                 linewidth = 0.25) +
  labs(x = "Total length (mm)",
       y = NULL,
       title = "(B)   Colony Run 2 (n = 72)") +
  scale_y_continuous(limits = c(0,20),
                     expand = c(0,0),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_continuous(limits = c(40,280),
                     expand = c(0,0),
                     breaks = seq(40, 280, by = 20)) +
  scale_fill_manual(values = c("grey65","grey40")) +
  theme_classic() +
  theme(legend.position = "top",
        plot.margin = margin(0, 10, 0, 0))

# Plot for offspring from W_2021 cohort
plot_2021 <- Offspring_2111 %>% 
  filter(Cohort_year == "W_2021") %>% 
  ggplot(aes(x = TotalLength, fill = Capture)) +
  geom_histogram(binwidth = 20,
                 center = 10,
                 color = "black",
                 linewidth = 0.25) +
  labs(x = NULL,
       y = NULL,
       title = "(C)   Colony Run 3 (n = 41)") +
  scale_y_continuous(limits = c(0,20),
                     expand = c(0,0),
                     breaks = seq(0, 20, by = 5)) +
  scale_x_continuous(limits = c(40,200),
                     expand = c(0,0),
                     breaks = seq(40, 200, by = 20)) +
  scale_fill_manual(values = c("grey65")) +
  theme_classic() +
  theme(legend.position = "top",
        plot.margin = margin(0, 10, 0, 0))

# Join the three plots together using the patchwork package
length_freqs <- plot_2019 + plot_2020 + plot_2021

ggsave(filename = "Wild_length_frequencies.pdf",
       plot = length_freqs,
       device = "pdf",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 4,
       width = 10,
       units = "in")

ggsave(filename = "Wild_length_frequencies.png",
       plot = length_freqs,
       device = "png",
       path = "X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Polished_plots_figures",
       height = 4,
       width = 10,
       units = "in")
