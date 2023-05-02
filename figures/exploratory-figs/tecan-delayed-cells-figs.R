#ATB
#Exploratory data analysis
#E or Smono in lactose with phage 11March2023

#load packages
library("tidyverse")
library("readxl")

#set date
date <- "28April2023"

#source tecan data cleaning script
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate paths
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#load data in
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#clean pfu data
wide_data <- pfus %>% 
  mutate(interaction = case_when(interaction == "None" & condition != "Start" ~ "no cells", TRUE ~ interaction)) %>%
  select(interaction, phage, plate, pfu, well) %>%
  pivot_wider(names_from = plate, values_from = pfu) %>%
  rename(pfu_E = E)

starting_phage <- wide_data %>%
  filter(interaction == "None") 

output <- wide_data %>%
  filter(interaction != "None") %>%
  mutate(phi_doublings = case_when(phage == "Phi" ~ log(pfu_E / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2))) %>%    
  mutate_all(~replace(., is.nan(.), 0)) %>%
  mutate_all(~replace(., is.infinite(.), log(1 / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2))) %>%
  mutate(interaction = case_when(interaction == "Smono" ~ "S Monoculture",
                                 interaction == "Emono" ~ "E Monoculture",
                                 interaction == "E Fac" ~ "E Facilitation",
                                 interaction == "Coop" ~ "Mutualism",
                                 interaction == "Comp" ~ "Competition",
                                 interaction == "Fac" ~ "Facilitation")) %>%
  pivot_longer(cols = ends_with("doublings"),
               names_to = "doubling_type",
               values_to = "doublings") %>%
  mutate(phage_type = case_when(doubling_type == "p22_doublings" ~ "Specialist phage",
                                doubling_type == "phi_doublings" ~ "Generalist phage",
                                TRUE ~ "NA")) %>%
  mutate(phage_interaction = case_when(phage == "P22" ~ "Specialist phage",
                                       phage == "Phi" ~ "Generalist phage",
                                       phage == "Phi + P22" ~ "Phage Competition")) %>%
  mutate(interaction = case_when(is.na(interaction) ~ "No Cells",
                                 TRUE ~ interaction))

#visualize pfus
output %>%
  inner_join(., plate_layout %>% select(well, timepoint), by = "well") %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(timepoint = factor(timepoint, levels = c(0, 2, 4, 8, 12, 24, 30))) %>%
  ggplot(aes(x = timepoint, y = doublings, color = interaction)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("growth rate (log2)")+
  xlab("hours incubating prior to addition of cells")+
  labs(color = "")