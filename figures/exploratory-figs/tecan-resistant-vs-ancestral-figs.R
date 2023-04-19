#ATB
#Exploratory data analysis
#Resistance ancestral S versus B2 high density colony isolates from 15March2023

#load packages
library("tidyverse")
library("readxl")

#set date
date <- "15March2023"

#source tecan data cleaning script
source(here::here("functions", "baranyi-helper-functions.R"))
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate 
OD_path <- here::here("experimental-data", "tecan-data", date, "resistant-isolates-with-ancestral", generate_paths("od", paste(date, "3", sep = "_"), "csv"))
pfu_path <- here::here("experimental-data", "tecan-data", date, "resistant-isolates-with-ancestral", generate_paths("pfu", paste(date, "3", sep = "_"), "xlsx"))

#load the data in
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "resistant-isolates-with-ancestral", "plate_layout_3.csv"))
OD_data<- import_tecan_data(OD_path, file_type = "OD") %>% inner_join(., plate_layout, by = "well")
resistance_data <- load_pfu_data(here::here("experimental-data", "tecan-data", "15March2023", "resistant-isolates-with-ancestral","resistance_15March2023_3.xlsx"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#clean data
cleaned_OD <- OD_data %>% filter(interaction != "none") 
growth_rates<- generate_baranyi_growth_data_from_OD(cleaned_OD) %>%
  inner_join(., plate_layout, by = "well")

wide_data <- pfus %>% 
  select(interaction, phage, plate, pfu, well) %>%
  pivot_wider(names_from = plate, values_from = pfu) %>%
  rename(pfu_S = S, pfu_E = E)

starting_phage <- wide_data %>%
  filter(interaction == "None" & is.na(well)) %>%
  mutate(ratio = pfu_E / (pfu_E + pfu_S))

cleaned_pfus <- wide_data %>%
  filter(phage == "Phi" & !is.na(well)) %>%
  mutate(phi_doublings = log(pfu_E / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2)) %>%    
  mutate_all(~replace(., is.infinite(.), log(1 / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2))) %>%
  mutate(media = case_when(well %in% c("B12", "C12", "D12") ~ "minimal_media",
                           well %in% c("E12", "F12", "G12") ~ "rich_media"))

OD_resistance <- resistance_data %>%
  inner_join(., OD_data %>% group_by(well) %>% slice_max(cycle) %>% ungroup(), by = "well")

#visualize
cleaned_OD %>% 
  unite("condition", interaction:phage, remove = FALSE, sep = " ") %>%
  ggplot(aes(x = cycle, y = OD, color = phage)) +
  geom_smooth(se = FALSE, span = 0.2) +
  facet_wrap(~interaction, ncol = 4)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

growth_rates %>% 
  filter(fit_variable != "y0") %>%
  ggplot(aes(x = interaction, y = model_fit, color = phage)) +
  geom_point() +
  facet_wrap(~fit_variable, ncol = 4)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

cleaned_pfus %>% ggplot(aes(x = media, y = phi_doublings)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

resistance_data %>% 
  unite("condition", Interaction:Phage, remove = FALSE, sep = " ") %>%
  pivot_longer(cols = ends_with("Resistant"), names_to = "phage_type", values_to = "percentage_resistant") %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x = Phage, y = percentage_resistant, color = phage_type)) +
  geom_boxplot() +
  facet_wrap(~Interaction)

OD_resistance %>% 
  ggplot(aes(x = OD, y = Phi_Resistant, color = Phage)) + 
  geom_point() + 
  facet_wrap(~Interaction)
  

