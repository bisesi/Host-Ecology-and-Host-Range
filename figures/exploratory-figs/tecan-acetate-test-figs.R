#ATB
#Exploratory data analysis
#Smono acetate test with high density gluc wells 15March2023

#load packages
library("tidyverse")
library("readxl")

#set date
date <- "15March2023"

#source tecan data cleaning script
source(here::here("functions", "baranyi-helper-functions.R"))
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate 
OD_path_first <- here::here("experimental-data", "tecan-data", date, generate_paths("od", paste(date, "1", sep = "_"), "csv"))
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
OD_path_second <- here::here("experimental-data", "tecan-data", date, generate_paths("od", paste(date, "2", sep = "_"), "csv"))

#load data in
plate_layout_first <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout_1.csv"))
OD_data_first <- import_tecan_data(OD_path_first, file_type = "OD") %>% inner_join(., plate_layout_first, by = "well")
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()
plate_layout_second <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout_2.csv"))
OD_data_second <- import_tecan_data(OD_path_second, file_type = "OD") %>% inner_join(., plate_layout_second, by = "well")

#clean data
cleaned_pfus <- clean_pfu_data(pfus)
cleaned_OD_first <- OD_data_first %>% filter(interaction != "none") 
pfus_and_final_density <- cleaned_pfus %>% inner_join(., cleaned_OD_first %>% slice_max(cycle) %>% select(cycle, well, OD, media), by = "well")
cleaned_OD_second <- OD_data_second %>% filter(interaction != "none") 
growth_rates_second <- generate_baranyi_growth_data_from_OD(cleaned_OD_second) %>%
  inner_join(., plate_layout_second, by = "well")

#visualize
cleaned_OD_first %>% 
  unite("condition", interaction:media, remove = FALSE, sep = " ") %>%
  ggplot(aes(x = cycle, y = OD, color = well)) +
  geom_smooth(se = FALSE, span = 0.2) +
  facet_wrap(~condition, ncol = 2)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())

pfus_and_final_density %>% 
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                              TRUE ~ doublings)) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_grid(media~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("number of doublings")+
  xlab("treatment condition")+
  labs(color = "")

cleaned_OD_second %>% 
  unite("condition", interaction:phage, remove = FALSE, sep = " ") %>%
  ggplot(aes(x = cycle, y = OD, color = phage)) +
  geom_smooth(se = FALSE, span = 0.2) +
  facet_wrap(~interaction, ncol = 3)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

growth_rates_second %>% 
  filter(fit_variable != "lag" & fit_variable != "y0") %>%
  ggplot(aes(x = interaction, y = model_fit, color = phage)) +
  geom_point() +
  facet_wrap(~fit_variable, ncol = 3)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())




