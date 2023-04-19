#ATB
#Exploratory data analysis
#E or Smono in lactose with phage 11March2023

#load packages
library("tidyverse")
library("readxl")

#set date
date <- "11April2023"

#source tecan data cleaning script
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate paths
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#load data in
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#clean pfu data
cleaned_pfus <- clean_pfu_data(pfus %>% 
                                 mutate(interaction = case_when(condition != "Start" ~ "no cells", 
                                                                TRUE ~ interaction))) %>%
  mutate(interaction = case_when(is.na(interaction) ~ "No Cells",
                                 TRUE ~ interaction))

#visualize pfus
cleaned_pfus %>%
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  filter(media != "none") %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~media) +
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
