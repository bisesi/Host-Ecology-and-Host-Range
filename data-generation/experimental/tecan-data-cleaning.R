#ATB
#Tecan and PFU analysis - clean and concatenate data
#Using csv files from tecan and xlsx of plaque assay results
#Generates final versions of data for visualization

#load packages
library("tidyverse")
library("readxl")

#source helper functions
source(here::here("functions", "baranyi-helper-functions.R"))
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate paths to import tecan and PFU data
OD_path <- here::here("experimental-data", "tecan-data", date, generate_paths("od", date, "csv"))
CFP_path <- here::here("experimental-data", "tecan-data", date, generate_paths("cfp", date, "csv"))
YFP_path <- here::here("experimental-data", "tecan-data", date, generate_paths("yfp", date, "csv"))
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#import tecan data, plate layout and pfus
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
OD_data <- import_tecan_data(OD_path, file_type = "OD") %>% inner_join(., plate_layout, by = "well")
CFP_data <- import_tecan_data(CFP_path, file_type = "CFP") %>% inner_join(., plate_layout, by = "well")
YFP_data <- import_tecan_data(YFP_path, file_type = "YFP") %>% inner_join(., plate_layout, by = "well")
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#bleedthrough correction
all_tecan_data <- OD_data %>% 
  inner_join(., CFP_data %>% select(c(CFP, cycle, well)), by = c("well", "cycle")) %>% 
  inner_join(., YFP_data %>% select(c(YFP, cycle, well)), by = c("well", "cycle"))
YFP_bleed <- all_tecan_data %>% adjust_FP_vs_FP("Emono", YFP, CFP, all_tecan_data %>% adjust_OD_vs_FP("Emono", CFP) %>% first())
CFP_bleed <- all_tecan_data %>% adjust_FP_vs_FP("Smono", CFP, YFP, all_tecan_data %>% adjust_OD_vs_FP("Smono", YFP) %>% first())
all_tecan_adjusted_FP <- adjust_FP_values(all_tecan_data, YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
  mutate(E_corrected_OD = adjusted_CFP * all_tecan_data %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
  mutate(S_corrected_OD = adjusted_YFP * all_tecan_data %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
  inner_join(., all_tecan_data %>% select(cycle, interaction, phage, well), by = c("well", "cycle"))

#clean pfu data and join with bleedthrough corrected OD data
cleaned_pfus <- clean_pfu_data(pfus)
pfus_and_final_density <- cleaned_pfus %>% inner_join(., all_tecan_adjusted_OD %>%
                                                       select(-c(interaction, phage)) %>% 
                                                      group_by(well) %>% slice_max(cycle), by = "well")

