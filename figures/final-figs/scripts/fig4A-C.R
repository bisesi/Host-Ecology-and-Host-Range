#ATB
#Tecan data figures
#Fig 4 - A, B, C

#load packages
library("tidyverse")
library("readxl")
library("patchwork")

#set date
date <- "13March2023"

#source tecan data cleaning script
source(here::here("functions", "tecan-data-helper-functions.R"))

#generate paths
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#load data in
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#clean pfu data
cleaned_pfus <- clean_pfu_data(pfus %>% 
                                 mutate(interaction = case_when(condition %in% c("C", "D") ~ "no cells", 
                                                                TRUE ~ interaction))) %>%
  mutate(interaction = case_when(is.na(interaction) ~ "No Cells",
                                 TRUE ~ interaction))

#part A
partA <- cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist",
                                                        phage_type == "Specialist phage" ~ "specialist"))%>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "e monoculture",
                                 interaction == "S Monoculture" ~ "s monoculture",
                                 interaction == "No Cells" ~ "no cells"))%>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(phage != "Phi + P22") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  ylab("growth rate (log2) ")+
  xlab("treatment condition")+
  labs(color = "phage type")

phage_dat <- cleaned_pfus

phage_model <- aov(doublings ~ interaction * phage * doubling_type, data = phage_dat)

phage_comparisons <- TukeyHSD(phage_model, conf.levels = 0.95)


#partB
date <- "15March2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, "resistant-isolates-with-ancestral", generate_paths("pfu", paste(date, "3", sep = "_"), "xlsx"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

wide_data <- pfus %>% 
  select(interaction, phage, plate, pfu, well) %>%
  pivot_wider(names_from = plate, values_from = pfu) %>%
  rename(pfu_S = S, pfu_E = E)

starting_phage <- wide_data %>%
  filter(interaction == "None" & is.na(well)) %>%
  mutate(ratio = pfu_E / (pfu_E + pfu_S))

cleaned_pfus_media <- wide_data %>%
  filter(phage == "Phi" & !is.na(well)) %>%
  mutate(phi_doublings = log(pfu_E / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2)) %>%    
  mutate_all(~replace(., is.infinite(.), log(1 / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull, 2))) %>%
  mutate(media = case_when(well %in% c("B12", "C12", "D12") ~ "minimal media",
                           well %in% c("E12", "F12", "G12") ~ "rich media"))

partB <- cleaned_pfus_media %>%
  ggplot(aes(x = media, y = phi_doublings)) +
  geom_boxplot(color = "#CA3542") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank()) +
  ylab("growth rate (log2)")+
  xlab("media")+
  labs(color = "")

media_dat <- cleaned_pfus_media

media_model <- aov(phi_doublings ~ media, data = media_dat)

media_comparisons <- TukeyHSD(media_model, conf.levels = 0.95)


#partC

#figure 4
partA + partB + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

