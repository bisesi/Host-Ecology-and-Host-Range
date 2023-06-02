#ATB
#Phage degradation assays 
#Supplemental Figure 2

#load packages
library("tidyverse")
library("readxl")
library("patchwork")

#part A - phage presence after 24 and 48 hours in minimal media without cells
date <- "3May2023"
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))
with_cost <- pfus_and_final_density %>% mutate(cost = "with cost")
with_cost_bacteria <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>% mutate(cost = "with cost")
no_cells <- clean_pfu_data(pfus %>% filter(interaction == "None")) %>%
  inner_join(., plate_layout %>% select(well, timepoint), by = "well")

partA <- no_cells %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(interaction = case_when(interaction == "No cells" ~ "no cells",
                                 TRUE ~ interaction))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(timepoint = case_when(timepoint == 24 ~ "hour 24",
                               timepoint == 48 ~ "hour 48",
                               TRUE ~ timepoint)) %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  #filter(phage != "both\nphage") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~timepoint) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "phage type")+
  ylim(-7, 2)

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

partB_data <- cleaned_pfus_media %>%
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (eh7)",
                                phage == "Specialist phage" ~ "specialist (p22vir)")) %>%
  mutate(media = case_when(media == "rich media" ~ "LB",
                           TRUE ~ media)) %>%
  mutate(date = "1")


#partC - minimal media components
date <- "11April2023"

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
partC_data <- cleaned_pfus %>%
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  filter(media != "none") %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                                   TRUE ~ doublings)) %>%
  mutate(media = case_when(media == "saline" ~ "0.85%\nsaline",
                           media == "H2O" ~ "water",
                           media == "MM-" ~ "metal\nfree",
                           media == "P-" ~ "phos\nfree",
                           media == "S-" ~ "sulf\nfree",
                           media == "P/S-" ~ "phos\n+\nsulf\nfree",
                           TRUE ~ media)) %>%
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (eh7)",
                                phage == "Specialist phage" ~ "specialist (p22vir)"))%>%
  filter(media != "0.85%\nsaline") %>%
  filter(phage != "P22") %>%
  mutate(date = "2") %>%
  select(interaction, phage, well, pfu_E, pfu_S, media, phage_type, date, doublings, doubling_type) %>%
  filter(doubling_type == "phi_doublings") %>%
  pivot_wider(names_from = doubling_type, values_from = doublings)
  

#part D - rich media components
date <- "19April2023"

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
partD_data <- cleaned_pfus %>%
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  filter(media != "none") %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(media = case_when(media == "glucose" ~ "minimal\nmedia",
                           media == "salt + yeast" ~ "salt\n+\nyeast",
                           media == "salt + tryptone" ~ "salt\n+\ntryp",
                           media == "yeast + tryptone" ~ "yeast\n+\ntryp",
                           media == "tryptone" ~ "tryp",
                           TRUE ~ media)) %>%
  mutate(media = factor(media, levels = c("yeast", "salt",
                                          "salt\n+\nyeast", "LB", "salt\n+\ntryp",
                                          "tryp", "yeast\n+\ntryp"))) %>%
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (eh7)",
                                phage == "Specialist phage" ~ "specialist (p22vir)"))%>%
  filter(media != "minimal\nmedia")%>%
  filter(media != "LB") %>%
  filter(phage != "P22") %>%
  mutate(date = "3") %>%
  select(interaction, phage, well, pfu_E, pfu_S, media, phage_type, date, doublings, doubling_type) %>%
  filter(doubling_type == "phi_doublings") %>%
  pivot_wider(names_from = doubling_type, values_from = doublings)
  

#part B visualization
all_partB_data <- rbind(partB_data, partC_data, partD_data)

partB <- all_partB_data %>%
  mutate(media = reorder(media, phi_doublings)) %>%
  group_by(date, media) %>%
  arrange(desc(phi_doublings)) %>%
  ungroup() %>%
  mutate(phi_doublings = case_when(phi_doublings == min(phi_doublings) ~ -5,
                               TRUE ~ phi_doublings)) %>%
  ggplot(aes(x = media, y = phi_doublings)) +
  geom_boxplot(color = eh7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  facet_wrap(~date, scales = "free_x", labeller = labeller(date = 
                                                             c("1" = "",
                                                               "2" = "",
                                                               "3" = "")))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank()) +
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "")+
  ylim(-7, 2)

legend <- get_legend(no_cells %>%
                       mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                                                    TRUE ~ doublings)) %>%
                       mutate(interaction = case_when(interaction == "No cells" ~ "no cells",
                                                      TRUE ~ interaction))%>%
                       mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "Generalist (EH7)",
                                                     phage_type == "Specialist phage" ~ "Specialist (P22*vir*)"))%>%
                       mutate(timepoint = case_when(timepoint == 24 ~ "hour 24",
                                                    timepoint == 48 ~ "hour 48",
                                                    TRUE ~ timepoint)) %>%
                       mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                                                phage == "Phi" ~ "generalist\nonly",
                                                phage == "Phi + P22" ~ "both\nphage")) %>%
                       mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
                       ggplot(aes(x = phage, y = doublings, color = phage_type)) +
                       facet_wrap(~timepoint) +
                       geom_boxplot() +
                       geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
                       theme_bw(base_size = 18)+
                       scale_color_manual(values = c(eh7, p22vir))+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             legend.position = "bottom",
                             legend.text = element_markdown(),
                             panel.grid.minor = element_blank(),
                             axis.title.x = element_blank(),
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("ln(final pfu / initial pfu)")+
                       labs(color = "species")+
                       ylim(-15, 12.5))

#figure 
supp_fig2 <- plot_grid(partA, partB, legend, labels = c("A", "B"), label_size = 26, rel_heights = c(1, 1, 0.1), ncol = 1)
                  

