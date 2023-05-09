#ATB
#Phage degradation assays 
#Supplemental Figure 3

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
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
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
  theme_bw()+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("treatment")+
  labs(color = "phage type")+
  ylim(-10, 12.5)

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
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (phi-C)",
                                phage == "Specialist phage" ~ "specialist (p22vir)"))%>%
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
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("media")+
  labs(color = "")

media_dat <- cleaned_pfus_media

media_model <- aov(phi_doublings ~ media, data = media_dat)

media_comparisons <- TukeyHSD(media_model, conf.levels = 0.95)

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
partC <- cleaned_pfus %>%
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  filter(media != "none") %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(media = case_when(media == "saline" ~ "0.85%\nsaline",
                           media == "H2O" ~ "water",
                           media == "MM-" ~ "metal\nfree",
                           media == "P-" ~ "phos\nfree",
                           media == "S-" ~ "sulf\nfree",
                           media == "P/S-" ~ "phos/sulf\nfree",
                           TRUE ~ media)) %>%
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (phi-C)",
                                phage == "Specialist phage" ~ "specialist (p22vir)"))%>%
  filter(media != "0.85%\nsaline") %>%
  filter(phage != "P22") %>%
  ggplot(aes(x = media, y = doublings)) +
  geom_boxplot(color = "#CA3542") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank()) +
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("media")+
  labs(color = "")

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
partD <- cleaned_pfus %>%
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
  mutate(media = factor(media, levels = c("minimal\nmedia", "yeast", "salt",
                                          "salt\n+\nyeast", "LB", "salt\n+\ntryp",
                                          "tryp", "yeast\n+\ntryp"))) %>%
  mutate(phage_type = case_when(phage == "Generalist phage" ~ "generalist (phi-C)",
                                phage == "Specialist phage" ~ "specialist (p22vir)"))%>%
  filter(media != "LB") %>%
  filter(phage != "P22") %>%
  ggplot(aes(x = media, y = doublings)) +
  geom_boxplot(color = "#CA3542") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank()) +
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("media")+
  labs(color = "")

legend <- get_legend(no_cells %>%
                       mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                                                    TRUE ~ doublings)) %>%
                       mutate(interaction = case_when(interaction == "No cells" ~ "no cells",
                                                      TRUE ~ interaction))%>%
                       mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                                     phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
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
                       theme_bw()+
                       scale_color_manual(values = c("#CA3542", "#27647B"))+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             legend.position = "bottom",
                             panel.grid.minor = element_blank(),
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("growth rate (ln(final pfu / initial pfu))")+
                       xlab("treatment")+
                       labs(color = "phage type")+
                       ylim(-10, 12.5))

#figure 
supp_fig3 <- plot_grid(plot_grid(partA, partB, labels = c("A", "B"), ncol = 1),
                       plot_grid(partC, partD, legend, labels = c("C", "D"), ncol = 1, rel_heights = c(1,1,0.2)))
                  

