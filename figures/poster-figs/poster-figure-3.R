#ATB
#poster figure 3
#all experimental results

#library
library("patchwork")
library("cowplot")

#load all data
date <- "8March2023"
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))
no_cost <- pfus_and_final_density %>% mutate(cost = "no cost")
no_cost_bacteria <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>% mutate(cost = "no cost")

date <- "3May2023"
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))
with_cost <- pfus_and_final_density %>% mutate(cost = "with cost")
with_cost_bacteria <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>% mutate(cost = "with cost")

date <- "13March2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus_starved <- load_pfu_data(pfu_path) %>% janitor::clean_names()
cleaned_pfus <- clean_pfu_data(pfus_starved %>% 
                                 mutate(interaction = case_when(condition %in% c("C", "D") ~ "no cells", 
                                                                TRUE ~ interaction))) %>%
  mutate(interaction = case_when(is.na(interaction) ~ "No Cells",
                                 TRUE ~ interaction))

all_data_phage <- rbind(no_cost, with_cost)
all_data_bacteria <- rbind(no_cost_bacteria, with_cost_bacteria)

#part A - starved cells or phage 24, 48 hour timepoints
partA <- cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                                        phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "starved\nE. coli\nmonoculture",
                                 interaction == "S Monoculture" ~ "starved\nS. enterica\nmonoculture",
                                 interaction == "No Cells" ~ "no\ncells",
                                 TRUE ~ interaction))%>%
  filter(interaction %in% c("starved\nE. coli\nmonoculture",
                            "starved\nS. enterica\nmonoculture",
                            "no\ncells")) %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(phage != "Phi + P22") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  ggplot(aes(x = interaction, y = doublings, color = phage)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  ylab("growth rate")+
  ylim(-10, 12.5)+
  labs(color = "phage type")

#part B - monoculture, competition, from original data and cost data
partB <- all_data_phage %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(!interaction %in% c("Facilitation", "Mutualism"))%>%
  mutate(interaction = case_when(interaction == "S Monoculture" ~ "S. enterica\nmonoculture",
                                 interaction == "Competition" ~ "competition")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)")) %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  filter(phage == "both\nphage") %>%
  ggplot(aes(x = cost, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("growth rate")+
  labs(color = "phage type")+
  ylim(-10, 12.5)

#legends
legend1 <- get_legend(cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                                                     phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
                        mutate(interaction = case_when(interaction == "E Monoculture" ~ "starved e monoculture",
                                                       interaction == "S Monoculture" ~ "starved s monoculture",
                                                       interaction == "No Cells" ~ "no cells",
                                                       TRUE ~ interaction))%>%
                        filter(interaction %in% c("starved e monoculture",
                                                  "starved s monoculture",
                                                  "no cells")) %>%
                        mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                                                     TRUE ~ doublings)) %>%
                        filter(phage != "Phi + P22") %>%
                        mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                                                 phage == "Phi" ~ "generalist\nonly",
                                                 phage == "Phi + P22" ~ "both\nphage")) %>%
                        mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
                        ggplot(aes(x = phage, y = doublings, color = phage_type)) +
                        facet_wrap(~interaction) +
                        geom_boxplot() +
                        geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
                        theme_bw(base_size = 18)+
                        theme(axis.title = element_text(), 
                              panel.background = element_rect(fill = "white"), 
                              plot.background = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "bottom",
                              legend.background = element_blank(),
                              strip.background = element_blank())+
                        scale_color_manual(values = c("#CA3542", "#27647B"))+
                        ylab("growth rate (ln(final pfu / initial pfu))")+
                        xlab("treatment")+
                        ylim(-10, 12.5)+
                        labs(color = "species"))


#partC
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
  mutate(media = case_when(media == "saline" ~ "0.85%\nsaline",
                           media == "H2O" ~ "water",
                           media == "MM-" ~ "metal\nfree",
                           media == "P-" ~ "phos\nfree",
                           media == "S-" ~ "sulf\nfree",
                           media == "P/S-" ~ "phos/sulf\nfree",
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
                           media == "salt + yeast" ~ "salt\nyeast",
                           media == "salt + tryptone" ~ "salt\ntryp",
                           media == "yeast + tryptone" ~ "yeast\ntryp",
                           media == "tryptone" ~ "tryp",
                           TRUE ~ media)) %>%
  mutate(media = factor(media, levels = c("yeast", "salt",
                                          "salt\nyeast", "LB", "salt\ntryp",
                                          "tryp", "yeast\ntryp"))) %>%
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
all_partC_data <- rbind(partB_data, partC_data, partD_data)

partC <- all_partC_data %>%
  mutate(media = reorder(media, phi_doublings)) %>%
  group_by(date, media) %>%
  arrange(desc(phi_doublings)) %>%
  ungroup() %>%
  ggplot(aes(x = media, y = phi_doublings)) +
  geom_boxplot(color = "#CA3542") +
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
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank()) +
  ylab("growth rate")+
  xlab("media")+
  labs(color = "")+
  ylim(-15, 12.5)

#full figure
top <- plot_grid(partA, partB, labels = c("A", "B"),  label_size = 32,
                 ncol = 1)
right <- plot_grid(partC, legend1, ncol = 1, labels = c("C"), label_size = 32, rel_heights = c(1,0.1))

deg_fig <- plot_grid(top, right, ncol = 1)

#create png
png(here::here("figures", "poster-figs", "top-fig3.png"), res = 200, width = 1200, height = 1300)
top
dev.off()

png(here::here("figures", "poster-figs", "right-fig3.png"), res = 200, width = 2000, height = 1300)
right
dev.off()
