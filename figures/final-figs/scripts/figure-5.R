#ATB
#starved cells assay, delayed cell addition in monoculture + comp vs same time addition
#Figure 5

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
partA <- cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
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
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("treatment")+
  ylim(-10, 12.5)+
  labs(color = "phage type")

#part B - monoculture, competition, from original data and cost data
partB <- all_data_phage %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(!interaction %in% c("Facilitation", "Mutualism"))%>%
  mutate(interaction = case_when(interaction == "S Monoculture" ~ "s monoculture",
                                 interaction == "Competition" ~ "competition")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)")) %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_grid(cost~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw()+
  scale_color_manual(values = c("#CA3542", "#27647B"))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("growth rate (ln(final pfu / initial pfu))")+
  xlab("treatment")+
  labs(color = "phage type")+
  ylim(-10, 12.5)

#part C bacterial dynamics without cost
partC <- all_data_bacteria %>%
  filter(cost == "no cost") %>%
  mutate(interaction = case_when(interaction == "Fac" ~ "facilitation",
                                 interaction == "Comp" ~ "competition",
                                 interaction == "Smono" ~ "s monoculture",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition" | interaction == "s monoculture")%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist only",
                           phage == "Phi" ~ "generalist only",
                           phage == "Phi + P22" ~ "both phage",
                           phage == "none" ~ "no phage")) %>%
  mutate(phage= factor(phage, levels = c("no phage", "specialist only", "generalist only", "both phage"))) %>%
  pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
  mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                             fluor == "S_corrected_OD" ~ "S. enterica")) %>%
  ggplot(aes(x = hours, y = OD, color = species))+
  geom_smooth(span = 0.2, se = FALSE)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))+
  facet_grid(interaction~phage)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  xlab("hours")+
  ylim(0, 0.4)+
  labs(color = "species")+
  xlim(0, 48)

#part D bacterial dynamics with cost
dummy <- all_data_bacteria %>% filter(cost == "with cost") %>% filter(hours == min(hours)) %>%
  mutate(hours = case_when(hours != 0 ~ -24,
                           TRUE ~ hours))

partD <- all_data_bacteria %>%
  filter(cost == "with cost") %>%
  rbind(., dummy) %>%
  mutate(interaction = case_when(interaction == "Fac" ~ "facilitation",
                                 interaction == "Comp" ~ "competition",
                                 interaction == "Smono" ~ "s monoculture",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition" | interaction == "s monoculture")%>%
  pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
  mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                             fluor == "S_corrected_OD" ~ "S. enterica")) %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist only",
                           phage == "Phi" ~ "generalist only",
                           phage == "Phi + P22" ~ "both phage",
                           phage == "none" ~ "no phage")) %>%
  mutate(phage= factor(phage, levels = c("no phage", "specialist only", "generalist only", "both phage"))) %>%
  ggplot(aes(x = hours + 24, y = OD, color = species))+
  geom_smooth(span = 0.2, se = FALSE)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))+
  facet_grid(interaction~phage)+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  xlab("hours")+
  labs(color = "species")+
  xlim(0, 48)+
  ylim(0, 0.4)+
  geom_vline(xintercept = 24, color = "red", linetype = "dashed")

#legends
legend1 <- get_legend(cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
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
                        theme_bw()+
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
                        labs(color = "phage type"))

legend2 <- get_legend(all_data_bacteria %>%
                        filter(cost == "with cost") %>%
                        rbind(., dummy) %>%
                        mutate(interaction = case_when(interaction == "Fac" ~ "facilitation",
                                                       interaction == "Comp" ~ "competition",
                                                       interaction == "Smono" ~ "s monoculture",
                                                       TRUE ~ interaction)) %>%
                        filter(interaction == "competition" | interaction == "s monoculture")%>%
                        pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
                        mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                                                   fluor == "S_corrected_OD" ~ "S. enterica")) %>%
                        mutate(phage = case_when(phage == "P22" ~ "specialist only",
                                                 phage == "Phi" ~ "generalist only",
                                                 phage == "Phi + P22" ~ "both phage",
                                                 phage == "none" ~ "no phage")) %>%
                        mutate(phage= factor(phage, levels = c("no phage", "specialist only", "generalist only", "both phage"))) %>%
                        ggplot(aes(x = hours + 24, y = OD, color = species))+
                        geom_smooth(span = 0.2, se = FALSE)+
                        scale_color_manual(values = c("#5ba300", "#e6308a"))+
                        facet_grid(interaction~phage)+
                        theme_bw()+
                        theme(axis.title = element_text(), 
                              panel.background = element_rect(fill = "white"), 
                              plot.background = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "bottom",
                              legend.background = element_blank(),
                              strip.background = element_blank())+
                        ylab("OD")+
                        xlab("hours")+
                        labs(color = "species")+
                        xlim(0, 48)+
                        ylim(0, 0.4)+
                        geom_vline(xintercept = 24, color = "red", linetype = "dashed"))

#complete fig
fig4 <- plot_grid(plot_grid(partA, partB, legend1, ncol = 1, rel_heights = c(1,1,0.1), labels = c("A", "B")),
                  plot_grid(partC, partD, legend2, ncol = 1, rel_heights = c(1,1,0.1), labels = c("C", "D")))

