#ATB
#Tecan data figures, competition and cooperation
#Figure 4

#library
library("patchwork")
library("cowplot")
library("ggpubr")
library("rstatix")

#set date
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#part A - growth of each phage on respective hosts
date <- "9May2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
pfus_partA <- load_pfu_data(pfu_path) %>% janitor::clean_names() %>% mutate(well = c(1:28))
cleaned_pfus <- clean_pfu_data(pfus_partA)

partA <- cleaned_pfus %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "e monoculture",
                                 interaction == "S Monoculture" ~ "s monoculture")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  filter(phage != "both\nphage") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
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
  labs(color = "phage type")


#part B - coop and comp phage densities growth rates
partB <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
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
  labs(color = "phage type")

legend <- get_legend(pfus_and_final_density %>%
                       mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                                                    TRUE ~ doublings)) %>%
                       filter(interaction == "Mutualism" | interaction == "Competition") %>%
                       mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                                      interaction == "Competition" ~ "competition"))%>%
                       mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                                     phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
                       mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                                                phage == "Phi" ~ "generalist\nonly",
                                                phage == "Phi + P22" ~ "both\nphage",
                                                phage == "none" ~ "no\nphage")) %>%
                       mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                                               "both\nphage"))) %>%
                       ggplot(aes(x = phage, y = doublings, color = phage_type)) +
                       facet_wrap(~interaction) +
                       geom_boxplot() +
                       geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
                       theme_bw()+
                       scale_color_manual(values = c("#CA3542", "#27647B"))+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom",
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("growth rate (ln(final pfu / initial pfu))")+
                       xlab("treatment")+
                       labs(color = "phage type"))

phage_dat <- pfus_and_final_density %>%
  filter(interaction == "mutualism" | interaction == "competition")

phage_model <- aov(doublings ~ interaction * phage * doubling_type, data = phage_dat)

phage_multiple_comparisons <- TukeyHSD(phage_model, conf.levels = 0.95)

#part C - coop and comp bacterial densities at final timepoint
partC <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (phi-C)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                interaction == "Comp" ~ "competition",
                                TRUE ~ interaction)) %>%
  filter(interaction == "mutualism" | interaction == "competition") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  ggplot(aes(x = phage, y = E_percentage, color = interaction))+
  geom_boxplot()+
  scale_color_manual(values = c("black", "grey60"))+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("E. coli fraction of final population")+
  xlab("treatment")+
  labs(color = "interaction")+
  ylim(0, 100)

E_percent_dat <- pfus_and_final_density %>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism" | interaction == "competition") %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100)

E_percent_model <- aov(E_percentage ~ interaction * phage, data = E_percent_dat)

E_percent_multiple_comparisons <- TukeyHSD(E_percent_model, conf.levels = 0.95)

#part D - OD600 of species coop or comp
partD <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism" | interaction == "competition")%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist only",
                           phage == "Phi" ~ "generalist only",
                           phage == "Phi + P22" ~ "both phage",
                           phage == "none" ~ "no phage")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                          "both phage"))) %>%
  filter(!well %in% c("B8", "B9")) %>%
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
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  xlab("hours")+
  labs(color = "species")+
  geom_vline(data = (data.frame(xint=20, interaction="competition", phage = "no phage") %>% mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                                                                                   "both phage")))), 
             aes(xintercept = xint), color = "red", linetype = "dashed")+
  geom_vline(data = (data.frame(xint=47, interaction="mutualism", phage = "no phage") %>% mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                                                                                    "both phage")))), 
             aes(xintercept = xint), color = "red", linetype = "dashed")
#fig 4
fig4 <- plot_grid(plot_grid(partA, partB, legend, labels = c("A", "B"), rel_heights = c(1,1, 0.15), ncol = 1),
                  plot_grid(partC, partD, labels = c("C", "D"), ncol = 1))



