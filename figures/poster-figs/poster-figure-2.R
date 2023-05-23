#ATB
#poster figure 
#experimental work

#library
library("patchwork")
library("cowplot")

#set date
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#get last of relevant data
date <- "9May2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
pfus_partA <- load_pfu_data(pfu_path) %>% janitor::clean_names() %>% mutate(well = c(1:28))
cleaned_pfus <- clean_pfu_data(pfus_partA)

#competition
competition_main_both <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-15, 15)+
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
  labs(color = "phage type")

competition_bacterial_OD <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition")%>%
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
  geom_smooth(span = 0.2, fill = "black")+
  scale_color_manual(values = c("#5ba300", "#e6308a"))+
  facet_wrap(~phage, ncol = 4)+
  theme_bw(base_size = 18)+
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
  labs(color = "species")

competition_bacterial_density <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  mutate(S_percentage = S_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  pivot_longer(E_percentage:S_percentage, names_to = "bacteria", values_to = "density")%>%
  group_by(bacteria, phage) %>%
  summarise(mean= mean(density), sd = sd(density)) %>%
  ggplot(aes(x = phage, y = mean, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  scale_fill_manual(values = c("#5ba300", "#e6308a"))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("percent of final population")+
  labs(color = "interaction")+
  ylim(0, 100)

bacterial_comp <- competition_bacterial_OD / competition_bacterial_density

competition_experimental <- competition_main_both + bacterial_comp

#mutualsim
mutualism_main_both <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-15, 15)+
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
  labs(color = "phage type")

mutualism_bacterial_OD <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism")%>%
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
  geom_smooth(span = 0.2, fill = "black")+
  scale_color_manual(values = c("#5ba300", "#e6308a"))+
  facet_wrap(~phage, ncol = 4)+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  ylim(0, 0.4)+
  xlab("hours")+
  labs(color = "species")

mutualism_bacterial_density <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  mutate(S_percentage = S_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  pivot_longer(E_percentage:S_percentage, names_to = "bacteria", values_to = "density")%>%
  group_by(bacteria, phage) %>%
  summarise(mean= mean(density), sd = sd(density)) %>%
  ggplot(aes(x = phage, y = mean, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  scale_fill_manual(values = c("#5ba300", "#e6308a"))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("percent of final population")+
  labs(color = "interaction")+
  ylim(0, 100)

bacterial_coop <- mutualism_bacterial_OD / mutualism_bacterial_density

mutualism_experimental <- mutualism_main_both + bacterial_coop

#full figure
legend_experimental <- get_legend(specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
                                    filter(time < 2000) %>%
                                    mutate(cost = gamma_sp / 20) %>%
                                    filter(cost == 5) %>%
                                    pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
                                    mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                                                               microbe == "S" ~ "S. enterica",
                                                               microbe == "gen" ~ "generalist (eh7)",
                                                               microbe == "sp" ~ "specialist (p22vir)")) %>%
                                    filter(microbe == "generalist (eh7)" | microbe == "specialist (p22vir)") %>%
                                    ggplot(aes(x = time / 10, y = biomass, color = microbe)) +
                                    geom_line(size = 2) +
                                    theme_bw(base_size = 18)+
                                    theme(axis.title = element_text(), 
                                          panel.background = element_rect(fill = "white"), 
                                          plot.background = element_blank(),
                                          legend.position = "bottom",
                                          panel.grid.minor = element_blank(),
                                          legend.background = element_blank(),
                                          strip.background = element_blank())+
                                    ylab("biomass")+
                                    xlab("time (a.u.)")+
                                    ylim(0, 250)+
                                    labs(color = "species")+
                                    scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                                                  "generalist (eh7)" = "#CA3542", "specialist (p22vir)" = "#27647B")))

experiment <- plot_grid(competition_experimental, mutualism_experimental, legend_experimental, ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1, 0.1),
                   label_size = 32)
#create png
png(here::here("figures", "poster-figs", "experiment.png"), res = 200, width = 2650, height = 2100)
experiment
dev.off()