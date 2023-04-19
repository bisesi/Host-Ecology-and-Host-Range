#ATB
#Tecan data figures
#Fig 3 - A, B, C, D

#library
library("patchwork")
library("ggpubr")
library("rstatix")

#set date
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#part A
partA <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition")%>%
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
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("growth rate (log2)")+
  xlab("treatment condition")+
  labs(color = "phage type")

phage_dat <- pfus_and_final_density %>%
  filter(interaction == "Mutualism" | interaction == "Competition")

phage_model <- aov(doublings ~ interaction * phage * doubling_type, data = phage_dat)

phage_multiple_comparisons <- TukeyHSD(phage_model, conf.levels = 0.95)

#part B
partB <- pfus_and_final_density %>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "Mutualism",
                                interaction == "Comp" ~ "Competition",
                                TRUE ~ interaction)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
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
  ylab("e coli fraction of final population")+
  xlab("treatment condition")+
  labs(color = "interaction")+
  ylim(0, 100)

E_percent_dat <- pfus_and_final_density %>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "Mutualism",
                                 interaction == "Comp" ~ "Competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100)

E_percent_model <- aov(E_percentage ~ interaction * phage, data = E_percent_dat)

E_percent_multiple_comparisons <- TukeyHSD(E_percent_model, conf.levels = 0.95)

#part C 
partC <- all_tecan_adjusted_OD %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "Mutualism",
                                 interaction == "Comp" ~ "Competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition")%>%
  pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
  mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                             fluor == "S_corrected_OD" ~ "S. enterica")) %>%
  ggplot(aes(x = cycle, y = OD, color = species))+
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
  xlab("cycle")+
  labs(color = "species")

#fig 3
fig3 <- (partA + partB) / (partC) + plot_annotation(tag_levels = "A")



