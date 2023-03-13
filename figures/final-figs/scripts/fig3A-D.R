#ATB
#Tecan data figures
#Fig 3 - A, B, C, D

#set date
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#part A
partA <- pfus_and_final_density %>%
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
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("e coli fraction of final population")+
  xlab("treatment condition")+
  labs(color = "bacterial interaction")+
  ylim(0, 100)


#part B
no_phage_controls <- generate_baranyi_growth_data_from_OD(OD_data) %>% ungroup() %>% 
  filter(fit_variable == "growth_rate") %>% inner_join(., plate_layout, by = "well") %>%
  filter(phage == "none" & interaction != "none")

partB <- pfus_and_final_density %>%
  inner_join(., generate_baranyi_growth_data_from_OD(OD_data) %>% ungroup() %>% 
               filter(fit_variable == "growth_rate"),
             by = "well") %>%
  select(fit_variable, interaction, phage, model_fit, well) %>%
  rbind(., no_phage_controls) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "Mutualism",
                                 interaction == "Comp" ~ "Competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>% 
  mutate(model_fit = case_when(is.na(model_fit) ~ 0,
                               TRUE ~ model_fit)) %>%
  ggplot(aes(x = phage, y = model_fit, color = interaction))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("baranyi growth rate")+
  xlab("treatment condition")+
  labs(color = "bacterial interaction")

#part C 
no_phage_controls_lag <- generate_baranyi_growth_data_from_OD(OD_data) %>% ungroup() %>% 
  filter(fit_variable == "lag") %>% inner_join(., plate_layout, by = "well") %>%
  filter(phage == "none" & interaction != "none")

partC <- pfus_and_final_density %>%
  inner_join(., generate_baranyi_growth_data_from_OD(OD_data) %>% ungroup() %>% 
               filter(fit_variable == "lag"),
             by = "well") %>%
  select(fit_variable, interaction, phage, model_fit, well) %>%
  rbind(., no_phage_controls_lag) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "Mutualism",
                                 interaction == "Comp" ~ "Competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition")%>%
  mutate(model_fit = case_when(is.na(model_fit) ~ 0,
                              TRUE ~ model_fit)) %>% 
  ggplot(aes(x = phage, y = model_fit, color = interaction))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("time to exponential phase")+
  xlab("treatment condition")+
  labs(color = "bacterial interaction")

#part D
partD <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition")%>%
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
  ylab("number of doublings")+
  xlab("treatment condition")+
  labs(color = "")

#fig 3
(partA + partB) / (partC + partD) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")



