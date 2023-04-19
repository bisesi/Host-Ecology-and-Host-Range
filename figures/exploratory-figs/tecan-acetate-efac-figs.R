#ATB
#E facilitation and acetate test
#22 March 2023

#set date of tecan run for analysis
date <- "22March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#adjusted OD graphs
E_OD <- all_tecan_adjusted_OD %>% 
  unite("condition", interaction:phage, sep = " ", remove = FALSE) %>%
  filter(interaction != "Emono") %>%
  ggplot(aes(x = cycle, y = E_corrected_OD, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition, ncol = 4) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

S_OD <- all_tecan_adjusted_OD %>% 
  unite("condition", interaction:phage, sep = " ", remove = FALSE) %>%
  filter(interaction != "Emono") %>%
  ggplot(aes(x = cycle, y = S_corrected_OD, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition, ncol = 4) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

#phage densities
pfus_and_final_density %>% 
  ggplot(aes(x = interaction, y = doublings, color = doubling_type)) + 
  geom_boxplot() + 
  facet_wrap(~phage_interaction) +
  theme_bw()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

pfus_and_final_density %>% 
  filter(phage_interaction == "Phage Competition") %>% 
  ggplot(aes(x = S_corrected_OD + E_corrected_OD, y = generalist_fitness, color = interaction)) + 
  geom_point() + 
  xlab("final OD")+
  ylab("generalist fitness")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()

allOD <- all_tecan_data %>% 
  unite("condition", interaction:phage, sep = " ", remove = FALSE) %>%
  filter(condition != "none none") %>%
  ggplot(aes(x = cycle, y = OD, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition, ncol = 4) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
