#ATB
#Tecan data figures
#Ac and succ figures

#set date of tecan run for analysis
date <- "10April2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#adjusted OD graphs
all_OD <- all_tecan_data %>%
  unite("condition", interaction:media, sep = " ", remove = FALSE) %>% 
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
  

E_OD <- all_tecan_adjusted_OD %>% 
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  unite("condition", interaction:media, sep = " ", remove = FALSE) %>%
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
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  unite("condition", interaction:media, sep = " ", remove = FALSE) %>%
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
  inner_join(., plate_layout %>% select(well, media), by = "well") %>%
  ggplot(aes(x = interaction, y = doublings, color = doubling_type)) + 
  geom_boxplot() + 
  facet_grid(media~phage_interaction) +
  theme_bw()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())