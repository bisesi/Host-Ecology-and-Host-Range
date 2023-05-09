#ATB
#Tecan data figures
#Fig 5 - A, B, C, D

#set date for glucose
date <- "8March2023"
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))
glucose <- pfus_and_final_density %>% mutate(media = "glucose")

#set date for acetate
date <- "22March2023"
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))
acetate <- pfus_and_final_density %>% mutate(media = "acetate")

#part A
partA <- glucose %>%
  rbind(., acetate) %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "S Monoculture")%>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~media) +
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
  labs(color = "")+
  scale_color_manual(values = c("#CA3542", "#27647B"))

media_dat <- glucose %>%
  rbind(., acetate)

smono <- media_dat %>% filter(media == "glucose" & interaction == "S Monoculture")

smono_model <- aov(doublings ~ doubling_type * phage, data = smono)

smono_comparisons <- TukeyHSD(smono_model, conf.levels = 0.95)


#part B
feb22OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "22Feb2023", "od_22Feb2023.csv"), file_type = "OD") 
mar8OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "8March2023", "od_8March2023.csv"), file_type = "OD")
mar15OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "15March2023", "initial-tecan", "od_15March2023_1.csv"), file_type = "OD")

feb22 <- load_pfu_data(here::here("experimental-data", "tecan-data", "22Feb2023", "pfu_22Feb2023.xlsx")) %>% 
  janitor::clean_names() %>% clean_pfu_data() %>% filter(phage == "Phi + P22" & interaction == "S Monoculture") %>%
  mutate(date = "feb22") %>% inner_join(., feb22OD %>% group_by(well) %>% slice_max(cycle), by = "well")

mar8 <- load_pfu_data(here::here("experimental-data", "tecan-data", "8March2023", "pfu_8March2023.xlsx")) %>% 
  janitor::clean_names() %>% clean_pfu_data() %>% filter(phage == "Phi + P22" & interaction == "S Monoculture") %>%
  mutate(date = "mar8") %>% inner_join(., mar8OD %>% group_by(well) %>% slice_max(cycle), by = "well")

mar15 <- load_pfu_data(here::here("experimental-data", "tecan-data", "15March2023", "initial-tecan", "pfu_15March2023.xlsx")) %>% 
  janitor::clean_names() %>% clean_pfu_data() %>% filter(phage == "Phi + P22" & interaction == "S Monoculture" & well %in% c("B7", "B8", "B9", "B10", "B11", "C2", "C3", "C4")) %>%
  mutate(date = "mar15") %>% inner_join(., mar15OD %>% group_by(well) %>% slice_max(cycle), by = "well")

partB <- feb22 %>%
  rbind(., mar8) %>%
  rbind(., mar15) %>%
  filter(doubling_type == "phi_doublings") %>%
  filter(well != "B9") %>%
  ggplot(aes(x = OD, y = generalist_fitness)) +
  geom_point() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("generalist fitness")+
  xlab("final OD")+
  labs(color = "")

od_corr <- feb22 %>%
  rbind(., mar8) %>%
  rbind(., mar15) %>%
  filter(doubling_type == "phi_doublings") %>%
  filter(well != "B9") 

od_model <- lm(generalist_fitness ~ OD, data = od_corr)

#partC
library("tidyverse")
library("readxl")
date <- "15March2023"
source(here::here("functions", "baranyi-helper-functions.R"))
source(here::here("functions", "tecan-data-helper-functions.R"))

resistance_data <- load_pfu_data(here::here("experimental-data", "tecan-data", "15March2023", "resistant-isolates-with-ancestral","resistance_15March2023_3.xlsx"))

high_density <- resistance_data %>% filter(Interaction %in% c("R2", "R4", "R8") & Phage == "None" & Condition != "Start") %>%
  mutate(well_type = "high density")

ancestral <- resistance_data %>% filter(Interaction == "S0" & Phage %in% c("None", "Phi + P22") & Condition != "Start" & Phi_Resistant != "NA") %>%
  mutate(well_type = case_when(Phage == "None" ~ "ancestral",
                               Phage == "Phi + P22" ~ "low density", 
                               TRUE ~ Phage))

partC <- high_density %>% 
  rbind(., ancestral) %>%
  mutate(interaction = "S Monoculture") %>%
  pivot_longer(cols = ends_with("Resistant"), names_to = "phage_type", values_to = "percentage_resistant") %>%
  mutate(phage_type = case_when(phage_type == "P22_Resistant" ~ "Specialist phage",
                                phage_type == "Phi_Resistant" ~ "Generalist phage")) %>%
  ggplot(aes(x = well_type, y = percentage_resistant, color = phage_type)) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("percent resistance by cross-streak")+
  xlab("condition")+
  labs(color = "")+
  scale_color_manual(values = c("#CA3542", "#27647B"))


#fig 5
fig5 <- partA + (partB / partC) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")



