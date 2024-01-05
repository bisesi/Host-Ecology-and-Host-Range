#ATB
#poster figure 
#experimental work

#library
library("patchwork")
library("cowplot")
library("ggtext")
library("tidyverse")
library("readxl")

ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"

#set date
date <- "22June2023"
#source tecan data cleaning script
source(here::here("functions", "tecan-data-helper-functions.R"))

#paths
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#import tecan data, plate layout and pfus
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

#clean pfu data
cleaned_pfus <- clean_pfu_data(pfus) %>% inner_join(., plate_layout %>% select(well, timepoint), by = "well")

#plot June 21 2023
plot <- cleaned_pfus %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(interaction = case_when(interaction == "S Monoculture" ~ "s monoculture",
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
  facet_grid(interaction~timepoint)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")


#over time
date <- "21June2023"

#paths
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))

#import tecan data, plate layout and pfus
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", date, "plate_layout.csv"))
pfus <- load_pfu_data(pfu_path) %>% janitor::clean_names()

wide_data <- pfus %>% 
  select(condition, interaction, phage, plate, pfu, well) %>%
  pivot_wider(names_from = plate, values_from = pfu) %>%
  rename(pfu_E = E)

output <- wide_data %>%
  filter(condition != "Start") %>%
  mutate(phi_doublings = case_when(phage == "Phi" ~ log(pfu_E / wide_data %>% filter(condition == "Start") %>% select(pfu_E) %>% pull()),
                                   TRUE ~ 0)) %>%  
  mutate_all(~replace(., is.infinite(.), log(1 / wide_data %>% filter(condition == "Start") %>% select(pfu_E) %>% pull()))) %>%
  drop_na() %>%
  mutate(interaction = case_when(interaction == "No Cells" ~ "no cells")) %>%
  pivot_longer(cols = ends_with("doublings"),
               names_to = "doubling_type",
               values_to = "doublings") %>%
  mutate(phage_type = case_when(doubling_type == "p22_doublings" ~ "Specialist phage",
                                doubling_type == "phi_doublings" ~ "Generalist phage",
                                TRUE ~ "NA")) %>%
  mutate(phage_interaction = case_when(phage == "P22" ~ "Specialist phage",
                                       phage == "Phi" ~ "Generalist phage",
                                       phage == "Phi + P22" ~ "Phage Competition")) %>%
  inner_join(., plate_layout %>% select(timepoint, media, well), by = "well")

#plot 21june2023
plot_time <- output %>% filter(media == "glucose") %>%
  group_by(timepoint) %>%
  summarize(mean = mean(pfu_E), sd = sd(pfu_E)) %>%
  mutate(timepoint= factor(timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "24", "48"))) %>%
  ggplot(aes(x = timepoint, y = mean)) +
  geom_point(size = 2)+
  #geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("pfu/mL")+
  labs(color = "species")

plot_LB <- output %>% filter(media == "LB") %>%
  ggplot(aes(x = media, y = doublings)) +
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-5, 5)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")



