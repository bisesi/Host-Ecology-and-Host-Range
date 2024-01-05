#ATB
#poster figure 3
#all experimental results

#library
library("patchwork")
library("cowplot")
library("ggtext")

ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"

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
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "starved<br>*E. coli*<br>monoculture",
                                 interaction == "S Monoculture" ~ "starved<br>*S. enterica*<br>monoculture",
                                 interaction == "No Cells" ~ "no cells",
                                 TRUE ~ interaction))%>%
  filter(interaction %in% c("starved<br>*E. coli*<br>monoculture",
                            "starved<br>*S. enterica*<br>monoculture",
                            "no cells")) %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(phage != "Phi + P22") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  ggplot(aes(x = interaction, y = doublings, color = phage)) +
  geom_boxplot(alpha = 0.5, varwidth = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_markdown(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  scale_color_manual(values = c(p22vir, eh7))+
  ylab("ln(final pfu / initial pfu)")+
  ylim(-5, 12.5)+
  labs(color = "phage type")

#part B - monoculture, competition, from original data and cost data
partB <- all_data_phage %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(!interaction %in% c("Facilitation", "Mutualism"))%>%
  mutate(interaction = case_when(interaction == "S Monoculture" ~ "*S. enterica*<br>monoculture",
                                 interaction == "Competition" ~ "competition")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)")) %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage")) %>%
  mutate(phage= factor(phage, levels = c("specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  filter(phage == "both\nphage") %>%
  mutate(cost = case_when(cost == "no cost" ~ "no\ncost",
                          cost == "with cost" ~ "with\ncost")) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  ggplot(aes(x = cost, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_markdown(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "phage type")+
  ylim(-5, 12.5)

#legends
legend1 <- get_legend(cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "Generalist (EH7)",
                                                                     phage_type == "Specialist phage" ~ "Specialist (P22*vir*)"))%>%
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
                              legend.text = element_markdown(),
                              legend.background = element_blank(),
                              strip.background = element_blank())+
                        scale_color_manual(values = c(eh7, p22vir))+
                        ylab("growth rate (ln(final pfu / initial pfu))")+
                        xlab("treatment")+
                        ylim(-10, 12.5)+
                        labs(color = "species"))


#full figure
top <- plot_grid(partA, partB, labels = c("A", "B"),  label_size = 28,
                 ncol = 2, rel_widths = c(0.75, 1))

full <- plot_grid(top, legend1, ncol = 1, rel_heights = c(1, 0.1))

#create png
png(here::here("figures", "poster-figs", "fig3.png"), res = 200, width = 2200, height = 900)
full
dev.off()

