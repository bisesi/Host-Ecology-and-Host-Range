#DDF Figure

#library
library("patchwork")
library("cowplot")
library("ggtext")

ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"

#set date
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#get last of relevant data
date <- "9May2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
pfus_partA <- load_pfu_data(pfu_path) %>% janitor::clean_names() %>% mutate(well = c(1:28))
cleaned_pfus <- clean_pfu_data(pfus_partA)


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
  filter(phage == "both\nphage") %>%
  group_by(phage_type) %>%
  summarize(mean = mean(doublings), sd = sd(doublings)) %>% 
  ggplot(aes(x = phage_type, y = mean, color = phage_type)) +
  geom_point(size = 7.5) +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-7, 15)+
  ggtitle("competition")+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")

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
  filter(phage == "both\nphage") %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  group_by(phage_type) %>%
  summarize(mean = mean(doublings), sd = sd(doublings)) %>% 
  ggplot(aes(x = phage_type, y = mean, color = phage_type)) +
  geom_point(size = 7.5) +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-7, 15)+
  ggtitle("mutualism")+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "phage type")

experimental <- competition_main_both / mutualism_main_both

#model
#generate all data for figures
source(here::here("data-generation", "model", "final-figs", "bacterial-dynamics-data-generation.R"))

both_phage_main_comp <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "Generalist (EH7)",
                             microbe == "sp" ~ "Specialist (P22vir)")) %>%
  filter(microbe == "Generalist (EH7)" | microbe == "Specialist (P22vir)") %>%
  filter(time == max(time)) %>%
  mutate(biomass = ifelse(microbe == "Specialist (P22vir)", 0, biomass)) %>%
  ggplot(aes(x = microbe, y = biomass, color = microbe)) +
  geom_point(size = 7.5)+
  theme_bw(base_size = 18)+
  ggtitle("competition")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("final phage biomass")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "Generalist (EH7)" = eh7, "Specialist (P22vir)" = p22vir))

both_phage_main_coop <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "Generalist (EH7)",
                             microbe == "sp" ~ "Specialist (P22vir)")) %>%
  filter(microbe == "Generalist (EH7)" | microbe == "Specialist (P22vir)") %>%
  filter(time == max(time)) %>%
  ggplot(aes(x = microbe, y = biomass, color = microbe)) +
  geom_point(size = 7.5)+
  theme_bw(base_size = 18)+
  ggtitle("mutualism")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("final phage biomass")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "Generalist (EH7)" = eh7, "Specialist (P22vir)" = p22vir))

model <- both_phage_main_comp / both_phage_main_coop


legend <- get_legend(specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
                                    filter(time < 2000) %>%
                                    mutate(cost = gamma_sp / 20) %>%
                                    filter(cost == 5) %>%
                                    pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
                                    mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                                                               microbe == "S" ~ "S. enterica",
                                                               microbe == "gen" ~ "generalist (eh7)",
                                                               microbe == "sp" ~ "specialist (p22vir)")) %>%
                                    filter(microbe == "generalist (eh7)" | microbe == "specialist (p22vir)") %>%
                                    ggplot(aes(x = time / 10, y = biomass, fill = microbe)) +
                                    geom_bar(stat = "identity")+
                                    theme_bw(base_size = 22)+
                                    theme(axis.title = element_text(), 
                                          panel.background = element_rect(fill = "white"), 
                                          plot.background = element_blank(),
                                          legend.position = "bottom",
                                          legend.text = element_markdown(),
                                          panel.grid.minor = element_blank(),
                                          legend.background = element_blank(),
                                          strip.background = element_blank())+
                                    ylab("biomass")+
                                    xlab("time (a.u.)")+
                                    ylim(0, 250)+
                                    labs(fill = "species")+
                                    scale_fill_manual(values = c("*E. coli*" = ecoli, "*S. enterica*" = senterica, 
                                                                 "Generalist (EH7)" = eh7, "Specialist (P22*vir*)" = p22vir)))



ddf_figure <- plot_grid(model, experimental, ncol = 2, rel_heights = c(1, 1),
                        label_size = 30)
#create png
png(here::here("figures", "poster-figs", "ddf-figure.png"), res = 200, width = 1500, height = 1500)
ddf_figure
dev.off()