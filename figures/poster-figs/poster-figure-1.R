#ATB
#poster figure
#modeling side

#load packages
library("patchwork")
library("cowplot")
library("ggtext")

#set color palette
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"

#generate all data for figures
source(here::here("data-generation", "model", "final-figs", "bacterial-dynamics-data-generation.R"))

#get dataset
single_phage <- specialist_gamma_comp_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "competition") %>%
  rbind(., specialist_gamma_comp_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_coop_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_comp_none %>% mutate(phage = "no phage") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_none %>% mutate(phage = "no phage") %>% mutate(interaction = "mutualism")) %>%
  filter((gamma_sp / 20) == 5) %>%
  filter(time < 2000)

#competition row
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
  ggplot(aes(x = microbe, y = biomass - 25, color = microbe)) +
  geom_point(size = 7.5)+
  theme_bw(base_size = 18)+
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

bacterial_dynamics_comp <- single_phage %>%
  rbind(., specialist_gamma_comp_both %>% mutate(phage = "both phage") %>% mutate(interaction = "competition")) %>%
  filter(time < 2000 & interaction == "competition") %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only", "both phage"))) %>%
  mutate(bacteria = case_when(phage == "no phage" & type == "E. coli" ~ bacteria + 0.05,
                              phage == "generalist only" & type == "E. coli" ~ bacteria + 0.05,
                              TRUE ~ bacteria)) %>%
  ggplot(aes(x = time / 100, y = bacteria, color = type)) +
  geom_line(size = 2) +
  facet_wrap(~phage, ncol = 4) +
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("bacterial\nbiomass")+
  xlab("time (a.u.)")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c(ecoli, senterica))

final_bacterial_dynamics_comp <- single_phage %>%
  rbind(., specialist_gamma_comp_both %>% mutate(phage = "both phage") %>% mutate(interaction = "competition")) %>%
  filter(time < 2000 & interaction == "competition") %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  filter(time == max(time)) %>%
  mutate(phage = case_when(phage == "both phage" ~ "both\nphage",
                           phage == "no phage" ~ "no\nphage",
                           phage == "generalist only" ~ "generalist\nonly",
                           phage == "specialist only" ~ "specialist\nonly",
                           TRUE ~ phage)) %>%
  mutate(E_percent = E / (E + S), S_percent = S / (E + S)) %>%
  pivot_longer(E_percent:S_percent, names_to = "bacteria", values_to = "density") %>%
  mutate(bacteria = case_when(bacteria == "E_percent" ~ "E. coli", 
                              bacteria == "S_percent" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  ggplot(aes(x = phage, y = density * 100, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  scale_fill_manual(values = c(ecoli, senterica))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("% of final\nco-culture")+
  labs(color = "interaction")+
  ylim(0, 100)

bacteria_comp <- bacterial_dynamics_comp / final_bacterial_dynamics_comp

competition <- plot_grid(both_phage_main_comp, bacteria_comp, rel_widths = c(0.5, 1.25))

#mutualism row
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

bacterial_dynamics_coop <- single_phage %>%
  rbind(., specialist_gamma_coop_both %>% mutate(phage = "both phage") %>% mutate(interaction = "mutualism")) %>%
  filter(time < 2000 & interaction == "mutualism") %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only", "both phage"))) %>%
  mutate(bacteria = case_when(phage == "no phage" & type == "E. coli" ~ bacteria + 0.05,
                              phage == "generalist only" & type == "E. coli" ~ bacteria + 0.05,
                              TRUE ~ bacteria)) %>%
  ggplot(aes(x = time / 100, y = bacteria, color = type)) +
  geom_line(size = 2) +
  facet_wrap(~phage, ncol = 4) +
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("bacterial\nbiomass")+
  xlab("time (a.u.)")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c(ecoli, senterica))

final_bacterial_dynamics_coop <- single_phage %>%
  rbind(., specialist_gamma_coop_both %>% mutate(phage = "both phage") %>% mutate(interaction = "mutualism")) %>%
  filter(time < 2000 & interaction == "mutualism") %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  filter(time == max(time)) %>%
  mutate(phage = case_when(phage == "both phage" ~ "both\nphage",
                           phage == "no phage" ~ "no\nphage",
                           phage == "generalist only" ~ "generalist\nonly",
                           phage == "specialist only" ~ "specialist\nonly",
                           TRUE ~ phage)) %>%
  mutate(E_percent = E / (E + S), S_percent = S / (E + S)) %>%
  pivot_longer(E_percent:S_percent, names_to = "bacteria", values_to = "density") %>%
  mutate(bacteria = case_when(bacteria == "E_percent" ~ "E. coli", 
                              bacteria == "S_percent" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly", "both\nphage"))) %>%
  ggplot(aes(x = phage, y = density * 100, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  scale_fill_manual(values = c(ecoli, senterica))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("% of final\nco-culture")+
  labs(color = "interaction")+
  ylim(0, 100)

bacteria_coop <- bacterial_dynamics_coop / final_bacterial_dynamics_coop

mutualism <- plot_grid(both_phage_main_coop, bacteria_coop, rel_widths = c(0.5, 1.25))

#final figure
legend <- get_legend(specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
                       filter(time < 2000) %>%
                       mutate(cost = gamma_sp / 20) %>%
                       filter(cost == 5) %>%
                       pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
                       mutate(microbe = case_when(microbe == "E" ~ "*E. coli*", 
                                                  microbe == "S" ~ "*S. enterica*",
                                                  microbe == "gen" ~ "Generalist (EH7)",
                                                  microbe == "sp" ~ "Specialist (P22vir)")) %>%
                       filter(microbe == "Generalist (EH7)" | microbe == "Specialist (P22vir)") %>%
                       ggplot(aes(x = time / 100, y = biomass, fill = microbe)) +
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

model <- plot_grid(competition, mutualism, legend, ncol = 1, labels = c("competition", "mutualism"), rel_heights = c(1, 1, 0.1),
          label_size = 30)

#create png
png(here::here("figures", "poster-figs", "model.png"), res = 200, width = 2650, height = 2100)
model
dev.off()
