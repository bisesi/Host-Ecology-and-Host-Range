#ATB
#Generate bacterial dynamics in a Lotka Volterra model figure
#Figure 2
#Uses bacterial dynamics data generation script

#load packages and data
library("patchwork")
library("cowplot")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "bacterial-dynamics-data-generation.R"))

#partA - application of individual phage
partA <- specialist_gamma_comp_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "competition") %>%
  rbind(., specialist_gamma_comp_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_coop_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_comp_none %>% mutate(phage = "no phage") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_none %>% mutate(phage = "no phage") %>% mutate(interaction = "mutualism")) %>%
  filter((gamma_sp / 20) == 5) %>%
  filter(time < 2000) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only"))) %>%
  mutate(bacteria = case_when(phage == "no phage" & type == "E. coli" ~ bacteria + 0.05,
                              phage == "generalist only" & type == "E. coli" ~ bacteria + 0.05,
                              TRUE ~ bacteria)) %>%
  ggplot(aes(x = time / 100, y = bacteria, color = type)) +
  geom_line(size = 1.5) +
  facet_grid(interaction~phage) +
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c(ecoli, senterica))

#partB - high cost in competition
partB_main <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
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
  geom_line(size = 1.5) +
  ggtitle("competition")+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "generalist (eh7)" = eh7, "specialist (p22vir)" = p22vir))

partB_inset <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "generalist (eh7)",
                             microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "E. coli" | microbe == "S. enterica") %>%
  ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
  geom_line(size = 1.5) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 2)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "generalist (eh7)" = eh7, "specialist (p22vir)" = p22vir))
partB <- ggdraw(partB_main) + 
  draw_plot(partB_inset, x = 0.17, y = 0.5, width = 0.35, height = 0.35)

#partC - cooperation at high cost
partC_main <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "generalist (eh7)",
                             microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "generalist (eh7)" | microbe == "specialist (p22vir)") %>%
  ggplot(aes(x = time / 100, y = biomass, color = microbe)) +
  geom_line(size = 1.5) +
  theme_bw(base_size = 18)+
  ggtitle("mutualism")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "generalist (eh7)" = eh7, "specialist (p22vir)" = p22vir))

partC_inset <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "Generalist (EH7)",
                             microbe == "sp" ~ "Specialist (P22vir)")) %>%
  filter(microbe == "E. coli" | microbe == "S. enterica") %>%
  ggplot(aes(x = time / 100, y = biomass, color = microbe)) +
  geom_line(size = 1.5) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 2)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = ecoli, "S. enterica" = senterica, 
                                "generalist (eh7)" = eh7, "specialist (p22vir)" = p22vir))
partC <- ggdraw(partC_main) + 
  draw_plot(partC_inset, x = 0.17, y = 0.5, width = 0.35, height = 0.35)

legend <- get_legend(specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
                       filter(time < 2000) %>%
                       mutate(cost = gamma_sp / 20) %>%
                       filter(cost == 5) %>%
                       pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
                       mutate(microbe = case_when(microbe == "E" ~ "*E. coli*", 
                                                  microbe == "S" ~ "*S. enterica*",
                                                  microbe == "gen" ~ "Generalist (EH7)",
                                                  microbe == "sp" ~ "Specialist (P22vir)")) %>%
                       filter(microbe == "E. coli" | microbe == "S. enterica") %>%
                       ggplot(aes(x = time / 100, y = biomass, color = microbe)) +
                       geom_line(size = 1.5) +
                       theme_bw(base_size = 18)+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom",
                             legend.text = element_markdown(),
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("biomass")+
                       xlab("time (a.u.)")+
                       ylim(0, 2)+
                       labs(color = "species")+
                       scale_color_manual(values = c("*E. coli*" = ecoli, "*S. enterica*" = senterica, 
                                                     "Generalist (EH7)" = eh7, "Specialist (P22vir)" = p22vir)))

#full plot
fig2 <- plot_grid(partA, 
                  plot_grid(partB, partC, labels = c("B", "C"), label_size = 26), 
                  ncol = 1,
                  labels = c("A"), label_size = 26,
                  legend, rel_heights = c(1, 1, .1))
