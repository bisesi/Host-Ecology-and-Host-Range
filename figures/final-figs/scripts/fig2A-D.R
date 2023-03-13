#ATB
#Paper Fig 2
#Generation and visualization code for parts A, B, C

#load visualization packages
library("patchwork")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "fig2A-D-data-generation.R"))

#part A
partA <- coop_sp_to_gen_pip_values %>% 
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "invasion")) %>%
  ggplot(aes(x = I_S_sp, y = I_E_gen, fill = pip)) +
  geom_tile()+
  scale_fill_manual(values = c("grey", "black"))+
  xlab("resident infectivity on shared prey")+
  ylab("mutant infectivity on alternative prey")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")

#part B
partB <- comp_sp_to_gen_pip_values %>% 
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = I_S_sp, y = I_E_gen, fill = pip)) +
  geom_tile()+
  scale_fill_manual(values = c("grey", "black"))+
  xlab("resident infectivity on shared prey")+
  ylab("mutant infectivity on alternative prey")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")

#part C
partC <- coop_gen_to_sp_pip_values %>% 
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = gamma_genS, y = gamma_sp, fill = pip)) +
  geom_tile()+
  scale_fill_manual(values = c("grey", "black"))+
  xlab("resident burst size on shared prey")+
  ylab("mutant burst size on shared prey")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  guides(fill="none")+
  theme_bw()+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")

#part D
partD <- comp_gen_to_sp_pip_values %>% 
  ungroup() %>%
  filter(mu2 == mu2_default) %>%
  filter(beta1 %in% c(0.3, 1.2, 1.5, 2.1)) %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = gamma_genS, y = gamma_sp, fill = pip)) +
  geom_tile()+
  facet_wrap(~beta1, labeller = labeller(beta1 = 
                                           c("0.3" = "shared prey coeff: 0.3",
                                             "1.2" = "shared prey coeff: 1.2",
                                             "1.5" = "shared prey coeff: 1.5",
                                             "2.1" = "shared prey coeff: 2.1")))+
  scale_fill_manual(values = c("grey", "black"))+
  xlab("resident burst size on shared prey")+
  ylab("mutant burst size on shared prey")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")

#all parts fig 2
fig2 <- (partA + partB) / (partC + partD) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

#load model generation code
setwd(here::here("figures", "final-figs", "imgs"))
ggsave("fig2.png")

