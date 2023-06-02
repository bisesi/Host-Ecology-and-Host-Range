#ATB
#ESS modeling figure
#Supplemental figure 2

#load visualization packages
library("patchwork")
library("cowplot")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "ess-data-generation.R"))

#part A
partA <- coop_sp_to_gen_pip_values %>% 
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "invasion")) %>%
  ggplot(aes(x = I_S_sp, y = I_E_gen)) +
  geom_tile(aes(fill = pip))+
  scale_fill_manual(values = c("black", "grey"))+
  xlab("resident infectivity on *S. enterica*")+
  ylab("mutant infectivity on *S. enterica*")+
  theme_bw(base_size = 18)+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank(),
        strip.background = element_blank())

#part B
partB <- comp_sp_to_gen_pip_values %>% 
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = I_S_sp, y = I_E_gen)) +
  geom_tile(aes(fill = pip))+
  xlab("resident infectivity on S. enterica")+
  ylab("mutant infectivity on S. enterica")+
  scale_fill_manual(values = c("black", "grey"))+
  theme_bw(base_size = 18)+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        legend.background = element_blank(),
        strip.background = element_blank())
  

#part C
partC <- coop_gen_to_sp_pip_values %>% 
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = gamma_genS, y = gamma_sp)) +
  geom_tile(aes(fill = pip))+
  scale_fill_manual(values = c("grey"))+
  xlab("resident burst size on *S. enterica*")+
  ylab("mutant burst size on *S. enterica*")+
  theme_bw(base_size = 18)+
  labs(fill = " ")+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank(),
        strip.background = element_blank())
  

#part D
partD <- comp_gen_to_sp_pip_values %>% 
  ungroup() %>%
  filter(mu2 == mu2_default) %>%
  filter(beta1 %in% c(0.3, 1.2, 1.5, 2.1)) %>%
  mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
  ggplot(aes(x = gamma_genS, y = gamma_sp)) +
  geom_tile(aes(fill = pip))+
  facet_wrap(~beta1, labeller = labeller(beta1 = 
                                           c("0.3" = "*S. enterica* competition coefficient: 0.3",
                                             "1.2" = "*S. enterica* competition coefficient: 1.2",
                                             "1.5" = "*S. enterica* competition coefficient: 1.5",
                                             "2.1" = "*S. enterica* competition coefficient: 2.1")))+
  scale_fill_manual(values = c("black", "grey"))+
  theme_bw(base_size = 18)+
  geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
  xlab("resident burst size on *S. enterica*")+
  ylab("mutant burst size on *S. enterica*")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        strip.text = element_markdown(),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank(),
        strip.background = element_blank())


#legend
legend <- get_legend(comp_sp_to_gen_pip_values %>% 
                       mutate(pip = case_when(eigen > 0 ~ "invasion", eigen < 0 ~ "resistance", TRUE ~ "resistance")) %>%
                       ggplot(aes(x = I_S_sp, y = I_E_gen)) +
                       geom_tile(aes(fill = pip))+
                       xlab("resident infectivity on S. enterica")+
                       ylab("mutant infectivity on S. enterica")+
                       scale_fill_manual(values = c("black", "grey"))+
                       theme_bw(base_size = 18)+
                       labs(fill = " ")+
                       geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom", 
                             legend.background = element_blank(),
                             strip.background = element_blank())
)

#all parts fig 2
supp_fig2 <- plot_grid(plot_grid(partA, partB, labels = c("A", "B"), label_size = 26), 
                      plot_grid(partC, partD, labels = c("C", "D"), label_size = 26), 
                      legend,
                       ncol = 1,
                      rel_heights = c(1, 1, .2))
