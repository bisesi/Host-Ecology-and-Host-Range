#ATB
#Paper Fig 6
#Generation and visualization code for figure using intrinsic mortality

#load visualization packages
library("patchwork")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "intrinsic-mortality-data-generation.R"))

#combine datasets
dilution_range <- dilution_gen
dilution_gen_comp <- dilution_gen_comp %>% mutate(interaction = "competition")
dilution_gen_coop <- dilution_gen_coop %>% mutate(interaction = "mutualism")
all_data <- rbind(dilution_gen_comp, dilution_gen_coop) %>% mutate(cost = dilution_gen / dilution) 

#partA
figure <- all_data %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  mutate(phage = case_when(phage == "gen" ~ "Generalist (EH7)",
                           phage == "sp" ~ "Specialist (P22*vir*)")) %>%
  filter(dilution_gen / dilution <= 5) %>%
  ggplot(aes(x = dilution_gen / dilution, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.25, size = 1.5)+
  facet_wrap(~interaction)+
  xlab("relative intrinsic mortality rate of specialist")+
  ylab("biomass")+
  theme_bw(base_size = 18) +
  labs(color = "species")+
  scale_color_manual(values = c("Generalist (EH7)" = eh7, "Specialist (P22*vir*)" = p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())

legend <- get_legend(all_data %>% ungroup() %>%
                       filter(time == max(time)) %>%
                       mutate(gen = round(gen, 3),
                              sp = round(sp, 3),
                              relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
                       pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
                       mutate(phage = case_when(phage == "gen" ~ "Generalist (EH7)",
                                                phage == "sp" ~ "Specialist (P22*vir*)")) %>%
                       filter(dilution_gen / dilution <= 5) %>%
                       ggplot(aes(x = dilution_gen / dilution, y = biomass, fill = phage))+
                       geom_bar(stat = "identity")+
                       facet_wrap(~interaction)+
                       xlab("relative intrinsic mortality rate")+
                       ylab("biomass")+
                       theme_bw(base_size = 18) +
                       labs(fill= "species")+
                       scale_fill_manual(values = c("Generalist (EH7)" = eh7, "Specialist (P22*vir*)" = p22vir))+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom",
                             legend.text = element_markdown(),
                             legend.background = element_blank(),
                             strip.background = element_blank()))

fig6 <- plot_grid(figure, legend, ncol = 1, rel_heights = c(1, 0.1))

