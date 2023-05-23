#ATB
#Paper Fig 6
#Generation and visualization code for parts A, B

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
fig6 <- all_data %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  mutate(phage = case_when(phage == "gen" ~ "generalist (eh7)",
                           phage == "sp" ~ "specialist (p22vir)")) %>%
  filter(dilution_gen / dilution <= 5) %>%
  ggplot(aes(x = dilution_gen / dilution, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.25)+
  facet_wrap(~interaction)+
  xlab("fitness cost of generalism")+
  ylab("biomass")+
  theme_bw(base_size = 18) +
  labs(color = "species")+
  scale_color_manual(values = c("generalist (eh7)" = "#CA3542", "specialist (p22vir)" = "#27647B"))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        strip.background = element_blank())

