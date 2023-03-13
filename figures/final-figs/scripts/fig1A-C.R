#ATB
#Paper Fig 1
#Generation and visualization code for parts A, B, C

#load visualization packages
library("patchwork")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "fig1A-C-data-generation.R"))

#part A
specialist_c_comp <- specialist_c_comp %>% mutate(cost = "attachment rate", interaction = "competition")
specialist_gamma_comp <- specialist_gamma_comp %>% mutate(cost = "burst size", interaction = "competition")
specialist_c_coop <- specialist_c_coop %>% mutate(cost = "attachment rate", interaction = "cooperation")
specialist_gamma_coop <- specialist_gamma_coop %>% mutate(cost = "burst size", interaction = "cooperation")
all_data_partA <- rbind(specialist_c_comp, specialist_gamma_comp, specialist_c_coop, specialist_gamma_coop)

partA <- all_data_partA %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  pivot_longer(cols = c(c_sp, gamma_sp), names_to = "cost_type", values_to = "cost_amount") %>%
  mutate(cost_amount = case_when(cost_type == "c_sp" ~ cost_amount / 0.001,
                                 cost_type == "gamma_sp" ~ cost_amount / 20)) %>%
  ggplot(aes(x = cost_amount, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.25)+
  theme_bw()+
  facet_grid(cost~interaction) +
  xlab("fitness cost of generalism")+
  ylab("biomass")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()

#part B
gamma_and_rate_coop <- gamma_and_rate_coop  %>% mutate(cost = "burst size", parameter = "rate", alpha2 = 1, c_sp = 1e-3)
c_and_rate_coop <- c_and_rate_coop  %>% mutate(cost = "attachment rate", parameter = "rate", alpha2 = 1, gamma_sp = 20)
gamma_and_alpha_coop <- gamma_and_alpha_coop  %>% mutate(cost = "burst size", parameter = "alpha", rate_e = 0.5, c_sp = 1e-3)
c_and_alpha_coop <- c_and_alpha_coop  %>% mutate(cost = "attachment rate", parameter = "alpha", rate_e = 0.5, gamma_sp = 20)
all_data_partB <- rbind(gamma_and_rate_coop, c_and_rate_coop, gamma_and_alpha_coop, c_and_alpha_coop) %>%
  ungroup() %>%
  filter(time == max(time)) %>%
  select(alpha2, gamma_sp, cost, parameter, rate_e, c_sp, gen, sp)

partB <- all_data_partB   %>%
  mutate(alpha_ratio = alpha2 / 1,
         rate_ratio = rate_e / 0.5, 
         gamma_ratio = gamma_sp / 20,
         c_ratio = c_sp / 1e-3) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  select(cost, parameter, alpha_ratio, rate_ratio, gamma_ratio, c_ratio, normalized) %>%
  pivot_longer(cols = c(c_ratio, gamma_ratio), names_to = "cost_type", values_to = "cost_amount") %>%
  filter((cost_type == "c_ratio" & cost == "attachment rate") | (cost_type == "gamma_ratio" & cost == "burst size")) %>%
  pivot_longer(cols = c(alpha_ratio, rate_ratio), names_to = "parameter_type", values_to = "parameter_value") %>%
  filter((parameter_type == "rate_ratio" & parameter == "rate") | (parameter_type == "alpha_ratio" & parameter == "alpha")) %>%
  ggplot(aes(x = cost_amount, y = parameter_value)) +
  geom_tile(aes(fill = normalized), width=1,height=1)+
  facet_grid(parameter~cost)+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  xlab("fitness cost of generalism")+
  ylab("relative growth advantage of alternative prey")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()+
  labs(fill = "specialist relative fitness")+
  xlim(0, 5)+
  ylim(0, 5)

#part C
gamma_and_rate_comp <- gamma_and_rate_comp  %>% mutate(cost = "burst size", parameter = "rate", beta2 = 0.9, c_sp = 1e-3)
c_and_rate_comp <- c_and_rate_comp  %>% mutate(cost = "attachment rate", parameter = "rate", beta2 = 0.9, gamma_sp = 20)
gamma_and_beta_comp <- gamma_and_beta_comp  %>% mutate(cost = "burst size", parameter = "beta", rate_e = 0.5, c_sp = 1e-3)
c_and_beta_comp <- c_and_beta_comp  %>% mutate(cost = "attachment rate", parameter = "beta", rate_e = 0.5, gamma_sp = 20)
all_data_partC <- rbind(gamma_and_rate_comp, c_and_rate_comp, gamma_and_beta_comp, c_and_beta_comp) %>%
  ungroup() %>%
  filter(time == max(time)) %>%
  select(beta2, gamma_sp, cost, parameter, rate_e, c_sp, gen, sp)

partC <- all_data_partC   %>%
  mutate(beta_ratio = beta2 / 0.9,
         rate_ratio = rate_e / 0.5, 
         gamma_ratio = gamma_sp / 20,
         c_ratio = c_sp / 1e-3) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  select(cost, parameter, beta_ratio, rate_ratio, gamma_ratio, c_ratio, normalized) %>%
  pivot_longer(cols = c(c_ratio, gamma_ratio), names_to = "cost_type", values_to = "cost_amount") %>%
  filter((cost_type == "c_ratio" & cost == "attachment rate") | (cost_type == "gamma_ratio" & cost == "burst size")) %>%
  pivot_longer(cols = c(beta_ratio, rate_ratio), names_to = "parameter_type", values_to = "parameter_value") %>%
  filter((parameter_type == "rate_ratio" & parameter == "rate") | (parameter_type == "beta_ratio" & parameter == "beta")) %>%
  ggplot(aes(x = cost_amount, y = parameter_value)) +
  geom_tile(aes(fill = normalized), width=1,height=1)+
  facet_grid(parameter~cost)+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  xlab("fitness cost of generalism")+
  ylab("relative growth advantage of alternative prey")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  theme_bw()+
  labs(fill = "specialist relative fitness")+
  xlim(0, 5)+
  ylim(0, 5)

#all parts fig 1
fig1 <- partA / (partB + partC) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

#load model generation code
setwd(here::here("figures", "final-figs", "imgs"))
ggsave("fig1.png")


