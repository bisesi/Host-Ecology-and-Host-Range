#ATB
#Modeling visualization, burst size and rate parameters
#Figure 3

#load visualization packages
library("patchwork")
library("cowplot")

#load model generation code
source(here::here("data-generation", "model", "final-figs", "cost-rate-data-generation.R"))

#part A
specialist_c_comp <- specialist_c_comp %>% mutate(cost = "attachment rate", interaction = "competition")
specialist_gamma_comp <- specialist_gamma_comp %>% mutate(cost = "burst size", interaction = "competition")
specialist_c_coop <- specialist_c_coop %>% mutate(cost = "attachment rate", interaction = "mutualism")
specialist_gamma_coop <- specialist_gamma_coop %>% mutate(cost = "burst size", interaction = "mutualism")
all_data_partA <- rbind(specialist_c_comp, specialist_c_coop, specialist_gamma_comp, specialist_gamma_coop)

partA <- all_data_partA %>% ungroup() %>%
  filter(time == max(time) & cost == "burst size") %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         abundance = sp / (sp + gen),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  pivot_longer(cols = c(c_sp, gamma_sp), names_to = "cost_type", values_to = "cost_amount") %>%
  mutate(cost_amount = case_when(cost_type == "c_sp" ~ cost_amount / 0.001,
                                 cost_type == "gamma_sp" ~ cost_amount / 20)) %>%
  filter(cost_type != "c_sp") %>%
  mutate(phage = case_when(phage == "gen" ~ "generalist (eh7)",
                           phage == "sp" ~ "specialist (p22vir)")) %>%
  ggplot(aes(x = cost_amount, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.25, size = 2)+
  theme_bw(base_size = 18)+
  facet_wrap(~interaction) +
  xlab("relative burst size")+
  labs(color = "phage type")+
  scale_color_manual(values = c("generalist (eh7)" = eh7, "specialist (p22vir)" = p22vir))+
  ylab("biomass")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())

legend1 <- get_legend(all_data_partA %>% ungroup() %>%
                        filter(time == max(time) & cost == "burst size") %>%
                        mutate(gen = round(gen, 3),
                               sp = round(sp, 3),
                               relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
                        pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
                        pivot_longer(cols = c(c_sp, gamma_sp), names_to = "cost_type", values_to = "cost_amount") %>%
                        mutate(cost_amount = case_when(cost_type == "c_sp" ~ cost_amount / 0.001,
                                                       cost_type == "gamma_sp" ~ cost_amount / 20)) %>%
                        filter(cost_type != "c_sp") %>%
                        mutate(phage = case_when(phage == "gen" ~ "generalist (eh7)",
                                                 phage == "sp" ~ "specialist (p22vir)")) %>%
                        ggplot(aes(x = cost_amount, y = biomass, fill = phage))+
                        geom_bar(stat = "identity") +
                        theme_bw()+
                        facet_wrap(~interaction) +
                        xlab("fitness cost of generalism")+
                        labs(fill= "species")+
                        scale_fill_manual(values = c("Generalist (EH7)" = eh7, "Specialist (P22*vir*)" = p22vir))+
                        ylab("biomass")+
                        theme_bw(base_size = 18)+
                        theme(axis.title = element_text(), 
                              panel.background = element_rect(fill = "white"), 
                              plot.background = element_blank(),
                              legend.position = "bottom",
                              panel.grid.minor = element_blank(),
                              legend.text = element_markdown(),
                              legend.background = element_blank(),
                              strip.background = element_blank()))

#part B
gamma_and_rate_coop <- gamma_and_rate_coop  %>% mutate(cost = "burst size", parameter = "rate", alpha2 = 1, c_sp = 1e-3)
c_and_rate_coop <- c_and_rate_coop  %>% mutate(cost = "attachment rate", parameter = "rate", alpha2 = 1, gamma_sp = 20)
gamma_and_alpha_coop <- gamma_and_alpha_coop  %>% mutate(cost = "burst size", parameter = "benefit coeff", rate_e = 0.5, c_sp = 1e-3)
c_and_alpha_coop <- c_and_alpha_coop  %>% mutate(cost = "attachment rate", parameter = "benefit coeff", rate_e = 0.5, gamma_sp = 20)
all_data_partB <- rbind(gamma_and_rate_coop, c_and_rate_coop, gamma_and_alpha_coop, c_and_alpha_coop) %>%
  ungroup() %>%
  filter(time == max(time)) %>%
  select(alpha2, gamma_sp, cost, parameter, rate_e, c_sp, gen, sp)

partB <- all_data_partB   %>%
  mutate(alpha_ratio = alpha2 / 1,
         rate_ratio = rate_e / 0.5, 
         gamma_ratio = gamma_sp / 20,
         c_ratio = c_sp / 1e-3) %>%
  filter(cost == "burst size") %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         abundance = sp / (gen + sp), 
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  select(cost, parameter, alpha_ratio, rate_ratio, gamma_ratio, c_ratio, normalized, abundance) %>%
  pivot_longer(cols = c(c_ratio, gamma_ratio), names_to = "cost_type", values_to = "cost_amount") %>%
  filter((cost_type == "c_ratio" & cost == "attachment rate") | (cost_type == "gamma_ratio" & cost == "burst size")) %>%
  pivot_longer(cols = c(alpha_ratio, rate_ratio), names_to = "parameter_type", values_to = "parameter_value") %>%
  filter((parameter_type == "rate_ratio" & parameter == "rate") | (parameter_type == "alpha_ratio" & parameter == "benefit coeff")) %>%
  mutate(parameter = case_when(parameter == "rate" ~ "growth\nrate",
                               parameter == "benefit coeff" ~ "mutualism\ncoefficient",
                               TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("mutualism\ncoefficient", "growth\nrate"))) %>%
  ggplot(aes(x = cost_amount, y = parameter_value)) +
  geom_tile(aes(fill = abundance))+
  facet_wrap(~parameter)+
  scale_fill_gradient2(low = eh7,
                       mid = "white",
                       high= p22vir,
                       midpoint = 0.5,
                       limits = c(0,1))+
  xlab("relative burst size")+
  ylab("relative advantage of *E. coli*")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_markdown(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  labs(fill = "specialist relative abundance")+
  xlim(0, 5)+
  ylim(0, 5)

#part C
gamma_and_rate_comp <- gamma_and_rate_comp  %>% mutate(cost = "burst size", parameter = "rate", beta2 = 0.9, c_sp = 1e-3)
c_and_rate_comp <- c_and_rate_comp  %>% mutate(cost = "attachment rate", parameter = "rate", beta2 = 0.9, gamma_sp = 20)
gamma_and_beta_comp <- gamma_and_beta_comp  %>% mutate(cost = "burst size", parameter = "comp coeff", rate_e = 0.5, c_sp = 1e-3)
c_and_beta_comp <- c_and_beta_comp  %>% mutate(cost = "attachment rate", parameter = "comp coeff", rate_e = 0.5, gamma_sp = 20)
all_data_partC <- rbind(gamma_and_rate_comp, c_and_rate_comp, gamma_and_beta_comp, c_and_beta_comp) %>%
  ungroup() %>%
  filter(time == max(time)) %>%
  select(beta2, gamma_sp, cost, parameter, rate_e, c_sp, gen, sp)

partC <- all_data_partC %>%
  filter(cost == "burst size") %>%
  mutate(beta_ratio = beta2 / 0.9,
         rate_ratio = rate_e / 0.5, 
         gamma_ratio = gamma_sp / 20,
         c_ratio = c_sp / 1e-3) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         abundance = sp / (gen + sp),
         relative_fitness = abs((((sp - start_density["sp"]) / start_density["gen"]) / ((gen - start_density["gen"]) / start_density["sp"])))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  select(cost, parameter, beta_ratio, rate_ratio, gamma_ratio, c_ratio, normalized, abundance) %>%
  pivot_longer(cols = c(c_ratio, gamma_ratio), names_to = "cost_type", values_to = "cost_amount") %>%
  filter((cost_type == "c_ratio" & cost == "attachment rate") | (cost_type == "gamma_ratio" & cost == "burst size")) %>%
  pivot_longer(cols = c(beta_ratio, rate_ratio), names_to = "parameter_type", values_to = "parameter_value") %>%
  filter((parameter_type == "rate_ratio" & parameter == "rate") | (parameter_type == "beta_ratio" & parameter == "comp coeff")) %>%
  mutate(parameter = case_when(parameter == "rate" ~ "growth\nrate",
                               parameter == "comp coeff" ~ "competition\ncoefficient",
                               TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("competition\ncoefficient", "growth\nrate"))) %>%
  ggplot(aes(x = cost_amount, y = parameter_value)) +
  geom_tile(aes(fill = abundance), height = 0.5, width = 0.5)+
  facet_wrap(~parameter)+
  scale_fill_gradient2(low = eh7,
                       mid = "white",
                       high= p22vir,
                       midpoint = 0.5,
                       limits = c(0,1))+
  xlab("relative burst size")+
  ylab("relative advantage of *E. coli*")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_markdown(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  labs(fill = "specialist relative abundance")+
  xlim(0, 5)+
  ylim(0, 5)

legend2 <- get_legend(all_data_partC %>%
                        filter(cost == "burst size") %>%
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
                        filter((parameter_type == "rate_ratio" & parameter == "rate") | (parameter_type == "beta_ratio" & parameter == "comp coeff")) %>%
                        ggplot(aes(x = cost_amount, y = parameter_value)) +
                        geom_tile(aes(fill = normalized), width=0.5,height=0.5)+
                        facet_wrap(~parameter)+
                        scale_fill_gradient2(low = eh7,
                                             mid = "white",
                                             high= p22vir,
                                             midpoint = 0.5,
                                             limits = c(0,1))+
                        xlab("fitness cost of generalism")+
                        ylab("relative advantage of E. coli")+
                        geom_vline(xintercept = 1, linetype = "dashed")+
                        geom_hline(yintercept = 1, linetype = "dashed")+
                        theme_bw()+
                        theme(axis.title = element_text(), 
                              panel.background = element_rect(fill = "white"), 
                              plot.background = element_blank(),
                              legend.position = "bottom",
                              legend.title = element_text(size=18),
                              legend.text = element_text(size=10), 
                              panel.grid.minor = element_blank(),
                              legend.background = element_blank(),
                              strip.background = element_blank())+
                        labs(fill = "specialist relative abundance")+
                        xlim(0, 5)+
                        ylim(0, 5))

#all parts fig 1
top <- plot_grid(partA, labels = c("A"), label_size = 26, ncol = 1)
bottom <- plot_grid(partC, partB, legend1, legend2, labels = c("B", "C"), label_size = 26, ncol = 2, rel_heights = c(1, 0.15))

fig3 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(0.9, 1))
