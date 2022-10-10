#ATB
#Fig 8
#resource changes and specialist benefit

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")
library("phaseR")

#blanket conditions
time = seq(from = 0 , to = 1e5, by = 1000)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

maxcost = 5

#comperation and gamma
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_comp_R['R'] * maxcost, by = 0.1)

cycled_gamma_comp <- data.frame()
for (i in 1:length(R)){
  parameters_comp_R['R']=R[i]
  output <- LVgeneralismcost(generalLV_comp_R, parameters_comp_R, "Competition", start_density, time, maxcost, "gamma")
  resources = R[i]
  output=data.frame(output) %>%
    mutate(resource_availability = resources) %>%
    mutate(cost_type = "Predator Efficiency (Gamma)")
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

comp_coexistence_gamma <- cycled_gamma_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = sp_burst, y = resource_availability))+
  geom_tile(aes(fill = coexistence))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Conversion Efficiency)")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_vline(xintercept = 20, linetype = "dashed")+
  labs(fill = " ")

comp_coexistence_gamma_fitness <- cycled_gamma_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))  %>%
  ggplot(aes(x = sp_burst, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Conversion Efficiency)")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  theme_bw()+
  geom_vline(xintercept = 20, linetype = "dashed")

#comperation and c
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_comp_R['R'] * maxcost, by = 0.1)

cycled_c_comp <- data.frame()
for (i in 1:length(R)){
  parameters_comp_R['R']=R[i]
  output <- LVgeneralismcost(generalLV_comp_R, parameters_comp_R, "Competition", start_density, time, maxcost, "consumption")
  resources = R[i]
  output=data.frame(output) %>%
    mutate(resource_availability = resources)  %>%
    mutate(cost_type = "Predator Attack Rate (c)")
  cycled_c_comp = rbind(cycled_c_comp, output)
}

comp_coexistence_c_fitness <- cycled_c_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = sp_c, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Attack rate)")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  theme_bw()+
  geom_vline(xintercept = 0.001, linetype = "dashed")
  
comp_coexistence_c <- cycled_c_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = sp_c, y = resource_availability))+
  geom_tile(aes(fill = coexistence))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Attack rate)")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  labs(fill = " ")


#cooperation and gamma
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_coop_R['R'] * maxcost, by = 0.1)

parameters_coop_R['k_s'] = 0.5
parameters_coop_R['k_e'] = 0.5

cycled_gamma_coop <- data.frame()
for (i in 1:length(R)){
  parameters_coop_R['R']=R[i]
  output <- LVgeneralismcost(generalLV_coop_R, parameters_coop_R, "Cooperation", start_density, time, maxcost, "gamma")
  resources = R[i]
  output=data.frame(output) %>%
    mutate(resource_availability = resources) %>%
    mutate(cost_type = "Predator Efficiency (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

coop_coexistence_gamma <- cycled_gamma_coop %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  filter(resource_availability <= 3) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = sp_burst, y = resource_availability))+
  geom_tile(aes(fill = coexistence))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Conversion Efficiency)")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_vline(xintercept = 20, linetype = "dashed")+
  labs(fill = " ")

coop_coexistence_gamma_fitness <- cycled_gamma_coop %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  filter(resource_availability <= 3) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = sp_burst, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  theme_bw()+
  geom_vline(xintercept = 20, linetype = "dashed")


#cooperation and c
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_coop_R['R'] * maxcost, by = 0.1)

cycled_c_coop <- data.frame()
for (i in 1:length(R)){
  parameters_coop_R['R']=R[i]
  output <- LVgeneralismcost(generalLV_coop_R, parameters_coop_R, "Cooperation", start_density, time, maxcost, "consumption")
  resources = R[i]
  output=data.frame(output) %>%
    mutate(resource_availability = resources)  %>%
    mutate(cost_type = "Predator Attack Rate (c)")
  cycled_c_coop = rbind(cycled_c_coop, output)
}

coop_coexistence_c <- cycled_c_coop %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  filter(resource_availability <= 3) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = sp_c, y = resource_availability))+
  geom_tile(aes(fill = coexistence))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Attack rate)")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  labs(fill = " ")

coop_coexistence_c_fitness <- cycled_c_coop %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  filter(resource_availability <= 3) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))  %>%
  ggplot(aes(x = sp_c, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  xlab("Cost of Generalism (Attack rate)")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  theme_bw()+
  geom_vline(xintercept = 0.001, linetype = "dashed")

#combined
fig8 <- (coop_coexistence_gamma + comp_coexistence_gamma) / (coop_coexistence_c + comp_coexistence_c) +
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")

fig8_fitness <- (coop_coexistence_gamma_fitness + comp_coexistence_gamma_fitness) / (coop_coexistence_c_fitness + comp_coexistence_c_fitness) +
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")
