#ATB
#Fig 5
#competition cost and ratio

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

#gamma and beta ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_beta2 <- seq(from = 0, to = parameters_comp['beta2'] * maxcost, by = 0.1)
gamma_sp <- seq(from = 0, to = parameters_comp['gamma_sp'] * maxcost, by = 1)

cycled_betagamma <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_beta2)){
  parameters = parameters_comp
  parameters['beta2']=s_beta2[i]
  for (j in 1:length(gamma_sp)){
    parameters['gamma_sp']=gamma_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(burst_sp = gamma_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), burst_sp = gamma_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(beta_ratio = s_beta2[i] / parameters["beta1"])
  cycled_betagamma = rbind(cycled_betagamma, output)
}

beta_gamma <- cycled_betagamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                                  final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                                  final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                                  final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = burst_sp, y = beta_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Betas (beta2 / beta1)")+
  xlab("Specialist Gamma")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

beta_gamma_fitness <- cycled_betagamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_sp, y = beta_ratio)) +
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Betas (beta2 / beta1)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()
  

#zeta and beta ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_beta2 <- seq(from = 0, to = parameters_comp['beta2'] * maxcost, by = 0.1)
c_sp <- seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = 0.00005)

cycled_betac <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_beta2)){
  parameters = parameters_comp
  parameters['beta2']=s_beta2[i]
  for (j in 1:length(c_sp)){
    parameters['c_sp']=c_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(consumption_sp = c_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), consumption_sp = c_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(beta_ratio = s_beta2[i] / parameters["beta1"])
  cycled_betac = rbind(cycled_betac, output)
}

beta_c <- cycled_betac %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                          final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                          final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                          final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = consumption_sp, y = beta_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Betas (beta2 / beta1)")+
  xlab("Specialist Zeta")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

beta_c_fitness <- cycled_betac %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = consumption_sp, y = beta_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Betas (beta2 / beta1)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

#gamma and mu ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0, to = parameters_comp['rate_e'] * maxcost, by = 0.1)
gamma_sp <- seq(from = 0, to = parameters_comp['gamma_sp'] * maxcost, by = 1)

mu_gamma <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(e_growthrange)){
  parameters = parameters_comp
  parameters['rate_e']=e_growthrange[i]
  for (j in 1:length(gamma_sp)){
    parameters['gamma_sp']=gamma_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(burst_sp = gamma_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), burst_sp = gamma_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(mu_ratio = e_growthrange[i] / parameters["rate_s"])
  mu_gamma = rbind(mu_gamma, output)
}

mus_gamma <- mu_gamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                         final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                         final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                         final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = burst_sp, y = mu_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Gamma")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

mus_gamma_fitness <- mu_gamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_sp, y = mu_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

#zeta and mu ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0, to = parameters_comp['rate_e'] * maxcost, by = 0.1)
c_sp <- seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = 0.00005)

mu_c <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(e_growthrange)){
  parameters = parameters_comp
  parameters['rate_e']=e_growthrange[i]
  for (j in 1:length(c_sp)){
    parameters['c_sp']=c_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(consumption_sp = c_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), consumption_sp = c_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(mu_ratio = e_growthrange[i] / parameters["rate_s"])
  mu_c = rbind(mu_c, output)
}

mus_c_fitness <- mu_c %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = consumption_sp, y = mu_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

mus_c <- mu_c %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = consumption_sp, y = mu_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Zeta")+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  labs(fill = " ")


#combined
fig5 <- (beta_gamma + beta_c) / (mus_gamma + mus_c) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

fig5_fitness <- (beta_gamma_fitness + beta_c_fitness) / (mus_gamma_fitness + mus_c_fitness) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


