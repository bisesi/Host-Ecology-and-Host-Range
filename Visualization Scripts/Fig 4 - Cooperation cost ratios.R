#ATB
#Fig 4
#Cooperation cost and ratio

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

#gamma and alpha ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_alpha2 <- seq(from = 0, to = parameters_coop['alpha2'] * maxcost, by = 0.1)
gamma_sp <- seq(from = 0, to = parameters_coop['gamma_sp'] * maxcost, by = 1)

cycled_alphagamma <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_alpha2)){
  parameters = parameters_coop
  parameters['alpha2']=s_alpha2[i]
  for (j in 1:length(gamma_sp)){
    parameters['gamma_sp']=gamma_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(burst_sp = gamma_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), burst_sp = gamma_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(alpha_ratio = s_alpha2[i] / parameters["alpha1"])
  cycled_alphagamma = rbind(cycled_alphagamma, output)
}

alpha_gamma <- cycled_alphagamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                                    final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                                    final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                                    final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = burst_sp, y = alpha_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Alphas (alpha2 / alpha1)")+
  xlab("Specialist Gamma")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

alpha_gamma_fitness <- cycled_alphagamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_sp, y = alpha_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Alphas (alpha2 / alpha1)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#zeta and alpha ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_alpha2 <- seq(from = 0, to = parameters_coop['alpha2'] * maxcost, by = 0.1)
c_sp <- seq(from = 0, to = parameters_coop['c_sp'] * maxcost, by = 0.00005)

cycled_alphac <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_alpha2)){
  parameters = parameters_coop
  parameters['alpha2']=s_alpha2[i]
  for (j in 1:length(c_sp)){
    parameters['c_sp']=c_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(consumption_sp = c_sp[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), consumption_sp = c_sp[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(alpha_ratio = s_alpha2[i] / parameters["alpha1"])
  cycled_alphac = rbind(cycled_alphac, output)
}

alpha_c <- cycled_alphac %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = consumption_sp, y = alpha_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Alphas (alpha2 / alpha1)")+
  xlab("Specialist Zeta")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

alpha_c_fitness <- cycled_alphac %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = consumption_sp, y = alpha_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Alphas (alpha2 / alpha1)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#gamma and mu ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0.3, to = parameters_coop['rate_e'] * maxcost, by = 0.1)
gamma_sp <- seq(from = 0, to = parameters_coop['gamma_sp'] * maxcost, by = 1)

mu_gamma <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(e_growthrange)){
  parameters = parameters_coop
  parameters['rate_e']=e_growthrange[i]
  for (j in 1:length(gamma_sp)){
    parameters['gamma_sp']=gamma_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
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
  mutate(coexistence = case_when(final_gen < 0.1 & final_sp < 0.1 ~ "Extinction",
                                                         final_gen > 0.1 & final_sp < 0.1 ~ "Generalist Only",
                                                         final_gen < 0.1 & final_sp > 0.1 ~ "Specialist Only",
                                                         final_gen > 0.1 & final_sp > 0.1 ~ "Coexistence")) %>%
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
  geom_vline(xintercept = 20, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#zeta and mu ratio
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0.3, to = parameters_coop['rate_e'] * maxcost, by = 0.1)
c_sp <- seq(from = 0, to = parameters_coop['c_sp'] * maxcost, by = 0.00005)

mu_c <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(e_growthrange)){
  parameters = parameters_coop
  parameters['rate_e']=e_growthrange[i]
  for (j in 1:length(c_sp)){
    parameters['c_sp']=c_sp[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
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

mus_c <- mu_c %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  select(final_gen, final_sp, mu_ratio, consumption_sp) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence"))%>%
  ggplot(aes(x = consumption_sp, y = mu_ratio))+
  geom_tile(aes(fill = coexistence))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Zeta")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")
  
mus_c_fitness <- mu_c %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  select(final_gen, final_sp, mu_ratio, consumption_sp) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))%>%
  ggplot(aes(x = consumption_sp, y = mu_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Mus (E / S)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_vline(xintercept = 0.001, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")


#combined
fig4 <- (alpha_gamma + alpha_c) / (mus_gamma + mus_c) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
  
fig4_fitness <- (alpha_gamma_fitness + alpha_c_fitness) / (mus_gamma_fitness + mus_c_fitness) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
  
  