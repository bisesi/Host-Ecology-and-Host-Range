#ATB
#Fig 3
#Coexistence 

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

#competition beta 
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_beta2 <- seq(from = 0, to = parameters_comp['beta2'] * maxcost, by = 0.1)
e_beta1 <- seq(from = 0, to = parameters_comp['beta1'] * maxcost, by = 0.1)

cycled_beta <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_beta2)){
  parameters = parameters_comp
  parameters['beta2']=s_beta2[i]
  for (j in 1:length(e_beta1)){
    parameters['beta1']=e_beta1[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(e_beta = e_beta1[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), e_beta = e_beta1[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(s_beta = s_beta2[i])
  cycled_beta = rbind(cycled_beta, output)
}

beta <- cycled_beta %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = e_beta, y = s_beta))+
  geom_tile(aes(fill = coexistence))+
  ylab("Beta2")+
  xlab("Beta1")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

beta_fitness <- cycled_beta %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = e_beta, y = s_beta))+
  geom_tile(aes(fill = normalized))+
  ylab("Beta2")+
  xlab("Beta1")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#competition mu
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_growthrange <- seq(from = 0, to = parameters_comp['rate_s'] * maxcost, by = 0.1)

cycled_mus_comp <- data.frame()
for (i in 1:length(s_growthrange)){
  parameters_comp['rate_s']=s_growthrange[i]
  output <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, maxcost, 0.1)
  s_mu = s_growthrange[i]
  output=data.frame(output) %>%
    mutate(s_mu = s_mu)
  cycled_mus_comp = rbind(cycled_mus_comp, output)
}

mus_comp <- cycled_mus_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                     final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                     final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                     final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = e_mu, y = s_mu))+
  geom_tile(aes(fill = coexistence))+
  ylab("S Mu")+
  xlab("E Mu")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

mus_comp_fitness <- cycled_mus_comp %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = e_mu, y = s_mu))+
  geom_tile(aes(fill = normalized))+
  ylab("S Mu")+
  xlab("E Mu")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#cooperation alpha
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_alpha2 <- seq(from = 0, to = parameters_coop['alpha2'] * maxcost, by = 0.1)
e_alpha1 <- seq(from = 0, to = parameters_coop['alpha1'] * maxcost, by = 0.1)

cycled_alpha <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(s_alpha2)){
  parameters = parameters_coop
  parameters['alpha2']=s_alpha2[i]
  for (j in 1:length(e_alpha1)){
    parameters['alpha1']=e_alpha1[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
            parms = parameters)
    out=data.frame(out) %>%
      mutate(e_alpha = e_alpha1[j])
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), final_sp = tail(out$sp, n = 1), e_alpha = e_alpha1[j]))
  }
  output=data.frame(total_fitness_e) %>%
    mutate(s_alpha = s_alpha2[i])
  cycled_alpha = rbind(cycled_alpha, output)
}

alpha <- cycled_alpha %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                         final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                         final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                         final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = e_alpha, y = s_alpha))+
  geom_tile(aes(fill = coexistence))+
  ylab("Alpha2")+
  xlab("Alpha1")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

alpha_fitness <- cycled_alpha %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = e_alpha, y = s_alpha))+
  geom_tile(aes(fill = normalized))+
  ylab("Alpha2")+
  xlab("Alpha1")+
  scale_fill_gradient2(low = "#CA3542",
                                      mid = "white",
                                      high= "#27647B",
                                      midpoint = 0.5,
                                      limits = c(0,1))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#cooperation mu
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_growthrange <- seq(from = 0.3, to = parameters_coop['rate_s'] * maxcost, by = 0.1)

cycled_mus <- data.frame()
for (i in 1:length(s_growthrange)){
  parameters_coop['rate_s']=s_growthrange[i]
  output <- LVgrowthrate(generalLV_coop, parameters_coop, "Cooperation", start_density, time, maxcost, 0.1)
  s_mu = s_growthrange[i]
  output=data.frame(output) %>%
    mutate(s_mu = s_mu)
  cycled_mus = rbind(cycled_mus, output)
}

mus_coop <- cycled_mus %>% 
  select(final_gen, final_sp, e_mu, s_mu) %>%
  add_row(e_mu = rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))), 
          s_mu = rep(seq(0, 2.5, by = 0.1), length(seq(0, 0.3, by = 0.1))),
          final_gen = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))))),
          final_sp = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1)))))) %>%
  add_row(s_mu = rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))), 
          e_mu = rep(seq(0, 2.5, by = 0.1), length(seq(0, 0.3, by = 0.1))),
          final_gen = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))))),
          final_sp = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1)))))) %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                                     final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                                     final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                                     final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = as.numeric(e_mu), y = as.numeric(s_mu)))+
  geom_tile(aes(fill = coexistence), height = 0.2, width = 0.2)+
  ylab("S Mu")+
  xlab("E Mu")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = " ")

mus_coop_fitness <- cycled_mus %>% 
  select(final_gen, final_sp, e_mu, s_mu) %>%
  add_row(e_mu = rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))), 
          s_mu = rep(seq(0, 2.5, by = 0.1), length(seq(0, 0.3, by = 0.1))),
          final_gen = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))))),
          final_sp = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1)))))) %>%
  add_row(s_mu = rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))), 
          e_mu = rep(seq(0, 2.5, by = 0.1), length(seq(0, 0.3, by = 0.1))),
          final_gen = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1))))),
          final_sp = rep(0, length(rep(seq(0, 0.3, by = 0.1), length(seq(0, 2.5, by = 0.1)))))) %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = as.numeric(e_mu), y = as.numeric(s_mu)))+
  geom_tile(aes(fill = normalized), height = 0.2, width = 0.2)+
  ylab("S Mu")+
  xlab("E Mu")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme_bw()+
  labs(fill = "Specialist Relative Fitness")

#patched
fig3 <- (beta + alpha) / (mus_comp + mus_coop) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
fig3_fitness <- (beta_fitness + alpha_fitness) / (mus_comp_fitness + mus_coop_fitness) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
