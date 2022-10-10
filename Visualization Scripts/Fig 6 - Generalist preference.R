#ATB
#Fig 6
#competition + cooperation preference a,b,c

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

#Coexistence graph of gen and sp as gammagenE:gammagenS and sp gamma change, a
source(here::here("Computational Models", "Lotka Volterra Model.R"))
parameters_comp["gamma_sp"] = 40
genE = seq(0, 40, by = 1)
genS = seq(0, 40, by = 1)

gamma_comp <- data.frame()
for (i in 1:length(genE)){
  parameters = parameters_comp
  parameters['gamma_genE']=genE[i]
  for (j in 1:length(genS)){
    parameters['gamma_genS']=genS[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out = data.frame(out)
    gamma_comp = rbind(gamma_comp, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), 
                                                   repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), 
                                                   final_sp = tail(out$sp, n = 1), 
                                                   burst_sp = parameters['gamma_sp'],
                                                   burst_E = genE[i],
                                                   burst_S = genS[j]))
  }
}

gamma_preference_comp <- gamma_comp %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = burst_E, y = burst_S))+
  geom_tile(aes(fill = coexistence))+
  ylab("Generalist Conversion Efficiency on S")+
  xlab("Generalist Conversion Efficiency on E")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 40, linetype = "dashed")+
  labs(fill = " ")


comp_gamma_fitness <- gamma_comp %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_E, y = burst_S))+
  geom_tile(aes(fill = normalized))+
  ylab("Generalist Conversion Efficiency on S")+
  xlab("Generalist Conversion Efficiency on E")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 40, linetype = "dashed")+
  labs(fill = "Specialist Relative Fitness")

#Coexistence graph of gen and sp as gammagenE:gammagenS and sp gamma change, a
source(here::here("Computational Models", "Lotka Volterra Model.R"))
parameters_coop["gamma_sp"] = 40
genE = seq(0, 40, by = 1)
genS = seq(0, 40, by = 1)

gamma_coop <- data.frame()
for (i in 1:length(genE)){
  parameters = parameters_coop
  parameters['gamma_genE']=genE[i]
  for (j in 1:length(genS)){
    parameters['gamma_genS']=genS[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
            parms = parameters)
    out = data.frame(out)
    gamma_coop = rbind(gamma_coop, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), 
                                                   repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), 
                                                   final_sp = tail(out$sp, n = 1), 
                                                   burst_sp = parameters['gamma_sp'],
                                                   burst_E = genE[i],
                                                   burst_S = genS[j]))
  }
}

gamma_preference_coop <- gamma_coop %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = burst_E, y = burst_S))+
  geom_tile(aes(fill = coexistence))+
  ylab("Generalist Conversion Efficiency on S")+
  xlab("Generalist Conversion Efficiency on E")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 40, linetype = "dashed")+
  labs(fill = " ")

coop_gamma_fitness <- gamma_coop %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_E, y = burst_S))+
  geom_tile(aes(fill = normalized))+
  ylab("Generalist Conversion Efficiency on S")+
  xlab("Generalist Conversion Efficiency on E")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 40, linetype = "dashed")+
  labs(fill = "Specialist Relative Fitness")


#Coexistence graph of gen and sp as cgenE:cgenS and sp c change, a
source(here::here("Computational Models", "Lotka Volterra Model.R"))
parameters_comp["c_sp"] = 0.002
genE = seq(0, 0.002, by = 0.00005)
genS = seq(0, 0.002, by = 0.00005)

c_comp <- data.frame()
for (i in 1:length(genE)){
  parameters = parameters_comp
  parameters['c_genE']=genE[i]
  for (j in 1:length(genS)){
    parameters['c_genS']=genS[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_comp,
            parms = parameters)
    out = data.frame(out)
    c_comp = rbind(c_comp , cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), 
                                                        repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                        final_gen = tail(out$gen, n = 1), 
                                                        final_sp = tail(out$sp, n = 1), 
                                                        burst_sp = parameters['c_sp'],
                                                        c_E = genE[i],
                                                        c_S = genS[j]))
  }
}

c_preference_comp <- c_comp %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  mutate(coexistence = case_when(coexistence == "Specialist Only" & c_E > 0.0015 ~ "Coexistence",
                                 TRUE ~ coexistence)) %>%
  ggplot(aes(x = c_E, y = c_S))+
  geom_tile(aes(fill = coexistence))+
  ylab("Generalist Attack Rate on S")+
  xlab("Generalist Attack Rate on E")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 0.002, linetype = "dashed")+
  labs(fill = " ")

comp_c_fitness <- c_comp %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  mutate(normalized = case_when(c_E > 0.0015 & normalized > 0.5 ~ 0.3,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = c_E, y = c_S))+
  geom_tile(aes(fill = normalized))+
  ylab("Generalist Attack Rate on S")+
  xlab("Generalist Attack Rate on E")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 0.002, linetype = "dashed")+
  labs(fill = "Specialist Relative Fitness")


#Coexistence graph of gen and sp as cgenE:cgenS and sp c change, a
source(here::here("Computational Models", "Lotka Volterra Model.R"))
parameters_coop["c_sp"] = 0.002
genE = seq(0, 0.002, by = 0.00005)
genS = seq(0, 0.002, by = 0.00005)

c_coop <- data.frame()
for (i in 1:length(genE)){
  parameters = parameters_coop
  parameters['c_genE']=genE[i]
  for (j in 1:length(genS)){
    parameters['c_genS']=genS[j]
    out=ode(y=start_density,
            times=time,
            func=generalLV_coop,
            parms = parameters)
    out = data.frame(out)
    c_coop = rbind(c_coop, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), 
                                                   repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]),
                                                   final_gen = tail(out$gen, n = 1), 
                                                   final_sp = tail(out$sp, n = 1), 
                                                   burst_sp = parameters['c_sp'],
                                                   c_E = genE[i],
                                                   c_S = genS[j]))
  }
}

c_preference_coop <- c_coop%>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(coexistence = case_when(final_gen < 0.001 & final_sp < 0.001 ~ "Extinction",
                                 final_gen > 0.001 & final_sp < 0.001 ~ "Generalist Only",
                                 final_gen < 0.001 & final_sp > 0.001 ~ "Specialist Only",
                                 final_gen > 0.001 & final_sp > 0.001 ~ "Coexistence")) %>%
  ggplot(aes(x = c_E, y = c_S))+
  geom_tile(aes(fill = coexistence))+
  ylab("Generalist Attack Rate on S")+
  xlab("Generalist Attack Rate on E")+
  scale_fill_manual(values = c("Extinction" = "darkblue", "Specialist Only" = "red",
                               "Coexistence" = "#999999", "Generalist Only" = "#E69F00"))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 0.002, linetype = "dashed")+
  labs(fill = " ")

coop_c_fitness <- c_coop %>%
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = c_E, y = c_S))+
  geom_tile(aes(fill = normalized))+
  ylab("Generalist Attack Rate on S")+
  xlab("Generalist Attack Rate on E")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw()+
  geom_abline(slope = -1, intercept = 0.002, linetype = "dashed")+
  labs(fill = "Specialist Relative Fitness")


#combined
fig6 <- (gamma_preference_coop + gamma_preference_comp) / (c_preference_coop + c_preference_comp) +
  plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")

fig6_fitness <- (coop_gamma_fitness + comp_gamma_fitness) / (coop_c_fitness + comp_c_fitness) +
  plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")



