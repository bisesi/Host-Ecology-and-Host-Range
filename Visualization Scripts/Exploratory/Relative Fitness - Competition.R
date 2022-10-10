#ATB
#relative fitness graphs - competition

#library and source scripts
library("tidyverse")
library("deSolve")
library("patchwork")

#blanket conditions
time = seq(from = 0 , to = 1e5, by = 1000)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

maxcost = 5

#Coexistence graph of gen and sp as gammas change (only sp, only gen, both but mostly one)
source(here::here("Computational Models", "Lotka Volterra Model.R"))
gamma_genE = seq(from = 0, to = parameters_comp['gamma_sp'] * maxcost, by = parameters_comp['gamma_sp'] / 20)
gamma_genS = seq(from = 0, to = parameters_comp['gamma_sp'] * maxcost, by = parameters_comp['gamma_sp'] / 20)

cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp['gamma_genS']=gamma_genS[i]
  parameters_comp['gamma_genE']=gamma_genE[i]
  output <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, maxcost, "gamma")
  gen_burst = gamma_genS[i]
  output=data.frame(output) %>%
    mutate(gen_burst = gen_burst) %>%
    mutate(cost_type = "Predator Efficiency (Gamma)")
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

burst_size <- cycled_gamma_comp %>% 
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = gen_burst, y = sp_burst))+
  geom_tile(aes(fill = normalized))+
  ylab("Specialist Gamma")+
  xlab("Generalist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#Coexistence graph of gen and sp as consumption change (only sp, only gen, both but mostly one) - comp + comp
source(here::here("Computational Models", "Lotka Volterra Model.R"))
c_genE = seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = 0.00005)
c_genS = seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = 0.00005)

cycled_c_comp <- data.frame()
for (i in 1:length(c_genS)){
  parameters_comp['c_genS']=c_genS[i]
  parameters_comp['c_genE']=c_genE[i]
  output <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, maxcost, "consumption")
  gen_c = c_genS[i]
  output=data.frame(output) %>%
    mutate(gen_c = gen_c) %>%
    mutate(cost_type = "Predator Attack Rate (c)")
  cycled_c_comp = rbind(cycled_c_comp, output)
}

consumption <- cycled_c_comp %>% 
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = gen_c, y = sp_c))+
  geom_tile(aes(fill = normalized))+
  ylab("Specialist Zeta")+
  xlab("Generalist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#Coexistence graph of gen and sp as both mus change (only sp, only gen, both but mostly one) - comp + comp
source(here::here("Computational Models", "Lotka Volterra Model.R"))
s_growthrange <- seq(from = 0, to = parameters_comp['rate_s'] * maxcost, by = 0.1)

cycled_mus <- data.frame()
for (i in 1:length(s_growthrange)){
  parameters_comp['rate_s']=s_growthrange[i]
  output <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, maxcost, 0.1)
  s_mu = s_growthrange[i]
  output=data.frame(output) %>%
    mutate(s_mu = s_mu)
  cycled_mus = rbind(cycled_mus, output)
}

mus <- cycled_mus %>% 
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
                       limits = c(0,1))

#Coexistence graph of gen and sp as both betas (only sp, only gen, both but mostly one) - comp
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
                       limits = c(0,1))

#Coexistence graph of gen and sp as mu1:mu2 and sp gamma change ~ fig 4\
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0, to = parameters_coop['rate_e'] * maxcost, by = 0.1)
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
                       limits = c(0,1))

#Coexistence graph of gen and sp as mu1:mu2 and sp consumption change
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0.3, to = parameters_comp['rate_e'] * maxcost, by = 0.1)
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

mus_c <- mu_c %>% 
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
                       limits = c(0,1))

#Coexistence graph of gen and sp as beta1:beta2 and sp gamma change
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
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_sp, y = beta_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Host Betas (beta2 / beta1)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#Coexistence graph of gen and sp as beta1:beta2 and sp consumption change
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
                       limits = c(0,1))

#Coexistence graph of gen and sp as gammagenE:gammagenS and sp gamma change
source(here::here("Computational Models", "Lotka Volterra Model.R"))
gamma_sp <- seq(from = 0, to = parameters_comp['gamma_sp'] * maxcost, by = 1)
gamma_genE <- seq(from = 0, to = parameters_comp['gamma_genE'] * maxcost, by = 1)

gamma_preference <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(gamma_genE)){
  parameters = parameters_comp
  parameters['gamma_genE']=gamma_genE[i]
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
    mutate(gen_ratio = gamma_genE[i] / parameters["gamma_genS"])
  gamma_preference = rbind(gamma_preference, output)
}

gamma_preference_gamma <- gamma_preference %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = burst_sp, y = gen_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Generalist Gamma (on E / on S)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#Coexistence graph of gen and sp as cgenE:cgenS and sp cv change
source(here::here("Computational Models", "Lotka Volterra Model.R"))
c_sp <- seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = 0.00005)
c_genE <- seq(from = 0, to = parameters_comp['c_genE'] * maxcost, by = 0.00005)

c_preference <- data.frame()
total_fitness_e <- data.frame()
for (i in 1:length(c_genE)){
  parameters = parameters_comp
  parameters['c_genE']=c_genE[i]
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
    mutate(gen_ratio = c_genE[i] / parameters["c_genS"])
  c_preference = rbind(c_preference, output)
}

c_preference_c <- c_preference %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = consumption_sp, y = gen_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab("Ratio of Generalist Zeta (on E / on S)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#cycle through resource availability vs gamma
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

burst_size_resources <- cycled_gamma_comp %>% 
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = sp_burst, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  xlab("Specialist Gamma")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))

#cycle through resource availability vs zeta
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

consumption_resources <- cycled_c_comp %>% 
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x = sp_c, y = resource_availability))+
  geom_tile(aes(fill = normalized))+
  ylab("Carrying Capacity (R)")+
  xlab("Specialist Zeta")+
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))


#look at all coexistence graphs
cost <- burst_size + consumption + plot_layout(guides = "collect")
growth <- beta + mus + plot_layout(guides = "collect")
beta_and_cost <- beta_gamma + beta_c + plot_layout(guides = "collect")
growth_and_cost <- mus_gamma + mus_c + plot_layout(guides = "collect")
preference <- gamma_preference_gamma + c_preference_c + plot_layout(guides = "collect")
resource <- burst_size_resources + consumption_resources + plot_layout(guides = "collect")

