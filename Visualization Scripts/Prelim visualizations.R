#prelim presentation figure

#benefit of specialization
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")
library("phaseR")

#GAMMA
time = seq(from = 0 , to = 1e5, by = 1000)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

maxcost = 5

#Coexistence graph of gen and sp as gammas change (only sp, only gen, both but mostly one) - coop + comp
source(here::here("Computational Models", "Lotka Volterra Model.R"))
gamma_genE = seq(from = 0, to = parameters_coop['gamma_sp'] * maxcost, by = parameters_coop['gamma_sp'] / 20)
gamma_genS = seq(from = 0, to = parameters_coop['gamma_sp'] * maxcost, by = parameters_coop['gamma_sp'] / 20)

cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop['gamma_genS']=gamma_genS[i]
  parameters_coop['gamma_genE']=gamma_genE[i]
  output <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, maxcost, "gamma")
  gen_burst = gamma_genS[i]
  output=data.frame(output) %>%
    mutate(gen_burst = gen_burst) %>%
    mutate(cost_type = "Predator Efficiency (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

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

total <- rbind(cycled_gamma_coop, cycled_gamma_comp) %>%
  filter(gen_burst == 20) %>%
  pivot_longer(cols = c(final_sp, final_gen), names_to = "Phage", values_to = "abundance") %>%
  select(label, cost, Phage, abundance)

plotA <- total %>% 
  mutate(rounded = round(abundance, 3)) %>% 
  mutate(phage = case_when(Phage == "final_sp" ~ "Specialist", Phage == "final_gen" ~ "Generalist", TRUE ~ Phage)) %>% 
  ggplot() + geom_line(aes(x = cost, y = rounded,color = phage), size = 2) + 
  theme_bw() + ylab("Equilibrium Abundance") + 
  xlab("Fitness Cost of Generalism (Burst Size)")+ 
  labs(color = " ")  +
  facet_wrap(~label)+
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.title=element_text(size=10))

#blanket conditions for benefit by interaction outcome
time = seq(from = 0 , to = 1e5, by = 1000)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

maxcost = 5

#mu by benefit comp
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

#mu by benefit coop
source(here::here("Computational Models", "Lotka Volterra Model.R"))
e_growthrange <- seq(from = 0, to = parameters_coop['rate_e'] * maxcost, by = 0.1)
gamma_sp <- seq(from = 0, to = parameters_coop['gamma_sp'] * maxcost, by = 1)

mu_gamma_coop <- data.frame()
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
  mu_gamma_coop = rbind(mu_gamma_coop, output)
}

#combine and make plot
mus_gamma_fitness <- mu_gamma %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  mutate(advantage = burst_sp / 20) %>%
  mutate(condition = "Competition")

plotB <- mu_gamma_coop %>% 
  mutate(final_gen = round(final_gen, 3)) %>%
  mutate(final_sp = round(final_sp, 3)) %>%
  mutate(relative_fitness = abs((((final_sp - 0.1) / 0.1) / ((final_gen - 0.1) / 0.1)))) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  mutate(advantage = burst_sp / 20) %>%
  mutate(condition = "Cooperation") %>%
  rbind(mus_gamma_fitness) %>%
  filter(advantage %in% seq(0,5,0.25)) %>%
  ggplot(aes(x = advantage, y = mu_ratio))+
  geom_tile(aes(fill = normalized))+
  ylab(" ")+
  xlab("Fitness Cost of Generalism (Burst Size)")+
  ylab("Relative Growth Rate of Specialist's Host")+
  scale_fill_gradient2(low = "#F8766D",
                       mid = "white",
                       high= "#00BFC4",
                       midpoint = 0.5,
                       limits = c(0,1))+
  labs(fill = "Specialist Relative Fitness")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()+
  facet_wrap(~condition)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.title=element_text(size=10))


