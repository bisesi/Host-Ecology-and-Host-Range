#ATB
#Fig 2
#Abundance as cost changes

#load packages and data
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

gamma <- total %>% 
  mutate(rounded = round(abundance, 3)) %>% 
  mutate(phage = case_when(Phage == "final_sp" ~ "Specialist Predator", Phage == "final_gen" ~ "Generalist Predator", TRUE ~ Phage)) %>% 
  ggplot() + geom_line(aes(x = cost, y = rounded,color = phage, linetype = label)) + 
  theme_bw() + ylab("Equilibrium Abundance") + 
  xlab("Fitness Cost of Generalism (Conversion Efficiency)") + 
  labs(color = " ", linetype = " ") + scale_color_manual(values = c("black", "grey60")) + 
  guides(linetype = guide_legend(override.aes = list(color = "black"))) + 
  geom_vline(xintercept = 2.08, color = "black", linetype = "dashed")

#CONSUMPTION
time = seq(from = 0 , to = 1e5, by = 1000)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

maxcost = 5

#Coexistence graph of gen and sp as cs change (only sp, only gen, both but mostly one) - coop + comp
source(here::here("Computational Models", "Lotka Volterra Model.R"))
c_genE = seq(from = 0, to = parameters_coop['c_sp'] * maxcost, by = parameters_coop['c_sp'] / 20)
c_genS = seq(from = 0, to = parameters_coop['c_sp'] * maxcost, by = parameters_coop['c_sp'] / 20)

cycled_c_coop <- data.frame()
for (i in 1:length(c_genS)){
  parameters_coop['c_genS']=c_genS[i]
  parameters_coop['c_genE']=c_genE[i]
  output <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, maxcost, "consumption")
  gen_burst = c_genS[i]
  output=data.frame(output) %>%
    mutate(gen_burst = gen_burst) %>%
    mutate(cost_type = "Predator Efficiency (c)")
  cycled_c_coop = rbind(cycled_c_coop, output)
}

c_genE = seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = parameters_comp['c_sp'] / 20)
c_genS = seq(from = 0, to = parameters_comp['c_sp'] * maxcost, by = parameters_comp['c_sp'] / 20)

cycled_c_comp <- data.frame()
for (i in 1:length(c_genS)){
  parameters_comp['c_genS']=c_genS[i]
  parameters_comp['c_genE']=c_genE[i]
  output <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, maxcost, "consumption")
  gen_burst = c_genS[i]
  output=data.frame(output) %>%
    mutate(gen_burst = gen_burst) %>%
    mutate(cost_type = "Predator Efficiency (c)")
  cycled_c_comp = rbind(cycled_c_comp, output)
}

total <- rbind(cycled_c_coop, cycled_c_comp) %>%
  filter(gen_burst == 0.001) %>%
  pivot_longer(cols = c(final_sp, final_gen), names_to = "Phage", values_to = "abundance") %>%
  select(label, cost, Phage, abundance)

consumption <- total %>% 
  mutate(rounded = round(abundance, 3)) %>% 
  mutate(phage = case_when(Phage == "final_sp" ~ "Specialist Predator", Phage == "final_gen" ~ "Generalist Predator", TRUE ~ Phage)) %>% 
  ggplot() + geom_line(aes(x = cost, y = rounded,color = phage, linetype = label)) + 
  theme_bw() + ylab("Equilibrium Abundance") + 
  xlab("Fitness Cost of Generalism (Attack Rate)") + 
  labs(color = " ", linetype = " ") + 
  scale_color_manual(values = c("black", "grey60")) + 
  guides(linetype = guide_legend(override.aes = list(color = "black"))) + 
  geom_vline(xintercept = 2.13, color = "black", linetype = "dashed")

#combined
fig2 <- gamma + consumption + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
