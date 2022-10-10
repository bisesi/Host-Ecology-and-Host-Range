#ATB
#Fig 7
#resource changes

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

#competition
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_comp_R['R'] * 5, by = 0.1)

R_range_comp <- data.frame()
for (i in 1:length(R)){
  parms = parameters_comp_R
  parms['R'] = R[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV_comp_R,
          parms = parms)
  out = data.frame(out) %>%
    mutate(carrying_capacity = R[i])
  R_range_comp <- rbind(R_range_comp, out)
}

comp <- R_range_comp %>% group_by(carrying_capacity) %>% slice_max(time) %>%
  ungroup() %>%
  mutate(generalist = round(gen, 3)) %>%
  mutate(specialist = round(sp, 3)) %>%
  mutate(alternative_prey = round(E, 3)) %>%
  mutate(shared_prey = round(S, 3)) %>%
  pivot_longer(c(generalist, specialist, alternative_prey, shared_prey), names_to = "Species") %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot() +
  geom_line(aes(x = carrying_capacity, y = normalized))+
  facet_wrap(~ Species)+
  theme_bw()+
  ylab("Normalized Equilibrium Abundance") +
  xlab("Carrying Capacity (R)")

#coperation
source(here::here("Computational Models", "Lotka Volterra Model.R"))
R = seq(from = 0, to = parameters_coop_R['R'] * 5, by = 0.1)

R_range_coop <- data.frame()
for (i in 1:length(R)){
  parms = parameters_coop_R
  parms['R'] = R[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV_coop_R,
          parms = parms)
  out = data.frame(out) %>%
    mutate(carrying_capacity = R[i])
  R_range_coop <- rbind(R_range_coop, out)
}

coop <- R_range_coop %>% group_by(carrying_capacity) %>% slice_max(time) %>%
  ungroup() %>%
  mutate(generalist = round(gen, 3)) %>%
  mutate(specialist = round(sp, 3)) %>%
  mutate(alternative_prey = round(E, 3)) %>%
  mutate(shared_prey = round(S, 3)) %>%
  pivot_longer(c(generalist, specialist, alternative_prey, shared_prey), names_to = "Species") %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot() +
  geom_line(aes(x = carrying_capacity, y = normalized))+
  facet_wrap(~ Species)+
  theme_bw()+
  ylab("Normalized Equilibrium Abundance") +
  xlab("Carrying Capacity (R)")

#combined
fig7 <- coop + comp + plot_annotation(tag_levels = "A")

