#load packages and data
library("tidyverse")
library("deSolve")
library("patchwork")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

#set some initial parameters - both phage
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
maxcost = 5

#competition, both phage
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_both <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#cooperation, both phage
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop_both <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#set some initial parameters - generalist only
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0)
maxcost = 5

#competition, both phage
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_gen <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#cooperation, generalist only
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop_gen <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#set some initial parameters - specialist only
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0, sp = 0.1)
maxcost = 5

#competition, specialist only
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_sp <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#cooperation, specialist only
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop_sp <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#set some initial parameters - no phage
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0, sp = 0)
maxcost = 5

#competition, no phage
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_none <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#cooperation, no phage
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop_none <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#plots
comp_both <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2500) %>%
  mutate(cost = gamma_sp / 20) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  ggplot(aes(x = time, y = bacteria, color = type)) +
  facet_wrap(~cost) +
  geom_line(size = 0.75) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("abundance")+
  xlab("time")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))

coop_both <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2500) %>%
  mutate(cost = gamma_sp / 20) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  ggplot(aes(x = time, y = bacteria, color = type)) +
  geom_line(size = 0.75) +
  facet_wrap(~cost) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("abundance")+
  xlab("time")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))

individual_phage <- specialist_gamma_comp_gen %>% mutate(phage = "Generalist only") %>% mutate(interaction = "Competition") %>%
  rbind(., specialist_gamma_comp_sp %>% mutate(phage = "Specialist only") %>% mutate(interaction = "Competition")) %>%
  rbind(., specialist_gamma_coop_gen %>% mutate(phage = "Generalist only") %>% mutate(interaction = "Mutualism")) %>%
  rbind(., specialist_gamma_coop_sp %>% mutate(phage = "Specialist only") %>% mutate(interaction = "Mutualism")) %>%
  rbind(., specialist_gamma_comp_none %>% mutate(phage = "No phage") %>% mutate(interaction = "Competition")) %>%
  rbind(., specialist_gamma_coop_none %>% mutate(phage = "No phage") %>% mutate(interaction = "Mutualism")) %>%
  filter((gamma_sp / 20) == 2) %>%
  filter(time < 2500) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  ggplot(aes(x = time, y = bacteria, color = type)) +
  geom_line(size = 0.75) +
  facet_grid(interaction~phage) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("abundance")+
  xlab("time")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))


(comp_both + coop_both) / (individual_phage) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

