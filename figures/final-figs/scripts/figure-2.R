#ATB
#Generate bacterial dynamics in a Lotka Volterra model figure
#Figure 2

#load packages and data
library("tidyverse")
library("deSolve")
library("patchwork")
library("cowplot")

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
maxcost = 15

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
maxcost = 15

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
maxcost = 15

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

#partA - application of individual phage
partA <- specialist_gamma_comp_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "competition") %>%
  rbind(., specialist_gamma_comp_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_gen %>% mutate(phage = "generalist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_coop_sp %>% mutate(phage = "specialist only") %>% mutate(interaction = "mutualism")) %>%
  rbind(., specialist_gamma_comp_none %>% mutate(phage = "no phage") %>% mutate(interaction = "competition")) %>%
  rbind(., specialist_gamma_coop_none %>% mutate(phage = "no phage") %>% mutate(interaction = "mutualism")) %>%
  filter((gamma_sp / 20) == 10) %>%
  filter(time < 2000) %>%
  pivot_longer(E:S, names_to = "type", values_to = "bacteria") %>%
  mutate(type = case_when(type == "E" ~ "E. coli", 
                          type == "S" ~ "S. enterica")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "generalist only", "specialist only"))) %>%
  ggplot(aes(x = time / 1000, y = bacteria, color = type)) +
  geom_line(size = 1) +
  facet_grid(interaction~phage) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  labs(color = "species")+
  ylim(0, 2)+
  scale_color_manual(values = c("#5ba300", "#e6308a"))

#partB - high cost in competition
partB_main <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                          microbe == "S" ~ "S. enterica",
                          microbe == "gen" ~ "generalist (phi-C)",
                          microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "generalist (phi-C)" | microbe == "specialist (p22vir)") %>%
  ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
  geom_line(size = 1) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                "generalist (phi-C)" = "#CA3542", "specialist (p22vir)" = "#27647B"))

partB_inset <- specialist_gamma_comp_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "generalist (phi-C)",
                             microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "E. coli" | microbe == "S. enterica") %>%
  ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
  geom_line(size = 1) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 2)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                "generalist (phi-C)" = "#CA3542", "specialist (p22vir)" = "#27647B"))

partB <- ggdraw(partB_main) + 
  draw_plot(partB_inset, x = 0.15, y = 0.6, width = 0.35, height = 0.35)

#partC - cooperation at high cost
partC_main <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "generalist (phi-C)",
                             microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "generalist (phi-C)" | microbe == "specialist (p22vir)") %>%
  ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
  geom_line(size = 1) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 250)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                "generalist (phi-C)" = "#CA3542", "specialist (p22vir)" = "#27647B"))

partC_inset <- specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
  filter(time < 2000) %>%
  mutate(cost = gamma_sp / 20) %>%
  filter(cost == 5) %>%
  pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
  mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                             microbe == "S" ~ "S. enterica",
                             microbe == "gen" ~ "generalist (phi-C)",
                             microbe == "sp" ~ "specialist (p22vir)")) %>%
  filter(microbe == "E. coli" | microbe == "S. enterica") %>%
  ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
  geom_line(size = 1) +
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("biomass")+
  xlab("time (a.u.)")+
  ylim(0, 2)+
  labs(color = "species")+
  scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                "generalist (phi-C)" = "#CA3542", "specialist (p22vir)" = "#27647B"))

partC <- ggdraw(partC_main) + 
  draw_plot(partC_inset, x = 0.15, y = 0.6, width = 0.35, height = 0.35)

legend <- get_legend(specialist_gamma_coop_both %>% filter((gamma_sp / 20) %in% seq(0, 5, by = 1)) %>%
                       filter(time < 2000) %>%
                       mutate(cost = gamma_sp / 20) %>%
                       filter(cost == 5) %>%
                       pivot_longer(E:sp, names_to = "microbe", values_to = "biomass") %>%
                       mutate(microbe = case_when(microbe == "E" ~ "E. coli", 
                                                  microbe == "S" ~ "S. enterica",
                                                  microbe == "gen" ~ "generalist (phi-C)",
                                                  microbe == "sp" ~ "specialist (p22vir)")) %>%
                       filter(microbe == "E. coli" | microbe == "S. enterica") %>%
                       ggplot(aes(x = time / 1000, y = biomass, color = microbe)) +
                       geom_line(size = 1) +
                       theme_bw()+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom",
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("biomass")+
                       xlab("time (a.u.)")+
                       ylim(0, 2)+
                       labs(color = "species")+
                       scale_color_manual(values = c("E. coli" = "#5ba300", "S. enterica" = "#e6308a", 
                                                     "generalist (phi-C)" = "#CA3542", "specialist (p22vir)" = "#27647B")))

#full plot
fig2 <- plot_grid(partA, 
                  plot_grid(partB, partC, labels = c("B", "C")), 
                  ncol = 1,
                  labels = c("A"), 
                  legend, rel_heights = c(1, 1, .1))
