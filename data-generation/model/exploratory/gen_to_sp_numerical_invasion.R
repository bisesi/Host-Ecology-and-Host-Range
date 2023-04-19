#test invasion analysis numerically
#competition and cooperation
#specialist invading generalist

#load packages and data
library("tidyverse")
library("deSolve")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

#set some initial parameters, 0.01% (2e-5), 1% (0.002), 10% (0.02), 50% (0.1)
time = seq(from = 0.1, to = 1e4, by = 10)
sp_0.001 = 2e-5
sp_1 = 0.002
sp_10 = 0.02
sp_50 = 0.1
start_density_0.001 <- c(E = 0.1, S = 0.1, gen = 0.1, sp = sp_0.001)
start_density_1 <- c(E = 0.1, S = 0.1, gen = 0.1, sp = sp_1)
start_density_10 <- c(E = 0.1, S = 0.1, gen = 0.1, sp = sp_10)
start_density_50 <- c(E = 0.1, S = 0.1, gen = 0.1, sp = sp_50)
maxcost = 5

#specialist gamma only
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_0.001 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density_0.001,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_comp_1 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density_1,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_comp_10 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density_10,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_comp_50 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density_50,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist gamma only, cooperation
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop_0.001 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density_0.001,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_coop_1 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density_1,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_coop_10 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density_10,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 2)) %>%
        as.data.frame()
    }
  )

specialist_gamma_coop_50 <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density_50,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 2)) %>%
        as.data.frame()
    }
  )

#figs
specialist_gamma_coop_0.001 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_0.001["sp"]) / start_density_0.001["gen"]) / ((gen - start_density_0.001["gen"]) / start_density_0.001["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp / 20, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("cooperation")

specialist_gamma_coop_1 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_1["sp"]) / start_density_1["gen"]) / ((gen - start_density_1["gen"]) / start_density_1["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp / 20, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("cooperation")

specialist_gamma_coop_10 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_10["sp"]) / start_density_10["gen"]) / ((gen - start_density_10["gen"]) / start_density_10["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp / 20, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("cooperation")

specialist_gamma_coop_50 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_50["sp"]) / start_density_50["gen"]) / ((gen - start_density_50["gen"]) / start_density_50["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp / 20, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("cooperation")

specialist_gamma_comp_0.001 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_0.001["sp"]) / start_density_0.001["gen"]) / ((gen - start_density_0.001["gen"]) / start_density_0.001["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("competition")

specialist_gamma_comp_1 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_1["sp"]) / start_density_1["gen"]) / ((gen - start_density_1["gen"]) / start_density_1["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("competition")

specialist_gamma_comp_10 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_10["sp"]) / start_density_10["gen"]) / ((gen - start_density_10["gen"]) / start_density_10["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("competition")

specialist_gamma_comp_50 %>% ungroup() %>%
  filter(time == max(time)) %>%
  mutate(gen = round(gen, 3),
         sp = round(sp, 3),
         relative_fitness = abs((((sp - start_density_50["sp"]) / start_density_50["gen"]) / ((gen - start_density_50["gen"]) / start_density_50["sp"])))) %>%
  pivot_longer(cols = c(gen, sp), names_to = "phage", values_to = "biomass") %>%
  ggplot(aes(x = gamma_sp, y = biomass, color = phage))+
  geom_smooth(se = FALSE, span = 0.2)+
  theme_bw()+
  xlab("benefit of specialism (gamma)")+
  ylab("biomass")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("competition")