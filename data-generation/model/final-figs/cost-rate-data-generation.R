#ATB
#Modeling results generation, altering both cost and interaction coefficients 
#Used for Figure 3 and Supplemental Figure 1

#load packages and data
library("tidyverse")
library("deSolve")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

#set some initial parameters
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
maxcost = 5

#PART A
#specialist gamma only - cooperation
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 20)

specialist_gamma_coop <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#specialist c only - cooperation
c_sp = seq(from = 0, to = parameters_coop_R['c_sp'] * maxcost, by = parameters_coop_R['c_sp'] / 20)

specialist_c_coop <- expand.grid(c_sp = c_sp)%>%
  group_by(c_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1,R = 1)) %>%
        as.data.frame()
    }
  )

#specialist gamma only - competition
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c only - competition
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 20)

specialist_c_comp <- expand.grid(c_sp = c_sp)%>%
  group_by(c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#PART B - cooperation parameters
#specialist gamma and ratio of rate s to rate e
rate_e <- seq(from = 0, to = parameters_coop_R['rate_s'] * maxcost, by = 0.1)
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 10)

gamma_and_rate_coop <- expand.grid(rate_e = rate_e, gamma_sp = gamma_sp) %>%
  group_by(rate_e, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = .$rate_e, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1,R = 1)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of rate s to rate e - cooperation
rate_e <- seq(from = 0, to = parameters_coop_R['rate_s'] * maxcost, by = 0.1)
c_sp = seq(from = 0, to = parameters_coop_R['c_sp'] * maxcost, by = parameters_coop_R['c_sp'] / 10)

c_and_rate_coop <- expand.grid(rate_e = rate_e, c_sp = c_sp) %>%
  group_by(rate_e, c_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = .$rate_e, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )

#specialist gamma and ratio of alpha1 to alpha2
s_alpha2 <- seq(from = 0, to = parameters_coop_R['alpha2'] * maxcost, by = 0.2)
gamma_sp = seq(from = 0, to = parameters_coop_R['gamma_sp'] * maxcost, by = parameters_coop_R['gamma_sp'] / 10)

gamma_and_alpha_coop <- expand.grid(alpha2 = s_alpha2, gamma_sp = gamma_sp) %>%
  group_by(alpha2, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = .$alpha2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1,R = 1)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of alpha1 to alpha2
s_alpha2 <- seq(from = 0, to = parameters_coop_R['alpha2'] * maxcost, by = 0.2)
c_sp = seq(from = 0, to = parameters_coop_R['c_sp'] * maxcost, by = parameters_coop_R['c_sp'] / 10)

c_and_alpha_coop <- expand.grid(alpha2 = s_alpha2, c_sp = c_sp) %>%
  group_by(alpha2, c_sp) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = .$alpha2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, k_e = 1, k_s = 1,R = 1)) %>%
        as.data.frame()
    }
  )

#PART C - competition parameters
#specialist gamma and ratio of rate s to rate e
rate_e <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 10)

gamma_and_rate_comp <- expand.grid(rate_e = rate_e, gamma_sp = gamma_sp) %>%
  group_by(rate_e, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = .$rate_e, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of rate s to rate e
rate_e <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 10)

c_and_rate_comp <- expand.grid(rate_e = rate_e, c_sp = c_sp) %>%
  group_by(rate_e, c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = .$rate_e, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist gamma and ratio of beta1 to beta2
s_beta2 <- seq(from = 0, to = parameters_comp_R['beta2'] * maxcost, by = 0.18)
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 10)

gamma_and_beta_comp <- expand.grid(beta2 = s_beta2, gamma_sp = gamma_sp) %>%
  group_by(beta2, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = .$beta2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of beta1 to beta2
s_beta2 <- seq(from = 0, to = parameters_comp_R['beta2'] * maxcost, by = 0.18)
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 10)

c_and_beta_comp <- expand.grid(beta2 = s_beta2, c_sp = c_sp) %>%
  group_by(beta2, c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = .$beta2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2,
                  dilution_sp = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )