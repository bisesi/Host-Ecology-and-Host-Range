#ATB
#Paper Fig 2
#Generation code for parts A, B and C of figure 2

#load packages
library("tidyverse")
library("deSolve")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

beta_comp = 1
R_coop = 1

#set some initial parameters - both phage
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
maxcost = 30

#competition, both phage
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_both <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = beta_comp, beta2 = beta_comp, rate_e = 0.5, rate_s = 0.5,
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
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = R_coop)) %>%
        as.data.frame()
    }
  )

#set some initial parameters - generalist only
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0)
maxcost = 5

#competition, generalist only
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma_comp_gen <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = beta_comp, beta2 = beta_comp, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 50,
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
                  gamma_genE = 20, gamma_genS = 50,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = R_coop)) %>%
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
          parms=c(beta1 = beta_comp, beta2 = beta_comp, rate_e = 0.5, rate_s = 0.5,
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
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = R_coop)) %>%
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
          parms=c(beta1 = beta_comp, beta2 = beta_comp, rate_e = 0.5, rate_s = 0.5,
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
                  dilution_E = 3e-2, dilution_S = 3e-2, dilution_gen = 3e-2, dilution_sp = 3e-2, k_e = 1, k_s = 1, R = R_coop)) %>%
        as.data.frame()
    }
  )