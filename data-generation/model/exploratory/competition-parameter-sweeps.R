#ATB
#Sensitivity Analyses
#Competition, Individual Parameter Sweeps

#load packages and data
library("tidyverse")
library("deSolve")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

#set some initial parameters
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
maxcost = 5

#specialist gamma only
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

specialist_gamma <- expand.grid(gamma_sp = gamma_sp)%>%
  group_by(gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c only
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 20)

specialist_c <- expand.grid(c_sp = c_sp)%>%
  group_by(c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#beta1 and beta2
s_beta2 <- seq(from = 0, to = parameters_comp_R['beta2'] * maxcost, by = 0.1)
e_beta1 <- seq(from = 0, to = parameters_comp_R['beta1'] * maxcost, by = 0.1)

betas <- expand.grid(beta2 = s_beta2, beta1 = e_beta1)%>%
  group_by(beta2, beta1) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = .$beta1, beta2 = .$beta2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#rate e and rate s
rate_s <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)
rate_e <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)

mus <- expand.grid(rate_s = rate_s, rate_e = rate_e) %>%
  group_by(rate_s, rate_e) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = .$rate_e, rate_s = .$rate_s,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist gamma and ratio of rate s to rate e
rate_s <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

gamma_and_rate <- expand.grid(rate_s = rate_s, gamma_sp = gamma_sp) %>%
  group_by(rate_s, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = .$rate_s,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of rate s to rate e
rate_s <- seq(from = 0, to = parameters_comp_R['rate_s'] * maxcost, by = 0.1)
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 20)

c_and_rate <- expand.grid(rate_s = rate_s, c_sp = c_sp) %>%
  group_by(rate_s, c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = .$rate_s,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist gamma and ratio of beta1 to beta2
s_beta2 <- seq(from = 0, to = parameters_comp_R['beta2'] * maxcost, by = 0.1)
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)

gamma_and_beta <- expand.grid(beta2 = s_beta2, gamma_sp = gamma_sp) %>%
  group_by(beta2, gamma_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = .$beta2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#specialist c and ratio of beta1 to beta2
s_beta2 <- seq(from = 0, to = parameters_comp_R['beta2'] * maxcost, by = 0.1)
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 20)

c_and_beta <- expand.grid(beta2 = s_beta2, c_sp = c_sp) %>%
  group_by(beta2, c_sp) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = .$beta2, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution = 3e-2, R = 2)) %>%
        as.data.frame()
    }
  )

#R and ratio of specialist gamma to generalist gamma
gamma_sp = seq(from = 0, to = parameters_comp_R['gamma_sp'] * maxcost, by = parameters_comp_R['gamma_sp'] / 20)
R = seq(from = 0, to = parameters_comp_R['R'] * maxcost, by = 0.1)

r_and_gamma <- expand.grid(gamma_sp = gamma_sp, R = R) %>%
  group_by(gamma_sp, R) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = .$gamma_sp, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2, R = .$R)) %>%
        as.data.frame()
    }
  )

#R and ratio of specialist c to generalist c
c_sp = seq(from = 0, to = parameters_comp_R['c_sp'] * maxcost, by = parameters_comp_R['c_sp'] / 20)
R = seq(from = 0, to = parameters_comp_R['R'] * maxcost, by = 0.1)

r_and_c <- expand.grid(c_sp = c_sp, R = R) %>%
  group_by(c_sp, R) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = .$c_sp,
                  dilution = 3e-2, R = .$R)) %>%
        as.data.frame()
    }
  )

#just R
R = seq(from = 0, to = parameters_comp_R['R'] * maxcost, by = 0.1)

r_alone <- expand.grid(R = R) %>%
  group_by(R) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution = 3e-2,R = .$R)) %>%
        as.data.frame()
    }
  )


