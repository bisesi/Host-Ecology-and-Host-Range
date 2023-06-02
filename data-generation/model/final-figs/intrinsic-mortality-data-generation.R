#ATB
#Paper Fig 6
#Generation code for parts A and B
#Alter intrinsic mortality as phage cost

#load packages and data
library("tidyverse")
library("deSolve")

#load models
source(here::here("ecological-models", "lotka-volterra-model.R"))

#set some initial parameters
time = seq(from = 0.1, to = 1e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
maxcost = 5
dilution = 0.0067

#PART A
dilution_gen = seq(from = 0, to = dilution * maxcost, by = dilution / 20)


dilution_gen_coop <- expand.grid(dilution_gen = dilution_gen)%>%
  group_by(dilution_gen) %>%
  do(
    {
      ode(func=generalLV_coop_R,y=start_density,times=time,
          parms=c(alpha1 = 1, alpha2 = 1, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = dilution, dilution_S = dilution, dilution_gen = .$dilution_gen, dilution_sp = dilution,
                  k_e = 1, k_s = 1, R = 1)) %>%
        as.data.frame()
    }
  )


#PART B
dilution_gen = seq(from = 0, to = dilution * maxcost, by = dilution / 20)

dilution_gen_comp <- expand.grid(dilution_gen = dilution_gen)%>%
  group_by(dilution_gen) %>%
  do(
    {
      ode(func=generalLV_comp_R,y=start_density,times=time,
          parms=c(beta1 = 0.9, beta2 = 0.9, rate_e = 0.5, rate_s = 0.5,
                  gamma_genE = 20, gamma_genS = 20,
                  gamma_sp = 20, c_genE = 1e-3, c_genS = 1e-3, c_sp = 1e-3,
                  dilution_E = dilution, dilution_S = dilution, dilution_gen = .$dilution_gen, dilution_sp = dilution, R = 2)) %>%
        as.data.frame()
    }
  )

