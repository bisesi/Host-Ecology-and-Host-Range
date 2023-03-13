#ATB
#Generate data for PIP, ESS
#Generalist to specialist equations 

#import packages and the equations
library("tidyverse")
source(here::here("ecological-models", "ess-equations.R"))

#set default parameters for gamma
gamma_Sgen_default = 20
beta1_default = 0.9
beta2_default = 0.9
mu1_default = 0.5
mu2_default = 0.5

#Resident = generalist, mutant = specialist
#initialize values to test
gamma_genS = seq(0, gamma_Sgen_default * 5, by = 10)
gamma_sp = seq(0, gamma_Sgen_default * 5, by = 10)
mu1 = seq(0, mu1_default * 5, by = 0.5)
mu2 = seq(0, mu2_default * 5, by = 0.5)
beta1 = seq(0, beta1_default * 5, by = 0.9)
beta2 = seq(0, beta2_default * 5, by = 0.9)

#Cooperation
coop_gen_to_sp_pip_values_gamma <- expand.grid(gamma_genS = gamma_genS, gamma_sp = gamma_sp, mu1 = mu1, mu2 = mu2) %>%
  group_by(gamma_genS, gamma_sp, mu1, mu2) %>%
  do(
    {
      coop_gen_to_sp(gamma_genS = .$gamma_genS, gamma_sp = .$gamma_sp, mu1 = .$mu1, mu2 = .$mu2) %>%
        as.data.frame()
    }
  )
coop_gen_to_sp_pip_values_gamma$eigen = coop_gen_to_sp_pip_values_gamma$.

#Competition
comp_gen_to_sp_pip_values_gamma <- expand.grid(gamma_genS = gamma_genS, gamma_sp = gamma_sp, mu1 = mu1, mu2 = mu2, beta1 = beta1, beta2 = beta2) %>%
  group_by(gamma_genS, gamma_sp, beta1, mu1, beta2, mu2) %>%
  do(
    {
      comp_gen_to_sp(gamma_genS = .$gamma_genS, gamma_sp = .$gamma_sp, beta1 = .$beta1, mu1 = .$mu1, beta2 = .$beta2, mu2 = .$mu2) %>%
        as.data.frame()
    }
  )
comp_gen_to_sp_pip_values_gamma$eigen = comp_gen_to_sp_pip_values_gamma$.

#set default parameters for c
c_Sgen_default = 1e-3
beta1_default = 0.9
beta2_default = 0.9
mu1_default = 0.5
mu2_default = 0.5

#Resident = generalist, mutant = specialist
#initialize values to test
c_genS = seq(0, c_Sgen_default * 5, by = 0.0025)
c_sp = seq(0, c_Sgen_default * 5, by = 0.0025)
mu1 = seq(0, mu1_default * 5, by = 0.5)
mu2 = seq(0, mu2_default * 5, by = 0.5)
beta1 = seq(0, beta1_default * 5, by = 0.9)
beta2 = seq(0, beta2_default * 5, by = 0.9)

#Cooperation
coop_gen_to_sp_pip_values_c <- expand.grid(c_genS = c_genS, c_sp = c_sp, mu1 = mu1, mu2 = mu2) %>%
  group_by(c_genS, c_sp, mu1, mu2) %>%
  do(
    {
      coop_gen_to_sp(c_genS = .$c_genS, c_sp = .$c_sp, mu1 = .$mu1, mu2 = .$mu2) %>%
        as.data.frame()
    }
  )
coop_gen_to_sp_pip_values_c$eigen = coop_gen_to_sp_pip_values_c$.

#Competition
comp_gen_to_sp_pip_values_c <- expand.grid(c_genS = c_genS, c_sp = c_sp, mu1 = mu1, mu2 = mu2, beta1 = beta1, beta2 = beta2) %>%
  group_by(c_genS, c_sp, beta1, mu1, beta2, mu2) %>%
  do(
    {
      comp_gen_to_sp(c_genS = .$c_genS, c_sp = .$c_sp, beta1 = .$beta1, mu1 = .$mu1, beta2 = .$beta2, mu2 = .$mu2) %>%
        as.data.frame()
    }
  )
comp_gen_to_sp_pip_values_c$eigen = comp_gen_to_sp_pip_values_c$.

