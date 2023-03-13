#ATB
#Generate data for PIP, ESS
#Specialist to generalist equations 

#import packages and the equations
library("tidyverse")
source(here::here("ecological-models", "ess-equations.R"))

#set default parameter
I_S_gen_default = 0.02
beta1_default = 0.9
mu1_default = 0.5

#Resident = specialist, mutant = generalist
#initialize values to test
I_E_gen = seq(0, I_S_gen_default * 5, by = 0.0025)
I_S_sp = seq(0, I_S_gen_default * 5, by = 0.0025)
mu1 = seq(0, mu1_default * 5, by = 0.25)
beta1 = seq(0, beta1_default * 5, by = 0.3)

#Cooperation, PIP with mu1, I_E_gen, and changing I_S_sp / I_S_gen ratio
coop_sp_to_gen_pip_values <- expand.grid(I_E_gen = I_E_gen, I_S_sp = I_S_sp, mu1 = mu1) %>%
  group_by(I_E_gen, I_S_sp, mu1) %>%
  do(
    {
      coop_sp_to_gen(I_E_gen = .$I_E_gen, I_S_sp = .$I_S_sp, mu1 = .$mu1) %>%
        as.data.frame()
    }
  )
coop_sp_to_gen_pip_values$eigen = coop_sp_to_gen_pip_values$.

#Competition
comp_sp_to_gen_pip_values <- expand.grid(I_E_gen = I_E_gen, I_S_sp = I_S_sp, mu1 = mu1, beta1 = beta1) %>%
  group_by(I_E_gen, I_S_sp, beta1, mu1) %>%
  do(
    {
      comp_sp_to_gen(I_E_gen = .$I_E_gen, I_S_sp = .$I_S_sp, beta1 = .$beta1, mu1 = .$mu1) %>%
        as.data.frame()
    }
  )
comp_sp_to_gen_pip_values$eigen = comp_sp_to_gen_pip_values$.


