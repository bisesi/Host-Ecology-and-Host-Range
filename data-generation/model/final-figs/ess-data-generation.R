#ATB
#Supplemental Figure 2
#ESS Modeling Data Generation

#import packages and the equations
library("tidyverse")
source(here::here("ecological-models", "ess-equations.R"))

#set default parameter for parts A and B
I_S_gen_default = 0.02
beta1_default = 0.9
mu1_default = 0.5

#initialize values to test
I_E_gen = seq(0, I_S_gen_default * 5, by = 0.0025)
I_S_sp = seq(0, I_S_gen_default * 5, by = 0.0025)
mu1 = seq(0, mu1_default * 5, by = 0.25)
beta1 = seq(0, beta1_default * 5, by = 0.3)

#Part A
coop_sp_to_gen_pip_values <- expand.grid(I_E_gen = I_E_gen, I_S_sp = I_S_sp) %>%
  group_by(I_E_gen, I_S_sp) %>%
  do(
    {
      coop_sp_to_gen(I_E_gen = .$I_E_gen, I_S_sp = .$I_S_sp) %>%
        as.data.frame()
    }
  )
coop_sp_to_gen_pip_values$eigen = coop_sp_to_gen_pip_values$.

#Part B
comp_sp_to_gen_pip_values <- expand.grid(I_E_gen = I_E_gen, I_S_sp = I_S_sp) %>%
  group_by(I_E_gen, I_S_sp) %>%
  do(
    {
      comp_sp_to_gen(I_E_gen = .$I_E_gen, I_S_sp = .$I_S_sp) %>%
        as.data.frame()
    }
  )
comp_sp_to_gen_pip_values$eigen = comp_sp_to_gen_pip_values$.

#initialize values for parts C and D
gamma_Sgen_default = 20
beta1_default = 0.9
beta2_default = 0.9
mu1_default = 0.5
mu2_default = 0.5

gamma_genS = seq(0, gamma_Sgen_default * 5, by = 10)
gamma_sp = seq(0, gamma_Sgen_default * 5, by = 10)
mu1 = seq(0, mu1_default * 5, by = 0.5)
mu2 = seq(0, mu2_default * 5, by = 0.5)
beta1 = seq(0, beta1_default * 5, by = 0.3)
beta2 = seq(0, beta2_default * 5, by = 0.3)

#part C
coop_gen_to_sp_pip_values <- expand.grid(gamma_genS = gamma_genS, gamma_sp = gamma_sp) %>%
  group_by(gamma_genS, gamma_sp) %>%
  do(
    {
      coop_gen_to_sp(gamma_genS = .$gamma_genS, gamma_sp = .$gamma_sp) %>%
        as.data.frame()
    }
  )
coop_gen_to_sp_pip_values$eigen = coop_gen_to_sp_pip_values$.

#part D
comp_gen_to_sp_pip_values <- expand.grid(gamma_genS = gamma_genS, gamma_sp = gamma_sp, mu2 = mu2, beta1 = beta1) %>%
  group_by(gamma_genS, gamma_sp, mu2, beta1) %>%
  do(
    {
      comp_gen_to_sp(gamma_genS = .$gamma_genS, gamma_sp = .$gamma_sp, beta1 = .$beta1, mu2 = .$mu2) %>%
        as.data.frame()
    }
  )
comp_gen_to_sp_pip_values$eigen = comp_gen_to_sp_pip_values$.
