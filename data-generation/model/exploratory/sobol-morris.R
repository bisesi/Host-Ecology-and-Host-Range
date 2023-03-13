#ATB
#Sensitivity analysis - Morris Screening and Sobol Variance

#load packages and data
library("deSolve")
library("ODEsensitivity")

#load models and helper functions
source(here::here("ecological-models", "lotka-volterra-model.R"))
source(here::here("functions", "sobol-morris-helper-functions.R"))

#set values for analyses
time = seq(from = 0.1, to = 1e3, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

params_comp = c("beta1", "beta2", "rate_e", "rate_s", "gamma_genE", "gamma_genS", "gamma_sp", "c_genE",
           "c_genS", "c_sp", "dilution", "R")

params_comp_min = c(0.1, 0.1, 0.1, 0.1, 15, 15, 15, 9e-04, 9e-04, 9e-04, 9e-04, 0)

params_comp_max = c(2.5, 2.5, 2.5, 2.5, 65, 65, 65, 0.01, 0.01, 0.01, 0.1, 5)

params_coop = c("alpha1", "alpha2", "rate_e", "rate_s", "gamma_genE", "gamma_genS", "gamma_sp", "c_genE",
                "c_genS", "c_sp", "dilution", "k_e", "k_s", "R")

params_coop_min = c(0.1, 0.1, 0.1, 0.1, 15, 15, 15, 9e-04, 9e-04, 9e-04, 9e-04, 0.1, 0.1, 0)

params_coop_max = c(2.5, 2.5, 2.5, 2.5, 65, 65, 65, 0.01, 0.01, 0.01, 0.1, 10, 10, 5)

#morris and sobol
morris_test_comp <- ODEmorris(mod = generalLV_comp_R, pars = params_comp,
                         state_init = start_density, times = time,
                         binf = params_comp_min, bsup = params_comp_max)

sobol_test_comp <- ODEsobol(mod = generalLV_comp_R, pars = params_comp,
                             state_init = start_density, times = time,
                            n = 2000, binf = params_comp_min, bsup = params_comp_max)

morris_test_coop <- ODEmorris(mod = generalLV_coop_R, pars = params_coop,
                              state_init = start_density, times = time,
                              binf = params_coop_min, bsup = params_coop_max)

sobol_test_coop <- ODEsobol(mod = generalLV_coop_R, pars = params_coop,
                            state_init = start_density, times = time,
                            n = 1000, binf = params_coop_min, bsup = params_coop_max)

#clean outputs
cleaned_morris_comp <- clean_morris_output(morris_test_comp)
gis_morris_comp <- calculate_morris_gi(cleaned_morris_comp, params_comp)
sobol_comp <- clean_sobol_output(sobol_test_comp)

cleaned_morris_coop <- clean_morris_output(morris_test_coop)
gis_morris_coop <- calculate_morris_gi(cleaned_morris_coop, params_coop)
sobol_coop <- clean_sobol_output(sobol_test_coop)




