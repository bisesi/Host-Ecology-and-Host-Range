#ATB
#Sensitivity analysis - Morris Screening and Sobol Variance
#completes morris screening and sobol variance analysis, cleans data using the sobol-morris-helper functions
#and then exports the results to csv files in the paper-supplement folder

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
           "c_genS", "c_sp", "dilution_E", "dilution_S", "dilution_gen", "dilution_sp", "R")

params_comp_min = c(0.1, 0.1, 0.1, 0.1, 15, 15, 15, 9e-04, 9e-04, 9e-04, 9e-04, 9e-04, 9e-04,9e-04, 0)

params_comp_max = c(2.5, 2.5, 2.5, 2.5, 65, 65, 65, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 5)

params_coop = c("alpha1", "alpha2", "rate_e", "rate_s", "gamma_genE", "gamma_genS", "gamma_sp", "c_genE",
                "c_genS", "c_sp", "dilution_E", "dilution_S", "dilution_gen", "dilution_sp", "k_e", "k_s", "R")

params_coop_min = c(0.1, 0.1, 0.1, 0.1, 15, 15, 15, 9e-04, 9e-04, 9e-04, 9e-04, 9e-04, 9e-04,9e-04, 0.1, 0.1, 0)

params_coop_max = c(2.5, 2.5, 2.5, 2.5, 65, 65, 65, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 10, 10, 5)

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
gis_morris_comp <- calculate_morris_gi(cleaned_morris_comp, params_comp) %>% 
  rename(parameter = type, gi_averaged_across_time = avg_gi, final_timepoint_gi = gi, focal_parameter = name) %>%
  select(-timepoint) %>% mutate(interaction_type = "competition")
sobol_comp <- clean_sobol_output(sobol_test_comp) %>% filter(name %in% c("gen", "sp")) %>%
  pivot_longer(beta1:R, names_to = "parameter", values_to = "value") %>%
  pivot_wider(names_from = index, values_from = value) %>% select(-time) %>% group_by(name, parameter) %>%
  mutate(avg_S = mean(S), avg_T = mean(T)) %>%
  rename(T_averaged_across_time = avg_T, S_averaged_across_time = avg_S, S_at_given_timepoint = S, T_at_given_timepoint = T) %>%
  ungroup() %>%
  group_by(name) %>%
  arrange(desc(T_averaged_across_time), .by_group = TRUE) %>%
  mutate(interaction_type = "competition")

cleaned_morris_coop <- clean_morris_output(morris_test_coop)
gis_morris_coop <- calculate_morris_gi(cleaned_morris_coop, params_coop) %>% 
  rename(parameter = type, gi_averaged_across_time = avg_gi, final_timepoint_gi = gi, focal_parameter = name) %>%
  select(-timepoint) %>% mutate(interaction_type = "mutualism")
sobol_coop <- clean_sobol_output(sobol_test_coop) %>% filter(name %in% c("gen", "sp")) %>%
  pivot_longer(alpha1:R, names_to = "parameter", values_to = "value") %>%
  pivot_wider(names_from = index, values_from = value) %>% select(-time) %>% group_by(name, parameter) %>%
  mutate(avg_S = mean(S), avg_T = mean(T)) %>%
  rename(T_averaged_across_time = avg_T, S_averaged_across_time = avg_S, S_at_given_timepoint = S, T_at_given_timepoint = T) %>%
  ungroup() %>%
  group_by(name) %>%
  arrange(desc(T_averaged_across_time), .by_group = TRUE) %>%
  mutate(interaction_type = "mutualism")

#save outputs to supplement
morris_all <- rbind(gis_morris_comp, gis_morris_coop) %>% group_by(interaction_type, focal_parameter) %>% 
  arrange(desc(gi_averaged_across_time), .by_group = TRUE)

sobol_all <- rbind(sobol_comp, sobol_coop) 

write_csv(morris_all, here::here("paper-supplement/morris-screening-global-indices.csv"))

write_csv(sobol_all, here::here("paper-supplement/sobol-analysis-first-and-total-indices.csv"))

