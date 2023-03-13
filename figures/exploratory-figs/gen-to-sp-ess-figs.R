#ATB
#Visualize PIP
#Generalist to specialist equations 

#source data
source(here::here("data-generation", "model", "exploratory", "gen-to-sp-ess-sweeps.R"))
library("patchwork")

#gamma
mu1_coop_gamma <- coop_gen_to_sp_pip_values_gamma %>% 
  filter(mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu1)

mu2_coop_gamma <- coop_gen_to_sp_pip_values_gamma %>% 
  filter(mu1 == mu1_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu2)

mu1_comp_gamma <- comp_gen_to_sp_pip_values_gamma %>% 
  filter(beta1 == beta1_default & beta2 == beta2_default & mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu1)

beta1_comp_gamma <- comp_gen_to_sp_pip_values_gamma %>% 
  filter(mu1 == mu1_default & beta2 == beta2_default & mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~beta1)

mu2_comp_gamma <- comp_gen_to_sp_pip_values_gamma %>% 
  filter(beta1 == beta1_default & beta2 == beta2_default & mu1 == mu1_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu2)

beta2_comp_gamma <- comp_gen_to_sp_pip_values_gamma %>% 
  filter(mu1 == mu1_default & beta1 == beta1_default & mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = gamma_sp, y = gamma_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~beta2)

#c
mu1_coop_c <- coop_gen_to_sp_pip_values_c %>% 
  filter(mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = c_sp , y = c_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu1)

mu2_coop_c <- coop_gen_to_sp_pip_values_c %>% 
  filter(mu1 == mu1_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = c_sp , y = c_genS, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu2)

betas_comp_c <- comp_gen_to_sp_pip_values_c %>% 
  filter(mu1 == mu1_default & mu2 == mu2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = c_sp, y = c_genS, fill = pip)) +
  geom_tile()+
  facet_grid(beta2~beta1)

mus_comp <- comp_gen_to_sp_pip_values_c %>% 
  filter(beta1 == beta1_default & beta2 == beta2_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = c_sp, y = c_genS, fill = pip)) +
  geom_tile()+
  facet_grid(mu1~mu2)