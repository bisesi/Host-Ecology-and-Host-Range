#ATB
#Visualize PIP
#Specialist to generalist equations 

#source data
source(here::here("data-generation", "model", "exploratory", "sp-to-gen-ess-sweeps.R"))
library("patchwork")

mu1_coop <- coop_sp_to_gen_pip_values %>% 
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = I_S_sp / I_S_gen_default, y = I_E_gen, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu1)

mu1_comp <- comp_sp_to_gen_pip_values %>% 
  filter(beta1 == seq(0, beta1_default * 5, by = 0.3)[4]) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = I_S_sp / I_S_gen_default, y = I_E_gen, fill = pip)) +
  geom_tile()+
  facet_wrap(~mu1)

beta1_comp <- comp_sp_to_gen_pip_values %>% 
  filter(mu1 == mu1_default) %>%
  ungroup() %>%
  mutate(pip = case_when(eigen > 0 ~ "positive", eigen < 0 ~ "negative", TRUE ~ "negative")) %>%
  ggplot(aes(x = I_S_sp / I_S_gen_default, y = I_E_gen, fill = pip)) +
  geom_tile()+
  facet_wrap(~beta1)
