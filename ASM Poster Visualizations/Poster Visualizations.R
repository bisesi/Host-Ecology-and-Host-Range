#ASM Poster Visualizations
#ATV
#May 30, 2022

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#intro graphs
#no phage
time = seq(from = 0 , to = 5e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0, sp = 0)

out_coop <- as.data.frame(ode(y=start_density,
                times=time,
                func=generalLV_coop,
                parms = parameters_coop,
                atol = 1e-14)) %>%
  mutate(label = "Cooperation")

out_comp <- as.data.frame(ode(y=start_density,
                times=time,
                func=generalLV_comp,
                parms = parameters_comp,
                atol = 1e-14)) %>%
  mutate(label = "Competition")

nophage <- rbind(out_coop, out_comp) %>%
  melt(id = c("time", "label")) %>%
  rename(Microbe = variable) %>%
  filter(Microbe == "E" | Microbe == "S") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(label ~ Microbe) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text()) + 
  ylab("Abundance")+
  xlab("Time")

#with phage, no cost
time = seq(from = 0 , to = 5e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

out_coop_with <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_coop,
                              parms = parameters_coop,
                              atol = 1e-14)) %>%
  mutate(label = "Cooperation")

out_comp_with <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_comp,
                              parms = parameters_comp,
                              atol = 1e-14)) %>%
  mutate(label = "Competition")

withphagehosts <- rbind(out_coop_with, out_comp_with) %>%
  melt(id = c("time", "label")) %>%
  rename(Microbe = variable) %>%
  filter(Microbe == "E" | Microbe == "S") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(label ~ Microbe) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text()) + 
  ylab("Abundance")+
  xlab("Time")

withphagephage <- rbind(out_coop_with, out_comp_with) %>%
  melt(id = c("time", "label")) %>%
  rename(Microbe = variable) %>%
  filter(Microbe == "gen" | Microbe == "sp") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(label ~ Microbe) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text()) + 
  ylab("Abundance")+
  xlab("Time")

#example plots
proof_of_concept <- nophage + (withphagehosts / withphagephage) +
  plot_annotation(tag_levels = 'A')

#ONLY COST OF GENERALISM
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

costofgeneralismcoop <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 10)
costofgeneralismcomp <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 10)

completedatacost <- rbind(costofgeneralismcoop, costofgeneralismcomp)

datasetbasegencost <- completedatacost %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

relativefitnessbasegencost <- datasetbasegencost %>%
  rbind(., c(10.0, "Cooperation", "Gamma", Inf)) %>%
  filter(type != "Consumption Rate") %>%
  mutate(log10 = log10(as.numeric(value))) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5,
                           TRUE ~ log10)) %>%
  ggplot(aes(x=as.numeric(cost), y = log10, color = Interaction)) +
  geom_line(size = 2)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text()) + 
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist")+
  xlab("Fold-Increase in Specialist Productivity")+
  ylim(-5, 5)

#COST AND GROWTH RATE
#growth rate and cost of generalism, cooperation
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_coop
gamma <- seq(from = 0, to = parameters['gamma_sp'] * 5, by = 0.01)

cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma)){
  parameters_coop['gamma_sp']=gamma[i]
  output <- LVgrowthrate(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, 0.1)
  generalism_cost = gamma[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

#growth rate and cost of generalism, competition
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_comp
gamma <- seq(from = 0, to = parameters['gamma_sp'] * 5, by = 0.01)

cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma)){
  parameters_comp['gamma_sp']=gamma[i]
  output <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, 0.1)
  generalism_cost = gamma[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

#heatmaps
cleaned_coop <- cycled_gamma_coop %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5.5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5.5,
                           TRUE ~ log10)) 

plot_coop <- cleaned_coop %>%
  filter(rate_cost %in% c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
  geom_tile(aes(fill = log10)) +
  scale_fill_gradient2(low = "red",
                         mid = "white",
                         high= "blue",
                         midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  labs(x = "Fold-Increase in E. Coli Growth Rate", y = "Fold-Increase in Cost of Phage Generalism (Productivity)",
       fill = "Specialist Relative Fitness") +
  theme(axis.title = element_text()) + 
  ylab("Fold-Increase in E. Coli Growth Rate")+
  xlab("Fold-Increase in Specialist Productivity") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

cleaned_comp <- cycled_gamma_comp %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5.5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5.5,
                           TRUE ~ log10))

plot_comp <- cleaned_comp %>%
  filter(rate_cost %in% c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
  geom_tile(aes(fill = log10)) +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "blue",
                       midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  labs(x = "Fold-Increase in E. Coli Growth Rate", y = "Fold-Increase in Cost of Phage Generalism (Productivity)",
       fill = "Specialist Relative Fitness") +
  theme(axis.title = element_text()) + 
  ylab("Fold-Increase in E. Coli Growth Rate")+
  xlab("Fold-Increase in Specialist Productivity") +
  geom_hline(yintercept=c(5), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

patched <- plot_coop / plot_comp + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom')


