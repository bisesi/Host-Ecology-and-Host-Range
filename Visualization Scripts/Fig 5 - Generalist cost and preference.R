#ATB
#Modeling host range paper
#Fig 5 - increasing specialist productivity and preference for E

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")

#load equations
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#growth rate and cost of generalism, cooperation and productivity
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

cost = 2.5

parameters <- parameters_coop
gamma_genS <- seq(from = 0, to = parameters['gamma_genS'] * cost, by = 0.001)

cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop['gamma_genS']=gamma_genS[i]
  output <- LVgeneralismcost(generalLV_coop, parameters, "Cooperation", start_density, time, 5, "gamma")
  generalism_preference = gamma_genS[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

#growth rate and cost of generalism, competition
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_comp
gamma_genS <- seq(from = 0, to = parameters['gamma_genS'] * cost, by = 0.001)

cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma)){
  parameters_comp['gamma_genS']=gamma_genS[i]
  output <- LVgeneralismcost(generalLV_comp, parameters, "Competition", start_density, time, 5, "gamma")
  generalism_preference = gamma_genS[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

#heatmaps
cleaned_coop <- cycled_gamma_coop %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("cost", "generalism_preference", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

plot_coop <- cleaned_coop %>%
  ggplot(aes(x = as.character(cost), y = as.character(round(generalism_preference,1)))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=24,face="bold"),
        axis.title.y  = element_text(size=24,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title=element_text(size=20))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Relative Generalist Conversion Eficiency on S. enterica")+
  xlab("Benefit of Specialism") +
  geom_hline(yintercept=c(6), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

cleaned_comp <- cycled_gamma_comp %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(normalized = exp(log2(value)) / (1 + exp(log2(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

plot_comp <- cleaned_comp %>%
  ggplot(aes(x = as.character(cost), y = as.character(round(generalism_preference,1)))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=24,face="bold"),
        axis.title.y  = element_text(size=24,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title=element_text(size=18))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Relative Generalist Conversion Eficiency on S. enterica")+
  xlab("Benefit of Specialism")+
  geom_hline(yintercept=c(6), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

fig5 <- plot_coop + plot_comp + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

