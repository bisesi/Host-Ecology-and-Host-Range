#ATB
#Modeling host range paper
#Fig 4b - heat map changing host growth rate and consumption rate

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

cost = 5

parameters <- parameters_coop
c <- seq(from = 0, to = parameters['c_sp'] * cost, by = 0.0001)

cycled_c_coop <- data.frame()
for (i in 1:length(c)){
  parameters_coop['c_sp']=c[i]
  output <- LVgrowthrate(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, 0.1)
  generalism_cost = c[i] / 1e-3
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Consumption Rate")
  cycled_c_coop = rbind(cycled_c_coop, output)
}

#growth rate and cost of generalism, competition
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_comp
c <- seq(from = 0, to = parameters['c_sp'] * cost, by = 0.0001)

cycled_c_comp <- data.frame()
for (i in 1:length(c)){
  parameters_comp['c_sp']=c[i]
  output <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, 0.1)
  generalism_cost = c[i] / 1e-3
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Consumption Rate")
  cycled_c_comp = rbind(cycled_c_comp, output)
}

#heatmaps
cleaned_coop <- cycled_c_coop %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

plot_coop <- cleaned_coop %>%
  add_row(rate_cost = rep(seq(0, 0.4, by = 0.2), length(seq(0, 5, by = 0.5))), generalism_cost = rep(seq(0, 5, by = 0.5), length(seq(0, 0.4, by = 0.2))),
          label = rep("Cooperation", length(rep(seq(0, 0.4, by = 0.2), length(seq(0, 5, by = 0.5))))),
          type = rep("Increasing E coli Growth Rate", length(rep(seq(0, 0.4, by = 0.2), length(seq(0, 5, by = 0.5))))),
          value = rep(NA, length(rep(seq(0, 0.4, by = 0.2), length(seq(0, 5, by = 0.5)))))) %>%
  filter(as.character(rate_cost) %in% c("0", "0.2", "0.4", "0.6", "0.8", "1", "2", "3", "4", "5")) %>%
  filter(as.character(generalism_cost) %in% c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5")) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
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
  ylab("Relative E. coli Growth Rate")+
  xlab("Benefit of Specialism") +
  geom_hline(yintercept=c(6), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

cleaned_comp <- cycled_c_comp %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

plot_comp <- cleaned_comp %>%
  filter(as.character(rate_cost) %in% c("0", "0.2", "0.4", "0.6", "0.8", "1", "2", "3", "4", "5")) %>%
  filter(as.character(generalism_cost) %in% c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5")) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
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
  ylab("Relative E. coli Growth Rate")+
  xlab("Benefit of Specialism")+
  geom_hline(yintercept=c(6), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

fig4b <- plot_coop + plot_comp + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))


