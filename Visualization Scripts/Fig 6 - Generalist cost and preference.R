#ATB
#Modeling host range paper
#Fig 6 - increasing specialist productivity and preference for E

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

gamma_genS = c(3.5e-2,3.5e-2,3.5e-2,1.75e-2,0)
gamma_genE = c(0,1.75e-2,3.5e-2,3.5e-2,3.5e-2)
c_genE = c(0, 0.001, 0.001, 0.001, 0.001)
c_genS = c(0.001, 0.001, 0.001, 0.001, 0)
labels = c("S specialist", "S preference", "True generalist", "E preference", "E specialist")

cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop["gamma_sp"] = 3.5e-2
  parameters_coop['gamma_genS']=gamma_genS[i]
  parameters_coop['gamma_genE']=gamma_genE[i]
  parameters_coop['c_genS']=c_genS[i]
  parameters_coop['c_genE']=c_genE[i]
  output <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, "gamma")
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

#growth rate and generalism benefit, competition
cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp["gamma_sp"] = 3.5e-2
  parameters_comp['gamma_genS']=gamma_genS[i]
  parameters_comp['gamma_genE']=gamma_genE[i]
  parameters_comp['c_genS']=c_genS[i]
  parameters_comp['c_genE']=c_genE[i]
  output <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, "gamma", 6)
  generalism_preference = labels[i]
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

cleaned_coop$generalism_preference <- factor(cleaned_coop$generalism_preference, levels = labels)

plot_coop <- cleaned_coop %>%
  ggplot(aes(x = generalism_preference, y = as.character(cost))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=14))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Benefit of Specialism")+
  xlab("Generalist Preference") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  coord_flip()


cleaned_comp <- cycled_gamma_comp %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("cost", "generalism_preference", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(normalized = exp(log2(value)) / (1 + exp(log2(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))
cleaned_comp$generalism_preference <- factor(cleaned_comp$generalism_preference, levels = labels)

plot_comp <- cleaned_comp %>%
  ggplot(aes(x = generalism_preference, y = as.character(cost))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=14))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Benefit of Specialism")+
  xlab("Generalist Preference") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  coord_flip()

fig6_noattach <- plot_coop + plot_comp + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

#reload
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#growth rate and cost of generalism, cooperation and productivity
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

gamma_genS = c(3.5e-2,3.5e-2,3.5e-2,1.75e-2,0)
gamma_genE = c(0,1.75e-2,3.5e-2,3.5e-2,3.5e-2)
labels = c("S specialist", "S preference", "True generalist", "E preference", "E specialist")

cycled_gamma_coop_with <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop["gamma_sp"] = 3.5e-2
  parameters_coop['gamma_genS']=gamma_genS[i]
  parameters_coop['gamma_genE']=gamma_genE[i]
  output <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, "gamma")
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_coop_with = rbind(cycled_gamma_coop_with, output)
}

#growth rate and generalism benefit, competition
cycled_gamma_comp_with <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp["gamma_sp"] = 3.5e-2
  parameters_comp['gamma_genS']=gamma_genS[i]
  parameters_comp['gamma_genE']=gamma_genE[i]
  output <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, "gamma", 6)
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_comp_with = rbind(cycled_gamma_comp_with, output)
}

#heatmaps
cleaned_coop_with <- cycled_gamma_coop_with %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("cost", "generalism_preference", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

cleaned_coop_with$generalism_preference <- factor(cleaned_coop_with$generalism_preference, levels = labels)

plot_coop_with <- cleaned_coop_with %>%
  ggplot(aes(x = generalism_preference, y = as.character(cost))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=14))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Benefit of Specialism")+
  xlab("Generalist Preference") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  coord_flip()


cleaned_comp_with <- cycled_gamma_comp_with %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("cost", "generalism_preference", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(normalized = exp(log2(value)) / (1 + exp(log2(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))
cleaned_comp_with$generalism_preference <- factor(cleaned_comp_with$generalism_preference, levels = labels)

plot_comp_with <- cleaned_comp_with %>%
  ggplot(aes(x = generalism_preference, y = as.character(cost))) +
  geom_tile(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0.5,
                       limits = c(0,1))+
  theme_bw() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=14))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Benefit of Specialism")+
  xlab("Generalist Preference") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  coord_flip()

fig6_withattach <- plot_coop_with + plot_comp_with + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

