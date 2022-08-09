#ATB
#Modeling host range paper
#Fig 6 - host preference no benefit

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")

#load equations
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#generate data for productivity
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

cost = 1.5

parameters <- parameters_coop
gamma_genS <- seq(from = 0.0, to = parameters['gamma_genS'] * cost, by = 0.001)

total_fitness_gamma_coop = data.frame()
cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop['gamma_genS']=gamma_genS[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_coop,
                parms = parameters_coop,
                atol = 1e-14)
    generalism_preference = gamma_genS[i] / parameters_coop["gamma_genE"]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_coop = rbind(total_fitness_gamma_coop, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

parameters <- parameters_comp
total_fitness_gamma_comp = data.frame()
cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp['gamma_genS']=gamma_genS[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_comp,
                parms = parameters_comp,
                atol = 1e-14)
  generalism_preference = gamma_genS[i] / parameters_comp["gamma_genE"]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_comp = rbind(total_fitness_gamma_comp, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

completedatacost_gamma <- rbind(cycled_gamma_coop %>% 
                                  mutate(label = "Cooperation") %>% 
                                  inner_join(., total_fitness_gamma_coop, by = "generalism_preference") %>% 
                                  group_by(generalism_preference) %>% 
                                  filter(time == max(time)), 
                                cycled_gamma_comp %>% 
                                  mutate(label = "Competition") %>% 
                                  inner_join(., total_fitness_gamma_comp, by = "generalism_preference") %>% 
                                  group_by(generalism_preference) %>% 
                                  filter(time == max(time)))

#generate data for consumption rate
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

cost = 2

parameters <- parameters_coop
c_genS <- seq(from = 0, to = parameters['c_genS'] * cost, by = 0.001)

total_fitness_c_coop = data.frame()
cycled_c_coop <- data.frame()
for (i in 1:length(c_genS)){
  parameters_coop['c_genS']=c_genS[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_coop,
                parms = parameters_coop,
                atol = 1e-14)
  generalism_preference = c_genS[i] / parameters_coop["c_genE"]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Consumption Rate")
  total_fitness_c_coop = rbind(total_fitness_c_coop, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_c_coop = rbind(cycled_c_coop, output)
}

parameters <- parameters_comp
total_fitness_c_comp = data.frame()
cycled_c_comp <- data.frame()
for (i in 1:length(c_genS)){
  parameters_comp['c_genS']=c_genS[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_comp,
                parms = parameters_comp,
                atol = 1e-14)
  generalism_preference = c_genS[i] / parameters_comp["c_genE"]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Consumption Rate")
  total_fitness_c_comp = rbind(total_fitness_c_comp, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_c_comp = rbind(cycled_c_comp, output)
}

completedatacost_c <- rbind(cycled_c_coop %>% 
                                  mutate(label = "Cooperation") %>% 
                                  inner_join(., total_fitness_c_coop, by = "generalism_preference") %>% 
                                  group_by(generalism_preference) %>% 
                                  filter(time == max(time)), 
                                cycled_c_comp %>% 
                                  mutate(label = "Competition") %>% 
                                  inner_join(., total_fitness_c_comp, by = "generalism_preference") %>% 
                                  group_by(generalism_preference) %>% 
                                  filter(time == max(time)))


#plots 
reprorate_gamma <- completedatacost_gamma %>%
  select(-c(time, E, S, gen, sp, repro_specialist)) %>%
  melt(id = c("generalism_preference", "label", "cost_type")) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x=as.numeric(generalism_preference), y = normalized)) +
  geom_line(aes(linetype = label), size = 2, color = "black")+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12),
        legend.position = "none")+
  #scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Reproductive Rate of Generalist Phage")+
  xlab("Relative Generalist Preference for S. enterica (Efficiency)")+
  ylim(0, 1)

reprorate_c <- completedatacost_gamma %>%
  select(-c(time, E, S, gen, sp, repro_specialist)) %>%
  melt(id = c("generalism_preference", "label", "cost_type")) %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x=as.numeric(generalism_preference), y = normalized)) +
  geom_line(size = 2, color = "black",
              aes(y = jitter(normalized, 0.75), x = jitter(as.numeric(generalism_preference), 0),
                  linetype = label))+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12),
        legend.position = "none")+
  #scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Reproductive Rate of Generalist Phage")+
  xlab("Relative Generalist Preference for S. enterica (Consumption Rate)")+
  ylim(0, 1)

fig6 <- reprorate_gamma + reprorate_c + plot_annotation(tag_levels = "A")
