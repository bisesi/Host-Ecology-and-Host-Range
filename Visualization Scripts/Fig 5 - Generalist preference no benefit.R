#ATB
#Modeling host range paper
#Fig 5 - host preference no benefit

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")

#load equations
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#generate loses attachment ability
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

gamma_genS = c(3.5e-2,3.5e-2,3.5e-2,1.75e-2,0)
gamma_genE = c(0,1.75e-2,3.5e-2,3.5e-2,3.5e-2)
c_genE = c(0, 0.001, 0.001, 0.001, 0.001)
c_genS = c(0.001, 0.001, 0.001, 0.001, 0)
labels = c("S specialist", "S preference", "True generalist", "E preference", "E specialist")

total_fitness_gamma_coop = data.frame()
cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop["gamma_sp"] = 3.5e-2
  parameters_coop['gamma_genS']=gamma_genS[i]
  parameters_coop['gamma_genE']=gamma_genE[i]
  parameters_coop['c_genS']=c_genS[i]
  parameters_coop['c_genE']=c_genE[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_coop,
                parms = parameters_coop,
                atol = 1e-14)
    generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_coop = rbind(total_fitness_gamma_coop, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

total_fitness_gamma_comp = data.frame()
cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp["gamma_sp"] = 3.5e-2
  parameters_comp['gamma_genS']=gamma_genS[i]
  parameters_comp['gamma_genE']=gamma_genE[i]
  parameters_comp['c_genS']=c_genS[i]
  parameters_comp['c_genE']=c_genE[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_comp,
                parms = parameters_comp,
                atol = 1e-14)
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_comp = rbind(total_fitness_gamma_comp, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

completedatacost_noattach <- rbind(total_fitness_gamma_coop %>%
                                  mutate(Interaction = "Cooperation") %>%
                                  mutate(relative_fitness = as.numeric(repro_specialist) / as.numeric(repro_generalist)),
                                total_fitness_gamma_comp %>%
                                  mutate(Interaction = "Competition") %>%
                                  mutate(relative_fitness = as.numeric(repro_specialist) / as.numeric(repro_generalist)))

#with ability to attach
#generate loses attachment ability
source(here::here("Computational Models", "Lotka Volterra Model.R"))
total_fitness_gamma_coop_attach = data.frame()
cycled_gamma_coop_attach <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_coop["gamma_sp"] = 3.5e-2
  parameters_coop['gamma_genS']=gamma_genS[i]
  parameters_coop['gamma_genE']=gamma_genE[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_coop,
                parms = parameters_coop,
                atol = 1e-14)
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_coop_attach = rbind(total_fitness_gamma_coop_attach, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_coop_attach = rbind(cycled_gamma_coop_attach, output)
}

total_fitness_gamma_comp_attach = data.frame()
cycled_gamma_comp_attach <- data.frame()
for (i in 1:length(gamma_genS)){
  parameters_comp["gamma_sp"] = 3.5e-2
  parameters_comp['gamma_genS']=gamma_genS[i]
  parameters_comp['gamma_genE']=gamma_genE[i]
  output <- ode(y=start_density,
                times=time,
                func=generalLV_comp,
                parms = parameters_comp,
                atol = 1e-14)
  generalism_preference = labels[i]
  output=data.frame(output) %>%
    mutate(generalism_preference = generalism_preference) %>%
    mutate(cost_type = "Productivity (Gamma)")
  total_fitness_gamma_comp_attach = rbind(total_fitness_gamma_comp_attach, cbind(repro_generalist = ((max(output$gen) - output$gen[1]) / output$gen[1]), repro_specialist = ((max(output$sp) - output$sp[1]) / output$sp[1]), generalism_preference = generalism_preference))
  cycled_gamma_comp_attach = rbind(cycled_gamma_comp_attach, output)
}

completedatacost_withattach <- rbind(total_fitness_gamma_coop_attach %>%
                                     mutate(Interaction = "Cooperation") %>%
                                     mutate(relative_fitness = as.numeric(repro_specialist) / as.numeric(repro_generalist)),
                                   total_fitness_gamma_comp_attach %>%
                                     mutate(Interaction = "Competition") %>%
                                     mutate(relative_fitness = as.numeric(repro_specialist) / as.numeric(repro_generalist)))


#plots 
reprorate_gamma <- completedatacost_withattach %>%
  mutate(type = "Attach") %>%
  rbind(., completedatacost_noattach %>% mutate(type = "No Attach")) %>%
  mutate(normalized = exp(log10(relative_fitness)) / (1 + exp(log10(relative_fitness)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized))

reprorate_gamma$generalism_preference <- factor(reprorate_gamma$generalism_preference, levels = labels)

fig5 <- reprorate_gamma %>%
  filter(type == "Attach") %>%
  ggplot(aes(x=generalism_preference, y = normalized, group = Interaction)) +
  geom_point(size = 3, aes(shape = Interaction))+
  theme_bw()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=14,face="bold"),
        axis.title.y  = element_text(size=14,face="bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.title=element_text(size=12),
        legend.position = "none")+
  #scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Relative Fitness of Specialist Phage")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)+
  xlab("Generalist Preference")+
  ylim(0, 1)
