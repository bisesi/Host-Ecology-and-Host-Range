#ATB
#Modeling host range paper
#Fig 3 - increasing specialist productivity 

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

parameters_coop["gamma_sp"] = 3.5e-2
parameters_coop["gamma_genS"] = 3.5e-2
parameters_coop["gamma_genE"] = 3.5e-2

parameters_comp["gamma_sp"] = 3.5e-2
parameters_comp["gamma_genS"] = 3.5e-2
parameters_comp["gamma_genE"] = 3.5e-2

costofgeneralismcoop_gamma <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, "gamma")
costofgeneralismcomp_gamma <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, "gamma")

completedatacost_gamma <- rbind(costofgeneralismcoop_gamma, costofgeneralismcomp_gamma)

datasetbasegencost_gamma <- completedatacost_gamma %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

#plots 
relativefitnessbasegencost_gamma <- datasetbasegencost_gamma %>%
  mutate(normalized = exp(log10(value)) / (1 + exp(log10(value)))) %>%
  mutate(normalized = case_when(is.nan(normalized) == TRUE ~ 1,
                                TRUE ~ normalized)) %>%
  ggplot(aes(x=as.numeric(cost), y = normalized, group = Interaction)) +
  geom_line(aes(linetype = Interaction), size = 1, color = "black")+
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
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12),
        legend.position = "none")+
  #scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist Phage")+
  xlab("Benefit of Specialism (Predator Efficiency)")+
  ylim(0, 1)


fig3 <- relativefitnessbasegencost_gamma 
