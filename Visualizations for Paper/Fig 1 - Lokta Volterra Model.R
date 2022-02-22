#ATB 
#Host ecology and host range paper
#Figure 1, Lotka Volterra Model

#load dependencies
library("ggplot2")
library("reshape")
library("tidyverse")
library("ggpubr")
library("deSolve")
library("gridExtra")
library("ggthemes")

source("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Visualizations for Paper/custom-theme.R")

#Cooperative parameters
parameters_coop = c(alpha1 = 5, 
                    alpha2 = 1.1, 
                    
                    beta1 = 1, 
                    beta2 = 1, 
                    
                    rate_e = 1.5,#intrinsic rate of growth of E
                    rate_s = 1.5, #intrinsic rate of growth of S
                    
                    gamma1 = 5e-5,
                    gamma2 = 5e-5,
                    gamma3 = 5e-5,
                    
                    dilution = 3e-7, #rate of dilution (emigration/death)
                    k_e = 10, 
                    k_s = 10, 
                    
                    c1 = 5e-6, 
                    c2 = 5e-6 
                    
)

#Competitive parameters
parameters_comp = c(alpha1 = 1, 
                    alpha2 = 1, 
                    
                    beta1 = 1,
                    beta2 = 0.9, 
                    
                    rate_e = 0.6089485,#intrinsic rate of growth of E
                    rate_s = 1.043844, #instrinsic rate of growth of S
                    
                    gamma1 = 5e-5,
                    gamma2 = 5e-5,
                    gamma3 = 5e-5,
                    
                    dilution = 3e-7, #rate of dilution (emigration/death)
                    k_e = 0, 
                    k_s = 0, 
                    
                    c1 = 5e-6, 
                    c2 = 5e-6 
                    
)

#Mathematical Model
generalLV <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    dE = rate_e * E * ((alpha1*S)/((alpha1*S) + k_e)) * (1-E-(beta1*S)) - c2*gen*E - dilution*E
    dS = rate_s * S * ((alpha2*E)/((alpha2*E) + k_s)) * (1-S-(beta2*E)) - c1*sp*S - c2*gen*S - dilution*S
    dgen = gamma1*gen*S + gamma2*gen*E - dilution*gen
    dsp = gamma3*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

#Get parameter Space
parm_range <- seq(from = 0.0001 , to = 15, by = 0.1)
time = seq(from = 0 , to = 5e5, by = 10)

#Cooperative loop - productivity
gamma3 <- 5e-5 * parm_range
gamma3_Coop = parameters_coop
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
gamma_coop = data.frame()
max_gamma_coop = data.frame()
for (i in 1:length(gamma3)){
  gamma3_Coop['gamma3']=gamma3[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV,
          parms = gamma3_Coop,
          rtol = 1e-14)
  cost_type = gamma3[i] / 5e-5
  out=data.frame(out) %>%
    mutate(biomass = E + S) %>%
    mutate(ratio = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(label = "Productivity") %>%
    mutate(type = "Cooperation")
  gamma_coop = rbind(gamma_coop, out)
  max_gamma_coop = rbind(max_gamma_coop, gamma_coop[gamma_coop$time == max(gamma_coop$time),])
  
}
#Competitive loop - productivity
gamma3 <- 5e-5 * parm_range
gamma3_Comp = parameters_comp
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)
gamma_comp = data.frame()
max_gamma_comp = data.frame()
for (i in 1:length(gamma3)){
  gamma3_Comp['gamma3']=gamma3[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV,
          parms = gamma3_Comp,
          rtol = 1e-14)
  cost_type = gamma3[i] / 5e-5
  out=data.frame(out) %>%
    mutate(biomass = E + S) %>%
    mutate(ratio = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(label = "Productivity") %>%
    mutate(type = "Competition")
  gamma_comp = rbind(gamma_comp, out)
  max_gamma_comp = rbind(max_gamma_comp, gamma_comp[gamma_comp$time == max(gamma_comp$time),])
}

#Productivity dataset
max_productivity <- rbind(max_gamma_coop %>% select(c(gen, sp, cost, label, type)), 
                          max_gamma_comp %>% select(c(gen, sp, cost, label, type)))
melted_productivity <- max_productivity %>%
  melt(id = c("label", "cost", "type")) %>%
  mutate(variable = case_when(variable == "gen" ~ "Generalist",
                           variable == "sp" ~ "Specialist")) %>%
  rename(c(Genotype = variable))

#Cooperative loop - rate of consumption
c1 <- 5e-6 * parm_range
c1_Coop = parameters_coop
c1_coop = data.frame()
max_c1_coop = data.frame()
for (i in 1:length(c1)){
  c1_Coop['c1']=c1[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV,
          parms = c1_Coop,
          rtol = 1e-14)
  cost_type = c1[i] / 5e-6
  out=data.frame(out) %>%
    mutate(biomass = E + S) %>%
    mutate(ratio = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(label = "Rate of Consumption") %>%
    mutate(type = "Cooperation")
  c1_coop = rbind(c1_coop, out)
  max_c1_coop = rbind(max_c1_coop, c1_coop[c1_coop$time == max(c1_coop$time),])
}

#Competitive loop - rate of consumption
c1 <- 5e-6 * parm_range
c1_Comp = parameters_comp
c1_comp = data.frame()
max_c1_comp = data.frame()
for (i in 1:length(c1)){
  c1_Comp['c1']=c1[i]
  out=ode(y=start_density,
          times=time,
          func=generalLV,
          parms = c1_Comp,
          rtol = 1e-14)
  cost_type = c1[i] / 5e-6
  out=data.frame(out) %>%
    mutate(biomass = E + S) %>%
    mutate(ratio = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(label = "Rate of Consumption") %>%
    mutate(type = "Competition")
  c1_comp = rbind(c1_comp, out)
  max_c1_comp = rbind(max_c1_comp, c1_comp[c1_comp$time == max(c1_comp$time),])
}

#Rate of consumption dataset
max_c1 <- rbind(max_c1_coop %>% select(c(gen, sp, cost, label, type)), 
                max_c1_comp %>% select(c(gen, sp, cost, label, type)))
melted_c1 <- max_c1 %>%
  melt(id = c("label", "cost", "type")) %>%
  mutate(variable = case_when(variable == "gen" ~ "Generalist",
                              variable == "sp" ~ "Specialist")) %>%
  rename(c(Genotype = variable))


#Plots
LVdataset <- rbind(melted_productivity, melted_c1)
plot_LV <- LVdataset %>%
  ggplot(aes(x=cost, y = log10(value), color = Genotype))+
  geom_line(size = 2)+
  facet_grid(label ~ type)+
  theme_bisesi()+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Equilibrium Phage Density (Log10)")+
  xlab("Fold-Increase in Specialist Infectivity")
