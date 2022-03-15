#ATB
#Resource explicit rates using relative fitness

#load packages
library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")
library("deSolve")

#load source
source(here::here("Visualizations for Paper", "custom-theme.R"))

#Set parameters with equal growth rates
parameters_Coop = c(mu_e = 0.5,
                    k_e_lactose = 7e-7,#speed of eating
                    k_e_met = 3e-7,#formerly 3e-7
                    c_e_lactose = 2e-12,#amount used to make one cell
                    p_e_acetate = 1e-13,#amount produced by one cell
                    c_s_acetate = 3e-13,#formerly 3e-13
                    p_s_met = 10e-14, #formerly 10e-14
                    c_e_met = 2e-14,#formerly 2e-14
                    mu_s = 0.5,   
                    k_s_acetate = 3e-7, #formerly 3e-7
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 =5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

Coop <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (Es < 1e-20){
      Es = 0
    }
    if (Er < 1e-20){
      Er = 0
    }
    if (Ss < 1e-20){
      Ss = 0
    }
    if (Sr < 1e-20){
      Sr = 0
    }
    if (lcts < 1e-50){
      lcts = 0
    }
    if (ac < 1e-50){
      ac = 0
    }
    if (met < 1e-50){
      met = 0
    }
    
    gEs = Es * mu_e * (lcts / (lcts + k_e_lactose)) * (met / (met + k_e_met))
    dEs = gEs - (Es*gen*adsorp_gen) - (Es*sp2*adsorp_sp2)
    dEdead = (Es*gen*adsorp_gen) + (Es*sp2*adsorp_sp2) 
    dEr = Er * mu_e * (lcts / (lcts + k_e_lactose)) * (met / (met + k_e_met)) 
    
    gSs = Ss * mu_s * (ac / (ac + k_s_acetate))
    dSs = gSs - (Ss*sp*adsorp_sp) - (Ss*gen*adsorp_gen) 
    dSdead = (Ss*sp*adsorp_sp) + (Ss*gen*adsorp_gen) 
    dSr = Sr * mu_s * (ac / (ac + k_s_acetate)) 
    
    dlcts =  (-(gEs+dEr) * c_e_lactose) 
    dmet = (-(gEs+dEr) * c_e_met)  + ((gSs+dSr) * p_s_met) 
    dac = (-(gSs+dSr) * c_s_acetate) + ((gEs+dEr) * p_e_acetate)
    
    dgen = (gen * burst_gen * Es * adsorp_gen) + (gen * burst_gen * Ss * adsorp_gen) 
    dsp = (sp * burst_sp * Ss * adsorp_sp) 
    dsp2 = (sp2 * burst_sp2 * Es * adsorp_sp2)
    
    return(list(c(dEs, dEr, dSs, dSr, dgen, dsp, dsp2, dEdead, dSdead, dlcts, dmet, dac)))
  })
}

parameters_Comp = c(mu_s = 0.05,#1.043844*0.05 based on experimental work
                    k_e_lactose = 7e-7,
                    c_e_lactose = 2e-12,
                    c_s_lactose = 2e-12,
                    mu_e = 0.05,  #0.6089485*0.05 based on experimental work 
                    k_s_lactose = 7e-7,
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

Comp <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    if (Es < 1e-20){
      Es = 0
    }
    if (Er < 1e-20){
      Er = 0
    }
    if (Ss < 1e-20){
      Ss = 0
    }
    if (Sr < 1e-20){
      Sr = 0
    }
    if (lcts < 1e-50){
      lcts = 0
    }
    
    gEs = Es * mu_e * (lcts / (lcts + k_e_lactose)) 
    dEs = gEs - (Es*gen*adsorp_gen) - (Es*sp2*adsorp_sp2) 
    dEdead = (Es*gen*adsorp_gen) + (Es*sp2*adsorp_sp2) 
    dEr = Er * mu_e * (lcts / (lcts + k_e_lactose))
    
    gSs = Ss * mu_s * (lcts / (lcts + k_s_lactose))
    dSs = gSs - (Ss*sp*adsorp_sp) - (Ss*gen*adsorp_gen)
    dSdead = (Ss*sp*adsorp_sp) + (Ss*gen*adsorp_gen)
    dSr = Sr * mu_s * (lcts / (lcts + k_s_lactose))
    
    dlcts = (-(gEs+dEr) * c_e_lactose) - ((gSs+dSr) * c_s_lactose)
    
    dgen = (gen * burst_gen * Es * adsorp_gen) + (gen * burst_gen * Ss * adsorp_gen)
    dsp = (sp * burst_sp * Ss * adsorp_sp)
    dsp2 = (sp2 * burst_sp2 * Es * adsorp_sp2)
    
    return(list(c(dEs, dEr, dSs, dSr, dgen, dsp, dsp2, dEdead, dSdead, dlcts)))
  })
}
parameters_Neutral = c(mu_e = 0.5,#0.8035*0.5, e rate of growth in lcts+met media
                       k_e_lactose = 7e-7,
                       c_e_lactose = 2e-12,
                       c_s_acetate = 2e-12,#formerly 3e-13
                       mu_s = 0.5, #1.071839*0.5, s rate of growth in gluc media
                       k_s_acetate = 7e-7, #formerly 3e-7
                       
                       burst_gen = 5,
                       burst_sp = 5,
                       burst_sp2 = 5,
                       adsorp_gen = 1e-8,
                       adsorp_sp = 1e-8,
                       adsorp_sp2 = 1e-8
)

Neutral <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (Es < 1e-20){
      Es = 0
    }
    if (Er < 1e-20){
      Er = 0
    }
    if (Ss < 1e-20){
      Ss = 0
    }
    if (Sr < 1e-20){
      Sr = 0
    }
    if (lcts < 1e-50){
      lcts = 0
    }
    if (ac < 1e-50){
      ac = 0
    }
    
    gEs = Es * mu_e * (lcts / (lcts + k_e_lactose)) 
    dEs = gEs - (Es*gen*adsorp_gen) - (Es*sp2*adsorp_sp2)
    dEdead = (Es*gen*adsorp_gen) + (Es*sp2*adsorp_sp2)
    dEr = Er * mu_e * (lcts / (lcts + k_e_lactose))
    
    gSs = Ss * mu_s * (ac / (ac + k_s_acetate))
    dSs = gSs - (Ss*sp*adsorp_sp) - (Ss*gen*adsorp_gen)
    dSdead = (Ss*sp*adsorp_sp) + (Ss*gen*adsorp_gen)
    dSr = Sr * mu_s * (ac / (ac + k_s_acetate))
    
    dlcts = (-(gEs+dEr) * c_e_lactose)
    dac = (-(gSs+dSr) * c_s_acetate)
    
    dgen = (gen * burst_gen * Es * adsorp_gen) + (gen * burst_gen * Ss * adsorp_gen)
    dsp = (sp * burst_sp * Ss * adsorp_sp)
    dsp2 = (sp2 * burst_sp2 * Es * adsorp_sp2)
    
    return(list(c(dEs, dEr, dSs, dSr, dgen, dsp, dsp2, dEdead, dSdead, dlcts, dac)))
  })
}

#Sweep burst size parameter space for cooperation
time = seq(from = 0 , to = 1e5, by = 10)
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0)
burst_sp <- c(0:50)
burst_Coop = parameters_Coop
burst_coop = data.frame()
total_fitness_burstcoop = data.frame()
for (i in 1:length(burst_sp)){
  burst_Coop['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = burst_Coop,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(cost = cost_type)
    total_fitness_burstcoop = rbind(total_fitness_burstcoop, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
   burst_coop = rbind(burst_coop, out)
}
total_fitness_burstcoop <- total_fitness_burstcoop %>% 
  mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                       ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                              ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA))))%>%
  mutate(label = "Cooperation")%>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Burst Size")

#Sweep burst size parameter space for competition
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3)
burst_Comp = parameters_Comp
burst_comp = data.frame()
total_fitness_burstcomp = data.frame()
for (i in 1:length(burst_sp)){
  burst_Comp['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = burst_Comp,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
   mutate(cost = cost_type)
  total_fitness_burstcomp = rbind(total_fitness_burstcomp, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
  burst_comp = rbind(burst_comp, out)
}
total_fitness_burstcomp <- total_fitness_burstcomp %>% 
  mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                       ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                              ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA))))%>%
  mutate(label = "Competition")%>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Burst Size")

#Sweep burst size parameter space for no interactions
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
burst_Neutral = parameters_Neutral
burst_neutral = data.frame()
total_fitness_burstneutral = data.frame()
for (i in 1:length(burst_sp)){
  burst_Neutral['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = burst_Neutral,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(cost = cost_type)
  total_fitness_burstneutral = rbind(total_fitness_burstneutral, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
  burst_neutral  = rbind(burst_neutral, out)
}
total_fitness_burstneutral <- total_fitness_burstneutral %>% 
  mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                       ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                              ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(label = "Neutral")%>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Burst Size")

#Sweep adsorption parameter space for cooperation
adsorp_sp <- burst_sp * 2e-9
adsorp_Coop = parameters_Coop #rename parameters or can't reuse after loop
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0) 
adsorp_coop = data.frame()
total_fitness_adsorpcoop = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Coop['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = adsorp_Coop,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(cost = cost_type)
  total_fitness_adsorpcoop = rbind(total_fitness_adsorpcoop, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
  adsorp_coop = rbind(adsorp_coop, out)
}
total_fitness_adsorpcoop <- total_fitness_adsorpcoop %>% 
  mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                       ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                              ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(label = "Cooperation")%>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Adsorption Rate")

#Sweep adsorption parameter space for competition
adsorp_Comp = parameters_Comp
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3) 
adsorp_comp = data.frame()
total_fitness_adsorpcomp = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Comp['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = adsorp_Comp,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(cost = cost_type)
  total_fitness_adsorpcomp = rbind(total_fitness_adsorpcomp, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
  adsorp_comp = rbind(adsorp_comp, out)
}
total_fitness_adsorpcomp <- total_fitness_adsorpcomp %>% 
  mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                       ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                              ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(label = "Competition") %>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Adsorption Rate")

#Sweep adsorption parameter space for no interactions
adsorp_Neutral = parameters_Neutral
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
adsorp_neutral = data.frame()
adsorp_neutral_fitness = data.frame()
total_fitness_adsorpneutral = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Neutral['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = adsorp_Neutral,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(cost = cost_type)
  total_fitness_adsorpneutral = rbind(total_fitness_adsorpneutral, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
  adsorp_neutral = rbind(adsorp_neutral, out)
}
total_fitness_adsorpneutral <- total_fitness_adsorpneutral %>% 
        mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                                        ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                               ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
        mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                                      ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                             ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
  mutate(label = "Neutral") %>%
  mutate(basefitness = repro_specialist / repro_generalist) %>%
  mutate(type = "Adsorption Rate")

#Combine data for visualization
all <- rbind(total_fitness_burstcomp, total_fitness_burstcoop, total_fitness_burstneutral,
             total_fitness_adsorpcomp, total_fitness_adsorpcoop, total_fitness_adsorpneutral)

#Visualizations

#Calculate relative fitness based on whichever has highest reproductive rate at time = 1e5
datasetswitch <- all %>%
  select(-c(repro_generalist, repro_specialist, basefitness)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Genotype = variable)
  
relativefitnessplot <- datasetswitch %>%
  ggplot(aes(x=cost, y = value, color = Genotype))+
  geom_line(size = 2)+
  facet_grid(type ~ label)+
  theme_bisesi()+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Relative Fitness (w)")+
  xlab("Fold-Increase in Specialist Infectivity")


#Calculate relative fitness as reproductive rate specialist / reproductive rate generalist
datasetbasegen <- all %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  select(-c(variable)) %>%
  filter(label == "Competition" | label == "Cooperation")

relativefitnessbasegen <- datasetbasegen %>%
  ggplot(aes(x=cost, y = value))+
  geom_line(size = 2)+
  facet_grid(type ~ label)+
  theme_bisesi()+
  geom_hline(yintercept=c(1), color = "red", linetype = "dashed", size = 0.5)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Relative Fitness of Specialist")+
  xlab("Fold-Increase in Specialist Infectivity")

#Final figures
fig3 <- (relativefitnessplot / relativefitnessbasegen) +
  plot_annotation(tag_levels = "A")
