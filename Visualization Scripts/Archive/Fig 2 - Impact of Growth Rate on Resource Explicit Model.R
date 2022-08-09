#ATB
#Resource Explicit with growth rate comparisons

#load packages
library("ggplot2")
library("reshape")
library("tidyverse")
library("deSolve")
library("patchwork")

#load source
source(here::here("Visualizations for Paper", "custom-theme.R"))

#Mathematical models
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

#Parameters when specialist host is faster
parameters_Coop = c(mu_e = 0.2,
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
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Comp = c(mu_s = 0.05,#1.043844*0.05 based on experimental work
                    k_e_lactose = 7e-7,
                    c_e_lactose = 2e-12,
                    c_s_lactose = 2e-12,
                    mu_e = 0.02,  #0.6089485*0.05 based on experimental work 
                    k_s_lactose = 7e-7,
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Neutral = c(mu_e = 0.2,#0.8035*0.5, e rate of growth in lcts+met media
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

#Sweep burst size parameter space - cooperation
time = seq(from = 0 , to = 1e5, by = 10)
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0)
burst_sp <- c(5:50) #set range of burst size parameters 
burst_Coop = parameters_Coop
burst_coop = data.frame()
max_burst_coop = data.frame()
for (i in 1:length(burst_sp)){
  burst_Coop['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = burst_Coop,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_coop = rbind(burst_coop, out)
  max_burst_coop = rbind(max_burst_coop, burst_coop[burst_coop$time == max(burst_coop$time),])
}

#Sweep burst size parameter space - competition
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3)
burst_Comp = parameters_Comp
burst_comp = data.frame()
max_burst_comp = data.frame()
for (i in 1:length(burst_sp)){
  burst_Comp['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = burst_Comp,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Competition")%>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_comp = rbind(burst_comp, out)
  max_burst_comp = rbind(max_burst_comp, burst_comp[burst_comp$time == max(burst_comp$time),])
}

#Sweep burst size parameter space - no interactions
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
burst_Neutral = parameters_Neutral
burst_neutral = data.frame()
max_burst_neutral = data.frame()
for (i in 1:length(burst_sp)){
  burst_Neutral['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = burst_Neutral,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_neutral  = rbind(burst_neutral, out)
  max_burst_neutral = rbind(max_burst_neutral, burst_neutral[burst_neutral$time == max(burst_neutral$time),])
}

#Sweep adsorption parameter space - cooperation
adsorp_sp <- 2e-9 * burst_sp
adsorp_Coop = parameters_Coop #rename parameters or can't reuse after loop
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0) 
adsorp_coop = data.frame()
max_adsorp_coop = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Coop['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = adsorp_Coop,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_coop = rbind(adsorp_coop, out)
  max_adsorp_coop= rbind(max_adsorp_coop, adsorp_coop[adsorp_coop$time == max(adsorp_coop$time),])
}

#Sweep adsorption parameter space - competition
adsorp_Comp = parameters_Comp
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3) 
adsorp_comp = data.frame()
max_adsorp_comp = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Comp['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = adsorp_Comp,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Competition") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_comp = rbind(adsorp_comp, out)
  max_adsorp_comp= rbind(max_adsorp_comp, adsorp_comp[adsorp_comp$time == max(adsorp_comp$time),])
}

#Sweep adsorption parameter space - no interactions
adsorp_Neutral = parameters_Neutral
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
adsorp_neutral = data.frame()
max_adsorp_neutral = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Neutral['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = adsorp_Neutral,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_neutral = rbind(adsorp_neutral, out)
  max_adsorp_neutral= rbind(max_adsorp_neutral, adsorp_neutral[adsorp_neutral$time == max(adsorp_neutral$time),])
}

max_fastS <- rbind(max_burst_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_neutral %>% select(gen, sp, cost, type, label, percentchange),
                   max_adsorp_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_neutral %>% select(gen, sp, cost, type, label, percentchange))

#Parameters when specialist host is slower
parameters_Coop = c(mu_e = 0.5,
                    k_e_lactose = 7e-7,#speed of eating
                    k_e_met = 3e-7,#formerly 3e-7
                    c_e_lactose = 2e-12,#amount used to make one cell
                    p_e_acetate = 1e-13,#amount produced by one cell
                    c_s_acetate = 3e-13,#formerly 3e-13
                    p_s_met = 10e-14, #formerly 10e-14
                    c_e_met = 2e-14,#formerly 2e-14
                    mu_s = 0.2,   
                    k_s_acetate = 3e-7, #formerly 3e-7
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Comp = c(mu_s = 0.02,#1.043844*0.05 based on experimental work
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

parameters_Neutral = c(mu_e = 0.5 ,#0.8035*0.5, e rate of growth in lcts+met media
                       k_e_lactose = 7e-7,
                       c_e_lactose = 2e-12,
                       c_s_acetate = 2e-12,#formerly 3e-13
                       mu_s = 0.2, #1.071839*0.5, s rate of growth in gluc media
                       k_s_acetate = 7e-7, #formerly 3e-7
                       
                       burst_gen = 5,
                       burst_sp = 5,
                       burst_sp2 = 5,
                       adsorp_gen = 1e-8,
                       adsorp_sp = 1e-8,
                       adsorp_sp2 = 1e-8
)

#Sweep burst size parameter space - cooperation
time = seq(from = 0 , to = 1e5, by = 10)
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0)
burst_sp <- c(5:50) #set range of burst size parameters 
burst_Coop = parameters_Coop
burst_coop = data.frame()
max_burst_coop = data.frame()
for (i in 1:length(burst_sp)){
  burst_Coop['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = burst_Coop,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_coop = rbind(burst_coop, out)
  max_burst_coop = rbind(max_burst_coop, burst_coop[burst_coop$time == max(burst_coop$time),])
}

#Sweep burst size parameter space - competition
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3)
burst_Comp = parameters_Comp
burst_comp = data.frame()
max_burst_comp = data.frame()
for (i in 1:length(burst_sp)){
  burst_Comp['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = burst_Comp,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Competition")%>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_comp = rbind(burst_comp, out)
  max_burst_comp = rbind(max_burst_comp, burst_comp[burst_comp$time == max(burst_comp$time),])
}

#Sweep burst size parameter space - no interactions
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
burst_Neutral = parameters_Neutral
burst_neutral = data.frame()
max_burst_neutral = data.frame()
for (i in 1:length(burst_sp)){
  burst_Neutral['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = burst_Neutral,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_neutral  = rbind(burst_neutral, out)
  max_burst_neutral = rbind(max_burst_neutral, burst_neutral[burst_neutral$time == max(burst_neutral$time),])
}

#Sweep adsorption parameter space - cooperation
adsorp_sp <- 2e-9 * burst_sp
adsorp_Coop = parameters_Coop #rename parameters or can't reuse after loop
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0) 
adsorp_coop = data.frame()
max_adsorp_coop = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Coop['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = adsorp_Coop,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_coop = rbind(adsorp_coop, out)
  max_adsorp_coop= rbind(max_adsorp_coop, adsorp_coop[adsorp_coop$time == max(adsorp_coop$time),])
}

#Sweep adsorption parameter space - competition
adsorp_Comp = parameters_Comp
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3) 
adsorp_comp = data.frame()
max_adsorp_comp = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Comp['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = adsorp_Comp,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Competition") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_comp = rbind(adsorp_comp, out)
  max_adsorp_comp= rbind(max_adsorp_comp, adsorp_comp[adsorp_comp$time == max(adsorp_comp$time),])
}

#Sweep adsorption parameter space - no interactions
adsorp_Neutral = parameters_Neutral
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
adsorp_neutral = data.frame()
max_adsorp_neutral = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Neutral['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = adsorp_Neutral,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_neutral = rbind(adsorp_neutral, out)
  max_adsorp_neutral= rbind(max_adsorp_neutral, adsorp_neutral[adsorp_neutral$time == max(adsorp_neutral$time),])
}

max_slowS <- rbind(max_burst_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_neutral %>% select(gen, sp, cost, type, label, percentchange),
                   max_adsorp_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_neutral %>% select(gen, sp, cost, type, label, percentchange))

#Parameters when specialist host and non-specialist host have same growth rate
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
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

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

parameters_Neutral = c(mu_e = 0.5 ,#0.8035*0.5, e rate of growth in lcts+met media
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

#Sweep burst size parameter space - cooperation
time = seq(from = 0 , to = 1e5, by = 10)
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0)
burst_sp <- c(5:50) #set range of burst size parameters 
burst_Coop = parameters_Coop
burst_coop = data.frame()
max_burst_coop = data.frame()
for (i in 1:length(burst_sp)){
  burst_Coop['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = burst_Coop,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_coop = rbind(burst_coop, out)
  max_burst_coop = rbind(max_burst_coop, burst_coop[burst_coop$time == max(burst_coop$time),])
}

#Sweep burst size parameter space - competition
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3)
burst_Comp = parameters_Comp
burst_comp = data.frame()
max_burst_comp = data.frame()
for (i in 1:length(burst_sp)){
  burst_Comp['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = burst_Comp,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "Competition")%>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_comp = rbind(burst_comp, out)
  max_burst_comp = rbind(max_burst_comp, burst_comp[burst_comp$time == max(burst_comp$time),])
}

#Sweep burst size parameter space - no interactions
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
burst_Neutral = parameters_Neutral
burst_neutral = data.frame()
max_burst_neutral = data.frame()
for (i in 1:length(burst_sp)){
  burst_Neutral['burst_sp']=burst_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = burst_Neutral,
          rtol = 1e-14)
  cost_type = burst_sp[i] / 5
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Burst Size") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  burst_neutral  = rbind(burst_neutral, out)
  max_burst_neutral = rbind(max_burst_neutral, burst_neutral[burst_neutral$time == max(burst_neutral$time),])
}

#Sweep adsorption parameter space - cooperation
adsorp_sp <- 2e-9 * burst_sp
adsorp_Coop = parameters_Coop #rename parameters or can't reuse after loop
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, met = 1e-12, ac = 0) 
adsorp_coop = data.frame()
max_adsorp_coop = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Coop['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Coop,
          parms = adsorp_Coop,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Cooperation") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_coop = rbind(adsorp_coop, out)
  max_adsorp_coop= rbind(max_adsorp_coop, adsorp_coop[adsorp_coop$time == max(adsorp_coop$time),])
}

#Sweep adsorption parameter space - competition
adsorp_Comp = parameters_Comp
start_density = c(Es = 1e4, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3) 
adsorp_comp = data.frame()
max_adsorp_comp = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Comp['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Comp,
          parms = adsorp_Comp,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "Competition") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_comp = rbind(adsorp_comp, out)
  max_adsorp_comp= rbind(max_adsorp_comp, adsorp_comp[adsorp_comp$time == max(adsorp_comp$time),])
}

#Sweep adsorption parameter space - no interactions
adsorp_Neutral = parameters_Neutral
start_density = c(Es = 0, Er = 0, Ss = 1e4, Sr = 0, gen = 1, sp = 1, sp2 = 0, Edead = 0, Sdead = 0,
                  lcts = 1e-3, ac = 1e-3)
adsorp_neutral = data.frame()
max_adsorp_neutral = data.frame()
for (i in 1:length(adsorp_sp)){
  adsorp_Neutral['adsorp_sp']=adsorp_sp[i]
  out=ode(y=start_density,
          times=time,
          func=Neutral,
          parms = adsorp_Neutral,
          rtol = 1e-14)
  cost_type = adsorp_sp[i] / 1e-8
  out=data.frame(out) %>%
    mutate(biomass = Es + Er + Ss + Sr) %>%
    mutate(phageratio_sp = gen / (sp+gen)) %>%
    mutate(cost = cost_type) %>%
    mutate(type = "Adsorption Rate") %>%
    mutate(label = "S0") %>%
    mutate(percentchange = phageratio_sp - 0.5)
  adsorp_neutral = rbind(adsorp_neutral, out)
  max_adsorp_neutral= rbind(max_adsorp_neutral, adsorp_neutral[adsorp_neutral$time == max(adsorp_neutral$time),])
}

max_equalES <- rbind(max_burst_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_burst_neutral %>% select(gen, sp, cost, type, label, percentchange),
                   max_adsorp_coop %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_comp %>% select(gen, sp, cost, type, label, percentchange), 
                   max_adsorp_neutral %>% select(gen, sp, cost, type, label, percentchange))

#Combine and clean datasets
max_slowS <- max_slowS %>% mutate(rate = "Slow Salmonella")
max_fastS <- max_fastS %>% mutate(rate = "Fast Salmonella")
max_equalES <- max_fastS %>% mutate(rate = "Equal Growth Rates")

all <- rbind(max_slowS, max_fastS, max_equalES)

meltedraw <- all %>%
  select(c(gen, sp, cost, type, label, rate)) %>%
  melt(id = c("label", "cost", "type", "rate")) %>%
  mutate(variable = case_when(variable == "gen" ~ "Generalist",
                              variable == "sp" ~ "Specialist")) %>%
  rename(c(Genotype = variable, Condition = label))

meltedpercentchange <- all %>%
  select(c(cost, type, label, percentchange, rate)) %>%
  melt(id = c("label", "cost", "type", "rate")) %>%
  rename(c(Condition = label))

#Visualize data
growthrateplotrawburst <- meltedraw %>%
  filter(type == "Burst Size" & Condition != "S0") %>%
  ggplot(aes(x=cost, y = log10(value), color = Genotype))+
  geom_line(size = 2)+
  facet_grid(Condition ~ rate)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  theme_bisesi()+
  ylab("Equilibrium Phage Density (Log10)")+
  xlab("Fold-Increase in Specialist Burst Size")

growthrateplotrawadsorp <- meltedraw %>%
  filter(type == "Adsorption Rate" & Condition != "S0") %>%
  ggplot(aes(x=cost, y = log10(value), color = Genotype))+
  geom_line(size = 2)+
  facet_grid(Condition ~ rate)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  theme_bisesi()+
  ylab("Equilibrium Phage Density (Log10)")+
  xlab("Fold-Increase in Specialist Adsorption Rate")

growthrateplotpercentchange <- meltedpercentchange %>%
  filter(Condition != "S0") %>%
  ggplot(aes(x=cost, y = value, color = Condition))+
  geom_line(size = 2)+
  facet_grid(type ~ rate)+
  theme_bisesi()+
  scale_color_manual(values=c("#000000", "#E66100"))+
  ylab("Net Change in Percent Generalists")+
  xlab("Fold-Increase in Specialist Infectivity")+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)

#Final figure
fig2 <- (growthrateplotpercentchange | (growthrateplotrawadsorp / growthrateplotrawburst + plot_layout(guides = "collect")))+ 
  plot_annotation(tag_levels= "A")

