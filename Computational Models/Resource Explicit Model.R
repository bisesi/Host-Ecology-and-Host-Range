#ATB
#Four species resource explicit model

#load packages
library("deSolve")

#Standard parameters
parameters_Coop = c(mu_e = 0.5,#rate of growth
                    k_e_lactose = 7e-7,#speed of eating
                    k_e_met = 3e-7,#formerly 3e-7
                    c_e_lactose = 2e-12,#amount used to make one cell
                    p_e_acetate = 1e-13,#amount produced by one cell
                    c_s_acetate = 3e-13,
                    p_s_met = 10e-14, 
                    c_e_met = 2e-14,
                    mu_s = 0.5,#rate of growth
                    k_s_acetate = 3e-7,
                    
                    burst_gen = 5,#burst size of the generalist
                    burst_sp = 5,
                    burst_sp2 =5,
                    adsorp_gen = 1e-8,#adsorption rate of the generalist
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Comp = c(mu_s = 0.05,#1.043844*0.05 based on experimental work, 0.0521922
                    k_e_lactose = 7e-7,
                    c_e_lactose = 2e-12,
                    c_s_lactose = 2e-12,
                    mu_e = 0.05,  #0.6089485*0.05 based on experimental work, 0.03044743
                    k_s_lactose = 7e-7,
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Facilitation = c(mu_e = 0.05,#rate of growth
                    k_e_glucose = 7e-7,#speed of eating
                    k_e_met = 3e-7,#formerly 3e-7
                    c_e_glucose = 2e-12,#amount used to make one cell
                    p_s_met = 10e-14, 
                    c_e_met = 2e-14,
                    mu_s = 0.05,#rate of growth
                    k_s_glucose = 3e-7,
                    c_s_glucose = 3e-13,
                    
                    burst_gen = 5,#burst size of the generalist
                    burst_sp = 5,
                    burst_sp2 =5,
                    adsorp_gen = 1e-8,#adsorption rate of the generalist
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Neutral = c(mu_e = 0.5,#0.8035*0.5, e rate of growth in lcts+met media
                       k_e_lactose = 7e-7,
                       c_e_lactose = 2e-12,
                       c_s_acetate = 2e-12,#formerly 3e-13
                       mu_s = 0.5, #1.071839*0.5, s rate of growth in gluc media
                       k_s_acetate = 7e-7, #formerly 3e-7
                       
                       burst_gen = 5,
                       burst_sp = 10,
                       burst_sp2 = 5,
                       adsorp_gen = 1e-8,
                       adsorp_sp = 1e-8,
                       adsorp_sp2 = 1e-8
)

#Mathematical Models
Coop <- function(t,n,parms){#host cooperation (obligate mutualism)
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

Comp <- function(t,n,parms){#host competition
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

Facilitation <- function(t,n,parms){#host cooperation (obligate mutualism)
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
    if (glucose < 1e-50){
      glucose = 0
    }
    if (met < 1e-50){
      met = 0
    }
    
    gEs = Es * mu_e * (glucose / (glucose + k_e_glucose)) * (met / (met + k_e_met))
    dEs = gEs - (Es*gen*adsorp_gen) - (Es*sp2*adsorp_sp2)
    dEdead = (Es*gen*adsorp_gen) + (Es*sp2*adsorp_sp2) 
    dEr = Er * mu_e * (glucose/ (glucose + k_e_glucose)) * (met / (met + k_e_met)) 
    
    gSs = Ss * mu_s * (glucose / (glucose  + k_s_glucose))
    dSs = gSs - (Ss*sp*adsorp_sp) - (Ss*gen*adsorp_gen)
    dSdead = (Ss*sp*adsorp_sp) + (Ss*gen*adsorp_gen)
    dSr = Sr * mu_s * (glucose  / (glucose  + k_s_glucose))
    
    dglucose =  (-(gSs+dSr) * c_s_glucose) 
    dmet = (-(gEs+dEr) * c_e_met)  + ((gSs+dSr) * p_s_met) 
    
    dgen = (gen * burst_gen * Es * adsorp_gen) + (gen * burst_gen * Ss * adsorp_gen) 
    dsp = (sp * burst_sp * Ss * adsorp_sp) 
    dsp2 = (sp2 * burst_sp2 * Es * adsorp_sp2)
    
    return(list(c(dEs, dEr, dSs, dSr, dgen, dsp, dsp2, dEdead, dSdead, dglucose, dmet)))
  })
}

Neutral <- function(t,n,parms){#no host interaction
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

#Function to cycle through changing cost of generalism with relative fitness
generalismcost <- function(model, parameters, modelname, startingdensity, timerange, maxcost){
  time = timerange
  start_density = startingdensity
  burst_range <- c(0:maxcost)
  parameters_burst = parameters
  parameters_adsorp = parameters
  burst = data.frame()
  total_fitness_burst = data.frame()
  for (i in 1:length(burst_range)){
    parameters_burst['burst_sp']=burst_range[i]
    out=ode(y=start_density,
            times=time,
            func=model,
            parms = parameters_burst,
            atol = 1e-14)
    cost_type = burst_range[i] / 5
    out=data.frame(out) %>%
      mutate(cost = cost_type)
    total_fitness_burst = rbind(total_fitness_burst, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
    burst = rbind(burst, out)
  }
  total_fitness_burst <- total_fitness_burst %>% 
    mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
    mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA))))%>%
    mutate(label = modelname)%>%
    mutate(basefitness = repro_specialist / repro_generalist) %>%
    mutate(type = "Burst Size")
  
  adsorp_range <- c(0:maxcost) * 2e-9
  start_density = startingdensity
  adsorp = data.frame()
  total_fitness_adsorp = data.frame()
  for (i in 1:length(adsorp_range)){
    parameters_adsorp['adsorp_sp']=adsorp_range[i]
    out=ode(y=start_density,
            times=time,
            func=model,
            parms = parameters_adsorp,
            atol = 1e-14)
    cost_type = adsorp_range[i] / 1e-8
    out=data.frame(out) %>%
      mutate(cost = cost_type)
    total_fitness_adsorp = rbind(total_fitness_adsorp, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
    adsorp = rbind(adsorp, out)
  }
  total_fitness_adsorp <- total_fitness_adsorp %>% 
    mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
    mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
    mutate(label = modelname)%>%
    mutate(basefitness = repro_specialist / repro_generalist) %>%
    mutate(type = "Adsorption Rate")
  
  alldata <- rbind(total_fitness_adsorp, total_fitness_burst)
  full <- list(alldata, burst, adsorp)
  return(full)
}

#Function to cycle through changing growth rates with relative fitness
growthrate <- function(model, parameters, modelname, startingdensity, timerange, maxcost, increment){
  time = timerange
  start_density = startingdensity
  e_growthrange <- seq(from = 0, to = parameters['mu_s'] * maxcost, by = increment)
  parameters_e = parameters
  parameters_s = parameters
  highe = data.frame()
  total_fitness_e = data.frame()
  for (i in 1:length(e_growthrange)){
    parameters_e['mu_e']=e_growthrange[i]
    out=ode(y=start_density,
            times=time,
            func=model,
            parms = parameters_e,
            atol = 1e-14)
    rate_cost = e_growthrange[i] / parameters_e['mu_s']
    out=data.frame(out) %>%
      mutate(rate_cost = rate_cost)
    total_fitness_e = rbind(total_fitness_e, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), rate_cost = rate_cost))
    highe = rbind(highe, out)
  }
  total_fitness_e <- total_fitness_e %>% 
    mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
    mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA))))%>%
    mutate(label = modelname)%>%
    mutate(basefitness = repro_specialist / repro_generalist) %>%
    mutate(type = "Increasing E coli Growth Rate")
  
  return(total_fitness_e)
}
