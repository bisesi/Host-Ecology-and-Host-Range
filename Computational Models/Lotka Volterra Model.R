#ATB
#Generalized Lotka Volterra model for 4 species system
#One generalist predation, one specialist predator

#load packages
library("deSolve")

#Cooperation parameters
parameters_coop <- c(
  #mutualism coefficients
  alpha1 = 1, 
  alpha2 = 1, 
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = .5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2, #rate of dilution (emigration/death) was 3e-2
  k_e = 1, #half-saturation
  k_s = 1 #half-saturation
)

#Competition parameters
parameters_comp <- c(
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 0.5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2 #rate of dilution (emigration/death) was 3e-2
)

#No interaction parameters
parameters_none <- c(
  #intrinsic growth rate
  rate_e = .5,
  rate_s = .5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2 #rate of dilution (emigration/death) was 3e-2
)

#Mathematical models
generalLV_coop <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
    
    #Coop    
    dE = rate_e * E * (alpha1*S/(alpha1*S + k_e)) * (1-E) - c_gen*gen*E - dilution*E
    dS = rate_s * S * (alpha2*E/(alpha2*E + k_s)) * (1-S) - c_sp*sp*S - c_gen*gen*S - dilution*S
    
    dgen = gamma_gen*gen*S + gamma_gen*gen*E - dilution*gen
    dsp = gamma_sp*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

generalLV_comp <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
    
    #Comp
    dE = rate_e * E * (2-E-(beta1*S)) - c_gen*gen*E - dilution*E
    dS = rate_s * S * (2-S-(beta2*E)) - c_sp*sp*S - c_gen*gen*S - dilution*S
    
    dgen = gamma_gen*gen*S + gamma_gen*gen*E - dilution*gen
    dsp = gamma_sp*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

generalLV_none <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
    
    #Independent
    dE = rate_e * E * (1.8-E) - c_gen*gen*E - dilution*E
    dS = rate_s * S * (1.8-S) - c_sp*sp*S - c_gen*gen*S - dilution*S
    
    dgen = gamma_gen*gen*S + gamma_gen*gen*E - dilution*gen
    dsp = gamma_sp*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

#Function to cycle through changing cost of generalism with relative fitness
LVgeneralismcost <- function(model, parameters, modelname, startingdensity, timerange, maxcost){
  time = timerange
  start_density = startingdensity
  if (modelname == "Cooperation"){
    gamma_range <- seq(from = 0, to = parameters['gamma_sp'] * maxcost, by = 0.01)
    gamma_range <- gamma_range[gamma_range < 0.19]
  }
  else {
    gamma_range <- seq(from = 0, to = parameters['gamma_sp'] * maxcost, by = 0.01)
  }
  parameters_gamma = parameters
  parameters_c = parameters
  gamma = data.frame()
  total_fitness_gamma = data.frame()
  for (i in 1:length(gamma_range)){
    parameters_gamma['gamma_sp']=gamma_range[i]
    out=ode(y=start_density,
            times=time,
            func=model,
            parms = parameters_gamma,
            atol = 1e-14)
    cost_type = gamma_range[i] / 2e-2
    out=data.frame(out) %>%
      mutate(cost = cost_type)
    total_fitness_gamma = rbind(total_fitness_gamma, cbind(repro_generalist = ((max(out$gen) - out$gen[1]) / out$gen[1]), repro_specialist = ((max(out$sp) - out$sp[1]) / out$sp[1]), cost = cost_type))
    gamma = rbind(gamma, out)
  }
  total_fitness_gamma <- total_fitness_gamma %>% 
    mutate(generalist = ifelse(repro_generalist > repro_specialist, repro_generalist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_generalist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA)))) %>%
    mutate(specialist = ifelse(repro_generalist > repro_specialist, repro_specialist / repro_generalist, 
                               ifelse(repro_specialist > repro_generalist, repro_specialist / repro_specialist, 
                                      ifelse(repro_specialist == repro_generalist, 1, NA))))%>%
    mutate(label = modelname)%>%
    mutate(basefitness = repro_specialist / repro_generalist) %>%
    mutate(type = "Gamma")
  return(total_fitness_gamma)
}

#Function to cycle through changing growth rates with relative fitness
LVgrowthrate <- function(model, parameters, modelname, startingdensity, timerange, maxcost, increment){
  time = timerange
  start_density = startingdensity
  if (modelname == "Cooperation"){
    e_growthrange <- seq(from = 0.3, to = parameters['rate_s'] * maxcost, by = increment)
  }
  else {
    e_growthrange <- seq(from = 0, to = parameters['rate_s'] * maxcost, by = increment)
  }
  parameters_e = parameters
  parameters_s = parameters
  highe = data.frame()
  total_fitness_e = data.frame()
  for (i in 1:length(e_growthrange)){
    parameters_e['rate_e']=e_growthrange[i]
    out=ode(y=start_density,
            times=time,
            func=model,
            parms = parameters_e,
            atol = 1e-14)
    rate_cost = e_growthrange[i] / parameters_e['rate_s']
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

  