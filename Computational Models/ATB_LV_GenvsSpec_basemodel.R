#ATB
#Generalized Lotka Volterra model for 4 species system
#One generalist predator, one specialist predator
#base models

#load packages
library("deSolve")

#Cooperation parameters
parameters_coop <- c(
  #mutualism coefficients
  alpha1 = 1, 
  alpha2 = 1, 
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 0.5,
  
  #phage productivity constants
  gamma_genE = 2e-2,
  gamma_genS = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_genE = 1e-3,
  c_genS = 1e-3,
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
  gamma_genE = 2e-2,
  gamma_genS = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_genE = 1e-3,
  c_genS = 1e-3,
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
    dE = rate_e * E * (alpha1*S/(alpha1*S + k_e)) * (1-E) - c_genE*gen*E - dilution*E
    dS = rate_s * S * (alpha2*E/(alpha2*E + k_s)) * (1-S) - c_sp*sp*S - c_genS*gen*S - dilution*S
    
    dgen = gamma_genS*gen*S + gamma_genE*gen*E - dilution*gen
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
    dE = rate_e * E * (2-E-(beta1*S)) - c_genE*gen*E - dilution*E
    dS = rate_s * S * (2-S-(beta2*E)) - c_sp*sp*S - c_genS*gen*S - dilution*S
    
    dgen = gamma_genS*gen*S + gamma_genE*gen*E - dilution*gen
    dsp = gamma_sp*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

