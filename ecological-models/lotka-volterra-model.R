#ATB
#Generalized Lotka Volterra model for 4 species system
#One generalist predator, one specialist predator

#load packages
library("deSolve")

#Cooperation parameters
parameters_coop_R <- c(
  #mutualism coefficients
  alpha1 = 1, 
  alpha2 = 1, 
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 0.5,
  
  #phage productivity constants
  gamma_genE = 20,
  gamma_genS = 20,
  gamma_sp = 20,
  
  #phage rate of consumption constants
  c_genE = 1e-3,
  c_genS = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution_E = 3e-2,
  dilution_S = 3e-2,
  dilution_gen = 3e-2,
  dilution_sp = 3e-2, #rate of dilution (emigration/death) was 3e-2
  k_e = 1, #half-saturation
  k_s = 1, #half-saturation
  R = 1
)

#Competition parameters
parameters_comp_R <- c(
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 0.5,
  
  #phage productivity constants
  gamma_genE = 20,
  gamma_genS = 20,
  gamma_sp = 20,
  
  #phage rate of consumption constants
  c_genE = 1e-3,
  c_genS = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution_E = 3e-2,
  dilution_S = 3e-2,
  dilution_gen = 3e-2,
  dilution_sp = 3e-2,
  R = 2 #rate of dilution (emigration/death) was 3e-2
)


#Mathematical models
generalLV_coop_R <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
    
    #Coop    
    dE = rate_e * E * (alpha1*S/(alpha1*S + k_e)) * (R-E) - c_genE*gen*E - dilution_E*E
    dS = rate_s * S * (alpha2*E/(alpha2*E + k_s)) * (R-S) - c_sp*sp*S - c_genS*gen*S - dilution_S*S
    
    dgen = gamma_genS*c_genS*gen*S + gamma_genE*c_genE*gen*E - dilution_gen*gen
    dsp = gamma_sp*c_sp*sp*S - dilution_sp*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

generalLV_comp_R <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
    
    #Comp
    dE = rate_e * E * (R-E-(beta1*S)) - c_genE*gen*E - dilution_E*E
    dS = rate_s * S * (R-S-(beta2*E)) - c_sp*sp*S - c_genS*gen*S - dilution_S*S
    
    dgen = gamma_genS*c_genS*gen*S + gamma_genE*c_genE*gen*E - dilution_gen*gen
    dsp = gamma_sp*c_sp*sp*S - dilution_sp*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

  