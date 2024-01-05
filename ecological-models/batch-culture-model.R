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
  k_e = 1, #half-saturation
  k_s = 1, #half-saturation
  R = 1
)

#Competition parameters
parameters_comp_R <- c(
  #competition coefficients
  beta1 = 1,
  beta2 = 1,
  
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
    dE = rate_e * E * (alpha1*S/(alpha1*S + k_e)) * (R-E) - c_genE*gen*E
    dS = rate_s * S * (alpha2*E/(alpha2*E + k_s)) * (R-S) - c_sp*sp*S - c_genS*gen*S 
    
    dgen = gamma_genS*c_genS*gen*S + gamma_genE*c_genE*gen*E
    dsp = gamma_sp*c_sp*sp*S
    
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
    dE = rate_e * E * (R-E-(beta1*S)) - c_genE*gen*E 
    dS = rate_s * S * (R-S-(beta2*E)) - c_sp*sp*S - c_genS*gen*S
    
    dgen = gamma_genS*c_genS*gen*S + gamma_genE*c_genE*gen*E 
    dsp = gamma_sp*c_sp*sp*S
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

