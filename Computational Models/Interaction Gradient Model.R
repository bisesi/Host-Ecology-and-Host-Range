#ATB
#Generalized Lotka Volterra model for 4 species system
#One generalist predator, one specialist predator
#Competition gradient 

#load packages
library("deSolve")

#Cooperation parameters
parameters <- c(
  #mutualism coefficients
  alpha1 = 10, 
  alpha2 = 5, 
  beta1 = 1,
  beta2 = 1,
  
  #interaction, large if competition, small if cooperation
  a = 1,
  R = 1, 
  
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
  dilution = 3e-2, #rate of dilution (emigration/death) was 3e-2
  k_e = 1, #half-saturation
  k_s = 1 #half-saturation
)

#Mathematical models
generalLV_gradient <- function(t,n,parms){
  with(as.list(c(t,n,parms)), {
    
    if (E < 1e-10){
      E = 0
    }
    if (S < 1e-10){
      S = 0
    }
      
    dE = rate_e * E * ((alpha1*S + a)/(alpha1*S + a + k_e)) * (R - E - beta1*S) - c_genE*gen*E - dilution*E
    dS = rate_s * S * ((alpha2*E + a)/(alpha2*E + a + k_s)) * (R - S - beta2*E) - c_sp*sp*S - c_genS*gen*S - dilution*S
    
    dgen = gamma_genS*c_genS*gen*S + gamma_genE*c_genE*gen*E - dilution*gen
    dsp = gamma_sp*c_sp*sp*S - dilution*sp
    
    return(list(c(dE, dS, dgen, dsp)))
  })
}

time = seq(from = 0 , to = 2e5, by = 10)

start_density <- c(E = 0, S = 0.1, gen = 0, sp = 0)

out = ode(y=start_density,
          times=time,
          func=generalLV_gradient,
          parms = parameters,
          atol = 1e-14)
