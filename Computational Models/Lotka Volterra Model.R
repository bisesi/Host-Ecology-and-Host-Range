#ATB
#Generalized Lotka Volterra model for 4 species system
#One generalist predation, one specialist predator

#load packages
library("deSolve")

#Parameter Desciptions
parameters = c(alpha1 = 1, #mutualism coefficient, set to 1 for competition
               alpha2 = 1, #mutualism coefficient, set to 1 for competition 
               
               beta1 = 1, #less than 1, competition coefficient of E on S, set to 1 for mutualism
               beta2 = 0.9, #less than 1, competition coefficient of E on S, set to 1 for mutualism
               
              rate_e = 0.6089485,#intrinsic rate of growth of E, 0.6089485
               rate_s = 1.043844, #instrinsic rate of growth of S, 1.043844
               
               gamma1 = 5e-5,#rate of generalist growth on S, combines adsorption/lysis/burst
               gamma2 = 5e-5,#rate of generalist growth on E, combines adsorption/lysis/burst
               gamma3 = 5e-5,#rate of specialist growth on S, combines adsorption/lysis/burst   
                    
               dilution = 3e-7, #rate of dilution (emigration/death)
               k_e = 0, # e half saturating constant, set to 0 for competition, set to 10 for mutualism
               k_s = 0, #s half saturating constant, set to 0 for competition, set to 10 for mutualism
                    
               c1 = 5e-6, #rate of specialist consumption of S, akin to adsorption (rate of finding)
               c2 = 5e-6 #rate of generalist consumption of E or S, assumed to be equivalent for both E and S

)

#Standard Parameters
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



  