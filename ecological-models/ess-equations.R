#ATB
#Based on eneralized Lotka Volterra model for 4 species system
#One generalist predator, one specialist predator
#ESS equations generated in Mathematica

#Resident type = specialist, mutant type = generalist
coop_sp_to_gen <- function(R = 2, I_E_gen = 0.02, dilution = 3e-2, I_S_gen = 0.02, I_S_sp = 0.02, mu1 = 0.5, k1 = 1) {
  return((R*I_E_gen) + dilution*(-1 + (I_S_gen / I_S_sp) - (I_E_gen / (mu1 + k1*mu1))))
}

comp_sp_to_gen <- function(R = 2, I_E_gen = 0.02, beta1 = 0.9, mu1 = 0.5, I_S_gen = 0.02, I_S_sp = 0.02, dilution = 3e-2) {
  return((R*I_E_gen) + (((beta1*mu1*I_E_gen - mu1*I_S_gen + (mu1 + I_E_gen)*I_S_sp)*dilution)/(mu1*I_S_sp)))
}

#Resident type = generalist, mutant type = specialist
coop_gen_to_sp <- function(R = 2, dilution = 3e-2, gamma_sp = 20, gamma_genE = 20, gamma_genS = 20, c_sp = 1e-3, c_genE = 1e-3, c_genS = 1e-3, mu1 = 0.5, mu2 = 0.5, k1 = 1, k2 = 1){
  num = -dilution + gamma_sp*c_sp*((1 + k2)*mu2*R*gamma_genE*c_genE^2 - (1 + k1)*mu1*c_genS*(R*gamma_genE*c_genE - dilution) + gamma_genE*c_genE*(-c_genE + c_genS)*dilution)
  denom = (1 + k2)*mu2*gamma_genE*c_genE^2 + (1 + k1)*mu1*gamma_genS*c_genS^2
  return(num / denom)
} 

comp_gen_to_sp <- function(R = 2, dilution = 3e-2, gamma_sp = 20, gamma_genE = 20, gamma_genS = 20, c_sp = 1e-3, c_genE = 1e-3, c_genS = 1e-3, mu1 = 0.5, mu2 = 0.5, beta1 = 0.9, beta2 = 0.9){
  num = -dilution + gamma_sp*c_sp*(gamma_genE*c_genE*(-c_genE + c_genS)*dilution + mu1*c_genS*(-R*c_genE*gamma_genE + dilution) + mu2*c_genE*(R*gamma_genE*c_genE - beta2*dilution))
  denom = mu1*c_genE*(-beta1*gamma_genE*c_genE + gamma_genS*c_genS) + mu2*c_genE*(gamma_genE*c_genE - beta2*c_genS*gamma_genS)
  return(num / denom)
} 