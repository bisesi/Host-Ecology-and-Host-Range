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
                    mu_s = 0.2,#rate of growth
                    k_s_acetate = 3e-7,
                    
                    burst_gen = 5,#burst size of the generalist
                    burst_sp = 5,
                    burst_sp2 =5,
                    adsorp_gen = 1e-8,#adsorption rate of the generalist
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Comp = c(mu_s = 0.0521922,#1.043844*0.05 based on experimental work
                    k_e_lactose = 7e-7,
                    c_e_lactose = 2e-12,
                    c_s_lactose = 2e-12,
                    mu_e = 0.03044743,  #0.6089485*0.05 based on experimental work 
                    k_s_lactose = 7e-7,
                    
                    burst_gen = 5,
                    burst_sp = 5,
                    burst_sp2 = 5,
                    adsorp_gen = 1e-8,
                    adsorp_sp = 1e-8,
                    adsorp_sp2 = 1e-8
)

parameters_Neutral = c(mu_e = 0.40175,#0.8035*0.5, e rate of growth in lcts+met media
                       k_e_lactose = 7e-7,
                       c_e_lactose = 2e-12,
                       c_s_acetate = 2e-12,#formerly 3e-13
                       mu_s = 0.5359195, #1.071839*0.5, s rate of growth in gluc media
                       k_s_acetate = 7e-7, #formerly 3e-7
                       
                       burst_gen = 5,
                       burst_sp = 5,
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