#ATB
#Fig 1
#Null isoclines

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")
library("phaseR")

#ZNGI without phage
model.coop <- function(t, y, parameters){
  with(as.list(parameters),{
    
    E<-y[1] 
    S<-y[2]
    dE = rate_e * E * (S/(S + k_e)) * (1-E) - dilution*E
    dS = rate_s * S * (E/(E + k_s)) * (1-S) - dilution*S
    
    list(c(dE,dS))
  })
}
params.coop<-c(rate_e = 0.5, k_e = 1, dilution = 3e-2, rate_s = 0.5, k_s = 1)

model.comp <- function(t, y, parameters){
  with(as.list(parameters),{
    
    E<-y[1] 
    S<-y[2]
    dE = rate_e * E * (2-E-(beta1*S)) -  dilution*E
    dS = rate_s * S * (2-S-(beta2*E)) - dilution*S
    
    list(c(dE,dS))
  })
}
params.comp<-c(rate_e = 0.5, beta1 = 0.9, dilution = 3e-2, rate_s = 0.5, beta2 = 0.9)

data.LV.coop <-as.data.frame(lsoda(c(E=0.1,S=0.1),seq(1,250,by=0.5), model.coop, params.coop))
data.LV.comp <-as.data.frame(lsoda(c(E=0.1,S=0.1),seq(1,250,by=0.5), model.comp, params.comp))

# plot the trajectories of the system
par(mfrow =c(1,2), mar = c(2.5,2.5,2.5,2.5))
plot(NULL, xlim=c(0,1),ylim=c(0,1), type = "l",xlab = "E. coli Density", ylab = "S. enterica Density", main = "Cooperation")
nullclines(model.coop, xlim=c(0,1),ylim=c(0,1), parameters=params.coop, system="two.dim",  col = c("black", "black"), add=TRUE, add.legend=FALSE)
text(x = 0.833, y = 0.275,                # Add text element
     "dE/dt")
text(x = 0.4, y = 0.875,                # Add text element
     "dS/dt")
plot(NULL, xlim=c(0,2),ylim=c(0,2), type = "l",xlab = "E. coli Density", ylab = "S. enterica Density", main = "Competition")
nullclines(model.comp, xlim=c(0,2),ylim=c(0,2), parameters=params.comp, system="two.dim",  col = c("black", "black"), add=TRUE, add.legend=FALSE)
text(x = 1.65, y = 0.02,                # Add text element
     "dE/dt")
text(x = 0.25, y = 1.5,                # Add text element
     "dS/dt")