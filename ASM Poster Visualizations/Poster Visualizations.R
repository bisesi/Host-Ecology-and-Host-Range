#ASM Poster Visualizations
#ATV
#May 30, 2022

#load packages and data
library("tidyverse")
library("deSolve")
library("reshape2")
library("patchwork")
library("ggthemes")
library("readxl")
library("rstatix")
library("ggpubr")
source(here::here("Computational Models", "Lotka Volterra Model.R"))

#################################################################################

#intro graphs
#no phage and isoclines, two-D
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
par(mfrow = c(1,4), mar = c(4.5,4.5,4.5,4.5))
plot(data.LV.coop$time,data.LV.coop$E, main = "E", xlab="Time", xlim = c(0,150), ylab="Abundance",type="l", cex.lab = 1.5)
mtext("A", side=2, line=1, at=1, col="black", las=1.15, cex = 1.5)
plot(data.LV.coop$time,data.LV.coop$S, main = "S", xlab="Time", xlim = c(0,150),ylab="",type="l", cex.lab = 1.5)
plot(data.LV.comp$time,data.LV.comp$E, main = "E", xlab="Time", xlim = c(0,150),ylab="Abundance",type="l", cex.lab = 1.5)
mtext("B", side=2, line=1, at=1.25, col="black", las=1, cex = 1.5)
plot(data.LV.comp$time,data.LV.comp$S, main = "S", xlab="Time", xlim = c(0,150),ylab="",type="l", cex.lab = 1.5)

#################################################################################

#cost of generalism with equivalent parameters
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

costofgeneralismcoop <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5)
costofgeneralismcomp <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 5)

completedatacost <- rbind(costofgeneralismcoop, costofgeneralismcomp)

datasetbasegencost <- completedatacost %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

relativefitnessbasegencost <- datasetbasegencost %>%
  #rbind(., c(10.0, "Cooperation", "Gamma", Inf)) %>%
  filter(type != "Consumption Rate") %>%
  mutate(log10 = log10(as.numeric(value))) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5,
                           TRUE ~ log10)) %>%
  ggplot(aes(x=as.numeric(cost), y = log10, color = Interaction)) +
  #geom_line()+
  geom_smooth(se=FALSE, size = 2)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist Phage")+
  xlab("Benefit of Specialism")+
  ylim(-5, 5)

#relative rate with equivalent parameters
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

ratecoop <- LVgrowthrate(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, 0.1)
ratecomp <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, 0.1)

ratecompletedatacost <- rbind(ratecoop, ratecomp)

datasetbasegencostrate <- ratecompletedatacost %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("rate_cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

relativefitnessbasegencostrate <- datasetbasegencostrate %>%
  #rbind(., c(10.0, "Cooperation", "Gamma", Inf)) %>%
  filter(type != "Consumption Rate") %>%
  mutate(log10 = log10(as.numeric(value))) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5,
                           TRUE ~ log10)) %>%
  ggplot(aes(x=as.numeric(rate_cost), y = log10, color = Interaction)) +
  geom_smooth(se = FALSE, size = 2)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist Phage")+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  xlab("Relative Growth Rate of Non-Specialist Host")+
  ylim(-5, 2)

part2 <- relativefitnessbasegencost / relativefitnessbasegencostrate + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

#################################################################################
#with phage

time = seq(from = 0 , to = 1e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

out_coop <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_coop,
                              parms = parameters_coop,
                              atol = 1e-14)) %>%
  mutate(Condition = "Cooperation")

out_comp <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_comp,
                              parms = parameters_comp,
                              atol = 1e-14))%>%
  mutate(Condition = "Competition")

par(mfrow = c(2,2), mar = c(4.5,4.5,4.5,4.5))
plot(out_coop$time,out_coop$E, main = "E", xlab="Time", xlim = c(0, 20000), ylim = c(0, 1.5), ylab="Abundance",type="l", cex.lab = 1.5)
plot(out_coop$time,out_coop$S, main = "S", xlab="Time", xlim = c(0, 20000), ylim = c(0, 1.5),ylab="",type="l", cex.lab = 1.5)
plot(out_coop$time,out_coop$gen, main = "Generalist", xlim = c(0, 20000),xlab="Time",  ylab="Abundance",type="l", cex.lab = 1.5)
plot(out_coop$time,out_coop$sp, main = "Specialist",xlim = c(0, 20000), xlab="Time",ylab="",type="l", cex.lab = 1.5)

par(mfrow = c(2,2), mar = c(4.5,4.5,4.5,4.5))
plot(out_comp$time,out_comp$E, main = "E", xlab="Time", xlim = c(0, 20000),ylim = c(0, 1.5),ylab="Abundance",type="l", cex.lab = 1.5)
plot(out_comp$time,out_comp$S, main = "S", xlab="Time", xlim = c(0, 20000),ylim = c(0, 1.5),ylab="",type="l", cex.lab = 1.5)
plot(out_comp$time,out_comp$gen, main = "Generalist", xlim = c(0, 20000),xlab="Time", ylab="Abundance",type="l", cex.lab = 1.5)
plot(out_comp$time,out_comp$sp, main = "Specialist", xlim = c(0, 20000),xlab="Time", ylab="",type="l", cex.lab = 1.5)

#################################################################################

#cost of generalism with in vitro parameters and bar chart
parameters_coop_invitro <- c(
  #mutualism coefficients
  alpha1 = 10, 
  alpha2 = 1, 
  
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 25,
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
parameters_comp_invitro <- c(
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 2.5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 2e-2,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2 #rate of dilution (emigration/death) was 3e-2
)
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

costofgeneralismcoop <- LVgeneralismcost(generalLV_coop, parameters_coop_invitro, "Cooperation", start_density, time, 5)
costofgeneralismcomp <- LVgeneralismcost(generalLV_comp, parameters_comp_invitro, "Competition", start_density, time, 5)

completedatacost <- rbind(costofgeneralismcoop, costofgeneralismcomp)

datasetbasegencost <- completedatacost %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

relativefitnessbasegencost <- datasetbasegencost %>%
  #rbind(., c(10.0, "Cooperation", "Gamma", Inf)) %>%
  filter(type != "Consumption Rate") %>%
  mutate(log10 = log10(as.numeric(value))) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5,
                           TRUE ~ log10)) %>%
  ggplot(aes(x=as.numeric(cost), y = log10, color = Interaction)) +
  #geom_line(size = 2)+
  geom_smooth(se = FALSE, size = 2)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist Phage")+
  xlab("Benefit of Specialism")+
  ylim(-5, 5)

#bar chart
parameters_coop_invitro <- c(
  #mutualism coefficients
  alpha1 = 10, 
  alpha2 = 1, 
  
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 25,
  rate_s = .5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 0.03,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2, #rate of dilution (emigration/death) was 3e-2
  k_e = 1, #half-saturation
  k_s = 1 #half-saturation
)

#Competition parameters
parameters_comp_invitro <- c(
  #competition coefficients
  beta1 = .9,
  beta2 = .9,
  
  #intrinsic growth rate
  rate_e = 0.5,
  rate_s = 2.5,
  
  #phage productivity constants
  gamma_gen = 2e-2,
  gamma_sp = 0.03,
  
  #phage rate of consumption constants
  c_gen = 1e-3,
  c_sp = 1e-3,
  
  #additional constants
  dilution = 3e-2 #rate of dilution (emigration/death) was 3e-2
)

time = seq(from = 0 , to = 5000, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

out_coop <- as.data.frame(ode(y=start_density,
                                     times=time,
                                     func=generalLV_coop,
                                     parms = parameters_coop_invitro,
                                     atol = 1e-14)) %>%
  mutate(Condition = "Cooperation")

out_comp <- as.data.frame(ode(y=start_density,
                                     times=time,
                                     func=generalLV_comp,
                                     parms = parameters_comp_invitro,
                                     atol = 1e-14))%>%
  mutate(Condition = "Competition")

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

out_none <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_comp,
                              parms = parameters_comp_invitro,
                              atol = 1e-14))%>%
  mutate(Condition = "S Monoculture")

alldata <- rbind(tail(out_coop, n = 1), tail(out_comp, n = 1), tail(out_none, n = 1)) %>%
  mutate(PercentGen = gen / (gen + sp)) %>%
  mutate(Change = PercentGen - 0.5)

barchartmodel <- alldata %>%
  ggplot(aes(x=Condition, y = Change)) +
  geom_bar(stat = "identity", position = position_dodge())+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 1)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12)) + 
  ylab("Change in Percent Generalist Phage")+
  xlab("Condition")+
  ylim(-0.5, 0.5)

patchedformodel <- relativefitnessbasegencost / barchartmodel 

#time series
timeseries <- rbind(out_coop, out_comp) %>%
  select(time, gen, sp, Condition) %>%
  pivot_longer(cols = c("gen", "sp")) %>%
  mutate(name = case_when(name == "gen" ~ "Generalist",
                          name == "sp" ~ "Specialist",
                          TRUE ~ name)) %>%
  rename(Phage = name) %>%
  filter(time <= 5000) %>%
  ggplot(aes(x = time, y = value, color = Phage)) +
  geom_line(size = 2)+
  facet_wrap(~Condition)+
  geom_line()+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12)) + 
  ylab("Abundance")+
  xlab("Time")

par(mfrow=c(2,2))
plot(testdf$time, testdf$E, type = "l", cex.lab = 1.5, main = "E. coli", col = "#E1BE6A", ylim = c(0,2), lwd = 2, ylab = "Abundance", xlab = "Time")
lines(test2df$time, test2df$E, type = "l", col = "#40B0A6", ylim = c(0,2), lwd = 2)
plot(testdf$time, testdf$S, type = "l", cex.lab = 1.5,main = "S. enterica", col = "#E1BE6A", ylim = c(0,2), lwd = 2, ylab = "", xlab = "Time")
lines(test2df$time, test2df$S, type = "l", col = "#40B0A6", ylim = c(0,2), lwd = 2)
plot(testdf$time, testdf$gen, cex.lab = 1.5,main = "Generalist", type = "l", col = "#E1BE6A", ylim = c(0,1200), lwd = 2, ylab = "Abundance", xlab = "Time")
lines(test2df$time, test2df$gen, type = "l", col = "#40B0A6", ylim = c(0,1200), lwd = 2)
plot(testdf$time, testdf$sp, cex.lab = 1.5,main = "Specialist", type = "l", col = "#E1BE6A", ylim = c(0,1200), lwd = 2, ylab = "", xlab = "Time")
lines(test2df$time, test2df$sp, type = "l", col = "#40B0A6", ylim = c(0,1200), lwd = 2)

realparams <- relativefitnessbasegencost + barchartmodel + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

#################################################################################

#COST AND GROWTH RATE
#growth rate and cost of generalism, cooperation
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_coop
gamma <- seq(from = 0, to = parameters['gamma_sp'] * 5, by = 0.01)

cycled_gamma_coop <- data.frame()
for (i in 1:length(gamma)){
  parameters_coop['gamma_sp']=gamma[i]
  output <- LVgrowthrate(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 5, 0.1)
  generalism_cost = gamma[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_coop = rbind(cycled_gamma_coop, output)
}

#growth rate and cost of generalism, competition
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

parameters <- parameters_comp
gamma <- seq(from = 0, to = parameters['gamma_sp'] * 5, by = 0.01)

cycled_gamma_comp <- data.frame()
for (i in 1:length(gamma)){
  parameters_comp['gamma_sp']=gamma[i]
  output <- LVgrowthrate(generalLV_comp, parameters_comp, "Competition", start_density, time, 5, 0.1)
  generalism_cost = gamma[i] / 2e-2
  output=data.frame(output) %>%
    mutate(generalism_cost = generalism_cost) %>%
    mutate(cost_type = "Productivity (Gamma)")
  cycled_gamma_comp = rbind(cycled_gamma_comp, output)
}

#heatmaps
cleaned_coop <- cycled_gamma_coop %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5.5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5.5,
                           TRUE ~ log10)) 

plot_coop <- cleaned_coop %>%
  filter(rate_cost %in% c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
  geom_tile(aes(fill = log10)) +
  scale_fill_gradient2(low = "#CA3542",
                         mid = "white",
                         high= "#27647B",
                         midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=24,face="bold"),
        axis.title.y  = element_text(size=24,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title=element_text(size=20))+
labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Relative E. coli Growth Rate")+
  xlab("Benefit of Specialism") +
  geom_hline(yintercept=c(3), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

cleaned_comp <- cycled_gamma_comp %>%
  dplyr::select(-c(repro_generalist, repro_specialist, generalist, specialist, cost_type)) %>%
  melt(id = c("rate_cost", "generalism_cost", "label", "type")) %>%
  dplyr::select(-c(variable)) %>%
  mutate(log10 = log10(value)) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5.5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5.5,
                           TRUE ~ log10))

plot_comp <- cleaned_comp %>%
  filter(rate_cost %in% c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)) %>%
  ggplot(aes(x = as.character(generalism_cost), y = as.character(round(rate_cost,1)))) +
  geom_tile(aes(fill = log10)) +
  scale_fill_gradient2(low = "#CA3542",
                       mid = "white",
                       high= "#27647B",
                       midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=24,face="bold"),
        axis.title.y  = element_text(size=24,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title=element_text(size=18))+
  labs(fill = "Specialist Phage Relative Fitness") +
  ylab("Relative E. coli Growth Rate")+
  xlab("Benefit of Specialism")+
  geom_hline(yintercept=c(5), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

patched <- plot_coop / plot_comp + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.title = element_text(size = 18))

#################################################################################

#net change in percent generalist
path = "25April2022 Phi vs P22 Flasks"
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("growthcurvecompassayMay13.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "growthcurvecompassayMay13.xlsx")
names(datalist) <- c("PFU_raw", "PFU_stat", "Percent_change", "Fullchange", "Allchange", "PFU_double")
list2env(datalist, .GlobalEnv)

#between replicate 
stat.test1 <- Fullchange %>%
  wilcox_test(NetChange ~ Condition) %>%
  add_significance()
stat.test1 <- stat.test1 %>% add_xy_position(fun = "mean_sd")

percentchange <- Fullchange %>%
  mutate(Condition = case_when(Condition == "Fac" ~ "Facilitation",
                               Condition == "Comp" ~ "Competition",
                               Condition == "Coop" ~ "Cooperation",
                               Condition == "S" ~ "S Monoculture",
                         TRUE ~ Condition)) %>%
  group_by(Condition) %>%
  summarize(Average = mean(NetChange), STDEV = sd(NetChange)) %>%
  ggplot(aes(x=Condition, y = Average)) +
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.2, position = position_dodge(0.01)) +
  geom_signif(comparisons = list(c("Competition", "Cooperation"), c("Competition", "Facilitation"), 
                                  c("Competition", "S Monoculture"), c("Cooperation", "Facilitation"), 
                                  c("Cooperation", "S Monoculture"), c("Facilitation", "S Monoculture")),
              test = "wilcox.test", step_increase = 0.075,
              map_signif_level = TRUE, tip_length = 0)+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 1)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title=element_text(size=18)) + 
  ylab("Change in Percent Generalist")+
  xlab("Condition")

#within and between replicates
change <- Allchange %>%
  mutate(Condition = case_when(Condition == "Fac1" ~ "Facilitation",
                               Condition == "Comp1" ~ "Competition",
                               Condition == "Coop1" ~ "Cooperation",
                               Condition == "S1" ~ "S Monoculture",
                               Condition == "Fac2" ~ "Facilitation",
                               Condition == "Comp2" ~ "Competition",
                               Condition == "Coop2" ~ "Cooperation",
                               Condition == "S2" ~ "S Monoculture",
                               Condition == "Fac3" ~ "Facilitation",
                               Condition == "Comp3" ~ "Competition",
                               Condition == "Coop3" ~ "Cooperation",
                               Condition == "S3" ~ "S Monoculture",
                               TRUE ~ Condition)) %>%
  filter(Condition != "Facilitation")

stat.test <- change %>%
  wilcox_test(PerChange ~ Condition) %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(fun = "mean_sd")

significance <- change %>%
  group_by(Condition) %>%
  summarize(Average = mean(PerChange), STDEV = sd(PerChange)) %>%
  ggplot(aes(x=Condition, y = Average)) +
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.2, position = position_dodge(0.01)) +
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 1)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        legend.title=element_text(size=18)) + 
  ylab("Change in Percent Generalist Phage")+
  xlab("Condition") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)

averages <- PFU_double %>%
  mutate(STDEV = case_when(Condition == "S1" & Timepoint == 24 & Type == "Specialist" ~ 0.145,
                         TRUE ~ STDEV)) %>%
  mutate(STDEV = case_when(Condition == "S3" & Timepoint == 24 & Type == "Specialist" ~ 1.04,
                           TRUE ~ STDEV)) %>%
  mutate(STDEV = case_when(Condition == "Coop3" & Timepoint == 48 & Type == "Generalist" ~ 0,
                           TRUE ~ STDEV)) %>%
  mutate(new = case_when(Condition == "Fac1" & Timepoint == 48 & Type == "Generalist" ~ "Comp1",
                         Condition == "Fac2" & Timepoint == 48 & Type == "Generalist" ~ "Comp2",
                         Condition == "Fac3" & Timepoint == 48 & Type == "Generalist" ~ "Comp3",
                         Condition == "Comp1" & Timepoint == 48 & Type == "Generalist" ~ "Fac1",
                         Condition == "Comp2" & Timepoint == 48 & Type == "Generalist" ~ "Fac2",
                         Condition == "Comp3" & Timepoint == 48 & Type == "Generalist" ~ "Fac3",
                         TRUE ~ Condition)) %>%
  mutate(new = case_when(new == "Fac1" ~ "Facilitation",
                         new == "Comp1" ~ "Competition",
                         new == "Coop1" ~ "Cooperation",
                         new == "S1" ~ "S Monoculture",
                         new == "Fac2" ~ "Facilitation",
                         new == "Comp2" ~ "Competition",
                         new == "Coop2" ~ "Cooperation",
                         new == "S2" ~ "S Monoculture",
                         new == "Fac3" ~ "Facilitation",
                         new == "Comp3" ~ "Competition",
                         new == "Coop3" ~ "Cooperation",
                         new == "S3" ~ "S Monoculture",
                         TRUE ~ Condition)) %>%
  filter(new != "Facilitation") %>%
  rename(Phage = Type) %>%
  group_by(new, Timepoint, Phage) %>%
  summarize(Average1 = mean(Average), STDEV1 = sd(STDEV)) %>%
  ungroup() %>%
  ggplot(aes(x=Timepoint, y = Average1, color = Phage))+
  facet_wrap(~new, ncol = 1, nrow = 3)+
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#CA3542", "#27647B")) +
  geom_errorbar(aes(ymin=Average1-STDEV1, ymax=Average1+STDEV1), width=0.2, position = position_dodge(0.01)) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12)) + 
  ylab("Number of Phage Doublings (log2)")+
  xlab("Timepoint (hours)")
  
exp <- significance + averages 

