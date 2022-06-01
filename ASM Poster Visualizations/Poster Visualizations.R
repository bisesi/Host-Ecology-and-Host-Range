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
source(here::here("Computational Models", "Lotka Volterra Model.R"))

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
par(mfrow = c(3,2), mar = c(4.5,4.5,4.5,4.5))
plot(data.LV.coop$time,data.LV.coop$E, main = "E", xlab="Time", ylab="Abundance",type="l", cex.lab = 1.5)
mtext("A", side=2, line=1, at=1, col="black", las=1.15, cex = 1.5)
plot(data.LV.coop$time,data.LV.coop$S, main = "S", xlab="Time", ylab="",type="l", cex.lab = 1.5)
plot(data.LV.comp$time,data.LV.comp$E, main = "E", xlab="Time", ylab="Abundance",type="l", cex.lab = 1.5)
mtext("B", side=2, line=1, at=1.25, col="black", las=1, cex = 1.5)
plot(data.LV.comp$time,data.LV.comp$S, main = "S", xlab="Time", ylab="",type="l", cex.lab = 1.5)
plot(data.LV.coop$E, data.LV.coop$S, col = "white", main = "Cooperation", xlab="dE/dt", ylab="dS/dt",type="l", cex.lab = 1.5, xlim=c(0,2.5), ylim=c(0,2.5))
nullclines(model.coop, xlim=c(0.1,2.5),ylim=c(0.1,2.5), col = c("black", "black"), parameters=params.coop, system="two.dim", add = TRUE, add.legend = FALSE)
mtext("C", side=2, line=1, at=3.15, col="black", las=1, cex = 1.5)
plot(data.LV.comp$E, data.LV.comp$S, col = "white", main = "Competition", xlab="dE/dt", ylab="",type="l", cex.lab = 1.5, xlim=c(0,2.5), ylim=c(0,2.5))
nullclines(model.comp, xlim=c(0.1,2.5),ylim=c(0.1,2.5), col = c("black", "black"), parameters=params.comp, system="two.dim", add = TRUE, add.legend = FALSE)

#with phage, no cost
time = seq(from = 0 , to = 5e4, by = 10)
start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

out_coop_with <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_coop,
                              parms = parameters_coop,
                              atol = 1e-14)) %>%
  mutate(label = "Cooperation")

out_comp_with <- as.data.frame(ode(y=start_density,
                              times=time,
                              func=generalLV_comp,
                              parms = parameters_comp,
                              atol = 1e-14)) %>%
  mutate(label = "Competition")

withphagehosts <- rbind(out_coop_with, out_comp_with) %>%
  melt(id = c("time", "label")) %>%
  rename(Microbe = variable) %>%
  filter(Microbe == "E" | Microbe == "S") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(label ~ Microbe) +
  theme_base() +
  theme(axis.title = element_text()) + 
  ylab("Abundance")+
  xlab("Time")

withphagephage <- rbind(out_coop_with, out_comp_with) %>%
  melt(id = c("time", "label")) %>%
  rename(Microbe = variable) %>%
  filter(Microbe == "gen" | Microbe == "sp") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_grid(label ~ Microbe) +
  theme_base() +
  theme(axis.title = element_text()) + 
  ylab("Abundance")+
  xlab("Time")

#example plots
proof_of_concept <- withphagehosts / withphagephage +
  plot_annotation(tag_levels = 'A')


#ONLY COST OF GENERALISM
time = seq(from = 0 , to = 5e5, by = 10)

start_density <- c(E = 0.1, S = 0.1, gen = 0.1, sp = 0.1)

costofgeneralismcoop <- LVgeneralismcost(generalLV_coop, parameters_coop, "Cooperation", start_density, time, 10)
costofgeneralismcomp <- LVgeneralismcost(generalLV_comp, parameters_comp, "Competition", start_density, time, 10)

completedatacost <- rbind(costofgeneralismcoop, costofgeneralismcomp)

datasetbasegencost <- completedatacost %>%
  select(-c(repro_generalist, repro_specialist, generalist, specialist)) %>%
  melt(id = c("cost", "label", "type")) %>%
  rename(Interaction = label) %>%
  select(-c(variable))

relativefitnessbasegencost <- datasetbasegencost %>%
  rbind(., c(10.0, "Cooperation", "Gamma", Inf)) %>%
  filter(type != "Consumption Rate") %>%
  mutate(log10 = log10(as.numeric(value))) %>%
  mutate(log10 = case_when(is.finite(log10) == FALSE & value == 0 ~ -5,
                           is.finite(log10) == FALSE & is.finite(value) == FALSE ~ 5,
                           TRUE ~ log10)) %>%
  ggplot(aes(x=as.numeric(cost), y = log10, color = Interaction)) +
  geom_line(size = 2)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text()) + 
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  ylab("Relative Fitness of Specialist")+
  xlab("Fold-Increase in Specialist Productivity")+
  ylim(-5, 5)

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
  scale_fill_gradient2(low = "red",
                         mid = "white",
                         high= "blue",
                         midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  labs(x = "Fold-Increase in E. Coli Growth Rate", y = "Fold-Increase in Cost of Phage Generalism (Productivity)",
       fill = "Specialist Relative Fitness") +
  theme(axis.title = element_text()) + 
  ylab("Fold-Increase in E. Coli Growth Rate")+
  xlab("Fold-Increase in Specialist Productivity") +
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
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "blue",
                       midpoint = 0,
                       limits = c(-5.5,5.5))+
  theme_fivethirtyeight() +
  labs(x = "Fold-Increase in E. Coli Growth Rate", y = "Fold-Increase in Cost of Phage Generalism (Productivity)",
       fill = "Specialist Relative Fitness") +
  theme(axis.title = element_text()) + 
  ylab("Fold-Increase in E. Coli Growth Rate")+
  xlab("Fold-Increase in Specialist Productivity") +
  geom_hline(yintercept=c(5), color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept=c(3), color = "black", linetype = "dashed", size = 1)

patched <- plot_coop / plot_comp + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

#net change in percent generalist
path = "25April2022 Phi vs P22 Flasks"
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("growthcurvecompassayMay13.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "growthcurvecompassayMay13.xlsx")
names(datalist) <- c("PFU_raw", "PFU_stat", "Percent_change", "Fullchange", "PFU_double")
list2env(datalist, .GlobalEnv)

percentchange <- Fullchange %>%
  mutate(Condition = case_when(Condition == "Fac" ~ "Facilitation",
                               Condition == "Comp" ~ "Competition",
                               Condition == "Coop" ~ "Cooperation",
                               Condition == "S" ~ "S Monoculture",
                         TRUE ~ Condition)) %>%
  ggplot(aes(x=Condition, y = NetChange)) +
  geom_bar(stat = "identity", position = "dodge")+
  stat_compare_means(comparisons = list( c("Competition", "Cooperation"), c("Competition", "Facilitation"), 
                                         c("Competition", "S Monoculture"), c("Cooperation", "Facilitation"), 
                                         c("Cooperation", "S Monoculture"), c("Facilitation", "S Monoculture")))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 2)+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 1)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  ylab("Change in Percent Generalist")+
  xlab("Condition")

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
  rename(Phage = Type) %>%
  group_by(new, Timepoint, Phage) %>%
  summarize(Average1 = mean(Average), STDEV1 = sd(STDEV)) %>%
  ungroup() %>%
  ggplot(aes(x=Timepoint, y = Average1, color = Phage))+
  facet_wrap(~new, ncol = 2, nrow = 2)+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=Average1-STDEV1, ymax=Average1+STDEV1), width=0.2, position = position_dodge(0.01)) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  ylab("Number of Doublings (log2)")+
  xlab("Timepoint (hours)")
  

