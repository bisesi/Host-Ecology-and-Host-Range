#ATB
#Bacterial growth curves without phage

#load packages
library("ggplot2")
library("reshape")
library("patchwork")
library("tidyverse")

#sources
source("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Computational Models/Baranyi-Functions.r")
setwd("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Experimental Data/17Sept2021 Tecan")
OD_path = "17sept.csv"

#Import and Clean Data
OD = read.table(OD_path, sep = ",", as.is = TRUE, row.names = 1)
OD = t(OD)
OD = data.frame(OD)
names(OD)[1:3] = c("cycle","seconds","temp")
OD = OD %>%
  gather(well, OD, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60)

growth_rate = OD %>%
  group_by(well) %>%
  summarize(model_fit = fit_baranyi(hour, log(OD), 5),
            fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
  ungroup()

#Analyze data for average growth rate plot
rates <- growth_rate %>% filter(fit_variable == "growth_rate")

Scomp <- rates %>% 
  filter(well == "G2" | well == "G3" | well == "G4") %>%
  mutate(average = mean(model_fit, na.rm = TRUE)) %>%
  mutate(stdev = sd(model_fit, na.rm = TRUE)) %>%
  mutate(label = "S in Glucose+Methionine")
Ecomp <- rates %>% 
  filter(well == "G5" | well == "G6" | well == "G7") %>%
  mutate(average = mean(model_fit, na.rm = TRUE)) %>%
  mutate(stdev = sd(model_fit, na.rm = TRUE)) %>%
  mutate(label = "E in Glucose+Methionine")
coop <- rates %>% 
  filter(well == "E8" | well == "E9" | well == "E10") %>%
  mutate(average = mean(model_fit, na.rm = TRUE)) %>%
  mutate(stdev = sd(model_fit, na.rm = TRUE)) %>%
  mutate(label = "ES in Lactose")
comp <- rates %>% 
  filter(well == "F2" | well == "F3" | well == "F4") %>%
  mutate(average = mean(model_fit, na.rm = TRUE)) %>%
  mutate(stdev = sd(model_fit, na.rm = TRUE)) %>%
  mutate(label = "ES in Glucose+Methionine")

cleanrates <- as.data.frame(rbind(Ecomp, Scomp, coop, comp))
order <- c("E in Glucose+Methionine", "S in Glucose+Methionine", "ES in Glucose+Methionine", "ES in Lactose")

#Analyze data for OD600 plots
OD_clean <- OD %>% 
  mutate(cycle = NULL, seconds = NULL, temp = NULL)

EcompOD <- OD_clean %>%
  filter(well == "G5" | well == "G6" | well == "G7")%>%
  group_by(hour)%>%
  summarize(average = mean(OD), stdev = sd(OD))%>%
  mutate(Condition = "E in Glucose+Methionine")
EcompODaverage <- OD_clean %>%
  filter(OD >= 0.1 & well == "G5" | OD >= 0.1 & well == "G6" | OD >= 0.1 & well == "G7") %>%
  group_by(well) %>%
  slice_min(n = 1, OD) %>%
  ungroup() %>%
  mutate(stdev = sd(hour)) %>%
  mutate(average = mean(hour)) %>%
  mutate(Condition = "E in Glucose+Methionine") %>%
  mutate(well = NULL, OD = NULL, hour = NULL) %>%
  head(n = 1)

ScompOD <- OD_clean %>%
  filter(well == "G2" | well == "G3" | well == "G4")%>%
  group_by(hour)%>%
  summarize(average = mean(OD), stdev = sd(OD))%>%
  mutate(Condition = "S in Glucose+Methionine")
ScompODaverage <- OD_clean %>%
  filter(OD >= 0.1 & well == "G2" | OD >= 0.1 & well == "G3" | OD >= 0.1 & well == "G4") %>%
  group_by(well) %>%
  slice_min(n = 1, OD) %>%
  ungroup() %>%
  mutate(stdev = sd(hour)) %>%
  mutate(average = mean(hour)) %>%
  mutate(Condition = "S in Glucose+Methionine") %>%
  mutate(well = NULL, OD = NULL, hour = NULL) %>%
  head(n = 1)

coopOD <- OD_clean %>%
  filter(well == "E8" | well == "E9" | well == "E10")%>%
  group_by(hour)%>%
  summarize(average = mean(OD), stdev = sd(OD))%>%
  mutate(Condition = "ES in Lactose")
coopODaverage <- OD_clean %>%
  filter(OD >= 0.1 & well == "E8" | OD >= 0.1 & well == "E9" | OD >= 0.1 & well == "E10") %>%
  group_by(well) %>%
  slice_min(n = 1, OD) %>%
  ungroup() %>%
  mutate(stdev = sd(hour)) %>%
  mutate(average = mean(hour)) %>%
  mutate(Condition = "ES in Lactose") %>%
  mutate(well = NULL, OD = NULL, hour = NULL) %>%
  head(n = 1)

compOD <- OD_clean %>%
  filter(well == "F2" | well == "F3" | well == "F4")%>%
  group_by(hour)%>%
  summarize(average = mean(OD), stdev = sd(OD))%>%
  mutate(Condition = "ES in Glucose+Methionine")
compODaverage <- OD_clean %>%
  filter(OD >= 0.1 & well == "F2" | OD >= 0.1 & well == "F3" | OD >= 0.1 & well == "F4") %>%
  group_by(well) %>%
  slice_min(n = 1, OD) %>%
  ungroup() %>%
  mutate(stdev = sd(hour)) %>%
  mutate(average = mean(hour)) %>%
  mutate(Condition = "ES in Glucose+Methionine") %>%
  mutate(well = NULL, OD = NULL, hour = NULL) %>%
  head(n = 1)

allODaverage <- rbind(EcompODaverage, ScompODaverage, coopODaverage, compODaverage) 
allOD <- rbind(EcompOD, ScompOD, coopOD, compOD)

#Plot 3 plots
#Average growth rate per hour
growthrateperhour <- cleanrates %>%
  ggplot(aes(x = factor(label, order), y = average, color = label))+
  geom_point()+
  theme_bisesi()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  geom_errorbar(aes(ymin=average-stdev, ymax=average+stdev), width=0.02)+
  labs(x = " ", y = "Average Growth Rate (per hour)")

fullODplot <- allOD %>%
  ggplot(aes(x = hour, y = average, colour = Condition))+
  geom_smooth(span = 0.1)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  theme_bisesi()+
  labs(x = "Time (hours)", y = "Average OD600 (Triplicate)")+
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8)) 

timetoOD0.1 <- allODaverage %>%
  ggplot(aes(x = factor(Condition, order), y = average, color = Condition))+
  geom_point()+
  theme_bisesi()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  geom_errorbar(aes(ymin=average-stdev, ymax=average+stdev), width=0.02)+
  labs(x = " ", y = "Average Time to 0.1 OD600 (hours)")

#arrange plots
figure7 <- ((growthrateperhour + timetoOD0.1) / fullODplot )+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")
