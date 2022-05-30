#ATB
#Bacterial densities following treatment with single and double phage

#load packages
library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")

#load source
source("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Visualizations for Paper/custom-theme.R")
setwd("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Experimental Data/9Feb2022 Flasks")
table <- read.csv("E-Sratios.csv")
tableS <- read.csv("rawCFUsS.csv")
tableE <- read.csv("rawCFUsE.csv")

#Clean and analyze data
#Ratio data
tableratio <- table %>%
  dplyr::rename(Phage = X) %>%
  melt(id = "Phage") %>%
  mutate(Timepoint = case_when(grepl("Start", variable) ~ "Time 0",
                               grepl("End", variable) ~ "Time 24")) %>%
  mutate(Condition = case_when(grepl("coop", variable) ~ "Cooperation",
                               grepl("comp", variable) ~ "Competition")) %>%
  group_by(Phage, Timepoint, Condition) %>%
  mutate(STDEV = sd(value),
         Average = mean(value)) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-c(variable, value))

#Net change data for ratios
competition <- tableratio %>%
  ungroup() %>%
  filter(Condition == "Competition") %>%
  mutate(netchange = round(Average - 0.978, 3)) %>%
  filter(Timepoint == "Time 24")

cooperation <- tableratio %>%
  ungroup() %>%
  filter(Condition == "Cooperation") %>%
  mutate(netchange = round(Average - 3.38, 3)) %>%
  filter(Timepoint == "Time 24")

complete <- full_join(competition, cooperation)

#Raw data
tableratioE <- tableE %>%
  dplyr::rename(Phage = X) %>%
  melt(id = "Phage") %>%
  mutate(Timepoint = case_when(grepl("Start", variable) ~ "Time 0",
                               grepl("End", variable) ~ "Time 24")) %>%
  mutate(Condition = case_when(grepl("coop", variable) ~ "Cooperation",
                               grepl("comp", variable) ~ "Competition",
                               grepl("S0", variable) ~ "S0"))%>%
  group_by(Phage, Timepoint, Condition) %>%
  mutate(STDEV = sd(value),
         Average = mean(value)) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-c(variable, value))

tableratioS <- tableS %>%
  dplyr::rename(Phage = X) %>%
  melt(id = "Phage") %>%
  mutate(Timepoint = case_when(grepl("Start", variable) ~ "Time 0",
                               grepl("End", variable) ~ "Time 24")) %>%
  mutate(Condition = case_when(grepl("coop", variable) ~ "Cooperation",
                               grepl("comp", variable) ~ "Competition",
                               grepl("S0", variable) ~ "S0"))%>%
  group_by(Phage, Timepoint, Condition) %>%
  mutate(STDEV = sd(value),
         Average = mean(value)) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-c(variable, value))

#Net change data for raw CFUs
competitionE <- tableratioE %>%
  ungroup() %>%
  filter(Condition == "Competition") %>%
  mutate(netchange = round(Average - 8391881.443, 3)) %>%
  filter(Timepoint == "Time 24")

cooperationE <- tableratioE %>%
  ungroup() %>%
  filter(Condition == "Cooperation") %>%
  mutate(netchange = round(Average - 1699763033, 3)) %>%
  filter(Timepoint == "Time 24")

competitionS <- tableratioS %>%
  ungroup() %>%
  filter(Condition == "Competition") %>%
  mutate(netchange = round(Average - 8584313.725, 3)) %>%
  filter(Timepoint == "Time 24")

cooperationS <- tableratioS %>%
  ungroup() %>%
  filter(Condition == "Cooperation") %>%
  mutate(netchange = round(Average - 503633491.3, 3)) %>%
  filter(Timepoint == "Time 24")

monoS <- tableratioS %>%
  ungroup() %>%
  filter(Condition == "S0") %>%
  mutate(netchange = round(Average - 17211764.71, 3)) %>%
  filter(Timepoint == "Time 24")

completeE <- full_join(competitionE, cooperationE)
completeSinter <- full_join(competitionS, cooperationS)
completeS <- full_join(completeSinter, monoS)

#Visualize data
#Ratio data
finalratios <- tableratio %>%
  filter(Timepoint == "Time 24") %>%
  ggplot(aes(x=Phage, y=Average, fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=.2,
                #position=position_dodge(.9))+
  ylab("Final E:S Ratio")+
  facet_wrap(~Condition)+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  xlab(" ")+
  theme_bisesi()+
  geom_hline(yintercept=c(1), color = "red", linetype = "dashed", size = 0.5)+
  theme(legend.position = "none")

startingratios <- tableratio %>%
  filter(Timepoint == "Time 0" & Phage == "P22") %>%
  ggplot(aes(x=Phage, y=Average, fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=.2,
                #position=position_dodge(.9))+
  ylab("Starting E:S Ratio")+
  facet_wrap(~Condition)+
  scale_fill_manual(values=c("#E1BE6A"))+
  xlab(" ")+
  theme_bisesi()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_hline(yintercept=c(1), color = "red", linetype = "dashed", size = 0.5)+
  theme(legend.position = "none")

netchange <- complete %>%
  ggplot(aes(x=Phage, y=netchange, fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=netchange-STDEV, ymax=netchange+STDEV), width=.2,
                #position=position_dodge(.9))+
  ylab("Net Change in Ratio of E:S")+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  xlab(" ")+
  facet_wrap(~Condition)+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  theme_bisesi()+
  theme(legend.position = "none")

#Raw CFUs
plotECFU <- tableratioE %>%
  filter(Timepoint == "Time 24" & Condition != "S0") %>%
  ggplot(aes(x=Phage, y=log10(Average), fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=log10(Average)-log10(STDEV), ymax=log10(Average)+log10(STDEV)), width=.2,
                #position=position_dodge(.9))+
  ylab("Final CFU/mL of E. coli (log10)")+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  xlab(" ")+
  facet_wrap(~Condition)+
  theme(axis.title.x=element_blank())+
  theme_bisesi()+
  theme(legend.position = "none")

plotSCFU <- tableratioS %>%
  filter(Timepoint == "Time 24") %>%
  ggplot(aes(x=Phage, y=log10(Average), fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=log10(Average)-log10(STDEV), ymax=log10(Average)+log10(STDEV)), width=.2,
  #position=position_dodge(.9))+
  ylab("Final CFU/mL of S. enterica (log10)")+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  xlab(" ")+
  facet_wrap(~Condition)+
  theme(axis.title.x=element_blank())+
  theme_bisesi()+
  theme(legend.position = "none")

#Raw net change
netchangeE <- completeE %>%
  filter(Condition != "S0" & Phage != "Control")%>%
  ggplot(aes(x=Phage, y=netchange, fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=netchange-STDEV, ymax=netchange+STDEV), width=.2,
                #position=position_dodge(.9))+
  ylab("Net Change in CFU/mL of E. coli")+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  facet_wrap(~Condition)+
  xlab(" ")+
  theme_bisesi()+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  theme(legend.position = "none")

netchangeS <- completeS %>%
  filter(Phage != "Control") %>%
  ggplot(aes(x=Phage, y=netchange, fill=Phage))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=netchange-STDEV, ymax=netchange+STDEV), width=.2,
                #position=position_dodge(.9))+
  ylab("Net Change in CFU/mL of S. enterica")+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  
  facet_wrap(~Condition)+
  xlab(" ")+
  theme_bisesi()+
  theme(legend.position = "none")

#Figures
#ratio with netchange plot
fig5a <- ((startingratios + finalratios) / netchange)+
  plot_annotation(tag_levels = "A")

#E raw with netchange
fig5b <- (plotECFU / netchangeE)+
  plot_annotation(tag_levels = "A")

#S raw with netchange
fig5c <- (plotSCFU / netchangeS)+
  plot_annotation(tag_levels = "A")

  


