#ATB
#Netchange in percent geneneralists

#load packages
library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")

#Set working directories and get source files
source("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Visualizations for Paper/custom-theme.R")
setwd("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Experimental Data/9Feb2022 Flasks")
table <- read.csv("phage-netchange.csv")

#Clean and analyze data
tablestart <- table %>%
  select(c(Generalists, Start, Start.STDEV)) %>%
  dplyr::rename(c(Time0 = Start, STDEV = Start.STDEV)) %>%
  melt(id = c("Generalists", "STDEV")) 

tableend <- table %>%
  select(c(Generalists, End, STDEV)) %>%
  dplyr::rename(c(Time24 = End, STDEV = STDEV)) %>%
  melt(id = c("Generalists", "STDEV")) %>%
  filter(Generalists == "Coop + P22 + phi" | Generalists == "Comp + P22 + phi" | Generalists == "S + P22 + phi")

melted <- rbind(tablestart, tableend)

melted <- melted %>%
  mutate(Generalists = ifelse(Generalists == "Coop + P22 + phi", "Cooperation",
                              ifelse(Generalists == "Comp + P22 + phi", "Competition",
                                     ifelse(Generalists == "S + P22 + phi", "S0", NA)))) %>%
  dplyr::rename(c(Hours = variable))

barchart <- table %>%
  select(c(Generalists, Netchange, STDEV)) %>%
  filter(Generalists == "Coop + P22 + phi" | Generalists == "Comp + P22 + phi" | Generalists == "S + P22 + phi") %>%
  mutate(Generalists = ifelse(Generalists == "Coop + P22 + phi", "Cooperation",
                              ifelse(Generalists == "Comp + P22 + phi", "Competition",
                                     ifelse(Generalists == "S + P22 + phi", "S0", NA))))

#Make plots
ratio <- melted %>%
  ggplot(aes(x = Generalists, y = value, color = Hours))+
  geom_point(size = 3, aes(shape = Hours))+
  theme_bisesi()+
  theme(axis.title.x=element_blank())+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6"))+
  ylab("Percent Generalists")+
  ylim(0, 1)+
  geom_errorbar(aes(ymin=value-STDEV, ymax=value+STDEV), width=.2)

barchartchange <- barchart %>%
  ggplot(aes(x=Generalists, y=Netchange, fill=Generalists))+ 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Netchange-STDEV, ymax=Netchange+STDEV), width=.2,
                position=position_dodge(.9))+
  ylab("Net Change in Percent Generalists")+
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6", "#000000"))+
  xlab(" ")+
  theme(axis.title.x=element_blank())+
  theme_bisesi()+
  ylim(-0.75, 0.75)+
  theme(legend.position = "none")

#Make figure
figure4 <- (barchartchange | ratio)+
  plot_annotation(tag_levels = "A")
