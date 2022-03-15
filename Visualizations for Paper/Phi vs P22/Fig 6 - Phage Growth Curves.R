#ATB
#Phage growth curves in isolation and competition assays

library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")

#Set directory and import data
source("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Visualizations for Paper/custom-theme.R")
setwd("/Users/abisesi/Desktop/PhD/Projects/Host-Ecology-and-Host-Range/Experimental Data/9Feb2022 Flasks")
table <- read.csv("24hourcurves.csv")

#Clean and analyze data for the first 2 hours of growth
dataset <- table %>%
  dplyr::rename(Time = Time..hours.) %>%
  melt(id = "Time") %>%
  mutate(Phage = case_when(grepl("Phi", variable) ~ "Phi",
                           TRUE ~ "P22")) %>%
  mutate(Treatment = case_when(grepl("Single", variable) ~ "Single",
                               TRUE ~ "Double")) %>%
  mutate(Condition = case_when(grepl("Coop", variable) ~ "Cooperation",
                               grepl("Comp", variable) ~ "Competition",
                               TRUE ~ "S0")) %>%
  mutate(Data = case_when(grepl("STDEV", variable) ~ "STDEV",
                               TRUE ~ "Average"))%>%
  dplyr::select(-c(variable)) %>%
  unite("Phage", Phage, Treatment, sep = " ") %>%
  split(.$Data)

merged <- cbind(dataset$Average, dataset$STDEV) 
colnames(merged) <- c("Time", "Average", "Phage", 
                      "Condition", "Data1", "Time2", "STDEV", 
                      "Phage2", "Condition2", "Data2")

#Clean and analyze data for full 24 hours of growth
dataset_24 <- table %>%
  dplyr::rename(Time = Time..hours.) %>%
  melt(id = "Time") %>%
  mutate(Phage = case_when(grepl("Phi", variable) ~ "Phi",
                           TRUE ~ "P22")) %>%
  mutate(Treatment = case_when(grepl("Single", variable) ~ "Single",
                               TRUE ~ "Double")) %>%
  mutate(Condition = case_when(grepl("Coop", variable) ~ "Cooperation",
                               grepl("Comp", variable) ~ "Competition",
                               TRUE ~ "S0")) %>%
  mutate(Data = case_when(grepl("STDEV", variable) ~ "STDEV",
                          TRUE ~ "Average"))%>%
  dplyr::select(-c(variable)) %>%
  unite("Phage", Phage, Treatment, sep = " ") %>%
  split(.$Data)

merged_24 <- cbind(dataset_24$Average, dataset_24$STDEV) 
colnames(merged_24) <- c("Time", "Average", "Phage", 
                      "Condition", "Data1", "Time2", "STDEV", 
                      "Phage2", "Condition2", "Data2")

#Plot 2 hour and 24 hour growth curves
plot2hour <- merged %>%
  dplyr::select(c(Time, Average, Phage, Condition, STDEV)) %>%
  filter(Time == 0 | Time == 1 | Time == 2) %>%
  ggplot(aes(x=Time*60, y = Average, color = Phage))+
  geom_line()+
  geom_point()+
  theme_bisesi()+
  facet_wrap(~Condition)+
  ylab("Number of Phage Doublings (log2)")+
  xlab("Time (minutes)")+
  ylim(-15, 15)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  geom_hline(yintercept=c(0), color = "black", linetype = "dashed", size = 0.5)+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=5, position = position_dodge(0.05))

plot24hour <- merged_24 %>%
  dplyr::select(c(Time, Average, Phage, Condition, STDEV)) %>%
  ggplot(aes(x=Time, y = Average, color = Phage))+
  geom_line()+
  geom_point()+
  theme_bisesi()+
  facet_wrap(~Condition)+
  ylab("Number of Phage Doublings (log2)")+
  xlab("Time (hours)")+
  ylim(-30, 30)+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  geom_hline(yintercept=c(0), color = "black", linetype = "dashed", size = 0.5)+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=1, position = position_dodge(0.05))

#Full Figure
figure6 <- (plot2hour / plot24hour)+
  plot_annotation(tag_levels = "A")+
  plot_layout(guide = "collect")
