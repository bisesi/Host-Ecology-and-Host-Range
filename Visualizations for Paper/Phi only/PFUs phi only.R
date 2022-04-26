#ATB
#Phage growth curves phi only

library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")
library("readxl")

#Set directory and import data
path = "25April2022 Phi vs P22 Flasks"
source(here::here("Visualizations for Paper", "custom-theme.R"))
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("24compassay.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "24compassay.xlsx")
names(datalist) <- c("PFU_raw", "PFU_stat", "CFU_raw", "CFU_stat")
list2env(datalist, .GlobalEnv)

#EDA stats
cfustatraw <- CFU_stat %>%
  filter(Condition != "Comp control" & Condition != "Coop control" & Condition != "Fac control" & Condition != "S control") %>%
  rename(PerE = 'Percent E') %>%
  ggplot(aes(x=Condition, y = Average))+
  facet_wrap(~Plate, nrow = 2)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.01))

cfustatcontrol <- CFU_stat %>%
  filter(Condition == "Comp control" | Condition == "Coop control" | Condition == "Fac control" | Condition == "S control") %>%
  rename(PerE = 'Percent E') %>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.01))

perEbar <- CFU_stat %>%
  filter(Condition != "S1" & Condition != "S2" & Condition != "S3" & Condition != "S control" )%>%
  rename(PerE = 'Percent E') %>%
  ggplot(aes(x=Condition, y = PerE))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)

pfustatraw <- PFU_stat %>%
  rename(PerG = 'Percent Gen') %>%
  ggplot(aes(x=Condition, y = Average))+
  facet_wrap(~Plate, nrow = 2)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.01))

perGbar <- PFU_stat %>%
  rename(PerG = 'Percent Gen')%>%
  ggplot(aes(x=Condition, y = PerG))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)

generalistchange <- PFU_stat %>%
  rename(PerG = 'Percent Gen') %>%
  filter(Condition != "Start") %>%
  select(Condition, PerG) %>%
  mutate(change = PerG - 0.5)

perchangebar <- generalistchange %>%
  ggplot(aes(x=Condition, y = change))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)





