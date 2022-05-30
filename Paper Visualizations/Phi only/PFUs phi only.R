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
sheetnames <- excel_sheets("spinvsfiltercompassayMay27.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "spinvsfiltercompassayMay27.xlsx")
names(datalist) <- c("PFU_raw_spin", "PFU_raw_filter", 
                     "PFU_stat_spin", "PFU_stat_filter",
                     "CFU_raw_spin", "CFU_raw_filter",
                     "CFU_stat_spin", "CFU_stat_filter")
list2env(datalist, .GlobalEnv)

#EDA stats
cfustatraw <- CFU_stat %>%
  filter(Condition != "Comp control" & Condition != "Coop control" & Condition != "Fac control" & Condition != "S control") %>%
  ggplot(aes(x=Condition, y = Average, fill= Plate))+
  facet_grid(Time~experiment)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.9))

cfustatcontrols <- CFU_stat_spin %>%
  filter(Condition == "Comp control" | Condition == "Coop control" | Condition == "Fac control" | Condition == "S control" | Condition == "Fac 5050 1" | Condition == "Fac 5050 2") %>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.9))

perEbar <- CFU_stat %>%
  filter(Condition != "S1" & Condition != "S2" & Condition != "S3" & Condition != "S control" )%>%
  rename(PerE = 'Per E') %>%
  ggplot(aes(x=Condition, y = PerE))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)

pfustatraw <- PFU_stat %>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  facet_grid(Time~experiment)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.02, position = position_dodge(0.9))

averages <- averages %>%
  ggplot(aes(x=Timepoint, y = log10(Average), color = Type))+
  facet_wrap(~Condition, nrow = 4)+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=log10(Average)-log10(STDEV), ymax=log10(Average)+log10(STDEV)), width=0.02, position = position_dodge(0.01))

doublings <- doublings %>%
  ggplot(aes(x=Timepoint, y = Average, color = Type))+
  facet_wrap(~Condition, nrow = 4)+
  geom_point()+
  geom_line()+
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
  mutate(change = PerG - 0.5) %>%
  mutate(category = case_when(grepl("S", Condition) ~ "S0",
                               grepl("Coop", Condition) ~ "Coop",
                               grepl("Comp", Condition) ~ "Comp",
                              grepl("Fac", Condition) ~ "Fac")) %>%
  group_by(category) %>%
  mutate(stdev = sd(change, na.rm = TRUE)) %>%
  mutate(avg = mean(change, na.rm = TRUE))

perchangebar <- generalistchange %>%
  ggplot(aes(x=category, y = avg))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), width=0.02, position = position_dodge(0.9))






