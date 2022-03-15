#ATB
#Phage growth curves phi only

library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")
library("readxl")

#Set directory and import data
source(here::here("Visualizations for Paper", "custom-theme.R"))
setwd(here::here("Experimental Data", "7March2022 Phi Only"))
sheetnames <- excel_sheets("tidydata.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "tidydata.xlsx")
names(datalist) <- c("Computations", "PFU_raw", "CFU_raw", "PFU_stat", "CFU_stat")
list2env(datalist, .GlobalEnv)
rm(Computations)

#EDA, raw data
Epfu <- PFU_raw %>%
  mutate(Rep = as.factor(Rep)) %>%
  filter(Plate == "E" & Condition != "Start1" & Condition != "Start2" & Condition != "Start3") %>%
  ggplot(aes(x=Day, y = log10(Value), color = Rep))+
  geom_point()+
  xlim(1, 3)+
  facet_wrap(~Condition, ncol = 3, nrow = 4)

Spfu <- PFU_raw %>%
  mutate(Rep = as.factor(Rep)) %>%
  filter(Plate == "S" & Condition != "Start1" & Condition != "Start2" & Condition != "Start3") %>%
  ggplot(aes(x=Day, y = log10(Value), color = Rep))+
  geom_point()+
  xlim(1, 3)+
  facet_wrap(~Condition, ncol = 3, nrow = 4)

Ecfu <- CFU_raw %>%
  mutate(Rep = as.factor(Rep)) %>%
  filter(Plate == "E" & Condition != "S1" & Condition != "S2" & Condition != "S3" & Condition != "S Control" & Value != 0) %>%
  ggplot(aes(x=Day, y = log10(Value), color = Rep))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition, ncol = 4, nrow = 3)

Scfu <- CFU_raw %>%
  mutate(Rep = as.factor(Rep)) %>%
  filter(Plate == "S" & Condition != "E1" & Condition != "E2" & Condition != "E3" & Condition != "E Control") %>%
  ggplot(aes(x=Day, y = log10(Value), color = Rep))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition, ncol = 4, nrow = 3)


#EDA stats
logPFU <- PFU_raw %>% 
  mutate(logtrans = log10(Value)) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>%
  group_by(Condition, Plate, Day) %>% 
  mutate(logavg = mean(logtrans)) %>% 
  mutate(logsd = sd(logtrans)) %>%
  select(Condition, Day, Plate, logavg, logsd) %>%
  distinct(.keep_all = TRUE)

pfustatplot <- logPFU %>%
  ungroup() %>%
  filter(Condition != "Start1" & Condition != "Start2" & Condition != "Start3") %>%
  ggplot(aes(x=Day, y = logavg, color = Plate))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition, nrow = 4, ncol = 3)+
  geom_errorbar(aes(ymin=logavg-logsd, ymax=logavg+logsd), width=25, position = position_dodge(0.01))+
  ylab("Triplicate Average PFU (log10)")

logCFU <- CFU_raw %>% 
  mutate(logtrans = log10(Value)) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>%
  group_by(Condition, Plate, Day) %>% 
  mutate(logavg = mean(logtrans)) %>% 
  mutate(logsd = sd(logtrans)) %>%
  select(Condition, Day, Plate, logavg, logsd) %>%
  distinct(.keep_all = TRUE)

cfustatplot <- logCFU %>%
  ggplot(aes(x=Day, y = logavg, color = Plate))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition)+
  geom_errorbar(aes(ymin=logavg-logsd, ymax=logavg+logsd), width=25, position = position_dodge(0.01))+
  ylab("Triplicate Average CFU (log10)")

#CFU PFU comparison
EOP <- PFU_stat %>% filter(!is.na(EOP)) %>%
  filter(Condition != "Start1" & Condition != "Start2" & Condition != "Start3") %>%
  select(Condition, Day, EOP)

comparison <- EOP %>%
  full_join(., logCFU, by = c("Condition", "Day"))

comparisonplot <- comparison %>%
  filter(Condition != "E control" & Condition != "S control" & Condition != "Comp control" & Condition != "Coop control")%>%
  filter(Condition != "E1" & Condition != "E2" & Condition != "E3")%>%
  ggplot(aes(x = Day))+
  geom_line(aes(y = EOP)) +
  geom_line(aes(y = logavg, color = Plate))+
  facet_wrap(~Condition, ncol =3)+
  ylab("Efficiency of Plating and Triplicate Average CFU (log10)")

compareE <- comparison %>%
  filter(Condition != "E control" & Condition != "S control" & Condition != "Comp control" & Condition != "Coop control")%>%
  filter(Condition == "E1" | Condition == "E2" | Condition == "E3")%>%
  ggplot(aes(x = Day))+
  geom_line(aes(y = log10(EOP))) +
  geom_line(aes(y = logavg, color = Plate))+
  facet_wrap(~Condition, ncol =3)+
  ylab("Efficiency of Plating and Triplicate Average CFU (log10)")

rawcompare <- logPFU %>% 
  mutate(label = "PFU") %>%
  mutate(type = case_when(Plate == "E" & label == "PFU" ~ "PFU on E",
                          Plate == "S" & label == "PFU" ~ "PFU on S",
                          Plate == "E" & label == "CFU" ~ "CFU on E",
                          Plate == "S" & label == "CFU" ~ "CFU on S",
                          TRUE ~ label)) %>%
  filter(Condition != "Start1" & Condition != "Start2" & Condition != "Start3") %>%
  rbind(., logCFU %>% 
          filter(Condition != "E control" & Condition != "S control" & Condition != "Comp control" & Condition != "Coop control") %>%
          mutate(label = "CFU") %>%
          mutate(type = case_when(Plate == "E" & label == "PFU" ~ "PFU on E",
                                Plate == "S" & label == "PFU" ~ "PFU on S",
                                Plate == "E" & label == "CFU" ~ "CFU on E",
                                Plate == "S" & label == "CFU" ~ "CFU on S",
                                TRUE ~ label)))

allraw <- rawcompare %>%
  ggplot(aes(x = Day, y = logavg, color = type))+
  geom_line()+
  facet_wrap(~Condition, ncol =3)+
  geom_errorbar(aes(ymin=logavg-logsd, ymax=logavg+logsd), width=0.05, position = position_dodge(0.01))+
  scale_color_manual(values=c("#E1BE6A", "#40B0A6", "#000000", "#E66100"))+
  ylab("Triplicate Average CFU or PFU (log10)")
