#ATB
#Phage growth curves phi only

library("ggplot2")
library("reshape")
library("tidyverse")
library("patchwork")
library("readxl")

#Set directory and import data
path = "28March2022 Phi Only"
source(here::here("Visualizations for Paper", "custom-theme.R"))
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("48hourtidy.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "48hourtidy.xlsx")
names(datalist) <- c("PFU_raw", "CFU_raw", "PFU_stat", "CFU_stat")
list2env(datalist, .GlobalEnv)
rm(Computations)

#EDA, raw data
Epfu <- PFU_raw %>%
  mutate(Rep = as.factor(Rep)) %>%
  filter(Plate == "E" & Condition != "Start1") %>%
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

#EDA stats, log10
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

#EDA stats, not log10
cfustatraw <- CFU_stat %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Day, y = Average, color = Plate))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition, nrow = 4)+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=25, position = position_dodge(0.01))

pfustatraw <- PFU_stat %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Day, y = Average, color = Plate))+
  geom_point()+
  geom_line()+
  xlim(1, 3)+
  facet_wrap(~Condition, nrow = 4)+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=25, position = position_dodge(0.01))

eopbar <- PFU_stat %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Day, y = log10(EOP)))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0), color = "red", linetype = "dashed", size = 0.5)+
  facet_wrap(~Condition, nrow = 4)

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


#specialist vs generalist data
rm(list = ls())

#Set directory and import data
source(here::here("Visualizations for Paper", "custom-theme.R"))
setwd(here::here("Experimental Data", "7March2022 Phi Only"))
sheetnames <- excel_sheets("evolved_fitness.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "evolved_fitness.xlsx")
names(datalist) <- c("Computations", "PFU_raw", "CFU_raw", "PFU_stat", "CFU_stat")
list2env(datalist, .GlobalEnv)
rm(Computations)

#condition ranges
cooprange <- c(LETTERS[1:9], LETTERS[19])
comprange <- c(LETTERS[10:18], LETTERS[20])

#cooperation               
cfustatrawcoop <- CFU_stat %>%
  filter(Condition %in% cooprange) %>%
  mutate(Condition = case_when(Condition == "A" ~ "Phi E 1",
                               Condition == "B" ~ "Phi E 2",
                               Condition == "C" ~ "Phi E 3",
                               Condition == "D" ~ "Phi S 1",
                               Condition == "E" ~ "Phi S 2",
                               Condition == "F" ~ "Phi S 3",
                               Condition == "G" ~ "Phi Gen 1",
                               Condition == "H" ~ "Phi Gen 2",
                               Condition == "I" ~ "Phi Gen 3",
                               Condition == "S" ~ "Control",
                               TRUE ~ Condition))%>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat="identity", color="black",
           position=position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.01, position = position_dodge(0.9))

pfustatrawcoop <- PFU_stat %>%
  filter(Condition %in% cooprange) %>%
  mutate(Condition = case_when(Condition == "A" ~ "Phi E 1",
                               Condition == "B" ~ "Phi E 2",
                               Condition == "C" ~ "Phi E 3",
                               Condition == "D" ~ "Phi S 1",
                               Condition == "E" ~ "Phi S 2",
                               Condition == "F" ~ "Phi S 3",
                               Condition == "G" ~ "Phi Gen 1",
                               Condition == "H" ~ "Phi Gen 2",
                               Condition == "I" ~ "Phi Gen 3",
                               TRUE ~ Condition))%>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat="identity", color="black",
           position=position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.01, position = position_dodge(0.9))

#competition             
cfustatrawcomp <- CFU_stat %>%
  filter(Condition %in% comprange) %>%
  mutate(Condition = case_when(Condition == "J" ~ "Phi E 1",
                               Condition == "K" ~ "Phi E 2",
                               Condition == "L" ~ "Phi E 3",
                               Condition == "M" ~ "Phi S 1",
                               Condition == "N" ~ "Phi S 2",
                               Condition == "O" ~ "Phi S 3",
                               Condition == "P" ~ "Phi Gen 1",
                               Condition == "Q" ~ "Phi Gen 2",
                               Condition == "R" ~ "Phi Gen 3",
                               Condition == "T" ~ "Control",
                               TRUE ~ Condition))%>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat="identity", color="black",
             position=position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.01, position = position_dodge(0.9))

pfustatrawcomp <- PFU_stat %>%
  filter(Condition %in% comprange) %>%
  mutate(Condition = case_when(Condition == "J" ~ "Phi E 1",
                               Condition == "K" ~ "Phi E 2",
                               Condition == "L" ~ "Phi E 3",
                               Condition == "M" ~ "Phi S 1",
                               Condition == "N" ~ "Phi S 2",
                               Condition == "O" ~ "Phi S 3",
                               Condition == "P" ~ "Phi Gen 1",
                               Condition == "Q" ~ "Phi Gen 2",
                               Condition == "R" ~ "Phi Gen 3",
                               TRUE ~ Condition)) %>%
  ggplot(aes(x=Condition, y = Average, fill = Plate))+
  geom_bar(stat="identity", color="black",
           position=position_dodge())+
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV), width=0.01, position = position_dodge(0.9))

#compare
compareCFUcoop <- CFU_stat %>%
  mutate(type = case_when(Plate == "E" ~ "CFU on E",
                          Plate == "S" ~ "CFU on S",
                          TRUE ~ Plate)) %>%
  filter(Condition %in% cooprange)

comparePFUcoop <- PFU_stat %>%
  mutate(type = case_when(Plate == "E" ~ "PFU on E",
                          Plate == "S" ~ "PFU on S",
                          TRUE ~ Plate)) %>%
  filter(Condition %in% cooprange)

comparecoop <- rbind(compareCFUcoop %>% select(Condition, Average, STDEV, type), comparePFUcoop %>% select(Condition, Average, STDEV, type))

plotcomparecoop <- comparecoop %>%
  ggplot(aes(x=Condition, y = Average,fill = type))+
  geom_bar(stat="identity", color="black",
           position=position_dodge())




