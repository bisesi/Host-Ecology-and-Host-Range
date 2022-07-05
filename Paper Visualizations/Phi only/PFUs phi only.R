#ATB
#Phage growth curves phi only

library("reshape")
library("tidyverse")
library("patchwork")
library("readxl")

#Set directory and import data
path = "25April2022 Phi vs P22 Flasks"
source(here::here("Paper Visualizations", "custom-theme.R"))
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("compassaybtubStrxA.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "compassaybtubStrxA.xlsx")
names(datalist) <- c("PFU_raw_17", "CFU_raw_17", "PFU_raw_23", "CFU_raw_23")
list2env(datalist, .GlobalEnv)

#get averages for 23
CFU_averaged_23 <- CFU_raw_23 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_CFU = mean(PFU, na.rm = TRUE), sd_CFU = sd(PFU, na.rm = TRUE))%>%
  mutate(date = "June23")

PFU_averaged_23 <- PFU_raw_23 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_PFU = mean(PFU, na.rm = TRUE), sd_PFU = sd(PFU, na.rm = TRUE))%>%
  mutate(date = "June23")

percentchange_23 <- PFU_averaged_23 %>%
  ungroup() %>%
  group_by(Condition) %>%
  pivot_wider(., names_from = Plate, values_from = Average_PFU:sd_PFU) %>%
  mutate(percentgen = Average_PFU_E / (Average_PFU_E + Average_PFU_S)) %>%
  mutate(changeinpercentgen = (percentgen - 0.528)*100)%>%
  mutate(date = "June23")

PFU_averaged_trans_23 <- PFU_raw_23 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_PFU = mean(Transformed, na.rm = TRUE), sd_PFU = sd(Transformed, na.rm = TRUE)) %>%
  mutate(date = "June23")

percentchange_trans_23 <- PFU_averaged_trans_23 %>%
  ungroup() %>%
  group_by(Condition) %>%
  pivot_wider(., names_from = Plate, values_from = Average_PFU:sd_PFU) %>%
  mutate(percentgen = Average_PFU_E / (Average_PFU_E + Average_PFU_S)) %>%
  mutate(changeinpercentgen = (percentgen - 0.234)*100) %>%
  mutate(date = "June23")


#get averages for 17
CFU_averaged_17 <- CFU_raw_17 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_CFU = mean(PFU, na.rm = TRUE), sd_CFU = sd(PFU, na.rm = TRUE))%>%
  mutate(date = "June17")

PFU_averaged_17 <- PFU_raw_17 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_PFU = mean(PFU, na.rm = TRUE), sd_PFU = sd(PFU, na.rm = TRUE))%>%
  mutate(date = "June17")

percentchange_17 <- PFU_averaged_17 %>%
  ungroup() %>%
  group_by(Condition) %>%
  pivot_wider(., names_from = Plate, values_from = Average_PFU:sd_PFU) %>%
  mutate(percentgen = Average_PFU_E / (Average_PFU_E + Average_PFU_S)) %>%
  mutate(changeinpercentgen = (percentgen - 0.901)*100)%>%
  mutate(date = "June17")

PFU_averaged_trans_17 <- PFU_raw_17 %>%
  group_by(Condition, Plate) %>%
  summarize(Average_PFU = mean(Transformed, na.rm = TRUE), sd_PFU = sd(Transformed, na.rm = TRUE))%>%
  mutate(date = "June17")

percentchange_trans_17 <- PFU_averaged_trans_17 %>%
  ungroup() %>%
  group_by(Condition) %>%
  pivot_wider(., names_from = Plate, values_from = Average_PFU:sd_PFU) %>%
  mutate(percentgen = Average_PFU_E / (Average_PFU_E + Average_PFU_S)) %>%
  mutate(changeinpercentgen = (percentgen - 0.414)*100)%>%
  mutate(date = "June17")

#EDA stats
CFU_together <- rbind(CFU_averaged_17, CFU_averaged_23)
PFU_together <- rbind(PFU_averaged_17, PFU_averaged_23)
PFU_together_trans <- rbind(PFU_averaged_trans_17, PFU_averaged_trans_23)
percentchange_together <- rbind(percentchange_17, percentchange_23)
percentchange_together_trans <- rbind(percentchange_trans_17, percentchange_trans_23)

cfuplotS <- CFU_together %>%
  filter(Condition != "Comp control" & Condition != "Coop control" & Condition != "Fac control" & Condition != "S control") %>%
  filter(Condition == "S1" | Condition == "S2" | Condition == "S3") %>%
  ggplot(aes(x=Condition, y = Average_CFU, fill= Plate))+
  facet_wrap(~date)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average_CFU-sd_CFU, ymax=Average_CFU+sd_CFU), width=0.02, position = position_dodge(0.9))

cfuplotnonS <- CFU_together %>%
  filter(Condition != "Comp control" & Condition != "Coop control" & Condition != "Fac control" & Condition != "S control") %>%
  filter(Condition != "S1" & Condition != "S2" & Condition != "S3") %>%
  filter(Condition != "Start" & Condition != "Start co" & Condition != "Start S") %>%
  ggplot(aes(x=Condition, y = Average_CFU, fill= Plate))+
  facet_wrap(~date)+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average_CFU-sd_CFU, ymax=Average_CFU+sd_CFU), width=0.02, position = position_dodge(0.9))

cfuplotstart <- CFU_together %>%
  filter(Condition == "Start" | Condition == "Start co" | Condition == "Start S") %>%
  ggplot(aes(x=Condition, y = Average_CFU, fill= Plate))+
  facet_wrap(~date, scales = "free_x")+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average_CFU-sd_CFU, ymax=Average_CFU+sd_CFU), width=0.02, position = position_dodge(0.9))

cfucontrols <- CFU_together %>%
  filter(Condition == "Comp control" | Condition == "Coop control" | Condition == "Fac control" | Condition == "S control") %>%
  ggplot(aes(x=Condition, y = Average_CFU, fill= Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~date)+
  geom_errorbar(aes(ymin=Average_CFU-sd_CFU, ymax=Average_CFU+sd_CFU), width=0.02, position = position_dodge(0.9))

pfuaveraged <- PFU_together %>%
  filter(Condition == "Start") %>%
  ggplot(aes(x=Condition, y = Average_PFU, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~date, ncol = 1)+
  geom_errorbar(aes(ymin=Average_PFU-sd_PFU, ymax=Average_PFU+sd_PFU), width=0.02, position = position_dodge(0.9))

pfuaveraged_trans <- PFU_together_trans %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Condition, y = Average_PFU, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~date, ncol = 1)+
  geom_errorbar(aes(ymin=Average_PFU-sd_PFU, ymax=Average_PFU+sd_PFU), width=0.02, position = position_dodge(0.9))

perGbar <- percentchange_together %>%
  filter(Condition != "Start") %>%
  select(Condition, changeinpercentgen, date) %>%
  ggplot(aes(x=Condition, y = changeinpercentgen))+
  facet_wrap(~date, ncol = 1) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)

perGbar_trans <- percentchange_together_trans %>%
  filter(Condition != "Start") %>%
  select(Condition, changeinpercentgen, date) %>%
  ggplot(aes(x=Condition, y = changeinpercentgen))+
  facet_wrap(~date) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)





