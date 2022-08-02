#ATB
#Phage growth curves phi only

library("reshape")
library("tidyverse")
library("patchwork")
library("readxl")

#Set directory and import data
path = "Phi vs P22 Flasks"
source(here::here("Paper Visualizations", "custom-theme.R"))
setwd(here::here("Experimental Data", path))
sheetnames <- excel_sheets("compassay17+23June2022.xlsx")
datalist <- lapply(sheetnames, read_excel, path = "compassay17+23June2022.xlsx")
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

PFU_raw_23_E <- PFU_raw_23 %>% filter(Plate == "E")
PFU_raw_23_S <- PFU_raw_23 %>% filter(Plate == "S")
sd_changeinpercentgen_23 <- PFU_raw_23_E %>%
  inner_join(., PFU_raw_23_S, by = c("Condition", "Rep")) %>%
  filter(Condition != "Start") %>%
  mutate(pergen = PFU.x / (PFU.y + PFU.x)) %>%
  mutate(changeinpercentgen = (pergen - 0.528)*100) %>%
  group_by(Condition) %>%
  summarize(sd_changeinpercentgen = sd(changeinpercentgen, na.rm = TRUE))

PFU_23_withdoubling <- PFU_raw_23 %>% 
  mutate(doublings = case_when(Plate == "E" ~ log2(PFU / 181), 
                               Plate == "S" ~ log2(PFU / 161), 
                               TRUE ~ PFU)) %>%
  group_by(Condition, Plate) %>%
  mutate(doubling_avg = mean(doublings, na.rm = TRUE),
         doubling_sd = sd(doublings, na.rm = TRUE)) %>%
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

PFU_raw_17_E <- PFU_raw_17 %>% filter(Plate == "E")
PFU_raw_17_S <- PFU_raw_17 %>% filter(Plate == "S")
sd_changeinpercentgen_17 <- PFU_raw_17_E %>%
  inner_join(., PFU_raw_17_S, by = c("Condition", "Rep")) %>%
  filter(Condition != "Start") %>%
  mutate(pergen = PFU.x / (PFU.y + PFU.x)) %>%
  mutate(changeinpercentgen = (pergen - 0.901)*100) %>%
  group_by(Condition) %>%
  summarize(sd_changeinpercentgen = sd(changeinpercentgen, na.rm = TRUE))

PFU_17_withdoubling <- PFU_raw_17 %>% 
  mutate(doublings = case_when(Plate == "E" ~ log2(PFU / 203), 
                               Plate == "S" ~ log2(PFU / 22.2), 
                               TRUE ~ PFU)) %>%
  group_by(Condition, Plate) %>%
  mutate(doubling_avg = mean(doublings, na.rm = TRUE),
         doubling_sd = sd(doublings, na.rm = TRUE)) %>%
  mutate(date = "June17")

#EDA stats
percentchange_17_new <- inner_join(percentchange_17, sd_changeinpercentgen_17, by = c("Condition"))
percentchange_23_new <- inner_join(percentchange_23, sd_changeinpercentgen_23, by = c("Condition"))
CFU_together <- rbind(CFU_averaged_17, CFU_averaged_23)
PFU_together <- rbind(PFU_averaged_17, PFU_averaged_23)
PFU_together_trans <- rbind(PFU_averaged_trans_17, PFU_averaged_trans_23)
percentchange_together <- rbind(percentchange_17_new, percentchange_23_new)
percentchange_together_trans <- rbind(percentchange_trans_17, percentchange_trans_23)
doublings <- rbind(PFU_17_withdoubling, PFU_23_withdoubling)

cfuplotS <- CFU_together %>%
  mutate(Interaction = case_when(grepl("Comp", Condition) ~ "Competition",
                                 grepl("Coop", Condition) ~ "Cooperation",
                                 grepl("Fac", Condition) ~ "Facilitation",
                                 grepl("S", Condition) ~ "S Monoculture",
                                 TRUE ~ Condition)) %>%
  mutate(Replicate = case_when(grepl("1", Condition) ~ "1",
                               grepl("2", Condition) ~ "2",
                               grepl("3", Condition) ~ "3",
                               TRUE ~ Condition)) %>%
  filter(Condition != "Comp control" & Condition != "Coop control" & Condition != "Fac control" & Condition != "S control") %>%
  filter(Condition == "S1" | Condition == "S2" | Condition == "S3") %>%
  ggplot(aes(x=Replicate, y = Average_CFU))+
  facet_grid(date~ Interaction, scales = "free_y")+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=Average_CFU-sd_CFU, ymax=Average_CFU+sd_CFU), width=0.02, position = position_dodge(0.9))+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  ylab("Final CFU/mL")+
  xlab("Replicate")+
  labs(fill = " ")
  

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
  filter(Condition != "Start") %>%
  ggplot(aes(x=Condition, y = Average_PFU, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~date, ncol = 1)+
  geom_errorbar(aes(ymin=Average_PFU-sd_PFU, ymax=Average_PFU+sd_PFU), width=0.02, position = position_dodge(0.9))+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  ylab("Final PFU/mL")+
  xlab("Replicate")

pfuaveraged_trans <- PFU_together_trans %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Condition, y = Average_PFU, fill = Plate))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~date, ncol = 1)+
  geom_errorbar(aes(ymin=Average_PFU-sd_PFU, ymax=Average_PFU+sd_PFU), width=0.02, position = position_dodge(0.9))

perGbar <- percentchange_together %>%
  mutate(Interaction = case_when(grepl("Comp", Condition) ~ "Competition",
                                 grepl("Coop", Condition) ~ "Cooperation",
                                 grepl("Fac", Condition) ~ "Facilitation",
                                 grepl("S", Condition) ~ "S Monoculture",
                                 TRUE ~ Condition)) %>%
  mutate(Replicate = case_when(grepl("1", Condition) ~ "1",
                                 grepl("2", Condition) ~ "2",
                                 grepl("3", Condition) ~ "3",
                                 TRUE ~ Condition)) %>%
  filter(Condition != "Start" & Interaction == "S Monoculture") %>%
  ggplot(aes(x=Replicate, y = changeinpercentgen))+
  facet_grid(date ~ Interaction) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  ylab("Change in Percent Generalist")+
  xlab("Replicate")+
  geom_errorbar(aes(ymin=changeinpercentgen-sd_changeinpercentgen, ymax=changeinpercentgen+sd_changeinpercentgen), width=0.02, position = position_dodge(0.9))


perGbar_trans <- percentchange_together_trans %>%
  filter(Condition != "Start") %>%
  select(Condition, changeinpercentgen, date) %>%
  ggplot(aes(x=Condition, y = changeinpercentgen))+
  facet_wrap(~date) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_hline(yintercept=c(0.5), color = "red", linetype = "dashed", size = 0.5)

doublingsplot <- doublings %>%
  mutate(Interaction = case_when(grepl("Comp", Condition) ~ "Competition",
                                 grepl("Coop", Condition) ~ "Cooperation",
                                 grepl("Fac", Condition) ~ "Facilitation",
                                 grepl("S", Condition) ~ "S Monoculture",
                                 TRUE ~ Condition)) %>%
  mutate(Phage = case_when(grepl("E", Plate) ~ "PhiC",
                                 grepl("S", Plate) ~ "P22",
                                 TRUE ~ Plate)) %>%
  mutate(Replicate = case_when(grepl("1", Condition) ~ "1",
                               grepl("2", Condition) ~ "2",
                               grepl("3", Condition) ~ "3",
                               TRUE ~ Condition)) %>%
  filter(Condition != "Start") %>%
  ggplot(aes(x=Replicate, y = doubling_avg, color = Interaction))+
  facet_grid(date ~ Phage) +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=doubling_avg-doubling_sd, ymax=doubling_avg+doubling_sd), width=0.02, position = position_dodge(0.01))+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  ylab("Number of Phage Doublings")+
  xlab("Replicate")


#mixed models / stats
PFU_raw_17_new <- PFU_raw_17 %>% 
  filter(Condition != "Start") %>%
  mutate(date = "June17") %>%
  rename(Techrep = Rep) %>%
  mutate(Interaction = case_when(grepl("Comp", Condition) ~ "Competition",
                                 grepl("Coop", Condition) ~ "Cooperation",
                                 grepl("Fac", Condition) ~ "Facilitation",
                                 grepl("S", Condition) ~ "S Monoculture",
                                 TRUE ~ Condition)) %>%
  mutate(Phage = case_when(grepl("E", Plate) ~ "PhiC",
                           grepl("S", Plate) ~ "P22",
                           TRUE ~ Plate)) %>%
  mutate(Replicate = case_when(grepl("1", Condition) ~ "1",
                               grepl("2", Condition) ~ "2",
                               grepl("3", Condition) ~ "3",
                               TRUE ~ Condition)) 

PFU_raw_23_new <- PFU_raw_23 %>% 
  filter(Condition != "Start") %>%
  mutate(date = "June23") %>%
  rename(Techrep = Rep) %>%
  mutate(Interaction = case_when(grepl("Comp", Condition) ~ "Competition",
                                 grepl("Coop", Condition) ~ "Cooperation",
                                 grepl("Fac", Condition) ~ "Facilitation",
                                 grepl("S", Condition) ~ "S Monoculture",
                                 TRUE ~ Condition)) %>%
  mutate(Phage = case_when(grepl("E", Plate) ~ "PhiC",
                           grepl("S", Plate) ~ "P22",
                           TRUE ~ Plate)) %>%
  mutate(Replicate = case_when(grepl("1", Condition) ~ "1",
                               grepl("2", Condition) ~ "2",
                               grepl("3", Condition) ~ "3",
                               TRUE ~ Condition)) 

PFU_raw_17_Enew <- PFU_raw_17_new %>% filter(Plate == "E")
PFU_raw_17_Snew <- PFU_raw_17_new %>% filter(Plate == "S")
sd_changeinpercentgen_17 <- PFU_raw_17_Enew %>%
  inner_join(., PFU_raw_17_Snew, by = c("Condition", "Techrep")) %>%
  mutate(pergen = PFU.x / (PFU.y + PFU.x)) %>%
  mutate(changeinpercentgen = ((pergen * (1 - 0.901)) / (0.901 * (1 - pergen)))) 

PFU_raw_23_Enew <- PFU_raw_23_new %>% filter(Plate == "E")
PFU_raw_23_Snew <- PFU_raw_23_new %>% filter(Plate == "S")
sd_changeinpercentgen_23 <- PFU_raw_23_Enew %>%
  inner_join(., PFU_raw_23_Snew, by = c("Condition", "Techrep")) %>%
  mutate(pergen = PFU.x / (PFU.y + PFU.x)) %>%
  mutate(changeinpercentgen = ((pergen * (1 - 0.528)) / (0.528 * (1 - pergen)))) 

stats <- rbind(sd_changeinpercentgen_17, sd_changeinpercentgen_23) 
stats_subset <- stats %>%
  select(c(changeinpercentgen, Interaction.y, date.y, Replicate.y, Techrep))
stats_subset$Techrep <- as.factor(stats_subset$Techrep)
stats_subset$Interaction.y <- as.factor(stats_subset$Interaction.y)
stats_subset$Replicate.y <- as.factor(stats_subset$Replicate.y)
stats_subset$date.y<- as.factor(stats_subset$date.y)

lmm <- lmer(changeinpercentgen ~ Interaction.y + (1 | date.y) + (1 | Replicate.y / Techrep), data = stats_subset,
            REML = FALSE)
lmm17 <- lmer(changeinpercentgen ~ Interaction.y + (1 | Replicate.y / Techrep), data = stats_subset %>% filter(date.y == "June17"),
            REML = FALSE)
lmm23 <- lmer(changeinpercentgen ~ Interaction.y + (1 | Replicate.y / Techrep), data = stats_subset %>% filter(date.y == "June23"),
              REML = FALSE)

stats_subset %>% 
  rename(Interaction = Interaction.y, Replicate = Replicate.y, Date = date.y) %>%
  ggplot(aes(x = Interaction, y = log10(changeinpercentgen), color = Replicate)) + 
  geom_jitter(height = 0.1, width = 0.1) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  facet_wrap(~Date)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x  = element_text(size=18,face="bold"),
        axis.title.y  = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title=element_text(size=12))+
  ylab("Change in Frequency of Generalist")+
  xlab("Interaction")


