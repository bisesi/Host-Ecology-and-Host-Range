#ATB 
#Tecan phi and S, comp, fac


#get data
source(here::here("Computational Models", "Baranyi-Functions.R"))
setwd(here::here("Experimental Data", "Tecan 2022"))

OD_path = "od600_23july2022.csv"

#conditions
smono = c("B2", "B3", "B4", "B5", "B6")
comp = c("B7", "B8", "B9", "B10", "B11")
fac = c("C2", "C3", "C4", "C5", "C6")
sphi = c("C7", "C8", "C9", "C10", "C11")
compphi = c("D2", "D3", "D4", "D5", "D6")
facphi = c("D7", "D8", "D9", "D10", "D11")
sp22 = c("E2", "E3", "E4", "E5", "E6")
compp22 = c("E7", "E8", "E9", "E10", "E11")
facp22 = c("F2", "F3", "F4", "F5", "F6")
sp22phi = c("F7", "F8", "F9", "F10", "F11")
compp22phi = c("G2", "G3", "G4", "G5", "G6")
facp22phi = c("G7", "G8", "G9", "G10", "G11")
none = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
         "B1", "B12", "C1", "C12", "D1", "D12", "E1", "E12", "F1", "F12", "G1", "G12",
         "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")

#Import and Clean Data
OD = read.table(OD_path, sep = ",", as.is = TRUE, row.names = 1)
OD = t(OD)
OD = data.frame(OD)
names(OD)[1:3] = c("cycle","seconds","temp")
OD = OD %>%
  gather(well, OD, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60) %>%
  mutate(condition = case_when(well %in% smono ~ "S Monoculture",
                               well %in% comp ~ "Competition",
                               well %in% fac ~ "Facilitation", 
                               well %in% sphi ~ "S Monoculture + Phi",
                               well %in% compphi ~ "Competition + Phi",
                               well %in% facphi ~ "Facilitation + Phi",
                               well %in% sp22 ~ "S Monoculture + P22",
                               well %in% compp22 ~ "Competition + P22",
                               well %in% facp22 ~ "Facilitation + P22",
                               well %in% sp22phi ~ "S Monoculture + P22 + Phi",
                               well %in% compp22phi ~ "Competition + P22 + Phi",
                               well %in% facp22phi ~ "Facilitation + P22 + Phi",
                               TRUE ~ well)) 

growth_rate = OD %>%
  group_by(well) %>%
  summarize(model_fit = fit_baranyi(hour, log(OD), 5),
            fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
  ungroup() %>%
  mutate(condition = case_when(well %in% smono ~ "S Monoculture",
                               well %in% comp ~ "Competition",
                               well %in% fac ~ "Facilitation", 
                               well %in% sphi ~ "S Monoculture + Phi",
                               well %in% compphi ~ "Competition + Phi",
                               well %in% facphi ~ "Facilitation + Phi",
                               well %in% sp22 ~ "S Monoculture + P22",
                               well %in% compp22 ~ "Competition + P22",
                               well %in% facp22 ~ "Facilitation + P22",
                               well %in% sp22phi ~ "S Monoculture + P22 + Phi",
                               well %in% compp22phi ~ "Competition + P22 + Phi",
                               well %in% facp22phi ~ "Facilitation + P22 + Phi",
                               TRUE ~ well)) 

#OD600 plots
allplot <- OD %>%
  filter(well %in% smono) %>%
  filter(!well %in% none) %>%
  drop_na() %>%
  ggplot(aes(x = hour, y = OD, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition) +
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
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
        legend.title=element_text(size=12))

growth <- growth_rate %>%
  filter(!well %in% none) %>%
  filter(fit_variable == "growth_rate") %>%
  ggplot(aes(x = well, y = model_fit)) +
  geom_point() +
  facet_wrap(~condition, scales = "free_x") +
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
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
  ylab("Growth Rate")+
  xlab("Replicate")

#CFP
CFP_path = "cfp_23july2022.csv"

#Import and Clean Data
CFP = read.table(CFP_path, sep = ",", as.is = TRUE, row.names = 1)
CFP= t(CFP)
CFP = data.frame(CFP)
names(CFP)[1:3] = c("cycle","seconds","temp")
CFP = CFP %>%
  gather(well, CFP, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60) %>%
  mutate(condition = case_when(well %in% smono ~ "S Monoculture",
                               well %in% comp ~ "Competition",
                               well %in% fac ~ "Facilitation", 
                               well %in% sphi ~ "S Monoculture + Phi",
                               well %in% compphi ~ "Competition + Phi",
                               well %in% facphi ~ "Facilitation + Phi",
                               well %in% sp22 ~ "S Monoculture + P22",
                               well %in% compp22 ~ "Competition + P22",
                               well %in% facp22 ~ "Facilitation + P22",
                               well %in% sp22phi ~ "S Monoculture + P22 + Phi",
                               well %in% compp22phi ~ "Competition + P22 + Phi",
                               well %in% facp22phi ~ "Facilitation + P22 + Phi",
                               TRUE ~ well)) 


#CFP plots
CFPplot <- CFP %>%
  filter(!well %in% none) %>%
  drop_na() %>%
  mutate(normalized_CFP = CFP - min(CFP)) %>%
  ggplot(aes(x = hour, y = CFP, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition) +
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
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
        legend.title=element_text(size=12))

#YFP
YFP_path = "yfp_23july2022.csv"

#Import and Clean Data
YFP = read.table(YFP_path, sep = ",", as.is = TRUE, row.names = 1)
YFP= t(YFP)
YFP = data.frame(YFP)
names(YFP)[1:3] = c("cycle","seconds","temp")
YFP = YFP %>%
  gather(well, YFP, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60) %>%
  mutate(condition = case_when(well %in% smono ~ "S Monoculture",
                               well %in% comp ~ "Competition",
                               well %in% fac ~ "Facilitation", 
                               well %in% sphi ~ "S Monoculture + Phi",
                               well %in% compphi ~ "Competition + Phi",
                               well %in% facphi ~ "Facilitation + Phi",
                               well %in% sp22 ~ "S Monoculture + P22",
                               well %in% compp22 ~ "Competition + P22",
                               well %in% facp22 ~ "Facilitation + P22",
                               well %in% sp22phi ~ "S Monoculture + P22 + Phi",
                               well %in% compp22phi ~ "Competition + P22 + Phi",
                               well %in% facp22phi ~ "Facilitation + P22 + Phi",
                               TRUE ~ well)) 


#CFP plots
YFPplot <- YFP %>%
  filter(!well %in% none) %>%
  drop_na() %>%
  mutate(normalized_YFP = YFP - min(YFP)) %>%
  ggplot(aes(x = hour, y = YFP, colour = well)) +
  geom_smooth(span = 0.1) +
  facet_wrap(~condition) +
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
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
        legend.title=element_text(size=12))

#CFP / YFP ratio
combined <- full_join(CFP, YFP, by = c("cycle", "seconds", "well", "hour", "condition")) %>%
  drop_na() %>%
  mutate(normalized_CFP = CFP - min(CFP)) %>%
  mutate(normalized_YFP = YFP - min(YFP))%>%
  mutate(ratio = normalized_CFP / normalized_YFP) %>%
  ggplot(aes(x = hour, y = ratio, colour = well)) +
  geom_line()+
  facet_wrap(~condition, scales = "free_y") +
  theme_fivethirtyeight()+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
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
  ylab("CFP / YFP") +
  xlab("Time (hours)")



