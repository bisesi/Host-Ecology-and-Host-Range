#ATB
#Tecan data figures, competition and cooperation
#Figure 4

#library
library("patchwork")
library("cowplot")

#set date for representative experimental run used in visualizations
date <- "8March2023"

#source tecan data cleaning script
source(here::here("data-generation", "experimental", "tecan-data-cleaning.R"))

#part A - growth of each phage on respective hosts
date <- "9May2023"
pfu_path <- here::here("experimental-data", "tecan-data", date, generate_paths("pfu", date, "xlsx"))
pfus_partA <- load_pfu_data(pfu_path) %>% janitor::clean_names() %>% mutate(well = c(1:28))
cleaned_pfus <- clean_pfu_data(pfus_partA)

partA <- cleaned_pfus %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "*E. coli*\nmonoculture",
                                 interaction == "S Monoculture" ~ "*S. enterica*\nmonoculture")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  filter(phage != "both\nphage") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-5,15)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_markdown(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")

comp_mut <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>% select(doublings, doubling_type, interaction, phage, phage_type)

alternative_partA <- cleaned_pfus %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "*E. coli*<br>monoculture",
                                 interaction == "S Monoculture" ~ "*S. enterica*<br>monoculture")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  filter(phage != "both\nphage") %>%
  select(doublings, doubling_type, interaction, phage, phage_type) %>% 
  rbind(., comp_mut) %>% 
  filter(phage != "both\nphage") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction, ncol = 4) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-7, 15)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text= element_markdown(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")

#part B - coop and comp phage densities growth rates
partB <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-7, 15)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")

alternative_partB <- pfus_and_final_density %>%
  mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                               TRUE ~ doublings)) %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(doublings = case_when(doublings == min(doublings) ~ -5,
                               TRUE ~ doublings)) %>%
  filter(phage == "both\nphage") %>%
  ggplot(aes(x = phage, y = doublings, color = phage_type)) +
  facet_wrap(~interaction) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  theme_bw(base_size = 18)+
  ylim(-7, 15)+
  scale_color_manual(values = c(eh7, p22vir))+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("ln(final pfu / initial pfu)")+
  labs(color = "species")

legend <- get_legend(pfus_and_final_density %>%
                       mutate(doublings = case_when(doublings == 0.0 ~ Inf,
                                                    TRUE ~ doublings)) %>%
                       filter(interaction == "Mutualism" | interaction == "Competition") %>%
                       mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                                      interaction == "Competition" ~ "competition"))%>%
                       mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "Generalist (EH7)",
                                                     phage_type == "Specialist phage" ~ "Specialist (P22vir)"))%>%
                       mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                                                phage == "Phi" ~ "generalist\nonly",
                                                phage == "Phi + P22" ~ "both\nphage",
                                                phage == "none" ~ "no\nphage")) %>%
                       mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                                               "both\nphage"))) %>%
                       ggplot(aes(x = phage, y = doublings, color = phage_type)) +
                       facet_wrap(~interaction) +
                       geom_boxplot() +
                       geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
                       theme_bw(base_size = 18)+
                       scale_color_manual(values = c(eh7, p22vir))+
                       theme(axis.title = element_text(), 
                             panel.background = element_rect(fill = "white"), 
                             plot.background = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "bottom",
                             legend.background = element_blank(),
                             strip.background = element_blank())+
                       ylab("ln(final pfu / initial pfu)")+
                       xlab("treatment")+
                       labs(color = "species"))


#part C - competition OD and final densities
partC_top <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition")%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist only",
                           phage == "Phi" ~ "generalist only",
                           phage == "Phi + P22" ~ "both phage",
                           phage == "none" ~ "no phage")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                          "both phage"))) %>%
  filter(!well %in% c("B8", "B9")) %>%
  pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
  mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                             fluor == "S_corrected_OD" ~ "S. enterica")) %>%
  ggplot(aes(x = hours, y = OD, color = species))+
  geom_smooth(span = 0.2, size = 1.5, 
              aes(ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y)))))+
  scale_color_manual(values = c(ecoli, senterica))+
  facet_wrap(~phage, ncol = 4)+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  xlab("hours")+
  ylim(0, 0.4)+
  labs(color = "species")
  #geom_vline(data = (data.frame(xint=20, interaction="competition", phage = "no phage") %>% mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                                                                                    #"both phage")))), 
             #aes(xintercept = xint), color = "red", linetype = "dashed")

partC_bottom <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  mutate(S_percentage = S_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  pivot_longer(E_percentage:S_percentage, names_to = "bacteria", values_to = "density")%>%
  group_by(bacteria, phage) %>%
  summarise(mean= mean(density), sd = sd(density)) %>%
  ggplot(aes(x = phage, y = mean, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  scale_fill_manual(values = c(ecoli, senterica))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("percent of final co-culture")+
  labs(color = "interaction")+
  ylim(0, 100)

partC <- partC_top / partC_bottom

#part D - OD600 of species coop or comp
partD_top <- all_tecan_adjusted_OD %>%
  mutate(time = cycle * 23, hours = (time / 60)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism")%>%
  mutate(phage = case_when(phage == "P22" ~ "specialist only",
                           phage == "Phi" ~ "generalist only",
                           phage == "Phi + P22" ~ "both phage",
                           phage == "none" ~ "no phage")) %>%
  mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                          "both phage"))) %>%
  filter(!well %in% c("B8", "B9")) %>%
  pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
  mutate(species = case_when(fluor == "E_corrected_OD" ~ "E. coli",
                             fluor == "S_corrected_OD" ~ "S. enterica")) %>%
  ggplot(aes(x = hours, y = OD, color = species))+
  geom_smooth(span = 0.2, size = 1.5, aes(ymax = after_stat(y + se * sqrt(length(y))),
                              ymin = after_stat(y - se * sqrt(length(y)))))+
  scale_color_manual(values = c(ecoli, senterica))+
  facet_wrap(~phage, ncol = 4)+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("OD")+
  ylim(0, 0.4)+
  xlab("hours")+
  labs(color = "species")
  
#geom_vline(data = (data.frame(xint=47, interaction="mutualism", phage = "no phage") %>% mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                                                                                    #"both phage")))), 
             #aes(xintercept = xint), color = "red", linetype = "dashed")

partD_bottom <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "mutualism") %>%
  mutate(phage = case_when(phage == "P22" ~ "specialist\nonly",
                           phage == "Phi" ~ "generalist\nonly",
                           phage == "Phi + P22" ~ "both\nphage",
                           phage == "none" ~ "no\nphage")) %>%
  mutate(phage = factor(phage, levels = c("no\nphage", "specialist\nonly", "generalist\nonly",
                                          "both\nphage"))) %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  mutate(S_percentage = S_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  pivot_longer(E_percentage:S_percentage, names_to = "bacteria", values_to = "density")%>%
  group_by(bacteria, phage) %>%
  summarise(mean= mean(density), sd = sd(density)) %>%
  ggplot(aes(x = phage, y = mean, fill = bacteria))+
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75))+
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), width=.2, position=position_dodge(width=0.75))+
  scale_fill_manual(values = c(ecoli, senterica))+
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(), 
        panel.background = element_rect(fill = "white"), 
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())+
  ylab("percent of final co-culture")+
  labs(color = "interaction")+
  ylim(0, 100)

partD <- partD_top / partD_bottom

legend2 <- get_legend(all_tecan_adjusted_OD %>%
                        mutate(time = cycle * 23, hours = (time / 60)) %>%
                        mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                                       interaction == "Comp" ~ "competition",
                                                       TRUE ~ interaction)) %>%
                        filter(interaction == "mutualism")%>%
                        mutate(phage = case_when(phage == "P22" ~ "specialist only",
                                                 phage == "Phi" ~ "generalist only",
                                                 phage == "Phi + P22" ~ "both phage",
                                                 phage == "none" ~ "no phage")) %>%
                        mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                "both phage"))) %>%
                        filter(!well %in% c("B8", "B9")) %>%
                        pivot_longer(cols = E_corrected_OD:S_corrected_OD, names_to = "fluor", values_to = "OD") %>%
                        mutate(species = case_when(fluor == "E_corrected_OD" ~ "*E. coli*",
                                                   fluor == "S_corrected_OD" ~ "*S. enterica*")) %>%
                        ggplot(aes(x = hours, y = OD, color = species))+
                        geom_smooth(span = 0.2, se = FALSE)+
                        scale_color_manual(values = c(ecoli, senterica))+
                        facet_wrap(~phage, ncol = 4)+
                        theme_bw(base_size = 18)+
                        theme(axis.title = element_text(), 
                              panel.background = element_rect(fill = "white"), 
                              plot.background = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "bottom",
                              legend.text = element_markdown(),
                              legend.background = element_blank(),
                              strip.background = element_blank())+
                        ylab("OD")+
                        xlab("hours")+
                        labs(color = "species"))
                        #geom_vline(data = (data.frame(xint=47, interaction="mutualism", phage = "no phage") %>% mutate(phage = factor(phage, levels = c("no phage", "specialist only", "generalist only",
                                                                                                                                                        #"both phage")))), 
                                   #aes(xintercept = xint), color = "red", linetype = "dashed"))

#fig 4
fig4 <- plot_grid(plot_grid(partA, partB, legend, labels = c("A", "B"), label_size = 26, rel_heights = c(0.25,0.25, 0.15), ncol = 1),
                  plot_grid(partC, partD, legend2, labels = c("C", "D"), label_size = 26, rel_heights = c(1,1, 0.15), ncol = 1))

alt_fig4 <- plot_grid(plot_grid(alternative_partA, alternative_partB, legend, labels = c("A", "B"), label_size = 26, rel_heights = c(0.3,0.3, 0.1), ncol = 1),
                  plot_grid(partC, partD, legend2, labels = c("C", "D"), label_size = 26, rel_heights = c(1,1, 0.15), ncol = 1), rel_widths = c(0.75, 0.9))



