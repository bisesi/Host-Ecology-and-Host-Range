#ATB
#Generate statistics
#Figures 4A, 4B, 4C, 4D, 5A, 5B, Supp 3A, 3B

#load packages
library("openxlsx")
library("ggpubr")
library("rstatix")

#generate all data for figure 4
source(here::here("figures", "final-figs", "scripts", "figure-4.R"))

partA_data_fig4 <- cleaned_pfus %>% select(interaction, phage, doublings) %>% 
  filter(phage != "Phi + P22") %>%
  filter(doublings != 0) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partB_data_fig4 <- pfus_and_final_density %>%
  filter(interaction == "Mutualism" | interaction == "Competition") %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (EH7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)")) %>%
  mutate(doubling_type = case_when(doubling_type == "phi_doublings" ~ "eh7_doublings",
                                   TRUE ~ doubling_type)) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partC_D_data_fig4 <- pfus_and_final_density %>%
  mutate(interaction = case_when(interaction == "Mutualism" ~ "mutualism",
                                 interaction == "Competition" ~ "competition"))%>%
  select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage) %>%
  rbind(., all_tecan_adjusted_OD %>% filter(interaction == "Coop" | interaction == "Comp") %>% 
          filter(phage == "none") %>% slice_max(cycle) %>%
          select(interaction, E_corrected_OD, S_corrected_OD, interaction, phage)) %>%
  mutate(interaction = case_when(interaction == "Coop" ~ "mutualism",
                                 interaction == "Comp" ~ "competition",
                                 TRUE ~ interaction)) %>%
  filter(interaction == "competition" | interaction == "mutualism") %>%
  mutate(E_percentage = E_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  mutate(S_percentage = S_corrected_OD / (E_corrected_OD + S_corrected_OD) * 100) %>%
  pivot_longer(E_percentage:S_percentage, names_to = "bacteria", values_to = "density") %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partA_anova_fig4 <- aov(doublings ~ interaction * phage, data = partA_data_fig4)
partA_multiple_comparisons_fig4 <- TukeyHSD(partA_anova_fig4, conf.levels = 0.95)

partB_anova_fig4 <- aov(doublings ~ interaction * phage * doubling_type, data = partB_data_fig4)
partB_multiple_comparisons_fig4 <- TukeyHSD(partB_anova_fig4, conf.levels = 0.95)

partC_D_anova_fig4 <- aov(density ~ interaction * phage * bacteria, data = partC_D_data_fig4)
partC_D_multiple_comparisons_fig4 <- TukeyHSD(partC_D_anova_fig4, conf.levels = 0.95)

supplemental_table_eight <- list('partA_anova' = partA_anova_fig4, 'partA_mult_comp' = tibble::rownames_to_column(data.frame(partA_multiple_comparisons_fig4$`interaction:phage`)),
                                 'partB_anova' = partB_anova_fig4, 'partB_mult_comp' = tibble::rownames_to_column(data.frame(partB_multiple_comparisons_fig4$`interaction:phage:doubling_type`)),
                                 'partC_D_anova' = partC_D_anova_fig4, 'partC_D_mult_comp' = tibble::rownames_to_column(data.frame(partC_D_multiple_comparisons_fig4$`interaction:phage:bacteria`)))
write.xlsx(supplemental_table_eight, file = here::here("paper-supplement", "tables", "supplemental-table-8.xlsx"))

#generate data for figure 5
rm(list = ls())
source(here::here("figures", "final-figs", "scripts", "figure-5.R"))

partA_data_fig5 <- cleaned_pfus %>% mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (EH7)",
                                                                  phage_type == "Specialist phage" ~ "specialist (p22vir)"))%>%
  mutate(interaction = case_when(interaction == "E Monoculture" ~ "starved e monoculture",
                                 interaction == "S Monoculture" ~ "starved s monoculture",
                                 interaction == "No Cells" ~ "no cells",
                                 TRUE ~ interaction))%>%
  filter(interaction %in% c("starved e monoculture",
                            "starved s monoculture",
                            "no cells")) %>%
  filter(phage != "Phi + P22") %>%
  mutate(doubling_type = case_when(doubling_type == "phi_doublings" ~ "eh7_doublings",
                                   TRUE ~ doubling_type)) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))


partB_data_fig5 <-all_data_phage %>%
  filter(!interaction %in% c("Facilitation", "Mutualism"))%>%
  mutate(interaction = case_when(interaction == "S Monoculture" ~ "s monoculture",
                                 interaction == "Competition" ~ "competition")) %>%
  mutate(phage_type = case_when(phage_type == "Generalist phage" ~ "generalist (eh7)",
                                phage_type == "Specialist phage" ~ "specialist (p22vir)")) %>%
  filter(phage == "Phi + P22") %>%
  mutate(doubling_type = case_when(doubling_type == "phi_doublings" ~ "eh7_doublings",
                                   TRUE ~ doubling_type)) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partA_anova_fig5 <- aov(doublings ~ interaction * phage, data = partA_data_fig5)
partA_multiple_comparisons_fig5 <- TukeyHSD(partA_anova_fig5, conf.levels = 0.95)

partB_anova_fig5 <- aov(doublings ~ interaction * cost * doubling_type, data = partB_data_fig5)
partB_multiple_comparisons_fig5 <- TukeyHSD(partB_anova_fig5, conf.levels = 0.95)

supplemental_table_nine <- list('partA_anova' = partA_anova_fig5, 'partA_mult_comp' = tibble::rownames_to_column(data.frame(partA_multiple_comparisons_fig5$`interaction:phage`)),
                                 'partB_anova' = partB_anova_fig5, 'partB_mult_comp' = tibble::rownames_to_column(data.frame(partB_multiple_comparisons_fig5$`interaction:cost:doubling_type`)))
write.xlsx(supplemental_table_nine, file = here::here("paper-supplement", "tables", "supplemental-table-9.xlsx"))

#generate stats for supplemental figure 3
rm(list = ls())
source(here::here("figures", "final-figs", "scripts", "supplemental-figure-3.R"))

partA_data_suppfig3 <- no_cells %>%
  mutate(interaction = case_when(interaction == "No cells" ~ "no cells",
                                 TRUE ~ interaction))%>%
  mutate(timepoint = case_when(timepoint == 24 ~ "hour 24",
                               timepoint == 48 ~ "hour 48",
                               TRUE ~ timepoint)) %>%
  mutate(doubling_type = case_when(doubling_type == "phi_doublings" ~ "eh7_doublings",
                                   TRUE ~ doubling_type)) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partB_panel1_data_suppfig3 <- partB_data %>%
  rename(eh7_doublings = phi_doublings) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partB_panel2_data_suppfig3 <- partC_data %>%
  rename(eh7_doublings = phi_doublings) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partB_panel3_data_suppfig3 <- partD_data %>% mutate(doublings = case_when(is.infinite(phi_doublings) ~ 0,
                                                                          TRUE ~ phi_doublings)) %>%
  rename(eh7_doublings = phi_doublings) %>%
  mutate(phage = case_when(phage == "Phi" ~ "EH7",
                           phage == "Phi + P22" ~ "EH7 + P22",
                           TRUE ~ phage))

partA_anova_suppfig3 <- aov(doublings ~ doubling_type * timepoint * phage, data = partA_data_suppfig3)
partA_multiple_comparisons_suppfig3 <- TukeyHSD(partA_anova_suppfig3, conf.levels = 0.95)

partB_panel1_ttest_suppfig3 <- t.test(eh7_doublings ~ media, data = partB_panel1_data_suppfig3)

partB_panel2_anova_suppfig3 <- aov(eh7_doublings ~ media, data = partB_panel2_data_suppfig3)
partB_panel2_multiple_comparisons_suppfig3 <- TukeyHSD(partB_panel2_anova_suppfig3, conf.levels = 0.95)

partB_panel3_anova_suppfig3 <- aov(doublings ~ media, data = partB_panel3_data_suppfig3)
partB_panel3_multiple_comparisons_suppfig3 <- TukeyHSD(partB_panel3_anova_suppfig3, conf.levels = 0.95)

supplemental_table_ten <- list('partA_anova' = partA_anova_suppfig3, 'partA_mult_comp' = tibble::rownames_to_column(data.frame(partA_multiple_comparisons_suppfig3$`doubling_type:timepoint:phage`)),
                                'partB_panel1' = data.frame(stat = partB_panel1_ttest_suppfig3$statistic, p.val = partB_panel1_ttest_suppfig3$p.value, 
                                                         parameter = partB_panel1_ttest_suppfig3$parameter, stderr = partB_panel1_ttest_suppfig3$stderr,
                                                         alt = partB_panel1_ttest_suppfig3$alternative, method = partB_panel1_ttest_suppfig3$method), 
                               'partB_panel2_anova' = partB_panel2_anova_suppfig3, 'partB_panel2_mult_comp' = tibble::rownames_to_column(data.frame(partB_panel2_multiple_comparisons_suppfig3$media)),
                               'partB_panel3_anova' = partB_panel3_anova_suppfig3, 'partB_panel3_mult_comp' = tibble::rownames_to_column(data.frame(partB_panel3_multiple_comparisons_suppfig3$media)))
write.xlsx(supplemental_table_ten, file = here::here("paper-supplement", "tables", "supplemental-table-10.xlsx"))

