#ATB
#Morris and Sobol Screening Visualizations

#load packages and data
library("reshape2")
library("patchwork")
library("ggcharts")
source(here::here("data-generation", "model", "exploratory", "sobol-morris.R"))

#competition
morris_final_comp <- bar_chart(gis_morris_comp, type, gi, facet = name) + 
  ylab("Global Morris Index at Final Timepoint") +
  xlab("Parameter")

morris_avg_comp <- bar_chart(gis_morris_comp, type, avg_gi, facet = name) + 
  ylab("Global Morris Index (Mean)") +
  xlab("Parameter")

sobol_first_order_comp <- sobol_comp %>%
  select(-c(time)) %>%
  filter(index == "S") %>%
  select(-c(index)) %>%
  filter(name %in% c("gen", "sp")) %>%
  melt(., id = c("name", "timepoint")) %>%
  ggplot(aes(x = as.numeric(timepoint), y = as.numeric(value), color = variable))+
  geom_line(size = 1.5)+
  facet_wrap(~name)+
  ylab("First Order") +
  xlab("Time")

sobol_total_comp <- sobol_comp %>%
  select(-c(time)) %>%
  filter(index == "T") %>%
  select(-c(index)) %>%
  filter(name %in% c("gen", "sp")) %>%
  melt(., id = c("name", "timepoint")) %>%
  ggplot(aes(x = as.numeric(timepoint), y = as.numeric(value), color = variable))+
  geom_line(size = 1.5)+
  facet_wrap(~name)+
  ylab("Total Order") +
  xlab("Time")

comp_figs <- (morris_final_comp + morris_avg_comp) / (sobol_first_order_comp + sobol_total_comp) + plot_layout(guides = "collect")

#cooperation
morris_final_coop <- bar_chart(gis_morris_coop, type, gi, facet = name) + 
  ylab("Global Morris Index at Final Timepoint") +
  xlab("Parameter")

morris_avg_coop <- bar_chart(gis_morris_coop, type, avg_gi, facet = name) + 
  ylab("Global Morris Index (Mean)") +
  xlab("Parameter")

sobol_first_order_coop <- sobol_coop %>%
  select(-c(time)) %>%
  filter(index == "S") %>%
  select(-c(index)) %>%
  filter(name %in% c("gen", "sp")) %>%
  melt(., id = c("name", "timepoint")) %>%
  ggplot(aes(x = as.numeric(timepoint), y = as.numeric(value), color = variable))+
  geom_line(size = 1.5)+
  facet_wrap(~name)+
  ylab("First Order") +
  xlab("Time")

sobol_total_coop <- sobol_coop %>%
  select(-c(time)) %>%
  filter(index == "T") %>%
  select(-c(index)) %>%
  filter(name %in% c("gen", "sp")) %>%
  melt(., id = c("name", "timepoint")) %>%
  ggplot(aes(x = as.numeric(timepoint), y = as.numeric(value), color = variable))+
  geom_line(size = 1.5)+
  facet_wrap(~name)+
  ylab("Total Order") +
  xlab("Time")

coop_figs <- (morris_final_coop + morris_avg_coop) / (sobol_first_order_coop + sobol_total_coop) + plot_layout(guides = "collect")

