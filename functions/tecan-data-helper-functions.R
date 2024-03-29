#ATB
#Tecan and plaque assay data analysis
#Helper functions
#Loads data, generates growth rate data, cleans PFU files, creates species-specific ODs

generate_paths <- function(tecan_file_type, date, ending){
  return(paste0(tecan_file_type, "_", date, ".", ending))
}

import_tecan_data <- function(path, file_type = "OD"){
  data <- read.table(path, sep = ",", as.is = TRUE, row.names = 1)
  data <- t(data)
  data <- data.frame(data)
  names(data)[1:3] = c("cycle","seconds","temp")
  data <- data %>%
    gather(well, OD, -cycle, -seconds, -temp) %>%
    mutate(hour = seconds / 60 / 60)
  if (file_type == "CFP"){
    data <- data %>% rename(CFP = OD)
  }
  if (file_type == "YFP"){
    data <- data %>% rename(YFP = OD)
  }
  if (file_type == "RFP"){
    data <- data %>% rename(RFP = OD)
  }
  return(data)
}

generate_baranyi_growth_data_from_OD <- function(file){
  output <- file %>%
    group_by(well) %>%
    summarize(model_fit = fit_baranyi(hour, log(OD), 5),
              fit_variable = c("growth_rate", "lag", "ymax", "y0"))
  return(output)
}

adjust_OD_vs_FP <- function(input, itx, fluorescence){
  output <- input %>% 
    filter(interaction == itx & phage == "none") %>% 
    group_by(well) %>% 
    mutate(adjusted = {{fluorescence}} - min({{fluorescence}})) %>% 
    ungroup() %>% 
    select(c(cycle, well, adjusted, OD))
  cycle_max <- input %>% select(cycle) %>% max()
  best_fit <- c()
  for (i in seq(5, cycle_max, by = 5)){
    best_fit[i] <- summary(lm(OD ~ adjusted, data = output %>% filter(cycle < i)))$r.squared
  }
  return(c(which.max(best_fit), summary(lm(OD ~ adjusted, data = output))$coefficients[2,1]))
}

adjust_FP_vs_FP <- function(input, itx, fluorescence1, fluorescence2, cycle_number){
 output <- input %>% 
    filter(interaction == itx & phage == "none") %>% 
    group_by(well) %>% 
    mutate(adjusted1 = {{fluorescence1}} - min({{fluorescence1}})) %>%
    ungroup() %>% 
    mutate(adjusted2 = {{fluorescence2}} - min({{fluorescence2}})) %>%
    ungroup() %>% 
    select(c(cycle, well, adjusted1, adjusted2)) %>%
    filter(cycle < cycle_number)
  return(summary(lm(adjusted1 ~ adjusted2, data = output))$coefficients[2,1])
}

adjust_FP_values <- function(input, YFP_bleed_value, CFP_bleed_value){
  wells <- c(input %>% filter(interaction != "none") %>% select(well) %>% unique())$well
  total_data <- input %>% select(c(cycle,well,YFP,CFP)) %>% pivot_wider(names_from = well, values_from = c("CFP","YFP"))
  store <- list()
  for (i in wells){
    j <- which(wells == i)
    new_matrix <- total_data %>% select(matches(i))
    dat <- apply(new_matrix, 1, function(x) solve(matrix(c(1,YFP_bleed_value, CFP_bleed_value,1), nrow=2, byrow = TRUE), c(x))) %>%
      as.data.frame() %>% 
      t() 
    colnames(dat) <- c("adjusted_CFP", "adjusted_YFP")
    store[[j]] <- dat
  }
  names(store) <- wells
  output <- map_df(store, ~as.data.frame(.x), .id="well") %>%
    group_by(well) %>%
    dplyr::mutate(cycle = c(1:n())) %>%
    ungroup()
  return(output)
}

load_pfu_data <- function(path){
  sheetnames <- excel_sheets(path)
  datalist <- lapply(sheetnames, read_excel, path = path)
  output <- datalist[[1]]
  return(output)
}

clean_pfu_data <- function(input){
  wide_data <- input %>% 
    select(condition, interaction, phage, plate, pfu, well) %>%
    pivot_wider(names_from = plate, values_from = pfu) %>%
    rename(pfu_S = S, pfu_E = E)
  starting_phage <- wide_data %>%
    filter(interaction == "None" & condition == "Start") %>%
    mutate(ratio = pfu_E / (pfu_E + pfu_S)) 
  output <- wide_data %>%
    filter(condition != "Start") %>%
    mutate(phi_doublings = case_when(phage == "Phi" ~ log(pfu_E / starting_phage %>% filter(phage == "Phi") %>% select(pfu_E) %>% pull),
                                     phage == "Phi + P22" ~ log(pfu_E / starting_phage %>% filter(phage == "Phi + P22") %>% select(pfu_E) %>% pull),
                                     phage == "P22" ~ 0,
                                     TRUE ~ 0)) %>%    
    mutate(p22_doublings = case_when(phage == "P22" ~ log(pfu_S / starting_phage %>% filter(phage == "P22") %>% select(pfu_S) %>% pull),
                                     phage == "Phi + P22" ~ log(pfu_S / starting_phage %>% filter(phage == "Phi + P22") %>% select(pfu_S) %>% pull),
                                     phage == "Phi" ~ 0,
                                     TRUE ~ 0)) %>%
    mutate(generalist_fitness = (((pfu_E) / (pfu_E + pfu_S)) / starting_phage %>% filter(phage == "Phi + P22") %>% select(ratio) %>% pull)) %>% 
    mutate(change_in_percent_generalist = (((pfu_E) / (pfu_E + pfu_S)) - starting_phage %>% filter(phage == "Phi + P22") %>% select(ratio) %>% pull)) %>%
    mutate_all(~replace(., is.nan(.), 0)) %>%
    mutate_all(~replace(., is.infinite(.), log(1 / starting_phage %>% filter(phage == "Phi + P22") %>% select(pfu_E) %>% pull))) %>%
    mutate(interaction = case_when(interaction == "Smono" ~ "S Monoculture",
                                   interaction == "Emono" ~ "E Monoculture",
                                   interaction == "E Fac" ~ "E Facilitation",
                                   interaction == "Coop" ~ "Mutualism",
                                   interaction == "Comp" ~ "Competition",
                                   interaction == "Fac" ~ "Facilitation",
                                   interaction == "None" ~ "No cells")) %>%
    pivot_longer(cols = ends_with("doublings"),
                 names_to = "doubling_type",
                 values_to = "doublings") %>%
    mutate(phage_type = case_when(doubling_type == "p22_doublings" ~ "Specialist phage",
                                  doubling_type == "phi_doublings" ~ "Generalist phage",
                                  TRUE ~ "NA")) %>%
    mutate(phage_interaction = case_when(phage == "P22" ~ "Specialist phage",
                             phage == "Phi" ~ "Generalist phage",
                             phage == "Phi + P22" ~ "Phage Competition"))
  return(output)
}
