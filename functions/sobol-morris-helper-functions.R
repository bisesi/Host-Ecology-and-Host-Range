#ATB
#Functions
#Sobol Morris Helper Functions
#Clean outputs

#required packages
library("tidyverse")

clean_morris_output <- function(input){
  new_list <- list()
  list_names <- names(input)
  for (i in 1:length(input)){
    new_list[[i]] <- input[[i]] %>% t() %>% data.frame() %>% 
      mutate(name = list_names[i]) %>%
      mutate(timepoint = c(1:nrow(input[[i]] %>% t() %>% data.frame())))
    rownames(new_list[[i]]) <- NULL
  }
  complete_morris <- bind_rows(new_list) %>% select(-c(time))
  return(complete_morris)
}

calculate_morris_gi <- function(input, parameters){
  cleaned_morris <- data.frame()
  for (i in parameters){
    new_dat <- input %>% select(c(grep(i, colnames(input)), name, timepoint)) %>%
      data.frame() %>% mutate(type = i) %>% rename_with(., ~gsub(paste0("_",i), "", .x, fixed = TRUE))
    cleaned_morris <- bind_rows(cleaned_morris, new_dat)
  }
  output<- cleaned_morris %>% 
    mutate(gi = sqrt(mu.star^2 + sigma^2)) %>% group_by(name, type) %>% mutate(avg_gi = mean(gi)) %>%
    slice_max(timepoint) %>%
    filter(name %in% c("gen", "sp"))
  return(output)
}

clean_sobol_output <- function(input) {
  new_list_S <- list()
  list_names <- names(input)
  for (i in 1:length(input)){
    new_list_S[[i]] <- input[[i]][1] %>% data.frame() %>% t() %>% 
      data.frame() %>% mutate(name = list_names[i]) %>% mutate(index = "S") %>% 
      mutate(timepoint = c(1:nrow(input[[i]][1] %>% data.frame() %>% t() %>% 
                                    data.frame())))
    rownames(new_list_S[[i]]) <- NULL
  }
  new_list_T <- list()
  for (i in 1:length(input)){
    new_list_T[[i]] <- input[[i]][2] %>% data.frame() %>% t() %>% 
      data.frame() %>% mutate(name = list_names[i]) %>% mutate(index = "T") %>% 
      mutate(timepoint = c(1:nrow(input[[i]][2] %>% data.frame() %>% t() %>% 
                                    data.frame())))
    rownames(new_list_T[[i]]) <- NULL
  }
  
  sobol_cleaned <- rbind(bind_rows(new_list_S), bind_rows(new_list_T))
  return(sobol_cleaned)
}


