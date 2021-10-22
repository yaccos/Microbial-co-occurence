library(micInt)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(Rcpp)
sourceCpp('generate_interaction_count.cpp')
significant_res <- readRDS('significant_res.rds')
q_cutoffs <- 10^(seq(from=log10(10^-100),to = log10(5*10^-2),by = 0.1))
interacton_count_per_table <- map(significant_res,~ interaction_count(.x,q_value_thresholds = q_cutoffs))
total_count_frame <- bind_rows(interacton_count_per_table,.id = "index")
total_count_frame$index %<>% as.integer()
info_frame <- readRDS(file = 'total_combined_combinations.rds')
annotated_count_frame <- total_count_frame %>% mutate(dataset = info_frame$dataset[index],filter_threshold =info_frame$threshold[index],
                             abundance_type = info_frame$abundance_type[index],
                             noise_level = info_frame$noise_level[index],
                             sim_measure_name = info_frame$sim_measure_name[index] %>% str_replace(fixed('_normal'),''),
                             positive = significant - negative) %>%
pivot_longer(cols = c("positive","negative"),names_to = "interaction_sign",values_to = "count_interactions")
annotated_count_frame$noise_level %<>% factor(levels = c("none","low","middle","high"))
saveRDS(object = annotated_count_frame,file = 'annotated_count_frame.rds')
significant_res <- readRDS('significant_res.rds')
top_res <- map(significant_res,~ .x[seq_len(min(10000,nrow(.x))),])
saveRDS(object = top_res,file = 'top_res.rds')
