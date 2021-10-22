library(micInt)
library(phyloseq)
library(tidyverse)
library(glue)
library(ggplot2)
library(magrittr)
library(matrixStats)
total_compute_frame <- readRDS('total_combined_combinations.rds')
res <- readRDS('total_combined_combinations_res.rds')
total_compute_frame$res <- res
total_compute_frame$taxonomy_string <- map(total_compute_frame$physeq,micInt::collapse_taxonomy)
total_compute_frame$significant_res <- pmap(total_compute_frame %$% list(res,taxonomy_string,sim_measures),
                                                                                 ~ micInt::create_interaction_table(data = ..1,taxonomy = ..2,
                                                                                                                    threshold.type = 'q',threshold.value = 0.01,
                                                                                                                    score_attributes = micInt::sim.measure.attributes(..3)))
saveRDS(object = total_compute_frame$significant_res,file = 'significant_res.rds',compress = 'xz')
saveRDS(object = total_compute_frame,file = 'total_added_results.rds',compress = FALSE)
