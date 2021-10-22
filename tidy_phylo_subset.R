# This is a neat function to prune phyloseq objects by samples by using tidyverse semantics
# and evalutation
library(magrittr)
library(phyloseq)
library(dplyr)
subset_samples_tidy <- function(physeq,...){
  physeq %>% `sample_data<-`(force(.) %>% sample_data() %>%  data.frame()  %>% filter(...))
}
