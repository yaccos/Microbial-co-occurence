library(micInt)
library(phyloseq)
library(tidyverse)
library(glue)
library(ggplot2)
library(magrittr)
library(matrixStats)
ss_physeq <- readRDS('selection_switch/ss_physeq.rds') %>% subset_samples(Reactor.name != 'Inoculum')
#  We prepare the base similarity measures to use
base_sim_measures <- similarity_measures(subset = c('pearson','spearman',
                                                    'euclidean','jaccard',
                                                    'gen_jaccard',
                                                    'cosine','bray_curtis'))
ss_relative_physeq <- microbiome::transform(ss_physeq,transform = 'compositional')
filter_thresholds <- c(low = 5e-4,middle=1e-3,high=5e-3)
microbiota_relative_list <- list(ss=ss_relative_physeq) 
analysis_frame_microbiota <- expand_grid(dataset = c('ss'),threshold = filter_thresholds)
analysis_frame_microbiota$data_filtered <- pmap(analysis_frame_microbiota,
                                                function(dataset,threshold,...) {
                                                  physeq <- microbiota_relative_list[[dataset]]
                                                  ds <- otu_table(physeq)
                                                  if(taxa_are_rows(physeq = physeq)){
                                                    OTU_max <- rowMaxs(ds)
                                                  }
                                                  else{
                                                    OTU_max <- colMaxs(ds)
                                                  }
                                                  phyloseq::prune_taxa(taxa = OTU_max > threshold,x=physeq)
                                                  })
analysis_frame_microbiota$ntaxa <- map_int(analysis_frame_microbiota$data_filtered, phyloseq::ntaxa)
analysis_frame_microbiota$proprotion_remaining_otus <- map(analysis_frame_microbiota$data_filtered,function(physeq){
  proportion_of_large_OTUs=phyloseq::sample_sums(physeq)
  ggplot()+geom_histogram(aes(proportion_of_large_OTUs),binwidth =0.01)+
    ylab('Number of samples')+xlab('Proportion of the sample made up by the kept OTUs')+
    scale_x_continuous(limits = c(0,1.05))  
}) 
analysis_frame_microbiota$proprotion_remaining_otus[[6]]

# Finds largest positive obsevation of the OTU abundances
min_dataset_FUN <- function(physeq){
  otu_table <- otu_table(physeq)
  matrix(otu_table)[otu_table > 0] %>% min()
}
analysis_frame_microbiota$min_dataset <- map_dbl(analysis_frame_microbiota$data_filtered,
                                             min_dataset_FUN)
magnitude_factors <- c(low=1,middle=10,high=100)
cell_count_column <- c(ss="events.ul")
analysis_frame_microbiota$noise_levels <- pmap(analysis_frame_microbiota, 
function(dataset,min_dataset,data_filtered,...){
  relative_physeq <- microbiome::transform(data_filtered,transform = 'compositional')
  absolute_physeq <- micInt::scale_by_column(data_filtered,cell_count_column[dataset])
  res <- tibble(abundance_type = c('relative','absolute'),physeq = list(relative_physeq,absolute_physeq))
  res$min_dataset <- map_dbl(res$physeq,min_dataset_FUN)
  res$sim_measures <- map(res$min_dataset,function(min_dataset){
    noise_levels <- magnitude_factors*min_dataset
    noised_similarity_measures <- map(noise_levels, ~micInt::noisify(sim.scores = base_sim_measures,
                                                                    magnitude = .x,noise = "normal"))
    
    return(c(list(none=base_sim_measures),noised_similarity_measures))
  })
  return(res)
}
)
microbiota_compute_frame <- analysis_frame_microbiota %>%
  select(dataset,threshold,noise_levels) %>%
  unnest(noise_levels) %>% 
  select(-min_dataset) %>% 
  mutate(noise_level = map(sim_measures,names)) %>% 
  unnest(sim_measures,noise_level) %>% 
  mutate(sim_measure_name = map(sim_measures,names)) %>%
  unnest(sim_measure_name,sim_measures)
total_compute_frame <- rbind(microbiota_compute_frame)
saveRDS(object = total_compute_frame, file = 'total_combined_combinations.rds',compress = 'xz')
