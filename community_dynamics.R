# This script shows PCoA plots of the dynamic nature of dynamics
library(micInt)
library(phyloseq)
library(tidyverse)
library(glue)
library(ggplot2)
library(magrittr)
library(matrixStats)
library(igraph)
library(rlang)
library(cowplot)
library(viridis)
library(vegan)
source("tidy_phylo_subset.R")
if(!dir.exists("figures")){
  dir.create("figures")
}
total_compute_frame <- readRDS('total_combined_combinations.rds')
selected_row <- total_compute_frame %>% filter(dataset == 'ss',noise_level == 'low', threshold == 5e-04, abundance_type == 'absolute',
                                               sim_measure_name == 'spearman_normal')
full_physeq <- selected_row$physeq %>%
  pluck(1)
full_sample_data <- sample_data(full_physeq) %>%
  data.frame() %>%
  mutate(replicate = str_extract(Reactor.name,"[:digit:]$"))
reactor_frame <- subdivide_by_environment(full_physeq,c("selection.group","Nutrient"))

reactor_frame$reactor_sorted <- purrr::map(reactor_frame$phyloseq,
                                          ~ subdivide_by_environment(.,variables = "Reactor.name"))

reactor_frame$time_series <- reactor_frame$reactor_sorted %>% 
  purrr::map(~ (magrittr::extract2(.,'phyloseq') %>%
                  purrr::map(~ micInt::OTU_time_series(table = .,
                                                       time_points = 'Day'))
  ))
all_time_plot <- reactor_frame$time_series %>%
  unlist(recursive = FALSE) %>% plot_trajectory(linetype = "Nutrient",label = FALSE,color = "replicate",label_size = 5)
all_time_plot$data$selection.group <- full_sample_data$selection.group[match(all_time_plot$data$time_series,
                                                                             full_sample_data$Reactor.name)]
all_time_plot$data$Nutrient <- full_sample_data$Nutrient[match(all_time_plot$data$time_series,
                                                                      full_sample_data$Reactor.name)]
all_time_plot$data$replicate <- full_sample_data$replicate[match(all_time_plot$data$time_series,
                                                               full_sample_data$Reactor.name)]
all_time_plot$data$Switch <- if_else(all_time_plot$data$time_points <= 28L,true = "Before switch",false = "After switch") %>% 
  as.factor() %>%
  relevel("Before switch")
all_time_plot$data$selection.regime.at.sampling <- all_time_plot$data %$% if_else(xor(selection.group == "RK",
                                                               Switch == "After switch"),
                                                           "r",
                                                           "K")
plot <- all_time_plot  +
  geom_path(size = 5,arrow = arrow(length = unit(0.2,'inches'))) +
  geom_point(size = 8,alpha = 0.5) +
  viridis::scale_color_viridis(discrete = TRUE) +
  geom_text(aes_string(label = "time_points"), size = 12, color = "black") +
  facet_grid(selection.regime.at.sampling ~ Switch) + theme_bw() + 
  theme(plot.title = element_text(size = 18),axis.title = element_text(size = 18),axis.text = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),legend.position = 'none')
ggsave(filename = "figures/ss_main_trajectories.pdf",
       plot = plot + theme_bw() + theme(legend.position = "none",axis.title = element_text(size = 55),
                                       axis.text = element_text(size = 55),
                                       strip.text = element_text(size = 55)),
       device = 'pdf',width = 40,
       height = 50, units = 'cm')

end_physeq <- full_physeq %>% subset_samples(Day %in% c(28L,50L))
adonis_res <- adonis(refine_data(end_physeq,renormalize = FALSE) %>% t() %>% cor(method = "spearman") %>%
         {subtract(1,.)} ~
         selection.group + selection.regime.at.sampling,
       data = end_physeq %>% sample_data() %>% data.frame(),permutations = 10^6L, parallel = 4L)
sink("PERMANOVA_res.txt")
print(adonis_res)
sink()
