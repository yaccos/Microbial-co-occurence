library(micInt)
library(phyloseq)
library(tidyverse)
library(glue)
library(ggplot2)
library(magrittr)
library(matrixStats)
set.seed(593205)
# options(warning.length = 8170L)
data <- readRDS(file = 'total_combined_combinations.rds')
ccrepe_jobs <- pmap(data, function(sim_measures,physeq,abundance_type,...){
  refined_table <- micInt::refine_data(OTU_table = physeq,renormalize = FALSE)
  list(x=refined_table,sim.score = sim_measures %>% sim_measure_function(),
                             renormalize = if_else(abundance_type == 'relative',
                                                   true = TRUE,false = FALSE))
}
)
ccrepe_res <- micInt::ccrepe_analysis(ccrepe_job = ccrepe_jobs,commonargs = list(errthresh = 1e-4,
                                                                                iterations = 1000L,memory.optimize = TRUE),
                                      parallel = TRUE, ncpus = parallel::detectCores())
saveRDS(object = ccrepe_res,file = 'total_combined_combinations_res.rds',compress = 'xz')
