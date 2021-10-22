library(micInt)
library(phyloseq)
library(tidyverse)
library(glue)
library(ggplot2)
library(magrittr)
library(matrixStats)
library(igraph)
library(viridis)
library(scales)
annotated_count_frame <- readRDS('annotated_count_frame.rds')

test_frame_ss_low_absolute <- annotated_count_frame %>% 
  filter(dataset == 'ss',abundance_type == 'absolute',
         sim_measure_name %in% c("spearman"),noise_level == "low",filter_threshold == 5e-4,
                                                               threshold > 10^-30, threshold < 0.05)
plotFUN <- function(frame){
  ggplot(frame) + geom_line(size = 5,mapping = aes(x=threshold,
                                          y=count_interactions,
                                          group =interaction_sign,
                                          linetype = factor(interaction_sign,levels = c("positive","negative")))
                            ) +
    scale_x_log10(limits = c(10^(-30),1), labels = trans_format("log10",math_format()), breaks = 10^ seq(-30,0,5)) +
    scale_y_continuous(labels = comma) +
    guides(linetype = guide_legend("Sign of interaction")) +
    xlab("q-value threshold") + ylab("Number of significant interactions")
}
plot <- plotFUN(test_frame_ss_low_absolute) + geom_vline(xintercept = 0.05,size = 2) + 
  annotate(geom = "text",
           label = "0.05",
           x = 0.1,
           y = test_frame_ss_low_absolute$count_interactions %>% max() %>% divide_by(2),
           angle = 90, 
           vjust = 1,
           size = 20)

print(plot)
ggsave(filename = "figures/ss_main_int_count.pdf",
       plot = plot + theme_bw() +  theme(legend.position = "none",axis.title = element_text(size = 75),
                                        axis.text = element_text(size = 55),
                                        strip.text = element_text(size = 55),
                                        axis.ticks = element_line(size = 3),
                                        axis.ticks.length = unit(1,"cm")),
       device = 'pdf',width = 60,
       height = 50, units = 'cm')
