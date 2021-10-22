library(micInt)
library(phyloseq)
library(tidyverse)
library(igraph)
library(magrittr)
library(glue)
library(viridis)
library(xtable)
library(ape)
library(kableExtra)
library(cowplot)
source("dynamic_networks_utils.R")
source("tidy_phylo_subset.R")
# The purpose of this script is to identify and explore interaction Modules (network communities) in the
# ReBoot interaction networks for the selection-switch dataset
options(stringsAsFactors = FALSE)
MAX_EDGES= 500
top_res <- readRDS('top_res.rds')
total_compute_frame <- readRDS('total_combined_combinations.rds')
total_compute_frame$top_res <- top_res
selected_row <- total_compute_frame %>% filter(dataset == 'ss',noise_level == 'low', threshold == 5e-04, abundance_type == 'absolute',
                                               sim_measure_name == 'spearman_normal')
reference_physeq <- total_compute_frame %>% filter(dataset == 'ss',noise_level == 'low', threshold == 5e-04, abundance_type == 'relative',
                                                   sim_measure_name == 'spearman_normal') %>% pluck("physeq") %>% pluck(1)
full_table <- selected_row %>% pluck('top_res') %>% pluck(1)
relevant_physeq <- selected_row %>% pluck("physeq") %>% pluck(1)
full_table %>% nrow()
full_table %>% filter(q.value < 1e-4) %>% nrow()
full_table %>% filter(q.value < 1e-10) %>% nrow()
full_table$q.value[MAX_EDGES+1]
table <- full_table %>% extract(seq_len(min(MAX_EDGES,nrow(.))),)
# We remove the negative interactions
table_positive = table %>%  subset(sim.score > 0)
graph=table_positive %>% micInt::as.edgelist() %>%
  igraph::graph_from_edgelist(directed = FALSE) %>% set_edge_attr('color',value='blue') %>% 
  set_edge_attr('lty',value = 1L)
# We are only interrested in connected components with a certain number of members
graph_components <- components(graph)
components_to_keep <- which(graph_components$csize > 10)
nodes_to_keep <- graph_components$membership %in% components_to_keep %>% which()
graph_reduced <- induced_subgraph(graph = graph,vids = nodes_to_keep)
# Finds the communities
# Keep this seed!
set.seed(2452)
components=graph_reduced %>% igraph::walktrap.community(steps = 20)
n_communities <- components$membership %>% unique() %>% length()
V(graph_reduced)$color <- viridis(n_communities)[components$membership]
V(graph_reduced)$membership <- components$membership
#  Adds edges of negative interactions
negative_graph = table %>% subset(sim.score < 0) %>% micInt::as.edgelist() %>% 
  igraph::graph_from_edgelist(directed = FALSE) %>% set_edge_attr('color',value='red') %>%
  set_edge_attr('lty',value = 2L)
new_graph <- union2(graph_reduced,negative_graph)
all_displayed_taxa <- V(new_graph) %>% names()
relevant_physeq <- selected_row %>% pluck('physeq') %>% pluck(1) %>% {phyloseq::prune_taxa(all_displayed_taxa,x=.)}
reduced_physeq <- relevant_physeq
t_table <- tax_table(reduced_physeq)
class_info <- t_table[,"Class"]
unique_classes <- class_info %>% unique()
pch_palette <- c(15:18,8)
# V(new_graph)$shape <- pch_palette[match(class_info,unique_classes)]
V(new_graph)$class_taxonomy <- class_info[V(new_graph)]
stat <- OTU_stats(reference_physeq)
stat_community <- merge(stat,data.frame(ID=components$names,community=components$membership %>% as.integer()),by = "ID")
node_mean_abundance <- stat[V(new_graph),"meanAbundance"]
layout <- layout_nicely(new_graph)
MIN_NODE_SIZE <- 5
MAX_NODE_SIZE <- 20
V(new_graph)$size <- 5
V(new_graph)$size <- MIN_NODE_SIZE +
  (MAX_NODE_SIZE - MIN_NODE_SIZE) *( -1 / log10(node_mean_abundance) + 1 / log10(node_mean_abundance %>% min()))
# shape_palette <- shapes() %>% {Filter(function(x) x %in% c('circle','square','rectangle','sphere','vreactangle','pie'),.)}
# V(new_graph)$shape <- shape_palette[match(class_info,unique_classes)]
# dev.off()
phylum_taxonomy <- tax_table(relevant_physeq)[all_displayed_taxa,"Phylum"] %>% as.factor()
# par(cex.main=1.2)
# layout[,1] %<>% add(100)
color_palette <- viridis(n_communities)
# createNetworkFromIgraph(new_graph,"ss_modular_graph")
message(glue("{V(new_graph) %>% length()} nodes in the graph"))
pdf('figures/ss_community_network.pdf',width = 20/cm(1),height = 22/cm(1))
par(cex = 2,mai = c(0,0,0,0))
new_graph %>% plot(layout= layout,
                   vertex.label='',mark.groups = communities(x = components),
                   mark.border = NA,
                   mark.col = viridis(n_communities, alpha = 0.4),
                   vertex.label.cex=1,vertex.label.color='black',
                   main='')
text(x = c(-.7,.2,.9,-.7),y= c(.95,1.2,.3,-0.65),labels = 1:4 %>% as.character(),cex = 3)
# legend("topright", legend = c(seq_len(n_communities)),
#        col = viridis(n_communities)
#        ,pch = 16L,
#        title = "Module membership")
# legend("bottomright", legend = V(new_graph)$class_taxonomy %>% unique(),
#        pch = shapes() %>% {Filter(function(x) x %in% c('circle','pie','square','sphere'),.)},
#        title = "Phylum taxonomy")
dev.off()
table <- stat_community %>%
  filter(community %in% c(3L,4L)) %>%
  select(community, meanAbundance,taxonomy) %>%
  rename(Module = community, `Mean abundance` = meanAbundance, Taxonomy = taxonomy) %>%
  arrange(Module, desc(`Mean abundance`)) %>%
  kableExtra::kable(caption = 'The OTUs found in interaction module 3 and 4 of the network together with their taxonomies and overall mean abundance',
                    label = 'modules_main',format = 'latex',booktabs = TRUE) %>%
  kable_styling(font_size = 6)
table
table %>% cat(file = "figures/modules_main.tex")
dynamic_graph <- new_graph %>% delete_edge_attr("lty") # For not showing dotted lines in the dynamic plots
all_days <- relevant_physeq %>% sample_data() %>% pluck("Day") %>% unique()
RK_H_plot_code <- ~ {
par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5,xaxs = "i")
dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'RK',
                  nutrient = 'H',type = 'rank',scale_factor = 3,
                  mark.groups = communities(x = components),
                  mark.col = viridis(n_communities, alpha = 0.4),
                  layout = layout)
}
KR_H_plot_code <- ~ {
  par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5, xaxs = "i")
  dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'KR',
                  nutrient = 'H',type = 'rank',scale_factor = 3,
                  mark.groups = communities(x = components),
                  mark.col = viridis(n_communities, alpha = 0.4),
                  layout = layout, xlim=c(0,0))
}
RK_L_plot_code <- ~ {
  par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5)
  dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'RK',
                  nutrient = 'L',type = 'rank',scale_factor = 3,
                  mark.groups = communities(x = components),
                  mark.col = viridis(n_communities, alpha = 0.4),
                  layout = layout)
}
KR_L_plot_code <- ~ {
  par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5)
  dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'KR',
                  nutrient = 'L',type = 'rank',scale_factor = 3,
                  mark.groups = communities(x = components),
                  mark.col = viridis(n_communities, alpha = 0.4),
                  layout = layout)
}
combined_plot_H <- plot_grid(RK_H_plot_code,KR_H_plot_code,align = "v",ncol = 1,labels = "auto",label_size = 50)
ggsave(filename = "ss_main_H_combined_dynamic_network.pdf",plot = combined_plot_H,path = 'figures/',
width = 36, height = 50, units = 'cm')
combined_plot_L <- plot_grid(RK_L_plot_code,KR_L_plot_code,align = "v",ncol = 1,labels = "auto",label_size = 50)
ggsave(filename = "ss_main_L_combined_dynamic_network.pdf",plot = combined_plot_H,path = 'figures/',
       width = 36, height = 50, units = 'cm')
interaction_trend(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,2,4,8,16,24,28,29,30,32,40,50),selection_group = 'KR',
                  nutrient = 'H',type = 'rank',subset = components[[3]])
relevant_physeq %>% subset_samples_tidy()
community_info_frame <- map(stat_community$community %>% unique(),
    ~ prune_taxa(stat_community %>% filter(community ==.x) %>% pluck("ID"),reference_physeq) %>%
      sample_sums() %>% 
      {tibble(sample = names(.),community = as.factor(.x),abundance = .)}
    ) %>% bind_rows() %>% left_join(sample_data(reference_physeq) %>% unclass() %>% `class<-`('data.frame') %>% rownames_to_column("sample"))
(ggplot(data=community_info_frame %>% filter(community %in% c(1L,4L)),
                                        mapping = aes(x=Day %>% as.numeric(),y=abundance,
                                                      color=selection.regime.before.switch.day.28,
                                                      linetype=community,
                                                      group=interaction(Reactor.name,community)))+geom_point(size = 5) + 
  geom_line(alpha=.8)+
  guides(color=guide_legend(title='Initial selection'),
         linetype=guide_legend(title='Network module'))+
  xlab('Day')+ylab('Relative abundance')+ theme_bw() +
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 30),axis.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),legend.position = c(.9,0.8))) %>% 
  {ggsave(plot=.,filename = "ss_module_count.pdf",path = 'figures/',
         width = 40, height = 50, units = 'cm')}


tablerelevant_OTUs <- as.character(stat_community$ID)
# plot_tree(physeq = reduced_physeq,shape = "membership")+guides(shape=guide_legend('Community membership'))
tree_dir <- 'figures'
# dir.create(tree_dir,recursive = TRUE,showWarnings = FALSE)
tree <- phy_tree(physeq = reduced_physeq)
pdf(file = '{tree_dir}/ss_community_phylogram.pdf' %>% glue(),width = 20/cm(1),height = 25/cm(1))
par(mar=c(0,0,0,0),xpd=TRUE,cex = 1)
plot.phylo(x = tree,show.node.label = FALSE,use.edge.length = TRUE,
           # direction = "upwards",
           show.tip.label = FALSE)
tiplabels(text = '',frame = 'none',pch = pch_palette[match(class_info,unique_classes)],cex = 1.25
          ,col = color_palette[stat_community$community[match(tree$tip.label,stat_community$ID)]])
legend(x =0.15 ,y=75, legend = c(seq_len(n_communities)), 
       col = color_palette
       ,pch = 16L,
       title = "Module membership",cex = 2)
legend(x=0.16,y=25, legend = unique_classes, 
       col = 'black'
       ,pch = pch_palette,
       title = "Class level taxonomy",cex = 1.5)
dev.off()

