library(micInt)
library(phyloseq)
library(tidyverse)
library(igraph)
library(magrittr)
library(glue)
library(viridis)
library(xtable)
library(ape)
library(rlang)
library(scales)
library(kableExtra)
library(cowplot)
source("dynamic_networks_utils.R")
# The purpose of this script is to identify and explore interaction Modules (network communities) in the
# ReBoot interaction networks for the selection-switch dataset, testing for robustness to filtering, q-cutoff value,
# random noise and choice of similarity measure
options(stringsAsFactors = FALSE)
MAX_EDGES = 500
top_res <- readRDS('top_res.rds')
if(!dir.exists("figures")){
  dir.create("figures")
}
total_compute_frame <- readRDS('total_combined_combinations.rds')
total_compute_frame$top_res <- top_res
annotated_count_frame <- readRDS('annotated_count_frame.rds')
main_filtering_criteria <-
  exprs(
    dataset = dataset == 'ss',
    noise = noise_level == 'low',
    filtering = threshold == 5e-04,
    abundance = abundance_type == 'absolute',
    sim_measure = sim_measure_name == 'spearman_normal',
    .named = TRUE
  )
custom_filtering_criteria <-
  exprs(
    noise = noise_level == 'middle',
    filtering = threshold == 1e-03,
    abundance = abundance_type == 'relative',
    sim_measure = sim_measure_name == 'pearson_normal',
    .named = TRUE
  )

int_count_cumulative_plot_function <- function(frame){
  ggplot(frame) + geom_line(size = 5, mapping = aes(x=threshold,
                                          y=count_interactions,
                                          group =interaction_sign,
                                          linetype = factor(interaction_sign,levels = c("positive","negative")))
  ) +
    scale_x_log10(limits = c(10^(-30),1), labels = trans_format("log10",math_format()), breaks = 10^ seq(-30,0,5)) +
    scale_y_continuous(labels = comma) +
    guides(linetype = guide_legend("Sign of interaction")) +
    xlab("q-value threshold") + ylab("Number of significant interactions") +
    geom_vline(xintercept = 0.05,size = 2) + 
    annotate(geom = "text",
             label = "0.05",
             x = 0.1,
             y = frame$count_interactions %>% max() %>% divide_by(2),
             angle = 90, 
             vjust = 1,
             size = 20)
}

for (parameter in names(custom_filtering_criteria)) {
  message(parameter)
  if(parameter == "sim_measure"){
    invisible()
  }

  this_filtering_criteria <-
    modifyList(main_filtering_criteria, custom_filtering_criteria[parameter])
  selected_row <-
    total_compute_frame %>% filter(!!!this_filtering_criteria %>% unname())
  reference_physeq <-
    total_compute_frame %>% filter(!!!(this_filtering_criteria %>%
                                         modifyList(list(
                                           abundance = expr(abundance_type == 'relative')
                                         )) %>%
                                         unname())) %>%
    pluck("physeq") %>%
    pluck(1)
  count_frame_filtering_criteria <- this_filtering_criteria
  # Yes, this is a very ugly hack to change threshold == Num to filter_threshold == Num
  count_frame_filtering_criteria$filtering[[2]] <- sym("filter_threshold")
  count_frame_filtering_criteria$sim_measure[[3]] %<>% as.character() %>% 
    str_extract("^[:alpha:]+")
  interaction_count_frame <- annotated_count_frame %>% 
    filter(!!! count_frame_filtering_criteria %>% unname(),
           threshold > 10^-30, threshold < 0.05)
  int_count_plot <- int_count_cumulative_plot_function(interaction_count_frame)
  ggsave(filename = "figures/ss_{parameter}_int_count.pdf" %>% glue(),
         plot = int_count_plot +theme_bw() + theme(legend.position = "none",axis.title = element_text(size = 75),
                                                   axis.text = element_text(size = 55),
                                                   strip.text = element_text(size = 55)),
         device = 'pdf',width = 60,
         height = 50, units = 'cm')
  
  full_table <- selected_row %>% pluck('top_res') %>% pluck(1)
  relevant_physeq <- selected_row %>% pluck("physeq") %>% pluck(1)
  full_table %>% nrow()
  full_table %>% filter(q.value < 1e-4) %>% nrow()
  full_table %>% filter(q.value < 1e-10) %>% nrow()
  full_table$q.value[MAX_EDGES + 1]
  table <- full_table %>% extract(seq_len(min(MAX_EDGES, nrow(.))), )
  # We remove the negative interactions
  table_positive = table %>%  subset(sim.score > 0)
  graph = table_positive %>% micInt::as.edgelist() %>%
    igraph::graph_from_edgelist(directed = FALSE) %>% set_edge_attr('color', value =
                                                                      'blue')  %>% 
    set_edge_attr('lty',value = 1L)
  # We are only interrested in connected components with a certain number of members
  graph_components <- components(graph)
  components_to_keep <- which(graph_components$csize > 10)
  nodes_to_keep <-
    graph_components$membership %in% components_to_keep %>% which()
  graph_reduced <-
    induced_subgraph(graph = graph, vids = nodes_to_keep)
  
  components = graph_reduced %>% igraph::walktrap.community(steps = 20)
  n_communities <- components$membership %>% unique() %>% length()
  V(graph_reduced)$color <-
    viridis(n_communities)[components$membership]
  V(graph_reduced)$membership <- components$membership
  #  Adds edges of negative interactions
  negative_graph = table %>% subset(sim.score < 0) %>% micInt::as.edgelist() %>%
    igraph::graph_from_edgelist(directed = FALSE) %>% set_edge_attr('color', value =
                                                                      'red') %>% 
    set_edge_attr('lty',value = 2L)
  new_graph = union2(graph_reduced, negative_graph)
  all_displayed_taxa <- V(new_graph) %>% names()
  relevant_physeq <-
    selected_row %>% pluck('physeq') %>% pluck(1) %>% {
      phyloseq::prune_taxa(all_displayed_taxa, x = .)
    }
  reduced_physeq <- relevant_physeq
  t_table <- tax_table(reduced_physeq)
  class_info <- t_table[, "Class"]
  unique_classes <- class_info %>% unique()
  pch_palette <- c(15:18, 8:14)
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
  phylum_taxonomy <-
    tax_table(relevant_physeq)[all_displayed_taxa, "Phylum"] %>% as.factor()
  par(cex.main = 1.2)
  color_palette <- viridis(n_communities)
  message(glue("{V(new_graph) %>% length()} nodes in the graph"))
  pdf(
    'figures/ss_{parameter}_community_network.pdf' %>% glue(),
    width = 20 / cm(1),
    height = 25 / cm(1)
  )
  # pdf('../Presentation-poster/Jakob/ss_community_network.pdf',width = 20/cm(1),height = 25/cm(1))
  new_graph %>% plot(
    layout = layout,
    vertex.label = '',
    mark.groups = communities(x = components),
    mark.border = NA,
    vertex.label.cex = 1,
    vertex.label.color = 'black',
    ylim=c(-0.5,1),
    main = ''
  )
  legend(
    "topright",
    legend = c(seq_len(n_communities)),
    col = viridis(n_communities)
    ,
    pch = 16L,
    cex = 2,
    title = "Module membership"
  )
  # legend("bottomright", legend = V(new_graph)$class_taxonomy %>% unique(),
  #        pch = shapes() %>% {Filter(function(x) x %in% c('circle','pie','square','sphere'),.)},
  #        title = "Phylum taxonomy")
  dev.off()
  sim_measure_type <- ifelse(parameter == 'sim_measure',"linear","rank")
  table <- stat_community %>%
    filter(community %in% c(3L, 4L)) %>%
    select(community, ID, meanAbundance, taxonomy) %>%
    rename(
      Module = community,
      Identifier = ID,
      `Mean abundance` = meanAbundance,
      Taxonomy = taxonomy
    ) %>%
    arrange(Module, desc(`Mean abundance`)) %>%
    kableExtra::kable(
      caption = 'The OTUs found in interaction module 3 and 4 of the network together with their taxonomies and overall mean abundance',
      label = 'modules_{parameter}' %>% glue(),
      format = 'latex',
      booktabs = TRUE
    ) %>%
    kable_styling(font_size = 6)
  table
  table %>% cat(file = "figures/modules_{parameter}.tex" %>% glue())
  dynamic_graph <- new_graph %>% delete_edge_attr("lty") # For not showing dotted lines in the dynamic plots
  all_days <-
    relevant_physeq %>% sample_data() %>% pluck("Day") %>% unique()
  # pdf(
  #   'figures/ss_{parameter}_RK_H_dynamic_network.pdf' %>% glue(),
  #   width = 20 / cm(1),
  #   height = 25 / cm(1)
  # )
  # par(mfrow = c(4, 3), mar = c(0, 0, 5, 0))
  RK_plot_code <- ~ {
    par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5)
    dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'RK',
                    nutrient = 'H',type = 'rank',scale_factor = 3,
                    mark.groups = communities(x = components),
                    mark.col = viridis(n_communities, alpha = 0.4),
                    layout = layout)
    # dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,2,4,8,16,24,28,29,30,32,40,50),selection_group = 'RK',
    #                 nutrient = 'H',type = 'rank',scale_factor = 3,
    #                 mark.groups = communities(x = components),layout = layout)
  }
  interaction_trend(
    graph = dynamic_graph,
    physeq = relevant_physeq,
    days = c(1, 2, 4, 8, 16, 24, 28, 29, 30, 32, 40, 50),
    selection_group = 'RK',
    nutrient = 'H',
    type = sim_measure_type,
    subset = components[[3]]
  )
  # dev.off()
  # pdf(
  #   'figures/ss_{parameter}_KR_H_dynamic_network.pdf' %>% glue(),
  #   width = 20 / cm(1),
  #   height = 25 / cm(1)
  # )
  # par(mfrow = c(4, 3), mar = c(0, 0, 5, 0))
  KR_plot_code <- ~ {
    par(mfrow = c(2,3),mar = c(0,0,1.5,0),cex = 1.5)
    dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,16,28,30,40,50),selection_group = 'KR',
                    nutrient = 'H',type = 'rank',scale_factor = 3,
                    mark.groups = communities(x = components),
                    mark.col = viridis(n_communities, alpha = 0.4),
                    layout = layout)
    # dynamic_z_score(graph = dynamic_graph,physeq = relevant_physeq,days = c(1,2,4,8,16,24,28,29,30,32,40,50),selection_group = 'KR',
    #                 nutrient = 'H',type = 'rank',scale_factor = 3,
    #                 mark.groups = communities(x = components),layout = layout)
  }
  # dev.off()
  interaction_trend(
    graph = dynamic_graph,
    physeq = relevant_physeq,
    days = c(1, 2, 4, 8, 16, 24, 28, 29, 30, 32, 40, 50),
    selection_group = 'KR',
    nutrient = 'H',
    type = sim_measure_type,
    subset = components[[3]]
  )
  combined_plot <- plot_grid(RK_plot_code,KR_plot_code,align = "v",ncol = 1,labels = "auto",label_size = 50)
  ggsave(filename = "ss_{parameter}_combined_dynamic_network.pdf" %>% glue(),
         plot = combined_plot,path = 'figures/',
         width = 35, height = 50, units = 'cm')
  
  relevant_physeq %>% subset_samples_tidy()
  community_info_frame <- map(
    stat_community$community %>% unique(),
    ~ prune_taxa(
      stat_community %>% filter(community == .x) %>% pluck("ID"),
      reference_physeq
    ) %>%
      sample_sums() %>%
      {
        tibble(
          sample = names(.),
          community = as.factor(.x),
          abundance = .
        )
      }
  ) %>% bind_rows() %>% left_join(
    sample_data(reference_physeq) %>% unclass() %>% `class<-`('data.frame') %>% rownames_to_column("sample")
  )
  (
    ggplot(
      data = community_info_frame %>% filter(community %in% c(1L, 4L)),
      mapping = aes(
        x = Day %>% as.numeric(),
        y = abundance,
        color = selection.regime.before.switch.day.28,
        linetype = community,
        group = interaction(Reactor.name, community)
      )
    ) + geom_point(size = 5) +
      geom_line(alpha = .8) +
      guides(
        color = guide_legend(title = 'Initial selection'),
        linetype = guide_legend(title = 'Network module')
      ) +
      xlab('Day') + ylab('Relative abundance')+ theme_bw() +
      theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = c(.9, 0.8)
      )
  ) %>%
    {
      ggsave(
        plot = .,
        filename = glue("ss_{parameter}_module_count.pdf"),
        path = 'figures/',
        width = 40,
        height = 50,
        units = 'cm'
      )
    }
  
  
  tablerelevant_OTUs <- as.character(stat_community$ID)
  # plot_tree(physeq = reduced_physeq,shape = "membership")+guides(shape=guide_legend('Community membership'))
  tree_dir <- 'figures/'
  # dir.create(tree_dir,recursive = TRUE,showWarnings = FALSE)
  tree <- phy_tree(physeq = reduced_physeq)
  pdf(
    file = '{tree_dir}/ss_{parameter}_community_phylogram.pdf' %>% glue(),
    width = 20 / cm(1),
    height = 25 / cm(1)
  )
  par(mar = c(0, 0, 0, 0), xpd = TRUE)
  plot.phylo(
    x = tree,
    show.node.label = FALSE,
    use.edge.length = TRUE,
    show.tip.label = FALSE
  )
  tiplabels(
    text = '',
    frame = 'none',
    pch = pch_palette[match(class_info, unique_classes)],
    cex = 1.25
    ,
    col = color_palette[stat_community$community[match(tree$tip.label, stat_community$ID)]]
  )
  legend(
    x = 0.2 ,
    y = 75,
    legend = c(seq_len(n_communities)),
    col = color_palette,
    cex = 1.5,
    pch = 16L,
    title = "Module membership"
  )
  legend(
    x = 0.17,
    y = 25,
    legend = unique_classes,
    col = 'black'
    ,
    pch = pch_palette,
    cex = 1.5,
    title = "Class level taxonomy"
  )
  dev.off()
  full_physeq <- selected_row$physeq %>%
    pluck(1)
  full_sample_data <- sample_data(full_physeq) %>%  
    data.frame() %>%
    mutate(replicate = str_extract(Reactor.name,"[:digit:]$"))
  reactor_frame <-
    subdivide_by_environment(full_physeq, c("selection.group", "Nutrient"))
  
  reactor_frame$reactor_sorted <- purrr::map(reactor_frame$phyloseq,
                                             ~ subdivide_by_environment(., variables = "Reactor.name"))
  
  reactor_frame$time_series <- reactor_frame$reactor_sorted %>%
    purrr::map( ~ (
      magrittr::extract2(., 'phyloseq') %>%
        purrr::map( ~ micInt::OTU_time_series(
          table = .,
          time_points = 'Day'
        ))
    ))
  all_time_plot <- reactor_frame$time_series %>%
    unlist(recursive = FALSE) %>% plot_trajectory(linetype = "time_series",
                                                  label = FALSE,color = "replicate",
                                                  label_size = 5)
  all_time_plot$data$selection.group <-
    full_sample_data$selection.group[match(all_time_plot$data$time_series,
                                           full_sample_data$Reactor.name)]
  all_time_plot$data$replicate <- full_sample_data$replicate[match(all_time_plot$data$time_series,
                                                                   full_sample_data$Reactor.name)]
  all_time_plot$data$Nutrient <-
    full_sample_data$Nutrient[match(all_time_plot$data$time_series,
                                    full_sample_data$Reactor.name)]
  all_time_plot$data$Switch <- if_else(all_time_plot$data$time_points <= 28L,true = "Before switch",false = "After switch") %>% 
    as.factor() %>%
    relevel("Before switch")
  all_time_plot$data$selection.regime.at.sampling <- all_time_plot$data %$% if_else(xor(selection.group == "RK",
                                                                                        Switch == "After switch"),
                                                                                    "r",
                                                                                    "K")
  
  trajectory_plot <-
    all_time_plot  +
    geom_path(size = 5,arrow = arrow(length = unit(0.2,'inches'))) +
    geom_point(size = 8,alpha = 0.5) +
    viridis::scale_color_viridis(discrete = TRUE) +
    geom_text(aes_string(label = "time_points"), size = 12, color = "black") +
    facet_grid(selection.regime.at.sampling ~ Switch) + theme_bw() + 
    theme(plot.title = element_text(size = 18),axis.title = element_text(size = 18),axis.text = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),legend.position = 'none')
  ggsave(
    filename = "figures/ss_{parameter}_trajectories.pdf" %>% glue(),
    plot = trajectory_plot + theme(legend.position = "none",axis.title = element_text(size = 55),
                                   axis.text = element_text(size = 55),
                                   strip.text = element_text(size = 55)),
    device = 'pdf',
    width = 40,
    height = 50,
    units = 'cm'
  )
}
