library(matrixStats)
library(igraph)
library(micInt)
source("tidy_phylo_subset.R")
extract_relevant_physeq <- function(physeq,selection_group,nutrient){
  subset_samples_tidy(physeq = physeq,
                      Nutrient == !! nutrient,selection.group == !! selection_group)
}
dynamic_z_score <- function(graph, physeq,days,selection_group,nutrient,
                            type = c('linear','rank')
                            ,scale_factor = 1,...){
  relevant_physeq <- extract_relevant_physeq(physeq,
                                             selection_group = selection_group,nutrient = nutrient)
  day_devided_physeq <- subdivide_by_environment(physeq = relevant_physeq %>% subset_samples_tidy(Day %in% days),
                                                 variables = "Day",
                                                 keep_variables = TRUE) 
  if(type == 'rank'){
    daily_averages <-  mutate(day_devided_physeq,abundance_value = map(phyloseq, . %>%
                                                                         otu_table %>%
                                                                         refine_data(renormalize = FALSE,abundance_cutoff = NULL) %>%
                                                                         as.matrix() %>%
                                                                         apply(MARGIN = 1,FUN = rank) %>%
                                                                         rowMeans()
    )
    )
    
  }
  else{
    daily_averages <-  mutate(day_devided_physeq,abundance_value = map(phyloseq, . %>%
                                                                         otu_table %>%
                                                                         refine_data(renormalize = FALSE,abundance_cutoff = NULL) %>%
                                                                         as.matrix() %>%
                                                                         colMeans()
    )
    )
  }
  overall_matrix <- exec(rbind,!!! daily_averages$abundance_value)
  overall_mean <- colMeans(overall_matrix)
  overall_sd <- colSds(overall_matrix)
  names(overall_sd) <- colnames(overall_matrix)
  NODE_MIN_SIZE <- 3
  for (i in seq_len(nrow(daily_averages))){
    day_mean <- daily_averages$abundance_value[[i]]
    Day <- daily_averages$Day[[i]]
    z_score <- (day_mean - overall_mean) / overall_sd
    # Sometimes the standard deviation can be zero due to OTUs present in a 
    # few samples only, so we need to account for this
    z_score[!is.finite(z_score)] <- 0
    z_score_graph <- z_score[V(graph) %>% names()]
    V(graph)$z_score <- z_score_graph
    #TODO linear function between max and min size
    V(graph)$size <- NODE_MIN_SIZE + scale_factor*abs(z_score_graph)
    node_color_cuts <- c(-Inf,-1,1,Inf)
    node_colors <- c("black","gray","orange")
    V(graph)$color <- node_colors[findInterval(V(graph)$z_score,node_color_cuts)]
    # V(graph)$color <- ifelse(z_score_graph > 0, "black","orange")
    E(graph)$prod_score <- graph %>%
      get.edgelist() %>% {z_score_graph[.[,1]]*z_score_graph[.[,2]]}
    edge_color_cuts <- c(-Inf,-0.3,0.3,Inf)
    edge_colors <- c("red","gray","blue")
    E(graph)$color <- edge_colors[findInterval(E(graph)$prod_score,edge_color_cuts)]
    plot(graph,vertex.label.cex=0.001,vertex.label.color='yellow',
         cex.main=0.15,main = glue("Day: {Day}"),...)
  }
}

interaction_trend <- function(graph, physeq,days,selection_group,nutrient,
                                type = c('linear','rank'),subset = NULL){
  if(!is.null(subset)){
    graph <- subgraph(graph = graph,subset)
  }
  relevant_physeq <- extract_relevant_physeq(physeq,
                                             selection_group = selection_group,nutrient = nutrient)
  day_devided_physeq <- subdivide_by_environment(physeq = relevant_physeq %>% subset_samples_tidy(Day %in% days),
                                                 variables = "Day",
                                                 keep_variables = TRUE) 
  if(type == 'rank'){
    daily_averages <-  mutate(day_devided_physeq,abundance_value = map(phyloseq, . %>%
                                                                         otu_table %>%
                                                                         refine_data(renormalize = FALSE,abundance_cutoff = NULL) %>%
                                                                         as.matrix() %>%
                                                                         apply(MARGIN = 1,FUN = rank) %>%
                                                                         rowMeans()
    )
    )
    
  }
  else{
    daily_averages <-  mutate(day_devided_physeq,abundance_value = map(phyloseq, . %>%
                                                                         otu_table %>%
                                                                         refine_data(renormalize = FALSE,abundance_cutoff = NULL) %>%
                                                                         as.matrix() %>%
                                                                         colMeans()
    )
    )
  }
  overall_matrix <- exec(rbind,!!! daily_averages$abundance_value)
  overall_mean <- colMeans(overall_matrix)
  overall_sd <- colSds(overall_matrix)
  names(overall_sd) <- colnames(overall_matrix)
  meanInteraction <- map_dbl(daily_averages$abundance_value,~ {
                           day_mean <- .x
                           z_score <- (day_mean - overall_mean) / overall_sd
                           z_score_graph <- z_score[V(graph) %>% names()]
                           prod_score <-  graph %>%
                           get.edgelist() %>% {z_score_graph[.[,1]]*z_score_graph[.[,2]]}
                           mean(prod_score)
  }
                           )
  tibble(Day = daily_averages$Day, meanInteraction)
}
dynamic_abundance <- function(graph, physeq,days,selection_group,nutrient
                              ,scale_factor = 1){
  relevant_physeq <- extract_relevant_physeq(physeq,
                                             selection_group = selection_group,nutrient = nutrient)
  specific_graph <- map(days, ~ {
    day_physeq <- relevant_physeq %>% subset_samples_tidy(Day == !! .x)
    stats <- OTU_stats(day_physeq)
    special_graph <- graph
    V(special_graph)$size <- scale_factor*stats[V(special_graph) %>% names(),"meanAbundance",drop=TRUE]
    special_graph
  })
  walk2(specific_graph,days,~ plot(.x,vertex.label.cex=0.001,vertex.label.color='black',
                                   cex.main=0.15,main = glue("Day: {.y}")))
  relevant_physeq %>% subdivide_by_environment(variables = "Reactor.name") %$% 
    map(phyloseq,~ OTU_time_series(.,time_points = "Day")) %>% plot_trajectory(label = TRUE)
}

union2<-function(g1, g2){
  
  #Internal function that cleans the names of a given attribute
  CleanNames <- function(g, target){
    #get target names
    gNames <- parse(text = (paste0(target,"_attr_names(g)"))) %>% eval 
    #find names that have a "_1" or "_2" at the end
    AttrNeedsCleaning <- grepl("(_\\d)$", gNames )
    #remove the _x ending
    StemName <- gsub("(_\\d)$", "", gNames)
    
    NewnNames <- unique(StemName[AttrNeedsCleaning])
    #replace attribute name for all attributes
    for( i in NewnNames){
      
      attr1 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_1"),"')"))) %>%
        eval
      attr2 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_2"),"')"))) %>%
        eval
      
      g <- parse(text = (paste0("set_",target,
                                "_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))"))) %>%
        eval
      
      g <- parse(text = (paste0("delete_",target,
                                "_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
      g <- parse(text = (paste0("delete_",target,
                                "_attr(g,'", paste0(i, "_2"),"')"))) %>% eval
      
    }
    
    return(g)
  }
  
  
  g <- igraph::union(g1, g2) 
  #loop through each attribute type in the graph and clean
  for(i in c("graph", "edge", "vertex")){
    g <- CleanNames(g, i)
  }
  
  return(g)
  
}