#function to get "all"one" shortest path edges from receptor_lig to intTFs based on weights
shortest_path_edges <- function(s,t,g){
  if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"))$res)
  if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
  edges=lapply(paths,edge_shortest,g)
  return(edges)
}

#function to get all shortest path edges from receptor_lig to intTFs
shortest_path_edges_all <- function(s,t,g){
  if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"),weight=NA)$res)
  if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
  edges=lapply(paths,edge_shortest,g)
  return(edges)
}

#function for returning all shortest paths
shortest.paths.all <- function(s,t,g){
    if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
    if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
    paths=(get.all.shortest.paths(g, s, t, mode = c("out"),weight=NA)$res)
    if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
    return(paths)
}

#funtion to get paths for each tf
paths.tf <- function(t,s,g) #x is a list of comp score
{
  paths=lapply(s,shortest.paths.all,t,g)
  #edges=lapply(paths,edge_shortest,g)
}

#function for retaining edges
edge_shortest <- function(path, graph)
{
  edges=E(graph, path=path)
}

# function to trim results
.trimResults <- function(results, active = TRUE) {

  res_trimmed = results[,1:2]
  if (active){
    res_trimmed <- head(res_trimmed, 10L)
    res_trimmed <- res_trimmed[res_trimmed[,2]>0.5,]
    colnames(res_trimmed) <- c('Active signaling hotspots', 'Compatibility score')
  } else
  {
    res_trimmed <- tail(res_trimmed, 10L)
    res_trimmed <- res_trimmed[res_trimmed[,2]<0.5,]
    colnames(res_trimmed) <- c('Inactive signaling hotspots', 'Compatibility score')
    res_trimmed <- res_trimmed[order(res_trimmed$'Compatibility score'),]
  }
  res_trimmed[,2] = round(res_trimmed[,2],4)
  rownames(res_trimmed) <- NULL
  return(res_trimmed)
}

#making the edge weight non-negative and taking the absolute
prun.int.g <- function(g){
  E(g)$weight=1-(abs(E(g)$weight))
  #function to delet edges from tf to dummy un the network
  del=incident(g, "Dummy", mode = c("in"))
  g <- delete.edges(g,del)
  g <- set.vertex.attribute(g,"name","Dummy","NICHE")
  return(g)
}

#making the edge weight abs
abs_weight <- function(g){
  E(g)$weight=(abs(E(g)$weight))
  return(g)
}



#function to generate the shortest path network given source, target and a network
shortest_path_network <- function (s,t,g, mst){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from TFs to dummy as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  if (mst==T){
    edges_g=lapply(s,shortest_path_edges,t,g)
    net <- subgraph.edges(g, unlist(edges_g), delete.vertices = T)
    net=mst(net)

  } else
  {edges_g=lapply(s,shortest_path_edges,t,g)
  net <- subgraph.edges(g, unlist(edges_g), delete.vertices = T)
  }
  return(net)
}

#' function for converting the input graph and the ints to a shortestpat network
#' @param a,i active and inactive source nodes or int
#' @param g input graph on which shortest path network must be inferred i.e. gintg
#' @param t terminal nodes
to_sp_net <- function(a,i,g,t){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  edges_a=lapply(a,shortest_path_edges,t,g)
  edges_i=lapply(i,shortest_path_edges,t,g)
  edges_d=lapply("Dummy",shortest_path_edges_dummy,c(a,i),g)
  #subnetwork from SP edges
  sp_sub_net=subgraph.edges(g, c(unlist(edges_a),unlist(edges_i),unlist(edges_d)), delete.vertices = T)
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(a), value="a_int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(i), value="i_int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(t), value="tf")
  return(sp_sub_net)
}

#' function for converting the input graph and the ints to a shortestpat network
#' @param a,i active and inactive source nodes or int
#' @param g input graph on which shortest path network must be inferred i.e. gintg
#' @param t terminal nodes
to_sp_net_int <- function(s,g,t,deg){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  edges_a=lapply(s,shortest_path_edges_all,t,g)
  edges_d=lapply("NICHE",shortest_path_edges_all,c(s),g)
  #classifying up and downregulated TFs
  up_t=nonterminal_DE_TFs(g,deg[deg$state==1,],non_interface_TFs)
  down_t=nonterminal_DE_TFs(g,deg[deg$state==-1,],non_interface_TFs)
  #subnetwork from SP edges
  sp_sub_net=subgraph.edges(g, c(unlist(edges_a),unlist(edges_d)), delete.vertices = T)
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(s), value="int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(up_t), value="upregulated")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(down_t), value="downregulated")
  return(sp_sub_net)
}

#' function for changing the visnetwork edge color
#' @param visg the visnetwork object
vis.edge.color <- function(visg){
  visg$edges$color[visg$edges$Effect==1]="green"
  visg$edges$color[visg$edges$Effect==-1]="red"
  return(visg)
}

#' Function to plot the Visnetowkr object
#' @param visg the visnetwork object
vis.net.plot <- function(visg){
  #hierarchy
  visNetwork(visg$nodes,visg$edges) %>% visNodes(visg, shape="box",color="grey") %>% visIgraphLayout(layout = "layout_as_tree",root="NICHE",flip.y = F) %>% visEdges(arrows = "to") %>%  visOptions(highlightNearest = list(enabled =TRUE, degree = 1, hover = T), nodesIdSelection = TRUE)  %>% visEdges(smooth = T) %>% visGroups(visg, groupname="int", shape="star", color="blue") %>% visGroups(visg, groupname="upregulated", color = "red",shape="triangle") %>% visGroups(visg, groupname="downregulated", color = "green", shape="triangle") %>% visPhysics(stabilization = FALSE) %>% visEdges(smooth = FALSE)
}

#' function for plotting a union network
#' @param visg is the visnetwork graph object
vis_plot_union <- function(visg){
  visge=lapply(visg,function(x) x$edges)
  visgn=lapply(visg,function(x) x$nodes)
  visgnj=join_all(visgn,by=c("id"),type="full")
  visgej=join_all(visge,by=c("from","to"),type="full")
  vis_union_A <- list(nodes=visgnj,edges=visgej)
}
