#' script to plot the network for each intermediate
#' in total there must be 40 individual networks (10active and 10 inactive for each phenotype)
#' The script stats from the current sighotspotter output
library(visNetwork)
#for running the sighotspotter command by command
species="MOUSE"
input_data="2i.txt"
cutoff=50
DE_Genes_data="DE_TFs_2iVSlif.txt"
percentile= 90
invert_DE = F
showprogress = F

#extracting top ten compatible nodes
toiintI=.trimResults(final_score,F)
toiintA=.trimResults(final_score,T)
toiintI=c(as.matrix(toiintI$`Inactive signaling hotspots`))
toiintA=c(as.matrix(toiintA$`Active signaling hotspots`))
twoi=c(toiintI,toiintA)

#pruning the integrated networks
gintg.p=prun.int.g(gintg)

#building networks for all intermediates for active signaling hotspots
sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes)

#building networks for inactive signaling hotspots
sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes)

#converting to visNetwork
#list for active intermediates
vis_net_A <- lapply(sp_int_A,toVisNetworkData)
#list for inactive intermediates
vis_net_I <- lapply(sp_int_I,toVisNetworkData)

#for edge color
vis_net_A <- lapply(vis_net_A,vis.edge.color)
vis_net_I <- lapply(vis_net_I,vis.edge.color)

#plotting, the function plots the network of molecule that corresponds with the rank
#to plot the intermediate with rank 2
vis.net.plot(vis_net_I[[2]])

#plot all the graphs and retain them in a variable
all_plots_A <- lapply(vis_net_A,vis.net.plot)
all_plots_I <- lapply(vis_net_I,vis.net.plot)

#display the network plot
all_plots_I[[1]] #displays the plot of rank 1 node in active signaling molecule

#to save network in Sif file for cytoscape
write.table(vis_net_A[[1]]$edges,paste(toiintA[1],".","sif",sep=""),row.names = F,quote=F)

#plotting the union of all active ints sp networks
vis.net.plot(vis_plot_union(vis_net_A))
vis.net.plot(vis_plot_union(vis_net_I))

--------------#ROUGH IGNORE______________________________________________________

#if mst needed
vis_sp_gintg <- toVisNetworkData(mst(sp_gintg))

#changing the edge attributes of the integrated network to non-negative weights
gintg_n=non_neg_weight(gintg)

#removing the edges from dummy to TFs as it affectes the shortest paths
del=as_adj_edge_list(gintg_n, mode = c("in"))$Dummy
gintg_n_d <- delete.edges(gintg_n,del)

#SP edges
edges_g=lapply(twoi,shortest_path_edges,nTF,gintg_n_d)
edges_d=lapply("Dummy",shortest_path_edges,twoi,gintg_n_d)
#subnetwork from SP edges
nd=subgraph.edges(gintg_n_d, c(unlist(edges_g),unlist(edges_d)), delete.vertices = T)

#SP network
ll_dummy=lapply("Dummy",shortest_path_network,toiintI,gintg_n,mst=F)
ll=lapply(toiintI,shortest_path_network,nTF,gintg_n,mst=F)
ll_gintg=lapply(twoi,shortest_path_network,nTF,gintg,mst=F)

#union of SP network
uni_sp_net<-do.call(igraph::union, c(ll,ll_dummy))
uni_sp_net_gintg<-do.call(igraph::union, ll_gintg)

#converting the SP network to Visnetwork
ll_net_visnet <- toVisNetworkData(gg)
ll_net_visnet_gintg <- toVisNetworkData(uni_sp_net_gintg)

#converting the list of graphs to list of visgraphs
oo=lapply(ll,toVisNetworkData)
oo_gintg=lapply(ll_gintg,toVisNetworkData)

#network to be show in the sighotspotter
edges_g=lapply("Dummy",shortest_path_net,toiintI,gintg)
shortest_subnetwork_d_i <- subgraph.edges(gintg, unlist(edges_g), delete.vertices = T)
edges_g=lapply("Gsk3b",shortest_path_net,nTF,gintg)
shortest_subnetwork_i_t <- subgraph.edges(gintg, unlist(edges_g), delete.vertices = T)
SSnet_D_int_nTF=(igraph::union(shortest_subnetwork_d_i,shortest_subnetwork_i_t))
SSnet_D_int_nTF_mst=mst(SSnet_D_int_nTF)

#network from dummy to nTF

#igraph to VisGraph
SSnet_D_int_nTF_data <- toVisNetworkData(SSnet_D_int_nTF)
SSnet_D_int_nTF_data_mst <- toVisNetworkData(SSnet_D_int_nTF_mst)

ll_net_visnet <- toVisNetworkData(ll_net)

#for different edge color
vis_sp_gintg$edges$color[vis_sp_gintg$edges$Effect==1]="green"
vis_sp_gintg$edges$color[vis_sp_gintg$edges$Effect==-1]="red"

#plotting the visnetwork
#optimized
visNetwork(vis_sp_gintg$nodes,vis_sp_gintg$edges) %>% visNodes(vis_sp_gintg, shape="box") %>% visIgraphLayout(layout = "layout_with_graphopt") %>% visEdges(arrows = "to") %>%  visOptions(highlightNearest = list(enabled =TRUE, degree = 4, hover = T), nodesIdSelection = TRUE)  %>% visEdges(smooth = T) %>% visGroups(vis_sp_gintg, groupname="tf", color = "orange") %>% visGroups(vis_sp_gintg, groupname="a_int", color = "green") %>% visGroups(vis_sp_gintg, groupname="i_int", color = "red")
#as tree
visNetwork(vis_sp_gintg$nodes,vis_sp_gintg$edges) %>% visNodes(vis_sp_gintg, shape="box") %>% visHierarchicalLayout() %>% visEdges(arrows = "to") %>%  visOptions(highlightNearest = list(enabled =TRUE, degree = 4, hover = T), nodesIdSelection = TRUE)  %>% visEdges(smooth = T) %>% visGroups(vis_sp_gintg, groupname="tf", color = "magenta") %>% visGroups(vis_sp_gintg, groupname="a_int", color = "green") %>% visGroups(vis_sp_gintg, groupname="i_int", color = "red")
