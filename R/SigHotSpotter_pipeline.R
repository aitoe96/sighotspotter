#' General pipeline for SigHotSpotter
#'
#' The function computes compatibility scores for signaling intermediates
#'
#' @param species Currently supported species: "HUMAN", "MOUSE"
#' @param input_data File name for input gene expression data
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param DE_Genes_data Differential expression dataset (1 for up-regulated, -1 for down-regulated genes)
#' @param percentile Predicted intermediates are taken into account above this threshold
#' @param invert_DE If the differential expression should be inverted, default = FALSE
#' @param showprogress shows progress bar in shiny app if set to TRUE, set it to FALSE in batch mode without GUI, default = TRUE
#' @return Compatibility scores
#' @export
SigHotSpotter_pipeline <- function(species, input_data, cutoff, DE_Genes_data, percentile, invert_DE = FALSE, showprogress = TRUE)
{
  ## Import necessary libraries
  #library(SigHotSpotter)
  library(plyr)
  library(igraph)
  library(Matrix)
  library(reshape2)
  library(parallel)
  library(RSpectra)
  library(dplyr)
  library(snow)
  library(visNetwork)

  ## Choose correct dataset according to species
  if(species == "MOUSE"){
    load(system.file("extdata", "MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))
  } else {
    if(species == "HUMAN"){
      load(system.file("extdata", "HUMAN_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))
    } else {
      stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
  }
  ## Checking if showprogress is needed
  #FIXME: is showprogress even necessary?
  if(showprogress && !shiny::isRunning()){
    showprogress = FALSE
    warning('Can not show progress without the shiny environment.')
  }

  ## Load the data
  if(showprogress){
    incProgress()
    setProgress(detail = "Loading data")
  }

  idata = read.table(input_data,header=TRUE,stringsAsFactors=FALSE,quote="")
  DE_Genes=read.table(DE_Genes_data,sep="\t",header=T)
  if(invert_DE)
  {
    DE_Genes[,2] = -DE_Genes[,2]
  }

  ## Preprocessing the data
  if(showprogress){
    incProgress()
    setProgress(detail = "Preprocessing")
  }

  subg=Data_preprocessing(idata,cutoff,species)

  ## Calculate stationary distribution of the MC
  if(showprogress){
    incProgress()
    setProgress(detail = "Steady state calculation")
  }

  Steady_state_true=Markov_chain_stationary_distribution(subg)

  # Uncomment the following for debug
  #out2_3=paste("SS_R_",args[1],"_",cutoff,".txt",sep="")
  #write.table(Steady_state_true,out2_3,sep="\t",row.names=FALSE)

  ## Retrieves high probability intermediates
  int=high_probability_intermediates(Steady_state_true, intermediates, percentile)

  if(showprogress){
    incProgress()
    setProgress(detail = "Finding signaling candidates")
  }

  gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )

  nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)

  ## Computing compatibility scores
  if(showprogress){
    incProgress()
    setProgress(detail = "Compatibility scores")
  }

  #comp_score for each TF individually
  score=lapply(nTF,comp_score_tf, int, gintg)
  #converting the nested list into a matrix whose row sum will give probability of each intermediate
  score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
  score_m_means=as.list(rowMeans(score_m))
  final_score=compatability_score(score_m_means,Steady_state_true,int)

  ## Computing networks for visualization
  if(showprogress){
    incProgress()
    setProgress(detail = "Computing networks")
  }

  trimmed_score_I=.trimResults(final_score,F)
  trimmed_score_A=.trimResults(final_score,T)
  toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
  toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
  twoi=c(toiintI,toiintA)

  #pruning the integrated networks
  gintg.p=prun.int.g(gintg)

  #building networks for all intermediates for active signaling hotspots
  sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)

  #building networks for inactive signaling hotspots
  sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)

  #converting to visNetwork
  vis_net_A <- lapply(sp_int_A,toVisNetworkData)
  vis_net_I <- lapply(sp_int_I,toVisNetworkData)

  #for edge color
  vis_net_A <- lapply(vis_net_A,vis.edge.color)
  vis_net_I <- lapply(vis_net_I,vis.edge.color)

  ret_value = list(final_score=final_score,
                   trimmed_score_A=trimmed_score_A,
                   trimmed_score_I=trimmed_score_I,
                   vis_net_A=vis_net_A,
                   vis_net_I=vis_net_I,
                   complete = gintg.p)
  return(ret_value)
}
