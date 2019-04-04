#' General pipeline for NicheSIG
#'
#' The function computes compatibility scores for signaling intermediates
#'
#' @param species Currently supported species: "HUMAN", "MOUSE"
#' @param input_data File name for input gene expression data
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param DE_Genes_data Differential expression dataset (1 for up-regulated, -1 for down-regulated genes)
#' @param percentile Predicted intermediates are taken into account above this threshold
#' @param invert_DE If the differential expression should be inverted, default = FALSE
#' @param showprogress shows progress bar if TRUE, set it to FALSE in batch mode, default = TRUE
#' @return Compatibility scores
#' @export
NicheSIG_pipeline <- function(species, input_data, cutoff, DE_Genes_data, percentile, invert_DE = FALSE, showprogress = TRUE)
{
  ## Import necessary libraries
  library(NicheSIG)
  library(plyr)
  library(igraph)
  library(Matrix)
  library(reshape2)
  library(parallel)
  library(RSpectra)
  library(dplyr)
  library(snow)

  ## Choose correct dataset according to species
  if(species == "MOUSE"){
    load(system.file("extdata", "MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "NicheSIG"))
  } else {
    if(species == "HUMAN"){
      load(system.file("extdata", "HUMAN_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "NicheSIG"))
    } else {
      stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
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
  #extract all paths from all sources and targets
  net_paths=lapply(int,all_path, nTF, gintg)
  #calculate each path weight
  net_weight=lapply(int,all_weight, nTF, gintg)
  #calculate compatibility score for each TF
  score_fast=lapply(c(nTF), score_pattern,net_paths,net_weight)
  score=(matrix(unlist(score_fast), ncol=length(score_fast), byrow=F))
  score_means=as.list(rowMeans(score))
  final_score=compatability_score(score_means,Steady_state_true,int)
  return(final_score)
}
