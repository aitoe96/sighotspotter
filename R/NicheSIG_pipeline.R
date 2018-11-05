#TODO: documentation
#!/usr/bin/env Rscript
##Script to extract the percent expressed cell and sum of product of expression betwen source and target nodes
# oarsub -l nodes=1/core=1,walltime=65:00:00 'source /etc/profile; module load lang/R/3.4.0-intel-2017a-X11-20170314-bare; Rscript MCSD_06022018.R OD_old.txt'
#USAGE: Rscript MCSD_06022018.R <input_data> <signaling_interactome> <cutoff for minimum no. of cell in which two gene must be expressed> <Number of random networks>

##  Read INPUT files
#args = commandArgs(trailing = TRUE)
#input_data = args[1]
#cutoff = args[2]
#DE_Genes_data = args[3]
#percentile = args[4]


NicheSIG_pipeline <- function(species, input_data, cutoff, DE_Genes_data, percentile, invert_DE = FALSE, showprogress = TRUE)
{
  ## Load necessary libraries
  library(NicheSIG)
  library(plyr)
  library(igraph)
  library(Matrix)
  library(reshape2)
  library(parallel)
  #load("HUMAN_Background_Network.RData")
  #TODO: load correct dataset
  if(species == "MOUSE"){
    load(system.file("extdata", "MOUSE_Background_Network.RData", package = "NicheSIG"))
  } else {
    if(species == "HUMAN"){
      load(system.file("extdata", "HUMAN_Background_Network.RData", package = "NicheSIG"))
    } else {
      stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
  }

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

  if(showprogress){
    incProgress()
    setProgress(detail = "Preprocessing")
  }

  subg=Data_preprocessing(idata,cutoff,species)

  Steady_state_true=Markov_chain_stationary_distribution(subg)
  #out2_3=paste("SS_R_",args[1],"_",cutoff,".txt",sep="")
  #write.table(Steady_state_true,out2_3,sep="\t",row.names=FALSE)

  if(showprogress){
    incProgress()
    setProgress(detail = "Steady state calculation")
  }

  int=high_probability_intermediates(Steady_state_true, intermediates, percentile)

  if(showprogress){
    incProgress()
    setProgress(detail = "Finding signaling candidates")
  }

  gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )

  nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)

  if(showprogress){
    incProgress()
    setProgress(detail = "Markov chain")
  }
  #comp score in serial
  comp_score=lapply(int,spcal_path_weight,nTF,gintg)
  #Compatability_score=unlist(comp_score)
  #comp_score1=cbind(as.data.frame(int),as.data.frame(unlist(comp_score)))
  #colnames(comp_score1)=c("Gene", "Compatability_Score")
  #comp_score1=join(comp_score1,Steady_state_true,by=c("Gene"),type="inner",match="first")
  #comp_score1=comp_score1[order(comp_score1$Compatability_Score,decreasing = TRUE),]

  if(showprogress){
    incProgress()
    setProgress(detail = "Compatibility scores")
  }
  comp_score1=compatability_score(comp_score,Steady_state_true,int)
#  out_CS=paste("compatability_score_",args[1],sep="")
#  out_SS=paste("Steady_state_",args[1],sep="")
#  names(comp_score1)=c("Gene", out_CS, out_SS)
  #out = paste("compatability_score_",input_data,"_",cutoff,"_",percentile,".txt",sep="")
  #write.table(comp_score1,out,sep="\t",quote=FALSE,row.names=FALSE)
  return(comp_score1)
}
