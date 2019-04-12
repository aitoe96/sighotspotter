#Script for generating subsamples of Background network for testing SigHotSpotter robustness
#load necessary libraries
library(SigHotSpotter)
library(plyr)
library(igraph)
library(Matrix)
library(reshape2)
library(parallel)
library(RSpectra)
library(dplyr)
library(snow)

#loading mouse background data
load(system.file("extdata", "MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))


##  Read INPUT files
input_data = system.file("extdata", "2i.txt", package = "SigHotSpotter")
#input_data2 = system.file("extdata", "lif.txt", package = "SigHotSpotter")
DE_Genes_data = system.file("extdata", "DE_BOOLEAN_2i_lif.txt", package = "SigHotSpotter")


idata = read.table(input_data,header=TRUE,stringsAsFactors=FALSE,quote="")
#idata = fread(input_data,header=TRUE,stringsAsFactors=FALSE,quote="")
DE_Genes=read.table(DE_Genes_data,sep="\t",header=T)
#invert Boolean gene expression
#DE_Genes[,2] = -DE_Genes[,2]
cutoff=50
percentile=90

subg=Data_preprocessing(idata,cutoff,"MOUSE")

Steady_state_true=Markov_chain_stationary_distribution(subg)

int=high_probability_intermediates(Steady_state_true,intermediates,percentile)

gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )

nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)

score=lapply(nTF,comp_score_tf, int, gintg)
#converting the nested list into a matrix whose row sum will give probability of each intermediate
score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
score_m_means=as.list(rowMeans(score_m))
#compatibility score for the real background network
final_score_true=compatability_score(score_m_means,Steady_state_true,int)

##iterative subsampling for estimating model robustness
i=100
DE_Genes_subsampled <- rep(list(DE_Genes), i)
non_interface_TFs_subsampled=rep(list(non_interface_TFs), i)
graph_subsampled=Data_preprocessing_subsample(idata,cutoff,i)
TF_TF_interactions_subsampled=lapply(rep(list(TF_TF_interactions), i), function(x) {x[-sample(1:nrow(x), ((nrow(x))*10)/100), ]})
Steady_state_true_subsampled=lapply(graph_subsampled,Markov_chain_stationary_distribution)
int_subsampled=lapply(Steady_state_true_subsampled, high_probability_intermediates,intermediates,percentile)
gintg_subsampled=mapply(integrate_sig_TF,graph_subsampled, Steady_state_true_subsampled,DE_Genes_subsampled,non_interface_TFs_subsampled, TF_TF_interactions_subsampled,SIMPLIFY = F)

nTF_subsampled=lapply(gintg_subsampled,nonterminal_DE_TFs,DE_Genes,non_interface_TFs)

final_score=list()
for (j in 1:i) {
  #score=lapply(nTF_subsampled[[j]],comp_score_tf, int_subsampled[[j]], gintg_subsampled[[j]])
  net_paths[[j]]=lapply(int_subsampled[[j]],all_path, nTF_subsampled[[j]], gintg_subsampled[[j]])
  net_weight[[j]]=lapply(int_subsampled[[j]],all_weight, nTF_subsampled[[j]], gintg_subsampled[[j]])
  score_fast=lapply(c(nTF_subsampled[[j]]), score_pattern,net_paths[[j]],net_weight[[j]])
  score_m=(matrix(unlist(score_fast), ncol=length(score_fast), byrow=F))
  score_m_means=as.list(rowMeans(score_m))
  #names()
  final_score[[j]]=compatability_score(score_m_means,Steady_state_true_subsampled[[j]],int_subsampled[[j]])
}
subsamples_scores=join_all(final_score,by="Gene",match="all")
subsampled_avg_activation_probability=(rowMeans(subsamples_scores[,colnames(subsamples_scores) == "Activation_probability"],na.rm=T))
subsampled_avg_values=cbind(subsamples_scores[1],subsampled_avg_activation_probability,rowMeans(subsamples_scores[,colnames(subsamples_scores) == "Steady_state"],na.rm=T))
names(subsampled_avg_values)=c("Gene","Avg_activation_probability","Avg_steady_state")
#final file with real and average of subsampled output
final_list=join(final_score_true,subsampled_avg_values,by=c("Gene"))
write.table(final_list,"subsampled.txt",row.names = F)
