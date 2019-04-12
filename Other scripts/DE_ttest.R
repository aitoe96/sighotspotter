#!/usr/bin/env Rscript
library (plyr)
args = commandArgs(trailing = TRUE)

pvalue_cal<- function(idx,a,b)
{
  l=t.test(a[idx,2:ncol(a)],b[idx,2:ncol(b)])
  mean=as.numeric(l$estimate)
  fc=log2(mean[1]/mean[2])
  list(as.data.frame(cbind(idx, mean[1], mean[2], fc, l$p.value)))
}

output_pval <- function(pval_list,idx)
{
  ll=join_all((pval_list),by=c("idx"),type="full")
  lol=cbind(t(t(names(idx))),ll[2:5])
  names(lol)=c("Gene","mean_population1","mean_population2","Log_Fold_change","p_value")
  adj_p=p.adjust(lol$`p_value`, method = "fdr")
  lol=cbind(lol,adj_p)
}

DEG_pvalue_singlecell <- function(a,b)
{
  idx <- seq_len(nrow(a))
  names(idx) <- (a[,1])
  ttest_p <- sapply(idx,pvalue_cal,a,b)
  pp=output_pval(ttest_p,idx)
  names(pp)=c("Gene","mean_population1","mean_population2","Log_Fold_change","p_value","adj_p")
  return(pp)
}


#read input data, first column is gene symbols
data1=read.table(args[1],sep="\t",quote="",header=T)
data2=read.table(args[2],sep="\t",quote="",header=T)
v=DEG_pvalue_singlecell(data1,data2)

out=paste("DEG_",args[1],"_",args[2],".txt",sep="")
out1=paste("DEG_SIG",args[1],"_",args[2],".txt",sep="")
out2=paste("DEG_BOOLEAN_",args[1],"_",args[2],sep="")

write.table(v,out,sep="\t",row.names=FALSE)

sub_v=subset(v, adj_p <0.05 & (Log_Fold_change > 1 | Log_Fold_change < -1))
write.table(sub_v,out1,sep="\t",row.names=FALSE)

bool=sub_v[,c("Gene","Log_Fold_change")]
bool$Log_Fold_change[bool$Log_Fold_change>0]<-1
bool$Log_Fold_change[bool$Log_Fold_change<0]<--1
names(bool)=c("Gene",out2)
write.table(bool,out2,sep="\t",row.names=FALSE,quote=F)

