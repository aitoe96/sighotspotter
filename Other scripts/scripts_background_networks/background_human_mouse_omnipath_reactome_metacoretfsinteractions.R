# Script to create the necessary background files for SigHotSpotter
# Inputs required are signaling interactome, list of TFs, list of ligands and receptors, list of kinases and phosphatases
# Outputs a markov chain based on the signaling interactome, interface-TFs, non-interface-TFs, signaling intermediates and other processed network files

library("biomaRt")
library("igraph")
library("plyr")

#read omnipath network
#omnipath network downloaded from http://archive.omnipathdb.org/
omni_net=read.table("omnipath_webservice_interactions__20180614-20181114.tsv",sep="\t",header=T)

#split the omnipath into mouse and human interaction
omnipath_rat=omni_net[omni_net$ncbi_tax_id_source==10116,]
omnipath_human=omni_net[omni_net$ncbi_tax_id_source==9606,]
omnipath_mouse=omni_net[omni_net$ncbi_tax_id_source==10090,]

#read reactome FI network
#data downloaded from http://reactomews.oicr.on.ca:8080/caBigR3WebApp2017/FIsInGene_071718_with_annotations.txt.zip
reactome_fi_real = read.table("FIsInGene_071718_with_annotations.txt",sep="\t",header=T)
#remove interactions without directionality
reactome_fi = reactome_fi_real[reactome_fi_real$Direction != "-",]
#swap columns of interactions with reverse edge directions
reactome_reverse1 = reactome_fi[reactome_fi$Direction=="<-",]
reactome_reverse2 = reactome_fi[reactome_fi$Direction=="|-",]
reactome_reverse3 = reactome_fi[reactome_fi$Direction=="|-|",] 
reactome_reverse4 = reactome_fi[reactome_fi$Direction=="<-|",] 
reactome_reverse5 = reactome_fi[reactome_fi$Direction=="<->",] 
reactome_reverse6 = reactome_fi[reactome_fi$Direction=="|->",] 
reactome_reverse = rbind(reactome_reverse1,reactome_reverse2,reactome_reverse3,reactome_reverse4,reactome_reverse5,reactome_reverse6)
reactome_re = reactome_reverse[,c(2,1,3,4,5)]
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("<->", 1, x, fixed = T)))
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("<-|", 1, x,fixed = T)))
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("|->", -1, x, fixed = T)))
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("|-|", -1, x,fixed = T)))
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("<-", 1, x, fixed = T)))
reactome_re = data.frame(lapply(reactome_re, function(x) gsub("|-", -1, x,fixed = T)))
reactome_re_f = reactome_re[,c(1,2,4)]
names(reactome_re_f) = c("Source","Target","Effect")
#removing <- and |- from reactomefi as they are in reactome_re but dont need to remove |-| in reactome_fi as the signs are changed
reactome_fi = reactome_fi[reactome_fi$Direction!="<-",]
reactome_fi = reactome_fi[reactome_fi$Direction!="|-",]

#changing -> to 1 and -| to -1 in reactome fi
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("<->", 1, x, fixed = T)))
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("<-|", -1, x,fixed = T)))
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("|->", 1, x, fixed = T)))
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("|-|", -1, x,fixed = T)))
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("->", 1, x, fixed = T)))
reactome_fi = data.frame(lapply(reactome_fi, function(x) gsub("-|", -1, x,fixed = T)))
reactome_fi_f = reactome_fi[,c(1,2,4)]
names(reactome_fi_f) = c("Source","Target","Effect")

#rbining reactome_fi_f and reactome_re_f
reactome_directed = rbind(reactome_re_f,reactome_fi_f)

#human gene games in the database
h_gene_names=unique(c(as.vector(omni_net$source_genesymbol),as.vector(omni_net$target_genesymbol)))

#convert human gene symbol to mouse and vice versa in the database
#download human and mosue gene symbols from biomart
human_mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human_gene_names=getBM(attributes = c("hgnc_symbol"),mart=human_mart)
mouse_gene_names=getBM(attributes = c("mgi_symbol"),mart=mouse_mart)

#get mouse gene symbol for a input list of human symbol obtained from the omnipath database
genes_human_mus = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = as.vector(c(human_gene_names))$hgnc_symbol, mart = human_mart, attributesL = c("mgi_symbol"), martL = mouse_mart, uniqueRows=T)
#write the human and mouse symbols
write.table(genes_human_mus,"human_mouse_gene_symbol.tsv",sep="\t",row.names = F,quote=F)

#function for generating mouse source and target columns for omnipath interaction database, x-human mouse gene names df, y-reactomeDB with only human interactions
#some spurious interactions in the mouse will be created since same human gene name is mapped to two mouse names.
convertdbhuman2mouse_reactome <- function(x,y){
  colnames(x)[1]<-"Source"
  colnames(x)[2]<-"source_mgi_genesymbol"
  y1=join(x,y,by=c("Source"),type="right")
  colnames(x)[1]<-"Target"
  colnames(x)[2]<-"target_mgi_genesymbol"
  y1=join(x,y1,by=c("Target"),type="right")
  y1=(y1[,c(3,1,4,2,5)])
  names(y1)[1]<-c("source_hgnc_genesymbol")
  names(y1)[2]<-c("target_hgnc_genesymbol")
  y1
}

#function for generating mouse source and target columns for omnipath interaction database, x-human mouse gene names df, y-reactomeDB with only human interactions
#some spurious interactions in the mouse will be created since same human gene name is mapped to two mouse names.
convertdbhuman2mouse_metacore <- function(x,y){
  colnames(x)[1]<-"Source"
  colnames(x)[2]<-"source_mgi_genesymbol"
  y1=join(x,y,by=c("Source"),type="right")
  colnames(x)[1]<-"Target"
  colnames(x)[2]<-"target_mgi_genesymbol"
  y1=join(x,y1,by=c("Target"),type="right")
  y1=(y1[,c(4,2,5)])
  names(y1)[1]<-c("Source")
  names(y1)[2]<-c("Target")
  y1=na.omit(y1)
  return(y1)
}

#function for generating mouse source and target columns for omnipath interaction database, x-human mouse gene names df, y-omnipathDB with only human interactions
#some spurious interactions in the mouse will be created since same human gene name is mapped to two mouse names.
convertdbhuman2mouse <- function(x,y){
  colnames(x)[1]<-"source_genesymbol"
  colnames(x)[2]<-"source_mgi_genesymbol"
  y1=join(x,y,by=c("source_genesymbol"),type="right")
  colnames(x)[1]<-"target_genesymbol"
  colnames(x)[2]<-"target_mgi_genesymbol"
  y1=join(x,y1,by=c("target_genesymbol"),type="right")
  y1=(y1[,c(3,4,1,2,5:24)])
  names(y1)[1]<-c("source_hgnc_genesymbol")
  names(y1)[3]<-c("target_hgnc_genesymbol")
  y1
}

#function for generating mouse source and target columns for omnipath interaction database, x-human mouse gene names df, y-omnipathDB with only human interactions
#some spurious interactions in the mouse will be created since same human gene name is mapped to two mouse names.
convertdbmouse2human <- function(x,y){
  colnames(x)[1]<-"source_hgnc_genesymbol"
  colnames(x)[2]<-"source_genesymbol"
  y1=join(x,y,by=c("source_genesymbol"),type="right")
  colnames(x)[1]<-"target_hgnc_genesymbol"
  colnames(x)[2]<-"target_genesymbol"
  y1=join(x,y1,by=c("target_genesymbol"),type="right")
  y1=(y1[,c(3,4,1,2,5:24)])
  names(y1)[1]<-c("source_mgi_genesymbol")
  names(y1)[3]<-c("target_mgi_genesymbol")
  y1=(y1[,c(2,1,4,3,5:24)])
  y1
}

#converting human symbols to mouse in reactome
mous_human_reactome_net=convertdbhuman2mouse_reactome(genes_human_mus,reactome_directed)

#generate mouse gene names for omnipathDB
mous_omni_net=convertdbhuman2mouse(genes_human_mus,omnipath_human)
#generate human gene names for mouse in omnipathdb
hgnc_omni_net=convertdbmouse2human(genes_human_mus,omnipath_mouse)
#generate human genes for the rat genes
rat_omni_net=convertdbmouse2human(genes_human_mus,omnipath_rat)
#merge the the three parts: !both mouse and rat have the same gene symbol
omnipath_DB_full=rbind(hgnc_omni_net,mous_omni_net, rat_omni_net)
write.table(omnipath_DB_full,"omnipath_DB_mgi_hgnc_rat.tsv",sep="\t",row.names = F,quote=F)



#extract only directed edges
omni_directed=omnipath_DB_full[omnipath_DB_full$is_directed==1,]
#extract only directed interactions i.e. either stimulated or inhibited 
omni_stmi_inhi=omnipath_DB_full[omnipath_DB_full$is_inhibition==1|omnipath_DB_full$is_stimulation==1,]
write.table(omni_stmi_inhi,"omnipath_db_directed_signed_with_contradictions.tsv",sep="\t",row.names = F,quote=F)
#interaction with both activation and inhibition !i.e contradicatory signs so to be removed
omnipath_signed <- omnipath_DB_full[(omnipath_DB_full$is_inhibition==1&omnipath_DB_full$is_stimulation==0)|(omnipath_DB_full$is_inhibition==0&omnipath_DB_full$is_stimulation==1),]
write.table(omnipath_signed,"omnipath_db_directed_signed.tsv",sep="\t",row.names = F,quote=F)

#extract a 3 column network from the omnipath directed network, discard those with contradicting signs
#changing the inhibitory interaction into -1 from the default 1 notation
net=(omnipath_signed[,c(1,2,3,4,8,9)])
net$is_inhibition<-net$is_inhibition*(-1)
net$Effect <- net$is_stimulation+net$is_inhibition
mouse_omnipath_net=(net[,c(2,4,7)])
human_omnipath_net=(net[,c(1,3,7)])
names(mouse_omnipath_net)=c("Source","Target","Effect") #this is omnipath signaling interactome !TF-TF interactions must be removed
names(human_omnipath_net)=c("Source","Target","Effect") #this is omnipath signaling interactome !TF-TF interactions must be removed
mouse_omnipath_net=na.omit(mouse_omnipath_net) #removing interaction with NA (un mapped mouse)
human_omnipath_net=na.omit(human_omnipath_net)
write.table(mouse_omnipath_net,"omnipath_interactome_signed_directed_3col.tsv",sep="\t",quote=F,row.names = F)

#extracting only mouse reactome interactions
mouse_reactome_net = mous_human_reactome_net[,c(3:5)]
mouse_reactome_net = na.omit(mouse_reactome_net) #removing interaction with NA (un mapped mouse)
names(mouse_reactome_net) = c("Source","Target","Effect")

#extracting only human reactome interactions
human_reactome_net = mous_human_reactome_net[,c(1,2,5)]
human_reactome_net = na.omit(human_reactome_net) 
names(mouse_reactome_net) = c("Source","Target","Effect")
names(human_reactome_net) = c("Source","Target","Effect")

#merging mouse_omnipath and reactome
mouse_global_net = rbind(mouse_omnipath_net,mouse_reactome_net)
#unique interactions
mouse_global_net = unique(mouse_global_net)
#removing conflicting interaction of omnipath and reactome
mouse_global_net=mouse_global_net[!duplicated(mouse_global_net[1:2]),]

#merging human omnipath and reactome
human_global_net = rbind(human_omnipath_net,human_reactome_net)
#unique interactions
human_global_net = unique(human_global_net)
#removing conflicting interaction of omnipath and reactome
human_global_net=human_global_net[!duplicated(human_global_net[1:2]),]

#reading TFs, mouse TFs downloaded from http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download used only Tfs and not Tf cofactors
tf_mouse=read.table("Mus_musculus_TF.txt",sep="\t",header=T)
tf_mouse=as.vector(tf_mouse$Symbol)
#reading TF, human TFs downloaded from http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download used only TFs and not TF cofactors
tf_human=read.table("Homo_sapiens_TF.txt",sep="\t",header=T)
tf_human=as.vector(tf_human$Symbol)  
  
#removing TF-TF interactions
remove_tf_tf_interactions <- function (net,tf_list){
net_g <- graph.data.frame(as.data.frame(net))
TF_int<-intersect(tf_list,V(net_g)$name)
TF_int_g <- induced.subgraph(net_g,TF_int)
net_wo_TFs <-  (net_g)-(TF_int_g) #subtracting the TF-TF interactions
el=as_edgelist(net_wo_TFs)
global_net_sig=as.data.frame(cbind(el,E(net_wo_TFs)$Effect))
names(global_net_sig)=c("Source","Target","Effect") #network without TF-TF interactions
return(global_net_sig)
}

subset_tf_tf_interactions <- function (net,tf_list){
  net_g <- graph.data.frame(as.data.frame(net))
  TF_int<-intersect(tf_list,V(net_g)$name)
  TF_int_g <- induced.subgraph(net_g,TF_int)
  el=as_edgelist(TF_int_g)
  TF_TF_interactions=as.data.frame(cbind(el,E(TF_int_g)$Effect))
  names(TF_TF_interactions)=c("Source","Target","Effect")
  return(TF_TF_interactions)
}

#removing mouse TF-TF interactions from omnipath and reactome combined network
mouse_global_net_sig = remove_tf_tf_interactions(mouse_global_net,tf_mouse)

#removing human TF-TF interactions interactions from omnipath and reactome combined network
human_global_net_sig = remove_tf_tf_interactions(human_global_net,tf_human)

#TF-TF interactions in omnipath and reactome
mouse_tf_tf_interactions =  subset_tf_tf_interactions(mouse_global_net,tf_mouse)
human_tf_tf_interactions = subset_tf_tf_interactions(human_global_net,tf_human)
write.table(mouse_tf_tf_interactions,"mouse_TF_TF_interactions_omnipath_reactome.sif",sep="\t",quote=F,row.names=F)
write.table(human_tf_tf_interactions,"human_TF_TF_interactions_omnipath_reactome.sif",sep="\t",quote=F,row.names=F)

#identifying interface TFs
subset_int_tfs <- function(net,tf_list){
global_net_sig_g <- graph.data.frame(as.data.frame(net))
int_tfs <- intersect(tf_list,V(global_net_sig_g)$name)
return(int_tfs)
}

#mouse int_TF
mouse_int_tfs = subset_int_tfs(mouse_global_net_sig,tf_mouse)
#human int_tfs
human_int_tfs = subset_int_tfs(human_global_net_sig,tf_human)

#Connect all int_TFs to Dummy and mark them as activaation
connect_int_tfs <- function(int_tfs){
dummy_vec=rep("Dummy",length(int_tfs))
TF_dummy_net=as.data.frame(cbind(int_tfs,dummy_vec,rep(1,length(int_tfs))))
names(TF_dummy_net)=c("Source","Target","Effect")
return(TF_dummy_net)
}

#Connect all int_TFs to Dummy and mark them as activaation for mouse and human
mouse_tf_dummy_net = connect_int_tfs(mouse_int_tfs)
human_tf_dummy_net = connect_int_tfs(human_int_tfs)

#connect Dummy node to all receptors and ligands, ligands and receptors are mouse orthlologs of human ligands and receptors from draft paper and plasma membrane gene ontology
#read the list of receptors and ligands
lig_receptors_mouse=read.table("ligand_receptor_mouse_090818.txt")
lig_receptors_human=read.table("ligand_receptor_human_090818.txt")

connect_lig_receptors <- function(lig_receptors_list,net){
lig_receptors_list <- intersect(as.vector((lig_receptors_list$V1)),V(graph.data.frame(as.data.frame(net)))$name)
dummy_vec=rep("Dummy",length(lig_receptors_mouse))
lig_rec_dummy_net=as.data.frame(cbind(dummy_vec,lig_receptors_mouse,rep(1,length(lig_receptors_mouse))))
names(lig_rec_dummy_net)=c("Source","Target","Effect")
return(lig_rec_dummy_net)
}

mouse_lig_rec_dummy_net = connect_lig_receptors(lig_receptors_mouse,mouse_global_net)
human_lig_rec_dummy_net = connect_lig_receptors(lig_receptors_human,human_global_net)

#Import kinases, phosphatases and merge with ligands and receptors for the set of signaling intermediates/molecules
#mouse kinases and phosphatases are downloaded from http://ekpd.biocuckoo.org/faq.php
kin_phos_mouse=read.table("mouse_kinases_phosphatases_07082018.txt")
kin_phos_human=read.table("human_kinases_phosphatases_07082018.txt")

#intermediates mouse, includes receptors, ligands, kinases and phosphatases present in the omnipath DB
subset_intermediates <- function(kin_phos_list,lig_receptors_list){ 
intermediates <- rbind(kin_phos_list,lig_receptors_list)
names(intermediates) <- "Gene"
return(intermediates)
}

intermediates_mouse=unique(subset_intermediates(kin_phos_mouse,lig_receptors_mouse))
intermediates_human=unique(subset_intermediates(kin_phos_human,lig_receptors_human))
#write the files
write.table(intermediates_mouse,"intermediates_mouse.txt",quote = F, row.names = F)

#non-interface TFs from metacore in human
#reading Metacore HUman-TF interactions
metacore_tf_interactions_human <- read.table("data_TransReg_TFHairball.txt")
#reading Metacore mouse-TF interactions
metacore_tf_interactions_mouse <- read.table("data_TransReg_TFHairball_mouse.txt")

#removing unspecified interactions and changing activation to 1 and inhibition to -1
process_metacore_intractions <- function(metacore_tf_interactions) {
metacore_tf_interactions <- metacore_tf_interactions[metacore_tf_interactions$V2!="Unspecified",]
metacore_tf_interactions=data.frame(lapply(metacore_tf_interactions, function(x) gsub("Activation", 1, x)))
metacore_tf_interactions=data.frame(lapply(metacore_tf_interactions, function(x) gsub("Inhibition", -1, x)))
metacore_tf_interactions=metacore_tf_interactions[c(1,3,2)]
names(metacore_tf_interactions)=c("Source","Target","Effect")
return(metacore_tf_interactions)
}

metacore_tf_interactions_human <- process_metacore_intractions(metacore_tf_interactions_human)
metacore_tf_interactions_mouse <- process_metacore_intractions(metacore_tf_interactions_mouse)

#merging metacore TF interactions and Omnipath/reactome TF interactions
subset_global_tf_tf_interaction <- function(tf_tf_ints,metacore_tf_ints){
global_TF_TF_interactions = rbind(tf_tf_ints, metacore_tf_ints)
#unique interactions
global_TF_TF_interactions = unique(global_TF_TF_interactions)
#removing conflicting interaction of omnipath and reactome
global_TF_TF_interactions=global_TF_TF_interactions[!duplicated(global_TF_TF_interactions[1:2]),]
return(global_TF_TF_interactions)
}

#global TF-interactions
human_global_tf_tf_interactions = subset_global_tf_tf_interaction(human_tf_tf_interactions,metacore_tf_interactions_human)

#convert human metacore to mouse metacore ##CONVERTING HUMAN TF NAMES TO MOUSE!!
#metacore_tf__tf_interactions_mouse=convertdbhuman2mouse_metacore(genes_human_mus,metacore_tf_interactions)

mouse_global_tf_tf_interactions = subset_global_tf_tf_interaction(mouse_tf_tf_interactions,metacore_tf_interactions_mouse)

#non-interface TFs in the global network
subset_non_int_tfs <- function(tf_list,global_tf_interactions,int_tfs){
TF_TF_interactions_g <- graph.data.frame(as.data.frame(global_tf_interactions))
non_int_tfs <- intersect(tf_list,V(TF_TF_interactions_g)$name)
non_int_tfs <- as.data.frame(setdiff(non_int_tfs,int_tfs))
names(non_int_tfs) <- "Gene"
return(non_int_tfs)
}

mouse_non_int_tfs = subset_non_int_tfs(tf_mouse,mouse_global_tf_tf_interactions,mouse_int_tfs)
human_non_int_tfs = subset_non_int_tfs(tf_human,human_global_tf_tf_interactions,human_int_tfs)
#write.table(non_int_tfs,"NON_interface_TF.txt", quote = F, row.names = F)

#Merge the signaling network with dummy node, dummy node is connected to receptors/ligands and int_TFs are connected to Dummy
global_net_sig_dummy_mouse <- rbind(mouse_global_net_sig,mouse_lig_rec_dummy_net,mouse_tf_dummy_net)

global_net_sig_dummy_human <- rbind(human_global_net_sig,human_lig_rec_dummy_net,human_tf_dummy_net)
#write.table(global_net_sig_dummy,"Mouse_signaling_interactome.sif", sep="\t", row.names = F,quote=F)

#renaming the elements for saving a .RData

Background_signaling_interactome <- global_net_sig_dummy_mouse
intermediates <- intermediates_mouse
non_interface_TFs <- mouse_non_int_tfs
TF_TF_interactions <- mouse_global_tf_tf_interactions


#Saving the four input files for NicheSIG as RData
save(Background_signaling_interactome, intermediates, non_interface_TFs, TF_TF_interactions, file = "MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData")

#renaming the elements for saving a .RData

Background_signaling_interactome <- global_net_sig_dummy_human
intermediates <- intermediates_human
non_interface_TFs <- human_non_int_tfs
TF_TF_interactions <- human_global_tf_tf_interactions

save(Background_signaling_interactome, intermediates, non_interface_TFs, TF_TF_interactions, file = "HUMAN_Background_Network_omnipath_reactome_metacore_01042019.RData")





