#library
.libPaths("/home/zc354/R/4.2")
library(dplyr)
library(tidyr)
library(tidyverse)
library(base)
library(do)
library(SparseM)

#parameters
args<-commandArgs(trailingOnly=TRUE)
data<-args[1]
gene_peak<-args[2]
model<-args[3]
start<-args[4]

savepath = paste("/home/zc354/GRN/output/trio/",data,"_",model,"_",gene_peak,sep="")
filename = paste("/dif_",start,".csv",sep="") 

#nonfilter <- read.table('/home/zc354/GRN/output/predicted_network/non_filter_sct_count_pbmc_500kb_motifmatchr.csv',sep=',',header=T)
#hvg_df = nonfilter[!duplicated(nonfilter$Gene2),]
#hvgs<-hvg_df[,'Gene2']
#saveRDS(hvgs,'/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/hvg_pbmc.Rds')

#lp_pbmc<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/LinkPeaks_0.05_0.05_SCT_matrix.Rds')
#saveRDS(lp_pbmc,'/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/pbmc_gene_peak_linkpeaks.Rds')

if(data=='pbmc'){
  #checked
  hvgs<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/hvg_pbmc.Rds')
  #2000
  #checked
  data_used<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/pbmc_rna.Rds')
  atac_used<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/pbmc_atac.Rds')
  #checked
  M1<- readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_TF_106120_1119.Rds")
  if(gene_peak =='500kb'){
    #checked
    hvg1<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/pbmc_gene_peak_500kb.Rds')
  }else{
    hvg1<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/pbmc_gene_peak_linkpeaks.Rds')
  }
}else{
  #checked
  hvgs<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/hvg_bmmc.Rds')
  #1608
  #checked
  data_used<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_rna.Rds')
  atac_used<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_atac.Rds')
  #checked
  M1<- readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_binary_matrix_bmmc.Rds")
  if(gene_peak=='500kb'){
    #checked
    hvg1<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_gene_peak_500kb.Rds')
  }
  else{
    #hvg1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_sctransform_8899_1870.Rds")
    #saveRDS(hvg1,'/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_gene_peak_linkpeaks.Rds')
    hvg1<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_gene_peak_linkpeaks.Rds')
  }
}

genes = base::union(hvgs,colnames(M1))
data_used = data_used[genes,]
hvgs <- base::intersect(hvgs,colnames(hvg1))
hvg1 <- hvg1[,hvgs]
w = which(rowSums(hvg1)==0)
if(length(w)!=0){
  hvg1<-hvg1[-w,]
}
atac = base::intersect(rownames(hvg1),rownames(M1))
hvg1<-hvg1[atac,]
M1<-M1[atac,]
atac_used = atac_used[atac,]

#gene_names = colnames(hvg1)
gene_names = readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/output/trio/bmmc_500kb_gene_diff.Rds')
length_gene = length(gene_names)
TF = colnames(M1)
start = as.numeric(start)
#if(length_gene%%2 ==0){
#  end=start+1
#}else{
#  if(start==(length_gene-length_gene%%2+1)){
#    end = length_gene
#  }else
#    end = start+1
#}
end = start
#####################################
library(GENIE3)
s=start
for(gene in gene_names[start:end]){
  print(paste(gene,s,sep="_"))
  y=0 #用于命名
  w=which(hvg1[,gene]!=0)
  if(length(w)==0){
    s=s+1
    next
  }
  atac=(rownames(hvg1))[w]
  atc=atac_used[atac,]
  wh=which(M1[atac[1],]!=0)
  TF=(colnames(M1))[wh]
  if(length(w)==1){
    if(length(wh)==1){
      s=s+1
      next
    }
  }
  
  if(gene %in% TF){
    w_TF = which(TF[]==gene)
    TF_used = TF[-w_TF]
  }else{
    TF_used = TF
  }
  
  if(length(w)==1){
    if(length(TF_used)==1){
      y=1
      rna=data_used[TF_used,]
      multiply=rna*atc
    }
    
    if(length(TF_used)>1){
      rna=data_used[TF_used,]
      for(k in 1:length(TF_used)) {
        rna[k,]=rna[k,]*atc
        rownames(rna)[k]=paste(TF_used[k],atac[1],sep="*")
      } #命名
      multiply=rna
    }
  }
  
  
  
  if(length(w)>1){
    if(length(TF_used)==1){
      y=1
      rna=data_used[TF_used,]
      multiply=rna*atc[1,]
    }
    
    if(length(TF_used)>1){
      rna=data_used[TF_used,]
      for(k in 1:length(TF_used)) {rna[k,]=rna[k,]*atc[1,]
      rownames(rna)[k]=paste(TF_used[k],atac[1],sep="*")} #命名
      multiply=rna
    }
  }
  
  if(length(w)>1){
    for(i in 2:length(atac)){
      wh=which(M1[atac[i],]!=0)
      TF=(colnames(M1))[wh]
      if(gene %in% TF){
        w_TF = which(TF[]==gene)
        TF_used = TF[-w_TF]
      }else{
        TF_used = TF
      }
      
      if(length(TF_used)==1){
        rna=data_used[TF_used,]
        multiply=rbind(multiply,rna*atc[i,])
        n=nrow(multiply) #命名
        rownames(multiply)[n]=paste(TF_used,atac[i],sep="*")
      }
      
      
      if(length(TF_used)>1){
        print("!")
        rna=data_used[TF_used,]
        for(k in 1:length(TF_used)) {rna[k,]=rna[k,]*atc[i,]
        rownames(rna)[k]=paste(TF_used[k],atac[i],sep="*") #命名
        }
        multiply=rbind(multiply,rna)
      }
    }
  }
  
  #命名
  if(y==1){
    wh=which(M1[atac[1],]!=0)
    TF_used=(colnames(M1))[wh]
    rownames(multiply)[1]=paste(TF_used,atac[1],sep="*")   
  }
  
  data=rbind(data_used[gene,],multiply)
  data=as.matrix(data)
  rownames(data)[1]=gene
  result=GENIE3(data,regulators = rownames(data)[2:nrow(data)],targets = gene,nCores = 4)
  result=getLinkList(result)
  result$regulatoryGene=as.character(result$regulatoryGene)
  result = separate(result,'regulatoryGene',c('peak','TF'),sep = "[*]",remove=TRUE)
  write.table(result,file = paste(savepath,filename,sep = ""),append = T,col.names = F,row.names = F,sep=",")
  s=s+1
}
