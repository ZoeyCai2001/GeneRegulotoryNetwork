library(dplyr)
library(tidyr)
library(tidyverse)
library(base)
#library(do)
library(SparseM)

args<-commandArgs(trailingOnly=TRUE)
data<-args[1]
gene_peak<-args[2]
model<-args[3]

dir = "/home/zc354/GRN/output/trio/"
writedir = '/home/zc354/GRN/output/trio/'
filename = paste(dir,data,"_",model,"_",gene_peak,".csv",sep="")

predicted_network = readRDS(paste(writedir,"predicted_network_eqtl/",data,"_",model,"_",gene_peak,'.Rds',sep=''))

predicted_network = predicted_network[order(-predicted_network[,'Edgeweight']),]

risk_ratio = data.frame(matrix(0,dim(predicted_network)[1],5))
colnames(risk_ratio) = c("positive","in significant eqtls","in background eqtls","thresholds","risk ratio")
risk_ratio[,"positive"]=c(1:dim(predicted_network)[1])
risk_ratio[,'thresholds'] = predicted_network[,'Edgeweight']

for(i in 2:dim(predicted_network)[1]){
  predicted_network[i,"A"]=sum(predicted_network[(i-1):i,'A'])
  predicted_network[i,"B"]=sum(predicted_network[(i-1):i,'B'])
  if(i%%1000==0){
    print(i)
  }
}

risk_ratio[,'in significant eqtls'] = predicted_network[,"A"]
risk_ratio[,'in background eqtls'] = predicted_network[,"B"]

significant_eqtls = 2214853
num_eqtl_variant = 1223408
num_eqtl_genes = 11902
risk_ratio[,'precision'] = risk_ratio[,'in significant eqtls']/risk_ratio[,'in background eqtls'] 
risk_ratio[,'risk_ratio'] = risk_ratio[,'precision']/significant_eqtls*num_eqtl_variant*num_eqtl_genes

saveRDS(risk_ratio,paste(writedir,"predicted_network_eqtl/risk_ratio_",data,"_",model,"_",gene_peak,'.Rds',sep=''))