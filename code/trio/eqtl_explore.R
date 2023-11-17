.libPaths("/home/zc354/R/4.2")
library(dplyr)
library(tidyr)
library(tidyverse)
library(base)
#library(do)
library(SparseM)
library(modEvA)

args<-commandArgs(trailingOnly=TRUE)
data<-args[1]
gene_peak<-args[2]
model<-args[3]

dir = "/home/zc354/GRN/output/trio/"
writedir = '/home/zc354/GRN/output/trio/'
filename = paste(dir,data,"_",model,"_",gene_peak,".csv",sep="")

#read eqtl
print("read eqtl")
eqtl = read.table(gzfile("/home/zc354/GRN/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"),sep='\t',header = T) 
eqtl = eqtl[,c('variant_id','gene_id','tss_distance','pval_nominal')]
eqtl = eqtl[which(abs(eqtl[,'tss_distance'])<=500000),]

#read egenes
print("read egenes")
egenes =read.table(gzfile("/home/zc354/GRN/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"),sep='\t',header = T)
egenes = egenes[,c('gene_id','gene_name','gene_chr','gene_start','gene_end','strand')]

#join egenes with eqtls
print("join egenes with eqtls")
eqtl<-eqtl %>%
  left_join(egenes, by = "gene_id")
eqtl<- na.omit(eqtl)

#genes and variants in eqtl
print('generate genes and variants in eqtl')
eqtl_genes = eqtl%>% distinct(gene_name, .keep_all=T)
eqtl_genes = eqtl_genes[,'gene_name']
eqtl_variant = eqtl%>% distinct(variant_id, .keep_all=T)
eqtl_variant = eqtl_variant[,'variant_id']

#read predicted network
print('read predicted network')
predicted_network = read.table(filename,header = T,sep=',')
colnames(predicted_network)<-c("TF","peak","Target_Gene",'Edgeweight')
w=which(predicted_network[,'Edgeweight']==0)
predicted_network = predicted_network[-w,]
#separate peaks
print('separate peaks')
predicted_network = separate(predicted_network,'peak',c('chr','start','end'),sep = "-",remove=F)
predicted_network['start']<-lapply(predicted_network['start'],as.numeric)
predicted_network['end']<-lapply(predicted_network['end'],as.numeric)

print('genes and peaks in predicted network')
hvgs_predicted_network = (predicted_network%>% distinct(Target_Gene, .keep_all = T))
hvgs_predicted_network = hvgs_predicted_network[,'Target_Gene']
peaks_predicted_network = predicted_network%>% distinct(peak, .keep_all = T)
peaks_predicted_network = peaks_predicted_network[,'peak']
TFs_predicted_network = predicted_network %>% distinct(TF,.keep_all = T)
TFs_predicted_network = TFs_predicted_network[,'peak']

print('filter eqtl')
w=which(eqtl[,'gene_name'] %in% hvgs_predicted_network)
eQTL_filtered = eqtl[w,]

print('variants and genes in filtered eqtl')
variant_filter <- eQTL_filtered %>% distinct(variant_id, .keep_all = T)
variant_filter <- variant_filter[,'variant_id']
#173368
gene_names_eqtl <- eQTL_filtered %>% distinct(gene_name, .keep_all = T)
gene_names_eqtl<-gene_names_eqtl[,'gene_name']

predicted_network[,'A']=0
predicted_network[,'B']=0

gene_onlyinpred <-dplyr::setdiff(hvgs_predicted_network,gene_names_eqtl)
w =which(predicted_network[,'Target_Gene'] %in% gene_onlyinpred)
#4884198
peak_matrix = predicted_network[-w,]
peak_matrix = peak_matrix%>% distinct(peak, .keep_all = T)
rownames(peak_matrix) = peak_matrix[,'peak']
#45525

eQTL_sep = separate(eQTL_filtered,'variant_id',c('chr','locus','else1','else2','else3'),sep = "_",remove=F)
eQTL_sep = eQTL_sep[,c('gene_name','variant_id','chr','locus','pval_nominal','tss_distance')]
eQTL_sep['locus']<-lapply(eQTL_sep['locus'],as.numeric)

variant_matrix = eQTL_sep %>% distinct(variant_id, .keep_all = T)
#173368
rownames(variant_matrix) = variant_matrix[,'variant_id']
chr = variant_matrix %>% distinct(chr, .keep_all =T)
chr = chr[,'chr']
t=list()
for(c in chr){
  w_variant=which(variant_matrix[,'chr']==c)
  w_peak = which(peak_matrix[,'chr']==c)
  t[[c]]=matrix(0,length(w_variant),length(w_peak))
  rownames(t[[c]])=variant_matrix[w_variant,'variant_id']
  colnames(t[[c]])=peak_matrix[w_peak,'peak']
  t[[c]] = cbind(variant_matrix[w_variant,'locus'],t[[c]])
  colnames(t[[c]])[1]='locus'
  peak_matrix_used = cbind(matrix(0,2,1),t(peak_matrix[w_peak,c('start','end')]))
  colnames(peak_matrix_used)[1] = 'locus'
  t[[c]] = rbind(peak_matrix_used,t[[c]])
  peaks = peak_matrix[w_peak,'peak']
  variants = variant_matrix[w_variant,'variant_id']
  s=0
  for(p in peaks){
    start = t[[c]]['start',p]
    end = t[[c]]['end',p]
    w=which(t[[c]][,'locus']>=start & t[[c]][,'locus']<=end)
    if(length(w)!=0){
      t[[c]][w,p]=1
    }
    s=s+1
    print(s)
  }
}

saveRDS(t,paste(writedir,"gene_variant/",data,"_",model,"_",gene_peak,".csv",sep=""))

w1=which(predicted_network[,'Target_Gene'] %in% gene_names_eqtl)
peak_1 = c()
for(c in chr){
  w=which((colSums(t[[c]])[2:dim(t[[c]])[2]]-t[[c]]['start',2:dim(t[[c]])[2]]-t[[c]]['end',2:dim(t[[c]])[2]])==0)
  peak_0 = colnames(t[[c]])[-w]
  peak_1 = c(peak_1,peak_0)
}
variant_1 = c()
for(c in chr){
  w=which((rowSums(t[[c]])[3:dim(t[[c]])[1]]-t[[c]][3:dim(t[[c]])[1],'locus'])==0)
  variant_0 = rownames(t[[c]])[-w]
  variant_1 = c(variant_0,variant_1)
}
w2 = which(predicted_network[,'peak'] %in% peak_1)
w=dplyr::intersect(w1,w2)
predicted_network_eqtl = predicted_network[w,]
predicted_network_0 = predicted_network[-w,]
predicted_network_eqtl[,'B']=1


for(i in 1:dim(predicted_network_eqtl)[1]){
  gene = predicted_network_eqtl[i,'Target_Gene']
  p = predicted_network_eqtl[i,'peak']
  c = predicted_network_eqtl[i,'chr']
  variants = eQTL_filtered[which(eQTL_filtered[,'gene_name']==gene),'variant_id']
  if(sum(t[[c]][variants,p])!=0){
    predicted_network_eqtl[i,'A']=1
  }
  if(i%%1000==0){
    print(i)
  }
}
predicted_network = rbind(predicted_network_eqtl,predicted_network_0)
output = data.frame('model' = model, 'data' = data,'gene_peak'=gene_peak,
                    '#of hvgs' = length(hvgs_predicted_network),'#of genes in eQTLs' = length(eqtl_genes),
                    '#of intersected genes' = length(dplyr::intersect(hvgs_predicted_network,gene_names_eqtl)),'#of peaks'=length(peaks_predicted_network),
                    '#of variants' = length(eqtl_variant), '#of peaks with variants' = length(peak_1), 
                    '#of variants on peaks' = length(variant_1),
                    '#of positive trio relations' = dim(predicted_network)[1],
                    '#of eQTLs' = dim(eqtl)[1],
                    '#of trio relations not in background eqtl' = dim(predicted_network_0)[1],
                    '#of trio relations in background eqtl' = dim(predicted_network_eqtl)[1],
                    '#of trio relations in significant eQTLs' = sum(predicted_network[,'A'])
)

saveRDS(predicted_network,paste(writedir,"predicted_network_eqtl/",data,"_",model,"_",gene_peak,'.Rds',sep=''))

write.table(output,paste(writedir,"predicted_network_eqtl/summary.csv",sep=''),row.names = FALSE,append = T,sep=',')
