.libPaths("/home/zc354/R/4.2")
library(dplyr)
library(tidyr)
library(tidyverse)
library(base)
#library(do)
library(SparseM)
library(modEvA)

eqtl = read.table(gzfile("/home/zc354/GRN/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"),sep='\t',header = T) 
eqtl = eqtl[,c('variant_id','gene_id','tss_distance','pval_nominal')]

print("read egenes")
egenes =read.table(gzfile("/home/zc354/GRN/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"),sep='\t',header = T)
egenes = egenes[,c('gene_id','gene_name','gene_chr','gene_start','gene_end','strand')]
egenes = egenes%>% distinct(gene_name, .keep_all=T)
colnames(egenes) = c('gene_id','gene_name','chr','start','end','strand')
gene_matrix = egenes
rownames(gene_matrix) = egenes[,'gene_name']

print('generate genes and variants in eqtl')
eqtl_variant = eqtl%>% distinct(variant_id, .keep_all=T)
eqtl_variant = eqtl_variant['variant_id']
eqtl_variant = separate(eqtl_variant,'variant_id',c('chr','locus','else1','else2','else3'),sep = "_",remove=F)
eqtl_variant = eqtl_variant[,c('variant_id','chr','locus')]
eqtl_variant['locus']<-lapply(eqtl_variant['locus'],as.numeric)
variant_matrix=eqtl_variant
chr = variant_matrix %>% distinct(chr, .keep_all =T)
chr = chr[,'chr']

t=list()
for(c in chr){
  print(c)
  w_variant=which(variant_matrix[,'chr']==c)
  w_gene = which(gene_matrix[,'chr']==c)
  variants_locus = variant_matrix[w_variant,'locus']
  gene_names = gene_matrix[w_gene,'gene_name']
  t[[c]]=matrix(rep(variants_locus,length(w_gene)),ncol = length(w_gene))
  rownames(t[[c]])=variant_matrix[w_variant,'variant_id']
  colnames(t[[c]])=gene_names
  s=0
  for(gene in gene_names){
    if(gene_matrix[gene,"strand"]=='-'){
      tss = gene_matrix[gene,"end"]
      t[[c]][,gene] = t[[c]][,gene]-tss
    }else{
      tss = gene_matrix[gene,"start"]
      t[[c]][,gene] = t[[c]][,gene]-tss
    }
    s=s+1
    print(s)
  }
}
saveRDS(t,'/gpfs/gibbs/pi/zhao/zc354/GRN/output/gene_SNP_distance.Rds')

t<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/output/gene_SNP_distance.Rds')
s=0
for(c in chr){
  s=s+sum(t[[c]]<500000)
}
print(s)
  