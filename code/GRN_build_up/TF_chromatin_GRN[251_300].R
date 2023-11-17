############################
#datacreation
############################
.libPaths("/home/zc354/R/4.2")
library(do)
library(plyr)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
counts <- Read10X_h5("/gpfs/gibbs/pi/zhao/zc354/GRN/data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "/gpfs/gibbs/pi/zhao/zc354/GRN/data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc.rna <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# create a Seurat object containing the ATAC adata
pbmc.atac <- CreateSeuratObject(
  counts = counts$Peaks,
  assay = "ATAC",
  fragments = fragpath,
  annotation = annotation
)

#change the sep ":" to "-"
dimnames<-pbmc.atac@assays[["ATAC"]]@counts@Dimnames[[1]]
dimnames<-Replace(dimnames,":","-")
pbmc.atac@assays[["ATAC"]]@counts@Dimnames[[1]]<-dimnames

library(GENIE3)

hvg1=readRDS("/gpfs/gibbs/pi/zhao/yw599/GRN/Keran/matrix_openchromt_hvg.Rds")
M1=readRDS("/gpfs/gibbs/pi/zhao/yw599/GRN/Keran/matrix_openchromt_TF.Rds")
gene_names=colnames(hvg1) #目标基因们
#gene_names chr[1:556]
#hvg1 dgCMatrix[23133*556]
#M1 dgCMatrix[23133*776]
for(gene in gene_names[251:300]){
  y=0 #用于命名
  w=which(hvg1[,gene]!=0)
  atac=(rownames(hvg1))[w] #"chr1-816881-817647"
  atc=pbmc.atac@assays$ATAC@counts[atac,]
  wh=which(M1[atac[1],]!=0)
  TF=(colnames(M1))[wh]
  
  
  if(length(w)==1){
    if(length(wh)==1){
      y=1
      rna=pbmc.rna@assays$RNA@counts[TF,]
      multiply=rna*atc
    }
    
    if(length(wh)>1){
      rna=pbmc.rna@assays$RNA@counts[TF,]
      for(k in 1:length(wh)) {
        rna[k,]=rna[k,]*atc
        rownames(rna)[k]=paste(TF[k],atac[1],sep="*")
      } #命名
      multiply=rna
    }
  }
  
  
  
  if(length(w)>1){
    if(length(wh)==1){
      y=1
      rna=pbmc.rna@assays$RNA@counts[TF,]
      multiply=rna*atc[1,]
    }
    
    if(length(wh)>1){
      rna=pbmc.rna@assays$RNA@counts[TF,]
      for(k in 1:length(wh)) {rna[k,]=rna[k,]*atc[1,]
      rownames(rna)[k]=paste(TF[k],atac[1],sep="*")} #命名
      multiply=rna
    }
  }
  
  if(length(w)>1){
    for(i in 2:length(atac)){
      wh=which(M1[atac[i],]!=0)
      TF=(colnames(M1))[wh]
      
      if(length(wh)==1){
        rna=pbmc.rna@assays$RNA@counts[TF,]
        multiply=rbind(multiply,rna*atc[i,])
        n=nrow(multiply) #命名
        rownames(multiply)[n]=paste(TF,atac[i],sep="*")
      }
      
      
      if(length(wh)>1){
        print("!")
        rna=pbmc.rna@assays$RNA@counts[TF,]
        for(k in 1:length(wh)) {rna[k,]=rna[k,]*atc[i,]
        rownames(rna)[k]=paste(TF[k],atac[i],sep="*") #命名
        }
        multiply=rbind(multiply,rna)
      }
      
      
    }
  }
  
  #命名
  if(y==1){
    wh=which(M1[atac[1],]!=0)
    TF=(colnames(M1))[wh]
    rownames(multiply)[1]=paste(TF,atac[1],sep="*")   
  }
  
  
  
  
  data=rbind(pbmc.rna@assays$RNA@data[gene,],multiply)
  data=NormalizeData(data,scale.factor = 10000)
  data=as.matrix(data)
  rownames(data)[1]=gene
  result=GENIE3(data,regulators = rownames(data)[2:nrow(data)],targets = gene,nCores = 4)
  result=getLinkList(result)
  result$regulatoryGene=as.character(result$regulatoryGene)
  a= strsplit(result$regulatoryGene[1],split="[*]")
  TFs=a[[1]][1]
  result$regulatoryGene[1]=TFs
  write.table(result[1,],file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",append = T,col.names = F,row.names = F,sep=",")
  
  for(i in 2:nrow(result)){
    #在TF*multiply变量中寻找不一样的TF，并加入TFs,因为要找的是哪些TF调控基因
    a= strsplit(result$regulatoryGene[i],split="[*]")
    if(!(a[[1]][1] %in% TFs)){ 
      
      result$regulatoryGene[i]=a[[1]][1] 
      TFs=c(TFs,a[[1]][1] )
      write.table(result[i,],file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",append = T,col.names = F,row.names = F,sep=",")  }
    
  }
  
}
t=read.table("/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",header = F,sep=',')
colnames(t) = c("Gene1","Gene2","Edgeweight")
write.table(t,file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",row.names = F,sep=",")





#TF 对应weight相加
t=read.table("/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",header = T,sep=',')
library(plyr)
t=aggregate(Edgeweight~Gene1+Gene2 , data=t , sum) #不能for循环， 慢死了
w=which(t$Edgeweight==0)
t=t[-w,]
t=t[order(-t$Edgeweight),] #weight降序排列
write.table(t,file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename_251_300.csv",row.names = F,sep=",")
