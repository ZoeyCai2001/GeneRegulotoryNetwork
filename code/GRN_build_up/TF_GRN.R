############################
#datacreation
############################
.libPaths("/home/zc354/R/4.2")
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

#####################################
library(GENIE3)
hvg1=readRDS("/gpfs/gibbs/pi/zhao/yw599/GRN/Keran/matrix_openchromt_hvg.Rds")
M1=readRDS("/gpfs/gibbs/pi/zhao/yw599/GRN/Keran/matrix_openchromt_TF.Rds")
gene_names=colnames(hvg1) #目标基因们
#gene_names chr[1:556]
#hvg1 dgCMatrix[23133*556]
#M1 dgCMatrix[23133*776]


for(gene in gene_names[1:40]){
  w=which(hvg1[,gene]!=0)
  atac=(rownames(hvg1))[w] #目标基因TSS周围的开放区间
  
  wh=which(M1[atac[1],]!=0)
  TF=(colnames(M1))[wh]
  if(length(atac)>1){
    for(i in 2:length(atac)){
      wh=which(M1[atac[i],]!=0)
      TF=union(TF,(colnames(M1))[wh])
    }
  }

  
  data=rbind(pbmc.rna@assays$RNA@data[gene,],pbmc.rna@assays$RNA@data[TF,])
  # data=NormalizeData(data,scale.factor = 10000)
  data=as.matrix(data) #是否normalization
  #  > class(data)
  #  [1] "dgCMatrix"
  #  attr(,"package")
  #  [1] "Matrix"
  rownames(data)[1]=gene
  result=GENIE3(data,regulators = TF,targets = gene,nCores = 4)
  #unable to find an inherited method for function 'GENIE3' for signature '"dgCMatrix"
  result=getLinkList(result)
  write.table(result,file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF__filename_1_40.csv",append = T,col.names = F,row.names = F,sep=",")
  #注意读取的时候header=F
}
t=read.table("/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF__filename_1_40.csv",header = F,sep=',')
colnames(t) = c("Gene1","Gene2","Edgeweight")  #与Beeline保持一致
write.table(t,file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF__filename_1_40.csv",row.names = F,sep=",")