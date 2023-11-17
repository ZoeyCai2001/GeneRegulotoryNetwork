.libPaths("/home/zc354/R/4.2")
library(Signac)
library(tidyr)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

savepath = "/gpfs/gibbs/pi/zhao/zc354/GRN/data/bmmc_linkpeaks.Rds"
p.value.cutoff = 0.05
score.cutoff = 0.05
expression.form = "SCT"

bmmc<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/sct_bmmc.Rds")
DefaultAssay(bmmc) <- "ATAC"

# first compute the GC content for each peak
bmmc <- RegionStats(bmmc, genome = BSgenome.Hsapiens.UCSC.hg38)

print("start LinkPeaks")
# link peaks to genes
bmmc <- LinkPeaks(
  object = bmmc,
  peak.assay = "ATAC",
  expression.assay = expression.form,
  expression.slot = "data",
  pvalue_cutoff = p.value.cutoff,
  score_cutoff = score.cutoff
)

saveRDS(bmmc,savepath)

bmmc<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/bmmc_linkpeaks.Rds")

link_peak_gene<-Links(bmmc)
link_peak_gene_df = as.data.frame(link_peak_gene)
link_peak_gene_df = link_peak_gene_df[,c("gene","peak","score")]

gene_df = link_peak_gene_df

gene_df_drop = gene_df[!duplicated(gene_df$gene),]['gene']
gene_names = gene_df_drop[,"gene"]
gene_length = length(gene_df_drop[,"gene"])
gene_df_drop[,'gene_number'] = c(1:gene_length)

atac_df = gene_df[!duplicated(gene_df$peak),]['peak']
atac_names = atac_df[,"peak"]
atac_length = length(atac_df[,"peak"])
atac_df[,'atac_number'] = c(1:atac_length)

gene_df_join= left_join(gene_df,gene_df_drop,by='gene')
gene_df_join= left_join(gene_df_join,atac_df,by='peak')

p<- gene_df_join[,'gene_number']
q<-gene_df_join[,'atac_number']
x<- seq(length=link_num_1, from=1, to=1)
TF_peak_matrix_dgc<- sparseMatrix(q,p,x=x)
rownames(TF_peak_matrix_dgc) = atac_names
colnames(TF_peak_matrix_dgc) = gene_names

saveRDS(TF_peak_matrix_dgc,"/gpfs/gibbs/pi/zhao/zc354/GRN/output/gene_peak_bmmc_linkpeaks.Rds")

