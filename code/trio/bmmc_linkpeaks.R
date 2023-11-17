library(Signac)
library(Seurat)

bmmc<-readRDS('/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc.Rds')
DefaultAssay(bmmc) <- "ATAC"

# first compute the GC content for each peak
bmmc <- RegionStats(bmmc, genome = BSgenome.Hsapiens.UCSC.hg38)

print("start LinkPeaks")
# link peaks to genes
bmmc_lp <- LinkPeaks(
  object = bmmc,
  peak.assay = "ATAC",
  expression.assay = expression.form,
  expression.slot = "data"
)

saveRDS(bmmc_lp,'/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_count/bmmc_linkpeaks.Rds')