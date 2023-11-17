library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)
# load the RNA and ATAC data
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
  sep = c("-", "-"),
  fragments = fragpath,
  annotation = annotation
)


# create ATAC assay and add it to the object
pbmc.atac[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c("-", "-"),
  fragments = fragpath,
  annotation = annotation
)
