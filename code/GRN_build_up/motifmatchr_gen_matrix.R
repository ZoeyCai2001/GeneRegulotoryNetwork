hvg1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_sctransform_55353_2000.Rds")
gene_names<-colnames(hvg1)

TF_peak_score = readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score_matrix.Rds")
TF_peak_binary = readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_binary_matrix.Rds")

w=which(colSums(TF_peak_score)[]==0)
#length(w)==0
w=which(rowSums(TF_peak_score)[]==0)
#length(w)==1527
TF_peak_score = TF_peak_score[-w,]

w=which(colSums(TF_peak_binary)[]==0)
#length(w)==0
w=which(rowSums(TF_peak_binary)[]==0)
#length(w)==1527
TF_peak_score = TF_peak_binary[-w,]

TF_names = colnames(TF_peak_binary)
atac_names = intersect(intersect(rownames(TF_peak_score),rownames(gene_peak_all)),rownames(pbmc@assays[["ATAC"]]@counts))

gene_peak = gene_peak_all[atac_names,gene_names]
TF_peak_score = TF_peak_score[atac_names,]
TF_peak_binary = TF_peak_binary[atac_names,]
atac_data = pbmc.atac@assays[["ATAC"]]@counts[atac_names,]
gene_data = pbmc.rna@assays[["RNA"]]@counts[union(gene_names,TF_names),]

saveRDS(gene_peak,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_normalizedata_106120_2000.Rds")
saveRDS(TF_peak_score,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_TF_score_106120_1119.Rds")
saveRDS(TF_peak_binary,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_TF_106120_1119.Rds")
saveRDS(atac_data,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/atac_count_106120.Rds")
saveRDS(gene_data,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/TF_hvg_sctransform_count_3001.Rds")
