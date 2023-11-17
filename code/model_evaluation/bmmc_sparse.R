pbmc_rna<-pbmc@assays[["RNA"]]@counts
#36601 11909
pbmc_atac<-pbmc@assays[["ATAC"]]@counts
#108344  11909
bmmc_rna<-bmmc@assays[["RNA"]]@counts
#13431 69249
bmmc_atac<-bmmc@assays[["ATAC"]]@counts
#116468  69249
pbmc_rna_1<-pbmc_rna
pbmc_rna_1@x<-rep(1, times=length(pbmc_rna_1@x))
pbmc_rna_0<-colSums(pbmc_rna_1)
pbmc_rna_sum<-colSums(pbmc_rna)

pbmc_atac_1<-pbmc_atac
pbmc_atac_1@x<-rep(1, times=length(pbmc_atac_1@x))
pbmc_atac_0<-colSums(pbmc_atac_1)
pbmc_atac_sum<-colSums(pbmc_atac)

bmmc_rna_1<-bmmc_rna
bmmc_rna_1@x<-rep(1, times=length(bmmc_rna_1@x))
bmmc_rna_0<-colSums(bmmc_rna_1)
bmmc_rna_sum<-colSums(bmmc_rna)

bmmc_atac_1<-bmmc_atac
bmmc_atac_1@x<-rep(1, times=length(bmmc_atac_1@x))
bmmc_atac_0<-colSums(bmmc_atac_1)
bmmc_atac_sum<-colSums(bmmc_atac)

count_pbmc_rna = data.frame(pbmc_rna_0,pbmc_rna_sum)
count_pbmc_atac = data.frame(pbmc_atac_0,pbmc_atac_sum)
count_bmmc_rna = data.frame(bmmc_rna_0,bmmc_rna_sum)
count_bmmc_atac = data.frame(bmmc_atac_0,bmmc_atac_sum)

write.table(count_pbmc_rna,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_pbmc_rna.csv",row.names = F,sep=",")
write.table(count_pbmc_atac,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_pbmc_atac.csv",row.names = F,sep=",")
write.table(count_bmmc_rna,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_bmmc_rna.csv",row.names = F,sep=",")
write.table(count_bmmc_atac,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_bmmc_atac.csv",row.names = F,sep=",")

hvg1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_normalizedata_106120_2000.Rds")
hvg_pbmc<-colnames(hvg1)
pbmc_rna<-pbmc@assays[["RNA"]]@counts[hvg_pbmc,]
pbmc_rna_1<-pbmc_rna
pbmc_rna_1@x<-rep(1, times=length(pbmc_rna_1@x))
pbmc_rna_0<-colSums(pbmc_rna_1)
pbmc_rna_sum<-colSums(pbmc_rna)

w=which(rowSums(hvg1)[]!=0)
atac_pbmc<-rownames(hvg1)[w]
pbmc_atac<-pbmc@assays[["ATAC"]]@counts[atac_pbmc,]
pbmc_atac_1<-pbmc_atac
pbmc_atac_1@x<-rep(1, times=length(pbmc_atac_1@x))
pbmc_atac_0<-colSums(pbmc_atac_1)
pbmc_atac_sum<-colSums(pbmc_atac)

hvg1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_sctransform_116468_1608.Rds")
hvg_bmmc<-colnames(hvg1)
bmmc_rna<-bmmc@assays[["RNA"]]@counts[hvg_bmmc,]
bmmc_rna_1<-bmmc_rna
bmmc_rna_1@x<-rep(1, times=length(bmmc_rna_1@x))
bmmc_rna_0<-colSums(bmmc_rna_1)
bmmc_rna_sum<-colSums(bmmc_rna)

w=which(rowSums(hvg1)[]!=0)
atac_bmmc<-rownames(hvg1)[w]
bmmc_atac<-bmmc@assays[["ATAC"]]@counts[atac_bmmc,]
bmmc_atac_1<-bmmc_atac
bmmc_atac_1@x<-rep(1, times=length(bmmc_atac_1@x))
bmmc_atac_0<-colSums(bmmc_atac_1)
bmmc_atac_sum<-colSums(bmmc_atac)

count_pbmc_rna_hvg = data.frame(pbmc_rna_0,pbmc_rna_sum)
count_bmmc_rna_hvg = data.frame(bmmc_rna_0,bmmc_rna_sum)
count_pbmc_atac_hvg = data.frame(pbmc_atac_0,pbmc_atac_sum)
count_bmmc_atac_hvg = data.frame(bmmc_atac_0,bmmc_atac_sum)


write.table(count_pbmc_rna_hvg,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_pbmc_rna_hvg.csv",row.names = F,sep=",")
write.table(count_bmmc_rna_hvg,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_bmmc_rna_hvg.csv",row.names = F,sep=",")
write.table(count_pbmc_atac_hvg,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_pbmc_atac_hvg.csv",row.names = F,sep=",")
write.table(count_bmmc_atac_hvg,file = "/gpfs/gibbs/pi/zhao/zc354/GRN/output/data_quality/count_bmmc_atac_hvg.csv",row.names = F,sep=",")
