.libPaths("/home/zc354/R/4.2")
library(TFBSTools)
library(motifmatchr)
library(GenomicRanges)
library(dplyr)
library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
p=readRDS("/gpfs/gibbs/pi/zhao/yw599/GRN/Keran/cisBP_human_pfms_2021.rds")
i=intersect(names(p),rownames(pbmc.rna@assays$RNA@counts))
p=p[i]
saveRDS(p,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/TF_PWM.Rds")
print("p has been saved.")
#p为出现在RNA-seq中的TF
#用FIMO扫描开放区间上可能结合的TF，注意此时p某些TF可能没有结合在任何区间
#length(p)=1134
#length(p)=1119

library(universalmotif)
mot=convert_motifs(p)
write_meme(mot,file="/gpfs/gibbs/pi/zhao/zc354/GRN/output/meme-5.3.3/humanmotif1.meme",version=4)

atac_names_all = rownames(gene_peak_all)

atac_df <- data.frame(atac_name = c(atac_names_all))
atac_sep = separate(atac_df, col = "atac_name", into = c("chr","start","end"), sep = "-" )

atac_sep["start"]<-lapply(atac_sep["start"],as.numeric)
atac_sep["end"]<-lapply(atac_sep["end"],as.numeric)

chr = c(atac_sep[,"chr"])
start = c(atac_sep[,"start"])
end = c(atac_sep[,"end"])

ir<-IRanges(start = start, end = end)
peaks<- GRanges(seqnames = chr, ranges = ir)

saveRDS(ir,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/atac_iRanges.Rds")
saveRDS(peaks,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/atac_iRanges.Rds")
print("peaks have been saved. MatchMotif start running.")

motif_ix_scores <- matchMotifs(p, peaks, genome = "hg38",out = "scores")

saveRDS(motif_ix_scores,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matchmotif_score_result.Rds")
print("result of matchmotifs has been saved.")

motif_match<-motifMatches(motif_ix_scores)
rownames(motif_match) = atac_names_all
print("motif_match has been generated")

motif_match_score<-motifScores(motif_ix_scores)
rownames(motif_match_score) = atac_names_all
print("motif_match_score has been generated")

motif_match_binary<-motifCounts(motif_ix_scores)
rownames(motif_match_binary) = atac_names_all
print("motif_match_binary has been generated")

saveRDS(motif_match,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score.Rds")
saveRDS(motif_match_score,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score_matrix.Rds")
saveRDS(motif_match_binary,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_binary_matrix.Rds")
#################################################################################################

motif_ix_scores_6e_5 <- matchMotifs(p, peaks, genome = "hg38",out = "scores",p.cutoff = 6e-05)

saveRDS(motif_ix_scores_6e_5 ,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matchmotif_score_6e_5_result.Rds")
print("result of matchmotifs has been saved.")

motif_match_6e_5<-motifMatches(motif_ix_scores_6e_5)
rownames(motif_match_6e_5) = atac_names_all
print("motif_match has been generated")

motif_match_score_6e_5<-motifScores(motif_ix_scores_6e_5)
rownames(motif_match_score_6e_5) = atac_names_all
print("motif_match_score has been generated")

motif_match_binary_6e_5<-motifCounts(motif_ix_scores_6e_5)
rownames(motif_match_binary_6e_5) = atac_names_all
print("motif_match_binary has been generated")

saveRDS(motif_match_6e_5,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score_6e_5.Rds")
saveRDS(motif_match_score_6e_5,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score_matrix_6e_5.Rds")
saveRDS(motif_match_binary_6e_5,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_binary_matrix_6e_5.Rds")

w=which(colSums(motif_match_score_6e_5)[]==0)
#0
w=which(rowSums(motif_match_score_6e_5)[]==0)
#1458
#######################################################################################################
TF_peak_FIMO<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_TF_score_55353_776.Rds")

atac_names_FIMO = rownames(TF_peak_FIMO)
TF_names=colnames(TF_peak_FIMO)
p_FIMO = p@listData[[TF_names]]

atac_df <- data.frame(atac_name = c(atac_names_all))
atac_sep = separate(atac_df, col = "atac_name", into = c("chr","start","end"), sep = "-" )

atac_sep["start"]<-lapply(atac_sep["start"],as.numeric)
atac_sep["end"]<-lapply(atac_sep["end"],as.numeric)

chr = c(atac_sep[,"chr"])
start = c(atac_sep[,"start"])
end = c(atac_sep[,"end"])

ir<-IRanges(start = start, end = end)
peaks<- GRanges(seqnames = chr, ranges = ir)

saveRDS(ir,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/atac_iRanges.Rds")
saveRDS(peaks,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/atac_iRanges.Rds")
print("peaks have been saved. MatchMotif start running.")

motif_ix_scores <- matchMotifs(p, peaks, genome = "hg38",out = "scores")

saveRDS(motif_ix_scores,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/matchmotif_score_result.Rds")
print("result of matchmotifs has been saved.")

motif_match<-motifMatches(motif_ix_scores)
rownames(motif_match) = atac_names_all
print("motif_match has been generated")

motif_match_score<-motifScores(motif_ix_scores)
rownames(motif_match_score) = atac_names_all
print("motif_match_score has been generated")

motif_match_binary<-motifCounts(motif_ix_scores)
rownames(motif_match_binary) = atac_names_all
print("motif_match_binary has been generated")

saveRDS(motif_match,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score.Rds")
saveRDS(motif_match_score,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_score_matrix.Rds")
saveRDS(motif_match_binary,"/gpfs/gibbs/pi/zhao/zc354/GRN/data/motif_match_binary_matrix.Rds")
