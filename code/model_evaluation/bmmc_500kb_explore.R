model1<-"non_filter"  #non_filter OR filter OR multiply
model2<-"filter"
hvg_version<-"sct"  #sct OR normdata
data_version<-"count"  #norm OR count
dir = "/home/zc354/GRN/output/network_2000/"
t_true <-read.table('/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv',header = T,sep=",")
t_result_1<-read.table(file =paste(writedir,model1,"_",hvg_version,"_",data_version,"_bmmc_500kb_motifmatchr_sampling",".csv",sep = ""),header = T,sep=",")
t_result_2<-read.table(file =paste(writedir,model2,"_",hvg_version,"_",data_version,"_bmmc_500kb_motifmatchr_sampling",".csv",sep = ""),header = T,sep=",")

t_true[,'Edge']=paste(t_true[,'Gene1'],'|',t_true[,'Gene2'],sep='')
t_result_1[,'Edge']=paste(t_result_1[,'Gene1'],'|',t_result_1[,'Gene2'],sep='')
t_result_2[,'Edge']=paste(t_result_2[,'Gene1'],'|',t_result_2[,'Gene2'],sep='')
gene_diff = setdiff(t_result_1[,'Gene2'],t_result_2[,'Gene2'])
w = which(t_result_1[,'Gene2'] %in% gene_diff)
t_result_1 = t_result_1[-w,]

true_tg = t_true[,'Gene2']
true_TF = t_true[,'Gene1']

t_result_1 <- t_result_1[order(-t_result_1$Edgeweight),]
t_result_2 <- t_result_2[order(-t_result_2$Edgeweight),]
trueedge = t_true[,'Edge']
prededge_1 = t_result_1[,'Edge']
prededge_2 = t_result_2[,'Edge']


pred_filter_1 = t_result_1[intersect(which(t_result_1[,'Gene2'] %in% true_tg),which(t_result_1[,'Gene1'] %in% true_TF)),]
pred_filter_2 = t_result_2[intersect(which(t_result_2[,'Gene2'] %in% true_tg),which(t_result_2[,'Gene1'] %in% true_TF)),]
gene_diff_filter = setdiff(pred_filter_1[,'Gene2'],pred_filter_2[,'Gene2'])

prededge_filter_1 = pred_filter_1[,'Edge']
prededge_filter_2 = pred_filter_2[,'Edge']

hvg1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_hvg_sctransform_116468_1608.Rds")
M1<-readRDS("/gpfs/gibbs/pi/zhao/zc354/GRN/data/matrix_openchromt_TF_116468_649.Rds")
gene_TF<-t(hvg1)%*%M1
gene_TF_pair = c()
for(j in 1:(length(gene_TF@p)-1)){
  num = gene_TF@p[j+1]-gene_TF@p[j]
  n_1 = gene_TF@p[j]+1
  n_2 = gene_TF@p[j+1]
  TF = (gene_TF@Dimnames[[2]])[j]
  gene = (gene_TF@Dimnames[[1]])[gene_TF@i[n_1:n_2]+1]
  gene_TF_pair[n_1:n_2] = paste(TF,"|",gene,sep = '')
  print(j)
}

num1=c()
num2=c()
num3=c()
num1_filtered = c()
num2_filtered = c()
num3_filtered = c()

for(i in 1:1000){
  pred_1 = intersect(prededge_1[1:(100*i)],trueedge)
  pred_2 = intersect(prededge_2[1:(100*i)],trueedge)
  num1[i] = length(pred_1)-length(pred_2)
  num2[i] = length(setdiff(setdiff(pred_1,pred_2),gene_TF_pair))
  w = which((pred_1 %in% pred_2) ==FALSE)
  num3[i] = length(which(t_result_1[w,'Gene2'] %in% gene_diff))
  
  pred_filtered_1 = intersect(prededge_filter_1[1:(100*i)],trueedge)
  pred_filtered_2 = intersect(prededge_filter_2[1:(100*i)],trueedge)
  num1_filtered[i] = length( pred_filtered_1) - length(pred_filtered_2)
  num2_filtered[i] = length(setdiff(setdiff(pred_filtered_1,pred_filtered_2),gene_TF_pair))
  w = which((pred_filtered_1[] %in% pred_filtered_2) ==FALSE)
  num3_filtered[i] = length(which(pred_filter_1[w,'Gene2'][] %in% gene_diff))
  print(i*100)
}

t = data.frame('diff' = num1, 'not_in_500kb' = num2, 'not_in_tg' = num3, 'diff_filter' = num1_filtered, 'not_in_500kb_filtered' = num2_filtered, 'not_in_tg_filtered' = num3_filtered)

write.table(t,file = '/gpfs/gibbs/pi/zhao/zc354/GRN/output/bmmc_result_explore_new.csv',row.names = F,sep=",")

