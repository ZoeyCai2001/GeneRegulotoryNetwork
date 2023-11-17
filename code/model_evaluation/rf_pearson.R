library(tidyverse)
model1="score"
model2="filter"
hvg_version="sct"
data_version="count"
#rf_path_1=paste('/home/zc354/network_2000/',model1,"_",hvg_version,"_",data_version,".csv",sep="")
#rf_path_2=paste('/home/zc354/network_2000/',model2,"_",hvg_version,"_",data_version,".csv",sep="")

rf_path_1=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group1.csv",sep="")
rf_path_2=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group2.csv",sep="")

pred_1<-read.table(file =rf_path_1,header = T,sep=",")
pred_2<-read.table(file =rf_path_2,header = T,sep=",")
pred_1=aggregate(Edgeweight~Gene1+Gene2 , data=pred_1 , sum) 
pred_2=aggregate(Edgeweight~Gene1+Gene2 , data=pred_2 , sum) 
gene_names_1 = pred_1['Gene2']
gene_names_1 <- gene_names_1 %>% distinct(Gene2, .keep_all = T)
gene_names_2 = pred_2['Gene2']
gene_names_2 <- gene_names_2 %>% distinct(Gene2, .keep_all = T)
gene_names = intersect(gene_names_1[,1],gene_names_2[,1])

r=vector()
p=vector()
for(i in 1:length(gene_names)){
  gene = gene_names[i]
  TF_1 = pred_1[which(pred_1['Gene2']==gene),]
  TF_2 = pred_2[which(pred_2['Gene2']==gene),]
  TF = rbind(TF_1,TF_2)
  TF <- TF %>% distinct(Gene1,.keep_all = T)
  TF$Edgeweight_1 <- TF$Edgeweight
  for(j in 1:dim(TF)[1]){
    tf = TF['Gene1'][j,]
    w_1 = which(TF_1['Gene1']==tf)
    w_2 = which(TF_2['Gene1']==tf)
    if(length(w_1)==0){
      TF[j,'Edgeweight']=0
    }else{
      TF[j,'Edgeweight']=TF_1[w_1,'Edgeweight']
    }
    
    if(length(w_2)==0){
      TF[j,'Edgeweight_1']=0
    }else{
      TF[j,'Edgeweight_1']=TF_2[w_2,'Edgeweight']
    }
  }
  r[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'],method = 'pearson')
  p[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'],method = 'spearman')
  print(i)
}
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_pearson.csv",sep=""),row.names = F,sep=",")
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_spearman.csv",sep=""),row.names = F,sep=",")

library(tidyverse)
model1="multiply"
model2="filter"
hvg_version="sct"
data_version="count"
#rf_path_1=paste('/home/zc354/network_2000/',model1,"_",hvg_version,"_",data_version,".csv",sep="")
#rf_path_2=paste('/home/zc354/network_2000/',model2,"_",hvg_version,"_",data_version,".csv",sep="")

rf_path_1=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group1.csv",sep="")
rf_path_2=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group2.csv",sep="")

pred_1<-read.table(file =rf_path_1,header = T,sep=",")
pred_2<-read.table(file =rf_path_2,header = T,sep=",")
pred_1=aggregate(Edgeweight~Gene1+Gene2 , data=pred_1 , sum) 
pred_2=aggregate(Edgeweight~Gene1+Gene2 , data=pred_2 , sum) 
gene_names_1 = pred_1['Gene2']
gene_names_1 <- gene_names_1 %>% distinct(Gene2, .keep_all = T)
gene_names_2 = pred_2['Gene2']
gene_names_2 <- gene_names_2 %>% distinct(Gene2, .keep_all = T)
gene_names = intersect(gene_names_1[,1],gene_names_2[,1])

r=vector()
p=vector()
for(i in 1:length(gene_names)){
  gene = gene_names[i]
  TF_1 = pred_1[which(pred_1['Gene2']==gene),]
  TF_2 = pred_2[which(pred_2['Gene2']==gene),]
  TF = rbind(TF_1,TF_2)
  TF <- TF %>% distinct(Gene1,.keep_all = T)
  TF$Edgeweight_1 <- TF$Edgeweight
  for(j in 1:dim(TF)[1]){
    tf = TF['Gene1'][j,]
    w_1 = which(TF_1['Gene1']==tf)
    w_2 = which(TF_2['Gene1']==tf)
    if(length(w_1)==0){
      TF[j,'Edgeweight']=0
    }else{
      TF[j,'Edgeweight']=TF_1[w_1,'Edgeweight']
    }
    
    if(length(w_2)==0){
      TF[j,'Edgeweight_1']=0
    }else{
      TF[j,'Edgeweight_1']=TF_2[w_2,'Edgeweight']
    }
  }
  r[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'])
  p[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'],method = 'spearman')
  print(i)
}
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_pearson.csv",sep=""),row.names = F,sep=",")
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_spearman.csv",sep=""),row.names = F,sep=",")

model1="score"
model2="filter"
hvg_version="sct"
data_version="count"
#rf_path_1=paste('/home/zc354/network_2000/',model1,"_",hvg_version,"_",data_version,".csv",sep="")
#rf_path_2=paste('/home/zc354/network_2000/',model2,"_",hvg_version,"_",data_version,".csv",sep="")

rf_path_1=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group1.csv",sep="")
rf_path_2=paste("/home/zc354/GRN/output/predicted_network/",model1,"_",hvg_version,"_",data_version,"_bmmc_linkpeaks_motifmatchr_group2.csv",sep="")

pred_1<-read.table(file =rf_path_1,header = T,sep=",")
pred_2<-read.table(file =rf_path_2,header = T,sep=",")
pred_1=aggregate(Edgeweight~Gene1+Gene2 , data=pred_1 , sum) 
pred_2=aggregate(Edgeweight~Gene1+Gene2 , data=pred_2 , sum) 
gene_names_1 = pred_1['Gene2']
gene_names_1 <- gene_names_1 %>% distinct(Gene2, .keep_all = T)
gene_names_2 = pred_2['Gene2']
gene_names_2 <- gene_names_2 %>% distinct(Gene2, .keep_all = T)
gene_names = intersect(gene_names_1[,1],gene_names_2[,1])

r=vector()
p=vector()
for(i in 1:length(gene_names)){
  gene = gene_names[i]
  TF_1 = pred_1[which(pred_1['Gene2']==gene),]
  TF_2 = pred_2[which(pred_2['Gene2']==gene),]
  TF = rbind(TF_1,TF_2)
  TF <- TF %>% distinct(Gene1,.keep_all = T)
  TF$Edgeweight_1 <- TF$Edgeweight
  for(j in 1:dim(TF)[1]){
    tf = TF['Gene1'][j,]
    w_1 = which(TF_1['Gene1']==tf)
    w_2 = which(TF_2['Gene1']==tf)
    if(length(w_1)==0){
      TF[j,'Edgeweight']=0
    }else{
      TF[j,'Edgeweight']=TF_1[w_1,'Edgeweight']
    }
    
    if(length(w_2)==0){
      TF[j,'Edgeweight_1']=0
    }else{
      TF[j,'Edgeweight_1']=TF_2[w_2,'Edgeweight']
    }
  }
  r[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'],method= 'pearson')
  p[i]=cor(TF['Edgeweight'],TF['Edgeweight_1'],method = 'spearman')
  print(i)
}
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_pearson.csv",sep=""),row.names = F,sep=",")
write.table(r,file=paste("/gpfs/gibbs/pi/zhao/zc354/GRN/output/",model1,"_donor_spearman.csv",sep=""),row.names = F,sep=",")

