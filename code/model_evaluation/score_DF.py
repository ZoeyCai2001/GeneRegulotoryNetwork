import numpy as np
import pandas as pd
from itertools import product, permutations, combinations, combinations_with_replacement
from sklearn.metrics import precision_recall_curve, roc_curve, auc

def outDF(model = "non_filter",hvg_version = "sct",data_version = "count",data = 'pbmc',gene_peak = '500kb', peak_TF = 'motifmatchr',
          directory = "/home/zc354/GRN/output/predicted_network/", true_path ='/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv'):
    
    pred_path=directory+model+"_"+hvg_version+"_"+data_version+"_"+data+"_"+gene_peak+"_"+peak_TF+'.csv'
    if(peak_TF == 'FIMO'):
        pred_path_1=directory+"filter_"+hvg_version+"_"+data_version+"_"+data+"_"+gene_peak+"_"+peak_TF+'.csv'
    else:
        pred_path_1=directory+"filter_sct_"+data_version+"_"+data+"_"+gene_peak+"_"+peak_TF+'.csv'
        
    
    trueEdgesDF = pd.read_csv(true_path,sep = ',', header = 0, index_col = None)
    predEdgesDF = pd.read_csv(pred_path, sep = ',', header =  0, index_col = None)
    predEdgesDF_1 = pd.read_csv(pred_path_1, sep = ',', header =  0, index_col = None)
    
    print('length of pred-network:'+str(len(predEdgesDF['Gene1'])))
    print('length of pred-network TF:'+str(len(set(predEdgesDF['Gene1']))))
    print('length of pred-network gene:'+str(len(set(predEdgesDF['Gene2']))))
    
    predEdgesDF.drop(predEdgesDF[predEdgesDF['Edgeweight'] == 0].index,inplace=True)

    trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
    predEdgesDF = predEdgesDF.loc[(predEdgesDF['Gene1'] != predEdgesDF['Gene2'])]

    trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
    predEdgesDF.drop_duplicates(keep = 'first',subset=['Gene1','Gene2'],inplace=True)

    print('length of true-network TF:'+str(len(set(trueEdgesDF['Gene1']))))
    print('length of true-network gene:'+str(len(set(trueEdgesDF['Gene2']))))
    
    print('length of pred-network TF_filtered:'+str(len(set(predEdgesDF['Gene1']))))
    print('length of pred-network gene_filtered:'+str(len(set(predEdgesDF['Gene2']))))
    
    true_TF = set(trueEdgesDF['Gene1'])
    true_tg = set(trueEdgesDF['Gene2'])

    pred_TF = set(predEdgesDF_1['Gene1'])
    pred_tg = set(predEdgesDF_1['Gene2'])
    
    print('number of TF in background:'+str(len(pred_TF)))
    print('number of GENE in background:'+str(len(pred_tg)))

    TFset = (true_TF & pred_TF) 
    tgset = (true_tg & pred_tg)
    
    print('length of TF_intersect:'+str(len(TFset)))
    print('length of tg_intersect:'+str(len(tgset)))

    trueEdgesDF = trueEdgesDF[(trueEdgesDF['Gene1'].isin(TFset))&(trueEdgesDF['Gene2'].isin(tgset))]
    predEdgesDF = predEdgesDF[(predEdgesDF['Gene1'].isin(TFset))&(predEdgesDF['Gene2'].isin(tgset))]
    print('length of TF in truenetwork:'+str(len(set(trueEdgesDF['Gene1'])))+'_'+'length of tg in truenetwork:'+str(len(set(trueEdgesDF['Gene2']))))
    print('length of TF in prednetwork:'+str(len(set(predEdgesDF['Gene1'])))+'_'+'length of tg in prednetwork:'+str(len(set(predEdgesDF['Gene2']))))
    
    #possibleEdges = set(product(set(trueEdgesDF['Gene1']),set(trueEdgesDF['Gene2'])))
    possibleEdges = set(product(TFset,tgset))
    possibleEdgesDF = pd.DataFrame(list(possibleEdges))
    possibleEdgesDF.columns=['Gene1','Gene2']
    
    TrueEdgesDF = pd.DataFrame(columns=['Edge'])
    TrueEdgesDF['Edge'] = trueEdgesDF['Gene1'].map(str) +"|"+ trueEdgesDF['Gene2'].map(str)
    PredEdgesDF = pd.DataFrame(columns=['Edge','Weight'])
    PredEdgesDF['Edge'] = predEdgesDF['Gene1'].map(str) +"|"+ predEdgesDF['Gene2'].map(str)
    PredEdgesDF['Weight'] = predEdgesDF['Edgeweight']
    
    PossibleEdgesDF = pd.DataFrame(columns=['Edge','Pred','True'])
    PossibleEdgesDF['Edge'] = possibleEdgesDF['Gene1'].map(str) +"|"+ possibleEdgesDF['Gene2'].map(str)
    
    #Generate the set of true/pred/possible edges
    TrueEdgesset=set(TrueEdgesDF['Edge'])
    PredEdgesset=set(PredEdgesDF['Edge'])
    PossibleEdgesset=set(PossibleEdgesDF['Edge'])
    
    print("length of trueEdges:"+str(len(TrueEdgesset)))
    print("length of predEdges:"+str(len(PredEdgesset)))
    print("length of backgroudEdges:"+str(len(PossibleEdgesset)))
    #Generate four possible conditions using set operation
    #inbothset represents edges in both true and pred network
    #onlyintrue set represents edges in true network but not in pred
    #onlyinpred set represents edges in pred but not in true network
    #inneitherset represent edges in neither of true and pred network but in possible network
    inbothset = (TrueEdgesset & PredEdgesset) & PossibleEdgesset
    onlyintrue = (TrueEdgesset - PredEdgesset) & PossibleEdgesset
    onlyinpred = (PredEdgesset - TrueEdgesset) & PossibleEdgesset
    inneitherset = PossibleEdgesset - (TrueEdgesset | PredEdgesset)
    print('length of inbothset:'+str(len(inbothset)))
    print('length of onlyintrue:'+str(len(onlyintrue)))
    print('length of onlyinpred:'+str(len(onlyinpred)))
    print('length of inneitherset:'+str(len(inneitherset)))
    #Set Edges as index in three DataFrames to facilitate indexing
    PredEdgesDF.set_index('Edge',inplace=True)
    TrueEdgesDF.set_index('Edge',inplace=True)
    PossibleEdgesDF.set_index('Edge',inplace=True)
    
    #Assign values to different types of edges
    inbothDF = PredEdgesDF.loc[list(inbothset),:]
    inbothDF['True'] = 1
    inbothDF.columns = ['Pred','True']

    onlyintrueDF = TrueEdgesDF.loc[list(onlyintrue),:]
    onlyintrueDF['Pred'] = 0
    onlyintrueDF['True'] = 1

    onlyinpredDF = PredEdgesDF.loc[list(onlyinpred),:]
    onlyinpredDF['True'] = 0
    onlyinpredDF.columns = ['Pred','True']

    inneitherDF = PossibleEdgesDF.loc[list(inneitherset),:]
    inneitherDF['Pred'] = 0
    inneitherDF['True'] = 0

    outDF = pd.concat([inbothDF,onlyintrueDF,onlyinpredDF,inneitherDF],axis=0)
    return PredEdgesDF,outDF

data_version = "count"
data = 'bmmc'
gene_peak = '500kb'
peak_TF = 'motifmatchr'

PredEdgesDF_non_filter,outDF_non_filter = outDF(model = "non_filter",data_version = data_version,data = data,gene_peak = gene_peak, peak_TF =  peak_TF)

PredEdgesDF_filter,outDF_filter = outDF(model = "filter",data_version = data_version,data = data,gene_peak = gene_peak, peak_TF =  peak_TF)

PredEdgesDF_multiply,outDF_multiply = outDF(model = "multiply",data_version = data_version,data = data,gene_peak = gene_peak, peak_TF =  peak_TF)

PredEdgesDF_score,outDF_score = outDF(model = "score",data_version = data_version,data = data,gene_peak = gene_peak, peak_TF =  peak_TF)

prec_non_filter, recall_non_filter, thresholds_non_filter = precision_recall_curve(y_true=outDF_non_filter ['True'],probas_pred=outDF_non_filter ['Pred'], pos_label=1)
prec_filter, recall_filter, thresholds_filter = precision_recall_curve(y_true=outDF_filter ['True'],probas_pred=outDF_filter ['Pred'], pos_label=1)
prec_multiply, recall_multiply, thresholds_multiply = precision_recall_curve(y_true=outDF_multiply['True'],probas_pred=outDF_multiply['Pred'], pos_label=1)
prec_score, recall_score, thresholds_score = precision_recall_curve(y_true=outDF_score['True'],probas_pred=outDF_score['Pred'], pos_label=1)

print('non_filter_auprc:'+str(auc(recall_non_filter,prec_non_filter)))
print('filter_auprc:'+str(auc(recall_filter,prec_filter)))
print('multiply_auprc:'+str(auc(recall_multiply,prec_multiply)))
print('score_auprc:'+str(auc(recall_score,prec_score)))

new_list = np.append(thresholds_score, 1)
DF = pd.DataFrame(prec_score)
DF['recall']= recall_score
DF['threshold']= new_list 
DF['positive']= 0
DF['true positive']= 0
DF.columns = ['precision','recall','thresholds','positive','true positive']

outDF_score  = outDF_score.sort_values(by="Pred",ascending=False)
trueDF = outDF_score.loc[outDF_score['True'] > 0]
true = set(trueDF.index.values)
print('score_begin')
for i in range(len(prec_score)):
    newDF = outDF_score.loc[outDF_score['Pred'] >= DF.loc[i,'thresholds']]
    DF.loc[i,'positive'] = newDF.shape[0]
    positives = set(newDF.index.values)
    DF.loc[i,'true positive'] = len(positives.intersection(true))
    if(i % 10 ==0):
        print(i)
        
DF.to_csv("/gpfs/gibbs/pi/zhao/zc354/GRN/output/topk_score_bmmc_500kb.csv",index=False,sep=',')
print('score end and saved as /gpfs/gibbs/pi/zhao/zc354/GRN/output/topk_score_bmmc_500kb.csv')
