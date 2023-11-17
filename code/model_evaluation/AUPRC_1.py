# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 09:37:01 2022

@author: Zhongyu_Cai
"""
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from itertools import product, permutations, combinations, combinations_with_replacement
from sklearn.metrics import precision_recall_curve, roc_curve, auc
import matplotlib.pyplot as plt

true_path='/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv'
pred_path='/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename_unnormalization.csv'
directed=False
save_path='/gpfs/gibbs/pi/zhao/zc354/GRN/output'
save_name='AUPRC_TF_unnorm_indir'


#Read true netword and predicted network
trueEdgesDF = pd.read_csv(true_path,sep = ',', header = 0, index_col = None)
predEdgesDF = pd.read_csv(pred_path, sep = ',', header =  0, index_col = None)
    
#Unification changes lowercase letters to uppercase letters
trueEdgesDF['Gene1'] = trueEdgesDF['Gene1'].str.upper()
trueEdgesDF['Gene2'] = trueEdgesDF['Gene2'].str.upper()
predEdgesDF['Gene1'] = predEdgesDF['Gene1'].str.upper()
predEdgesDF['Gene2'] = predEdgesDF['Gene2'].str.upper()
    
#Drop index in predEdgesDF that Edgeweight=0
predEdgesDF.drop(predEdgesDF[predEdgesDF['Edgeweight'] == 0].index,inplace=True)
    
#Drop selfEdges
trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
predEdgesDF = predEdgesDF.loc[(predEdgesDF['Gene1'] != predEdgesDF['Gene2'])]
    
#Drop duplicates in trueEdgesDF and predEdgesDF
trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
predEdgesDF.drop_duplicates(keep = 'first', subset=['Gene1','Gene2'], inplace=True)
    
#Generate all possibleEdges using Cartesian product of trueEdgesDF['Gene1'] & trueEdgesDF['Gene2']
possibleEdges = set(product(set(trueEdgesDF['Gene1']),set(trueEdgesDF['Gene2'])))
possibleEdgesDF = pd.DataFrame(list(possibleEdges))
possibleEdgesDF.columns=['Gene1','Gene2']
    
#drop selfEdges
possibleEdgesDF = possibleEdgesDF.loc[(possibleEdgesDF['Gene1'] != possibleEdgesDF['Gene2'])]
    
if directed:
   #Represent an edge as a string: Gene1|Gene2 and stored in DataFrame
   TrueEdgesDF = pd.DataFrame(columns=['Edge'])
   TrueEdgesDF['Edge'] = trueEdgesDF['Gene1'].map(str) +"|"+ trueEdgesDF['Gene2'].map(str)
   PredEdgesDF = pd.DataFrame(columns=['Edge','Weight'])
   PredEdgesDF['Edge'] = predEdgesDF['Gene1'].map(str) +"|"+ predEdgesDF['Gene2'].map(str)
   PredEdgesDF['Weight']=predEdgesDF['Edgeweight']
   PossibleEdgesDF = pd.DataFrame(columns=['Edge','Pred','True'])
   PossibleEdgesDF['Edge'] = possibleEdgesDF['Gene1'].map(str) +"|"+ possibleEdgesDF['Gene2'].map(str)  
        
   #Generate the set of true/pred/possible edges
   TrueEdgesset=set(TrueEdgesDF['Edge'])
   PredEdgesset=set(PredEdgesDF['Edge'])
   PossibleEdgesset=set(PossibleEdgesDF['Edge'])
        
   #Generate four possible conditions using set operation
   #inbothset represents edges in both true and pred network
   #onlyintrue set represents edges in true network but not in pred
   #onlyinpred set represents edges in pred but not in true network
   #inneitherset represent edges in neither of true and pred network but in possible network
   inbothset = (TrueEdgesset & PredEdgesset) & PossibleEdgesset
   onlyintrue = (TrueEdgesset - PredEdgesset) & PossibleEdgesset
   onlyinpred = (PredEdgesset - TrueEdgesset) & PossibleEdgesset
   inneitherset = PossibleEdgesset - (TrueEdgesset | PredEdgesset)
        
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
    
else:
    #Generate outDF copied from possibleEdgesDF
    outDF = pd.DataFrame(columns=['Edge','Pred','True'])
    outDF['Edge'] = possibleEdgesDF['Gene1'].map(str) +"|"+ possibleEdgesDF['Gene2'].map(str)
        
    #Generate a funtion of the dictionary tree for trueEdgesDF to facilitate searching
    def gen_dict():
         gene2_dict={}
         for gene2,df in trueEdgesDF.groupby('Gene2'):
             gene2_dict[gene2]=set(df['Gene1'])
         return gene2_dict
        
    gene2_dict=gen_dict()
        
    #Generate a function of searching for trueEdgesDF
    def main_fun(row):
        global gene2_dict
        t1,t2=row['Gene1'],row['Gene2']
        # First search for (gene1==t1 & gene2==t2)
        if gene2_dict.get(t2):
            temp=gene2_dict[t2]
            if t1 in temp:
                return 1
          # Then search for (gene1==t2 & gene2==t1)
        if gene2_dict.get(t1):
            temp=gene2_dict[t1]
            if t2 in temp:
                return 1
        return 0
        
    import swifter
    swifter.register_modin()
        
    #Using apply function to speed up assigning values
    Value=possibleEdgesDF.swifter.apply(main_fun,axis=1)
        
    #Assigning values for TrueEdges
    outDF['True'] = Value
        
    #Generate a funtion of the dictionary tree for predEdgesDF to facilitate searching
    #gene2_dict_2 is a dictionary including set('Gene2')
    #gene2_dict_2['a gene name in Gene2'] is a dictionary of gene1 which means [Gene1,Gene2] is in PredEdgesDF
    def gen_dict_2():
         gene2_dict_2={}
         for gene2,df in predEdgesDF.groupby('Gene2'):
             gene2_dict_2[gene2]=dict(zip(list(df['Gene1']),list(df['Edgeweight'])))
         return gene2_dict_2
        
    gene2_dict_2 = gen_dict_2()
        
        #Generate a function of searching for trueEdgesDF
    def main_fun_2(row):
        global gene2_dict_2
        #res to record EdgeWeights in predEdgesDF
        res=[]
        t1 ,t2=row['Gene1'],row['Gene2']
        # First search for (gene1==t1 & gene2==t2)
        if gene2_dict_2.get(t2):
            temp=gene2_dict_2[t2]
            if t1 in temp:
              res.append(temp[t1])
            
        # Then search for (gene1==t2 & gene2==t1)
        if gene2_dict_2.get(t1):
            temp=gene2_dict_2[t1]
            if t2 in temp:
                res.append(temp[t2])
        if res:
            return max(res)
        else:
            return 0
        
    #Using apply function to speed up assigning values
    swifter.register_modin()
    Weight = possibleEdgesDF.swifter.apply(main_fun_2,axis=1)
        
    #Assigning values for PredEdges
    outDF['Pred'] = Weight
        
    #Generate outDF to use precision_recall_curve
    outDF.set_index('Edge',inplace = True)
    
#Calculate precision_recall curve
prec, recall, thresholds = precision_recall_curve(y_true=outDF['True'],probas_pred=outDF['Pred'], pos_label=1)
    
#Calculate Area Under Precision_Recall Curve
auprc = auc(recall, prec)
    
# Make PR curves
plt.xlim(0,1)    
plt.ylim(0,1)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.plot(prec, recall, "blue")  
plt.savefig(save_path+save_name+'.pdf')
plt.savefig(save_path+save_name+'.png')
plt.clf()
    
print(save_name+':'+str(auprc))