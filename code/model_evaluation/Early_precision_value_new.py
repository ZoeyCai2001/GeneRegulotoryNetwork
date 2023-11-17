# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:48:23 2022

@author: Zhongyu_Cai
"""

import numpy as np
import pandas as pd
from itertools import product,permutations

true_path = '/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv'
pred_path = '/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename_unnormalization.csv'
directed = False
TFEdges = True
save_path = '/gpfs/gibbs/pi/zhao/zc354/GRN/output'
save_name = 'EPV'

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

if TFEdges:
    
    # Consider only edges going out of TFs
    # Get a list of all possible TF to gene interactions 
    uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
    possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))
    
    # Get a list of all possible interactions 
    possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))     
    possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)

else:
    #Consider all possible edges generated from trueEdgesDF
    possibleEdges = set(product(set(trueEdgesDF.Gene1),set(trueEdgesDF.Gene2)))
    possibleEdgesDF = pd.DataFrame(list(possibleEdges))
    possibleEdgesDF.columns=['Gene1','Gene2']
    #drop selfEdges
    possibleEdgesDF = possibleEdgesDF.loc[(possibleEdgesDF['Gene1'] != possibleEdgesDF['Gene2'])]
            
#Generate the possible edges dictionary
TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

#Generate the true edges dataframe
trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]

#Generate the possible edges dataframe         
predEdgesDF['Edges'] = predEdgesDF['Gene1'] + "|" + predEdgesDF['Gene2']
# limit the predicted edges to the genes that are in the ground truth
predEdgesDF['Edges'] = predEdgesDF[predEdgesDF['Edges'].isin(TrueEdgeDict)]
predEdgesDF.sort_values(by=['Edgeweight'],ascending=False)    
    
# we want to ensure that we do not include
# edges without any edge weight
# so check if the non-zero minimum is
# greater than the edge weight of the top-kth
# node, else use the non-zero minimum value.
predEdgesDF.Edgeweight = predEdgesDF.Edgeweight.round(6)
predEdgesDF.Edgeweight = predEdgesDF.Edgeweight.abs()

def EarlyPrec(k=10000):
    
    # Use num True edges or the number of
    # edges in the dataframe, which ever is lower
    numEdges = len(trueEdges)
    maxk = min(predEdgesDF.shape[0], numEdges)
    if k>maxk:
        print('k should be smaller than'+str(maxk))
    edgeweightTopk = predEdgesDF.iloc[k-1].Edgeweight

    nonZeroMin = np.nanmin(predEdgesDF.Edgeweight.replace(0, np.nan).values)
    bestVal = max(nonZeroMin, edgeweightTopk)
    dataset = {}
    newDF = predEdgesDF.loc[(predEdgesDF['Edgeweight'] >= bestVal)]
    dataset = set(newDF['Gene1'] + "|" + newDF['Gene2'])
   
    #Calculate Early prcision value
    Erec = {}    
    intersectionSet = dataset.intersection(trueEdges)
    Erec = len(intersectionSet)/k  

    return Erec
        
a = np.zeros(10)
for i in range(0,10,1):
      a[i] = (EarlyPrec(k = 100*(i+1)))
      
np.savetxt('/gpfs/gibbs/pi/zhao/zc354/GRN/output/model_evaluation/TF_filename_unnormalization_EPV.txt',a)