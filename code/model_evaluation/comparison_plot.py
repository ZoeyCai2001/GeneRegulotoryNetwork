# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 17:13:48 2022

@author: Zhongyu_Cai
"""

import numpy as np
import pandas as pd
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
import matplotlib.pyplot as plt

def gen_df(path='/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename.csv'):
    #Read true netword and predicted network
    EdgesDF = pd.read_csv(path,sep = ',', header = 0, index_col = None)

    #Unification changes lowercase letters to uppercase letters
    EdgesDF['Gene1'] = EdgesDF['Gene1'].str.upper()

    #Drop index in predEdgesDF that Edgeweight=0
    if 'Edgeweight' in EdgesDF.columns.values.tolist() :
        EdgesDF.drop(EdgesDF[EdgesDF['Edgeweight'] == 0].index,inplace=True)
        EdgesDF = EdgesDF.sort_values(by=['Edgeweight'],ascending=False)

    #Drop selfEdges
    EdgesDF = EdgesDF.loc[(EdgesDF['Gene1'] != EdgesDF['Gene2'])]

    #Drop duplicates in trueEdgesDF and predEdgesDF
    EdgesDF.drop_duplicates(keep = 'first', subset=['Gene1','Gene2'], inplace=True)
    
    #Generate a column to represent Edges
    EdgesDF['Edges'] = EdgesDF['Gene1'] + "|" + EdgesDF['Gene2']
    
    return EdgesDF

#Define a function to calculate number of top k weight edges in predEdgesDF 
#that is in true network
def num_topk(predEdgesDF = predEdgesDF_1, k = 10000):
    maxk = min(predEdgesDF.shape[0], trueEdgesDF.shape[0])
    if k>maxk:
        print('k should be smaller than '+str(maxk))
        return 0
    edgeweightTopk = predEdgesDF.iloc[k-1].Edgeweight
    
    nonZeroMin = np.nanmin(predEdgesDF.Edgeweight.replace(0, np.nan).values)
    bestVal = max(nonZeroMin, edgeweightTopk)
    dataset = {}
    newDF = predEdgesDF.loc[(predEdgesDF['Edgeweight'] >= bestVal)]
    dataset = set(newDF['Edges'])
    
    num=len(dataset & trueEdgesset)
    return num

#
def plot(pred_1 = predEdgesDF_1, pred_2 = predEdgesDF_2, maxk = 40000, step_len = 200):
    #Generate a array to restore difference_value
    n = maxk/step_len
    difference_value = np.zeros(int(n))
    for i in range(0,int(n)):
        difference_value[i] = num_topk(predEdgesDF = pred_1, k = (i+1)*step_len)-num_topk(predEdgesDF = pred_2, k = (i+1)*step_len)
    
    #plot
    x = np.arange(step_len,maxk+step_len,step_len)
    plt.figure(figsize=(20,12))
    plt.xticks([x for x in range(maxk + 1) if x % 10000 == 0])   
    plt.yticks([x for x in range(181) if x % 20 == 0])   
    plt.plot( x,difference_value,marker = '.',markersize = 3,linewidth=1,color = 'black')  
    plt.savefig(save_path+save_name+'.pdf')
    plt.savefig(save_path+save_name+'.png')
    plt.clf()
    

def comparison_plot(true_path='/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv',pred_path_1='/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename_unnormalization.csv',pred_path_2='/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename.csv',directed=False,save_path='/gpfs/gibbs/pi/zhao/zc354/GRN/output',save_name='/TF-TF_chromation_directed',maxk = 40000, step_len = 200)
   
    trueEdgesDF = gen_df(path = true_path)
    predEdgesDF_1 = gen_df(path = pred_path_1)
    predEdgesDF_2 = gen_df(path = pred_path_2)
    if(directed):
        trueEdgesset = set(trueEdgesDF['Edges'])
    else:
        trueEdgesset = set(trueEdgesDF['Edges'])|set(trueEdgesDF['reverseEdges']
    plot(pred_1 = predEdgesDF_1, pred_2 = predEdgesDF_2, maxk = maxk, step_len = step_len)   
    return 0
comparison_plot(true_path='/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv',pred_path_1='/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename_unnormalization.csv',pred_path_2='/gpfs/gibbs/pi/zhao/zc354/GRN/output/multiply_filename.csv',directed=False,save_path='/gpfs/gibbs/pi/zhao/zc354/GRN/output',save_name='/TF-TF_chromation_directed',maxk = 40000, step_len = 200)