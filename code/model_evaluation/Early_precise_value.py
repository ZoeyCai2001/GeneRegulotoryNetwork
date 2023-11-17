#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 13:20:21 2022

@author: chenkeran
"""

import os
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
def EarlyPrec(TFEdges = False,k=3670):
    '''
    Computes early precision for a given algorithm for each dataset.
    We define early precision as the fraction of true 
    positives in the top-k edges, where k is the number of
    edges in the ground truth network (excluding self loops).
    
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the early precision is computed.
    :type algorithmName: str
      
            
    :returns:
        A dataframe containing early precision values
        for a given algorithm for each dataset.
    '''
    rankDict = {}
    trueEdgesDF = pd.read_csv('/gpfs/gibbs/pi/zhao/yw599/GRN/BEELINE-Networks/Networks/human/Non-specific-ChIP-seq-network.csv', sep = ',',  header = 0, index_col = None)
  #  trueEdgesDF = pd.read_csv('/Users/chenkeran/Desktop/服务器上传/STRING9-GRN.csv', sep = '\t',  header = 0, index_col = None)
    trueEdgesDF = trueEdgesDF.loc[(trueEdgesDF['Gene1'] != trueEdgesDF['Gene2'])]
    trueEdgesDF.drop_duplicates(keep = 'first', inplace=True)
    trueEdgesDF.reset_index(drop=True, inplace=True)
    

     

    predDF = pd.read_csv('/gpfs/gibbs/pi/zhao/zc354/GRN/output/TF_filename_unnormalization.csv', sep=",", header=0, index_col=None)
    predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
    predDF.drop_duplicates(keep = 'first', inplace=True)
    predDF.reset_index(drop=True, inplace=True)  
        
        
    if TFEdges:
        
            # Consider only edges going out of TFs
            
            # Get a list of all possible TF to gene interactions 
        uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
        possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))

            # Get a list of all possible interactions    当然也可以self调控，把下两行去掉就行了
        possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))     
        possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
            
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
        trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]
    
        
        predDF['Edges'] = predDF['Gene1'] + "|" + predDF['Gene2']
            # limit the predicted edges to the genes that are in the ground truth
        predDF = predDF[predDF['Edges'].isin(TrueEdgeDict)]

      

            # we want to ensure that we do not include
            # edges without any edge weight
            # so check if the non-zero minimum is
            # greater than the edge weight of the top-kth
            # node, else use the non-zero minimum value.
        predDF.Edgeweight = predDF.Edgeweight.round(6)
        predDF.Edgeweight = predDF.Edgeweight.abs()

            # Use num True edges or the number of
            # edges in the dataframe, which ever is lower
        numEdges = len(trueEdges)    
        maxk = min(predDF.shape[0], numEdges)
        print(maxk)
      #  print("\nEdges considered ", maxk)
        edgeweightTopk = predDF.iloc[k-1].Edgeweight  #maxk-1， k-1
      #  print(edgeweightTopk)

        nonZeroMin = np.nanmin(predDF.Edgeweight.replace(0, np.nan).values)
        bestVal = max(nonZeroMin, edgeweightTopk)
        dataset = {}
        newDF = predDF.loc[(predDF['Edgeweight'] >= bestVal)]
        dataset = set(newDF['Gene1'] + "|" + newDF['Gene2'])
   
        Eprec = {}
        Erec = {}    
        intersectionSet = dataset.intersection(trueEdges)
        Eprec = len(intersectionSet)/len(dataset)
        Erec = len(intersectionSet)/k  #603*2/294/293=0.014, 603/152/142 =0.02793736,  3670*2/718/717=0.01425,3670/225/493=0.03308542
        print(Eprec)
        print(Erec)
        print(len(trueEdges))
        
    return(Erec)

a = np.zeros(400)
for i in range(0,3,1):
      a[i] = (EarlyPrec(TFEdges = True,k = 100*(i+1)))
      print(i)

np.savetxt('/gpfs/gibbs/pi/zhao/zc354/GRN/output/model_evaluation/TF_filename_unnormalization_EPV.txt',a)
