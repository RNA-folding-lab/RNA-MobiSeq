########################
# Author：YZ Shi
# Date:   2023.8
# Loc.:   WTU
# Description：calculate the F1 score for the score map from CODA
########################

import sys,os,math
from collections import defaultdict
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import argparse,gzip
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
from math import sqrt
import time
t0 = time.time()
basepair = ['AU','UA','GC','CG','GU','UG']
######### Read sequence ##########
def read_seq(file):
    from itertools import islice
    seq = []
    with open(file,'r') as f:
        for line in islice(f,1,None):
            seq = line.strip('\n')
    seqLen = len(seq)
    print(' sequence:',seq,'\n','seqLen:',seqLen)
    return seq, seqLen
## Read the final_score and preprocessing 
def Read_score(score_path,L):
    score_matrix = np.zeros((L,L),dtype=float)      
    with open(score_path,'r') as df:
        for row, line in enumerate(df):
            if row >= L:
                break 
            score_list = list(map(float, line.strip().split('\t')))  
            score_matrix[row, :len(score_list)] = score_list[:L] 
    return score_matrix
def Score2pair(score,threshold=0.0):
    pair = np.where(score > threshold, 1, 0)
    return pair
###based on pair info, get the stems
def pair2stem(pairs):    
    Stems = []
    edge = np.stack(pairs.nonzero())
    if edge.size > 0:
        i_consistancy = (edge[0][1:] - edge[0][:-1]) == 1
        j_consistancy = (edge[1][1:] - edge[1][:-1]) == -1
        consistancy = np.pad(i_consistancy & j_consistancy, (1, 0), constant_values=True)
        segments = np.cumsum(~consistancy)
        for segment in np.unique(segments):
            Stems.append(edge[:, segments == segment].T)
    return Stems    
## Calculate the energy for the scoring matrix
def Energy_calculation(score,seq,wt=1.0):
    delt_H = {'AAUU':-6.82,'UUAA':-6.82,'AUUA':-9.38,'UAAU':-7.69,'CUGA':-10.48,'AGUC':-10.48,
        'CAGU':-10.44,'UGAC':-10.44,'GUCA':-11.40,'ACUG':-11.40,'GACU':-12.44,'UCAG':-12.44,
        'CGGC':-10.64,'GGCC':-13.39,'CCGG':-13.39,'GCCG':-14.88,'AGUU':-3.96,'UUGA':-3.96,
        'AUUG':-7.39,'GUUA':-7.39,'CGGU':-5.56,'UGGC':-5.56,'CUGG':-9.44,'GGUC':-9.44,
        'GGCU':-7.03,'UCGG':-7.03,'GUCG':-11.09,'GCUG':-11.09,'GAUU':-10.38,'UUAG':-10.38,
        'GGUU':-17.82,'UUGG':-17.82,'GUUG':-13.83,'UGAU':-0.96,'UAGU':-0.96,'UGGU':-12.64}
    delt_S = {'AAUU':-19.0,'UUAA':-19.0,'AUUA':-26.7,'UAAU':-20.5,'CUGA':-27.1,'AGUC':-27.1,
        'CAGU':-26.9,'UGAC':-26.9,'GUCA':-29.5,'ACUG':-29.5,'GACU':-32.5,'UCAG':-32.5,
        'CGGC':-26.7,'GGCC':-32.7,'CCGG':-32.7,'GCCG':-36.9,'AGUU':-11.6,'UUGA':-11.6,
        'AUUG':-21.0,'GUUA':-21.0,'CGGU':-13.9,'UGGC':-13.9,'CUGG':-24.7,'GGUC':-24.7,
        'GGCU':-16.8,'UCGG':-16.8,'GUCG':-28.8,'GCUG':-28.8,'GAUU':-31.8,'UUAG':-31.8,
        'GGUU':-56.7,'UUGG':-56.7,'GUUG':-46.9,'UGAU':-1.8,'UAGU':-1.8,'UGGU':-38.9}
    #basepair = set(delt_H.keys())
    pair = Score2pair(score,0.5)    
    stems = pair2stem(pair)
    count = np.count_nonzero(pair > 0)
    u_bp, u_st = 0.0, 0.0
    for stem in stems:
        L = len(stem)
        for j in range(L):
            pair1, pair2 = stem[j]
            seq_pair = str(seq[pair1]) + str(seq[pair2])
            if seq_pair in basepair:
                u_bp -= score[pair1, pair2]
                if j != L - 1:
                    pair3, pair4 = stem[j + 1]
                    seq_pair_next = str(seq[pair3]) + str(seq[pair4])
                    if seq_pair_next in basepair:
                        s1 = str(seq[pair1])+str(seq[pair3])+str(seq[pair2])+str(seq[pair4])
                        U_stack = delt_H.get(s1, 0) - 298 * 0.001 * delt_S.get(s1, 0)
                        u_st += wt * min(U_stack, 0.0)      # Avoid positive stacking energy
    
    U = (u_bp / count) + u_st if count > 0 else 0.0
    print(u_bp / count, u_st, U)
    return U                   
####Read base pair from DSSR####
def native_bp(bp_file,L,val=1):
    native = open(bp_file,"r")
    lines = native.readlines()
    bplist = []
    for line in lines:   ##
        line = line.rstrip("\n") 
        sublist = line.split(" ")
        if int(sublist[0])<int(sublist[2]) and int(sublist[2])>0:
            bplist.append(sublist)
    #print(bplist)
    pairs = np.zeros((L, L))
    signal = []
    for i in range(1,L):
        for j in range(i):
            for row in bplist:
                if int(row[2])>0 and (int(row[0])==j+1 and int(row[2])==i+1):
                    if (j,i) not in signal:
                        signal.append((j,i))   #
                        pairs[i,j] = val
    return pairs

####Score evaluation####
def Score_evaluation(mat1,mat2,L,val_max,val_min):
    rightuptri = []
    for i in range(L):
       for j in range(L):
          if i<j:
             rightuptri.append((i,j))            
    thlst = np.arange(val_min,val_max,0.1).tolist()
    tprlst, fprlst, pcslst, F1lst, MCClst = [],[],[],[],[]
    for threshold in thlst:
        TP = len([i for i in rightuptri if (mat1[i[0],i[1]] >= threshold and mat2[i[1],i[0]]==val_max)])
        FN = len([i for i in rightuptri if (mat1[i[0],i[1]] <= threshold and mat2[i[1],i[0]]==val_max)])
        FP = len([i for i in rightuptri if (mat1[i[0],i[1]] >= threshold and mat2[i[1],i[0]]!=val_max)])
        TN = len([i for i in rightuptri if (mat1[i[0],i[1]] <= threshold and mat2[i[1],i[0]]!=val_max)])
        TPR = TP/(TP+FN)
        FPR = FP/(FP+TN)
        precision = TP/(TP+FP)
        F1 = 2*TP/(2*TP+FN+FP)
        MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        tprlst.append(TPR)
        fprlst.append(FPR)
        pcslst.append(precision)
        F1lst.append(F1)
        MCClst.append(MCC)
    evadf = pd.DataFrame({'threshold':thlst,'TPR(recall)':tprlst,'FPR':fprlst,'precision':pcslst,'F1':F1lst,'MCC':MCClst})
    Max_F1 = max(evadf['F1'])
    Max_MCC = max(evadf['MCC']) 
    return Max_F1
if __name__ == '__main__':        
    
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=str, help='The score data path (e.g., 1.0/0.1/2.0/score)')
    parser.add_argument('sequence', type=str, help='The fasta sequence file (e.g., data/5TPY.fasta)')
    parser.add_argument('-o', '--out_path', type=str, default='./', help='The output path ')
    parser.add_argument('-n', '--native_secondary', type=str, default=None, 
        help='The path to the native secondary structure file (optional,e.g., 5TPY-native.txt) for result evaluation (e.g., F1 score  or score map)')

    args = parser.parse_args()

    score_path = args.data
    score_file = score_path + 'score/final_score.txt'
    sequence_file = args.sequence
    output_path = args.out_path
    sequence, seqLen = read_seq(sequence_file)    
    score = Read_score(score_file,seqLen)
    U = Energy_calculation(score,sequence,1.0)
    output_file_path = f"{output_path}U.txt"
    if args.native_secondary:
        native_bp_file = args.native_secondary   
        pairs = native_bp(native_bp_file,seqLen,1) 
        F1 = Score_evaluation(score,pairs,seqLen,1,0)
        
        with open(output_file_path, 'a') as output_file:
            output_file.write(f"{score_path}\t{round(U, 3)}\t{round(F1, 3)}\n")

        print(F1,'\t',U)
    else:
        with open(output_file_path, 'a') as output_file:
            output_file.write(f"{score_path}\t{round(U, 3)}\n")
    run_time = time.time() - t0
    print("The run time: ",run_time,"s")

