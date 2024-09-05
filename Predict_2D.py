########################
# Author：YZ Shi
# Date:   2022.8
# Loc.:   SZBL
# Description：Do MC simulation to predict RNA contact map 
#              using CODA score as constraint
# version: 
#       1.0-Weaken the single bp & GU_END  
########################
import sys, os, argparse
from collections import defaultdict  
from numpy import exp
import numpy  as np       
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import random
from math import sqrt
from itertools import islice
import time
t0 = time.time()
#######hyper-parameter######
wt = 0.76         # The ratio of MC energy
U_lone_m = 0.8    # penalty of lone base pair
U_AU_end_m = 0.8  # penalty of AU/GU at end
######### Read sequence ##########
def read_seq(file):
    with open(file,'r') as f:
        seq = next(islice(f, 1, None)).strip()
    seqLen = len(seq)
    print("Read sequence successfully (seq length): ",seqLen)
    return seq,seqLen
## Read the final_score calculated from CODA 
def Read_score(score_path,L):
    score_matrix = np.zeros((L,L),dtype=float) 
    row = 0           
    with open(score_path,'r') as df:
        for line in df:
            score_list = list(map(float, line.strip().split('\t')))
            if len(score_list) != L:
                raise ValueError(f"Row {row} length {len(score_list)} does not match expected length {L}")
            score_matrix[row:] = score_list
            row+=1  
    val_max, val_min = np.nanmax(score_matrix), np.nanmin(score_matrix)
    print("Score max/mean/min: ",val_max,np.nanmean(score_matrix),val_min)  
    return score_matrix   

###based on pair info, get the stems
def pair2stem(pairs):    
    Stems = []
    edge = np.stack(pairs.nonzero())
    if edge.size > 0:        
        i_consistancy = (edge[0][1:] - edge[0][:-1])==1
        j_consistancy = (edge[1][1:] - edge[1][:-1])==-1
        consistancy = np.pad(i_consistancy & j_consistancy, (1,0), constant_values=True)
        segment = np.cumsum(~consistancy)
        N = len(segment) 
        for i in range(N):
            stem = (edge[:, segment==i].T)
            Stems.append(stem) 
    return Stems
    
## Monte Carlo simulating Anneaning
def exponential_decay(T, alpha):
    return alpha * T

def linear_decay(T, beta):
    return T - beta

def MC_SA(score_matrix,seq,high_T,low_T,alpha,step,Out):   
    seqLen = len(seq)
    ## Thermodynamics parameters (Refs.: Turner et al. Biochemistry 1998, 37, 14719-14735 & J. Mol. Biol. (1999) 288, 911-940 & Biochemisty 2012, 3508-3522)
    '''
    stack = {'AAUU':-0.93,'UUAA':-0.93,'AUUA':-1.1,'UAAU':-1.33,'CUGA':-2.08,'AGUC':-2.08,
             'CAGU':-2.11,'UGAC':-2.11,'GUCA':-2.24,'ACUG':-2.24,'GACU':-2.35,'UCAG':-2.35,
             'CGGC':-2.36,'GGCC':-3.26,'CCGG':-3.26,'GCCG':-3.42,'AGUU':-0.35,'UUGA':-0.35,
             'AUUG':-0.90,'GUUA':-0.90,'CGGU':-1.25,'UGGC':-1.25,'CUGG':-1.77,'GGUC':-1.77,
             'GGCU':-1.80,'UCGG':-1.80,'GUCG':-2.15,'GCUG':-2.15,'GAUU':-0.51,'UUAG':-0.51,
             'GGUU':-0.25,'UUGG':-0.25,'GUUG':+0.72,'UGAU':-0.39,'UAGU':-0.39,'UGGU':-0.57}
    '''
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
    f1 = open(Out + 'mc_energy.txt', 'w')
    f2 = open(Out + 'secondary_structure.txt', 'w')
    f3 = open(Out + 'secondary_structure_minE.txt', 'w')
    pair = np.zeros([seqLen,seqLen],dtype=np.int32)
    base = [0 for _ in range(seqLen)]
    U_min = 1000.0    
    U_before = 0.0    ###Energy Before Changes    
    tt = 0
    random.seed(seed) ##set random_seed for multi-MC & do MC simulations
    candidate_pairs = [(i1, i2) for i1 in range(seqLen - 4) for i2 in range(i1 + 4, seqLen) if seq[i1] + seq[i2] in basepair]
    T = high_T
    MC_step = step * seqLen
    print('random_seed:',seed,'MC_steps: ', MC_step)
    while T > low_T:
        TK = T + 273.15
        D = 0.002*TK 
        if T <= 25.0:
            MC_step = step * seqLen * 2
        #U_lone = max(U_lone_m - 0.008 * T, 0)
        #U_AU_end = max(U_AU_end_m - 0.008 * T, 0)
        #U_lone = U_lone_m / (1+exp(0.08*(T-50.0))) 
        #U_AU_end = U_AU_end_m / (1+exp(0.08*(T-50.0)))
        U_lone = U_lone_m
        U_AU_end = U_AU_end_m
        for t in range(MC_step):          
            i1, i2 = random.choice(candidate_pairs)  ###Randomly choose two possible pairing bases 
            ###Changes the pairing bases
            pair_a,base_a = update_pairing(i1,i2,base,pair,seqLen,T)
            ###Energy after Changes based on Stems
            U_after = 0.0
            stems_a = pair2stem(pair_a)
            for i in range(len(stems_a)):
                L=len(stems_a[i])            
                for j in range(L):
                    pair1=stems_a[i][j][0]
                    pair2=stems_a[i][j][1]
                    U_after += -score_matrix[pair1][pair2]
                    if L==1:
                        U_after += wt*U_lone
                    else:
                        if j!=L-1:
                            pair3=stems_a[i][j+1][0] 
                            pair4=stems_a[i][j+1][1]
                            s1 = str(seq[pair1])+str(seq[pair3])+str(seq[pair2])+str(seq[pair4])
                            U_stack = delt_H[s1]-TK*0.001*delt_S[s1]
                            U_after += wt * (U_stack if U_stack < 0 else 0)
                        if (j==0 or j==L-1):
                            if ((seq[pair1]=='U' or seq[pair2]=='U')):
                                U_after += wt*U_AU_end
            ###Metroplis
            delt = U_after - U_before
            if delt <= 0 or exp(-delt / D) > random.random():
                base, pair, stems = base_a, pair_a, stems_a
                U_before=U_after
            if U_before < U_min:
                U_min, pair_minE = U_before, pair.copy();
            if (t+1) % (MC_step*0.01) == 0:
                tt += 1
                f1.write(f"{tt} {T} {t} {round(U_before, 3)}\n")
                f1.flush()
            if (t+1)%(MC_step*0.5)==0:    
                print(f"{T:.1f} {t} {U_before:.2f}")
    
        for i in range(seqLen):
            for j in range(seqLen):
                if pair[i,j]>0:
                   f2.write(f"{T} {t} {i} {j} {seq[i]} {seq[j]}\n")
                   f2.flush()            
        T = exponential_decay(T, alpha)
    f1.close()
    f2.close()
    for i in range(seqLen):
        for j in range(seqLen):            
            if pair_minE[i,j]>0:
               f3.write(f"{i} {j} {seq[i]} {seq[j]}\n")
               f3.flush()
    f3.close()
    ### Save the contact_map of conf with minE
    np.save(Out+"contact_minE.npy",pair_minE)
    ###Plot the energy vs MC_step
    energy = pd.read_table(Out + 'mc_energy.txt',sep =' ',header=None)
    energy.columns = ['t0','T','t','U']
    fig = plt.figure(figsize = (10,5),dpi = 80)
    plt.plot(energy['t0'],energy['U'])
    plt.xlabel('MC steps')
    plt.ylabel('Energy')
    plt.savefig(Out + 'MC_energy.png',dpi=100)
    #plt.show()
    plt.close()
    dbn = pair2dbn(pair_minE, seq,Out)
    plot_heatmap(Out+"MC_minE_contact_map",pair_minE,seqLen,1,0,1) 
    return pair, pair_minE, U_min, dbn
###Change/update the base-pairing
def update_pairing(i1,i2,base,pair,seqLen,T):
    pair_a = pair.copy()
    base_a = base.copy()

    def break_pair(index):
        for j in range(seqLen):
            if pair[index, j] == 1 or pair[j, index] == 1:
                pair_a[index, j] = 0
                pair_a[j, index] = 0
                base_a[j] = 0
                break
    if base[i1] == 0 and base[i2] == 0:
        pair_a[i1, i2] = 1
        base_a[i1] = 1
        base_a[i2] = 1
    else:
        if base[i1] == 1:
            break_pair(i1)
        if base[i2] == 1:
            break_pair(i2)
        prob_pair = 1 - exp(-0.02 * T)
        if random.random() <= prob_pair:
            pair_a[i1, i2] = 1
            base_a[i1] = 1
            base_a[i2] = 1

    return pair_a, base_a

###add line to one figure
def abline(slope, intercept):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, ls='--', c='black')
##Score_Heatmap plot
def plot_heatmap(name,mat,L,val_max,val_min=-10.0,ismask=0):
    mask = np.zeros((L, L))
    if ismask==1:
        for i in range(1,L):
            for j in range(i):
                if mat[i, j] != val_max:                
                    mask[i, j] = True
                if mat[j, i] != val_max:                
                    mask[j, i] = True
    #mask[np.tril_indices_from(mask)] = True  #set the below as 1
    plt.figure(figsize=(7, 5))
    #sns.set(font_scale=1.5)
    ax = sns.heatmap(mat,mask=mask,cmap="YlGn",vmin=val_min,vmax=val_max)   
    abline(1,0) 
    plt.xlim(0,L)
    plt.ylim(L,0)
    plt.xticks(rotation=90)
    #plt.yticks(range(0,L,5))
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    plt.savefig(name +'.png', dpi=300)   

############ output the dbn format of 2D structure with minimum energy ##########
def pair2dbn(pair_matrix, seq,Out):
    index = np.argwhere(pair_matrix == 1)
    pair_final = index[index[:,0]<index[:,1]]
    dbn = ['.'] * len(seq)
    open_brackets = ['(', '[', '{', '<']
    close_brackets = [')', ']', '}', '>']
    for pair in pair_final:
        i, j = pair 
        if dbn[i] == '.' and dbn[j] == '.':
            existing_open = [b for b in dbn[i:j] if b in open_brackets]
            existing_close = [b for b in dbn[i:j] if b in close_brackets]
            for open_b, close_b in zip(open_brackets, close_brackets):
                if open_b not in existing_open and close_b not in existing_close:
                    dbn[i] = open_b
                    dbn[j] = close_b
                    break               
    dbn_str = "".join([x for x in dbn])
    with open(Out + '2D.dot', 'w') as out:
        out.write(dbn_str)
    return dbn_str
## Read the native base-base interaction derived from DSSR
## only used to evaluate the results
def native_bp(file,L):
    with open(file, "r") as native:
        lines = native.readlines()
    bplist = [line.split() for line in lines if int(line.split()[0]) < int(line.split()[2]) and int(line.split()[2]) > 0]
    pairs = np.zeros((L, L))
    for row in bplist:
        j, i = int(row[0]) - 1, int(row[2]) - 1
        pairs[i, j] = 1
        pairs[j, i] = 1
    return pairs
def Integrate_map(map1, map2, L):
    final_map = np.zeros((L, L))
    # Copy the lower triangle from map1 to final_map
    final_map[np.tril_indices(L, -1)] = map1[np.tril_indices(L, -1)]
    # Copy the upper triangle from map2 to final_map
    final_map[np.triu_indices(L, 1)] = map2[np.triu_indices(L, 1)]
    return final_map
## Evaluate the results
def F1_calculation(score_matrix,pair,rightuptri):
    TP = sum(score_matrix[i, j] ==1 and pair[j, i] == 1 for i, j in rightuptri)
    FN = sum(score_matrix[i, j] ==0 and pair[j, i] == 1 for i, j in rightuptri)
    FP = sum(score_matrix[i, j] ==1 and pair[j, i] == 0 for i, j in rightuptri)
    TN = sum(score_matrix[i, j] ==0 and pair[j, i] == 0 for i, j in rightuptri)

    TPR = TP / (TP + FN) if TP + FN > 0 else 0
    FPR = FP / (FP + TN) if FP + TN > 0 else 0
    precision = TP / (TP + FP) if TP + FP > 0 else 0
    F1 = 2 * TP / (2 * TP + FN + FP) if 2 * TP + FN + FP > 0 else 0
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 0
    print("F1:",F1,"MCC",MCC)
    return F1, MCC, precision, TPR, FPR, TP, FP, FN, TN
def Evaluation_score(pair,pair_minE,seq,native_file,U_min,Out):
    L = len(seq)
    pair_native = native_bp(native_file,L)
    pair_matrix_minE = Integrate_map(pair_native, pair_minE, L)
    pair_matrix = Integrate_map(pair_native, pair, L)
    rightuptri = [(i, j) for i in range(L) for j in range(i+1, L)]
    with open(Out + 'MC_result.txt', 'w') as output_result:
        header = '\tF1\tMCC\tprecision\tTPR\tFPR\tTP\tFP\tFN\tTN\n'
        output_result.write(header)
        ## Evaluation contact map from low T
        print("Result for Lowest T:")     
        F1, MCC, precision, TPR, FPR, TP, FP, FN, TN = F1_calculation(pair_matrix,pair_native,rightuptri)
        lowT_result = f'Lowest T:\t{F1:.3f}\t{MCC:.3f}\t{precision:.3f}\t{TPR:.3f}\t{FPR:.3f}\t{TP}\t{FP}\t{FN}\t{TN}\n'
        output_result.write(lowT_result)

        ## Evaluation for contact map from minE
        print("Result for Lowest U")
        F1, MCC, precision, TPR, FPR, TP, FP, FN, TN = F1_calculation(pair_matrix_minE,pair_native,rightuptri)
        lowU_result = f'Lowest U:\t{F1:.3f}\t{MCC:.3f}\t{precision:.3f}\t{TPR:.3f}\t{FPR:.3f}\t{TP}\t{FP}\t{FN}\t{TN}\n'
        output_result.write(lowU_result)

    plot_heatmap(Out+"MC_contact_map_native",pair_matrix_minE,seqLen,1,0,1) 
    #plot_heatmap(Out+"MC_END_contact_map",pair_matrix,seqLen,1,0,1) 


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=str, help='The score data file (e.g., score/final_score.txt)')
    parser.add_argument('sequence', type=str, help='The fasta sequence file (e.g., data/5TPY.fasta)')
    parser.add_argument('-o', '--out_path', type=str, default='./', help='The output path ')
    parser.add_argument('-n', '--native_secondary', type=str, default=None, 
        help='The path to the native secondary structure file (optional,e.g., 5TPY-native.txt) for result evaluation (e.g., F1 score  or score map)')
    parser.add_argument('-ht', '--high_T', type=float, default=120.0, 
        help='The initial high temperature for MC annealing')
    parser.add_argument('-lt', '--low_T', type=float, default=22.0, 
        help='The end low temperature in MC annealing')
    parser.add_argument('-a', '--alpha', type=float, default=0.88, 
        help='The annealing rate (0.8-1.0)')
    parser.add_argument('-s', '--step', type=int, default=3000, 
        help='The step in each temperature in MC annealing (total_step = step * Len(seq))')
    parser.add_argument('-seed', '--seed', type=int, default=0, 
        help='The random seed in current simulation')
    
    args = parser.parse_args()
    
    score_path = args.data
    sequence_file = args.sequence
    
    seed = args.seed
    if args.seed:
        seed = args.seed
        Outpath = args.out_path + '2D/' + str(seed) + '/'
    else:
        Outpath = args.out_path + '2D/'
    os.makedirs(Outpath, exist_ok=True)   ##mkdir Outpath
    high_T, low_T, alpha = args.high_T, args.low_T, args.alpha
    MC_step = args.step
    basepair = ['AU','UA','GC','CG','GU','UG']

    sequence, seqLen = read_seq(sequence_file)
    score = Read_score(score_path,seqLen)
    pair, pair_minE, U_min, dbn = MC_SA(score,sequence,high_T,low_T,alpha,MC_step,Outpath)
    print(f'The_lowest_U:\t{U_min:.3f}\n2D: {dbn}')
    if args.native_secondary:
        native_file = args.native_secondary
        Evaluation_score(pair,pair_minE,sequence,native_file,U_min,Outpath)
    
    run_time = time.time() - t0
    print("The run time: ",run_time,"s")




