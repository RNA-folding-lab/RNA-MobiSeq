########################
# Author：YZ Shi, Z Zhang, D Luo, YQ Huang
# Date:   2022.8
# Loc.:   SZBL & WTU
# Description: Analyze the processed high-throughput sequencing data to infer base-pairing
# Version: 0.1 (Code optimization @ Griffith University by YZ in 2024)
########################
import sys, os, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
import seaborn as sns
from math import sqrt,exp,log
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR
import time
t0 = time.time()
######### Read sequence ##########
def read_seq(file):
    seq = []
    with open(file,'r') as f:
        next(f)
        seq = f.readline().strip()
    seqLen = len(seq)
    print("Read sequence successfully (seq length): ",seqLen)
    return seq,seqLen
##Get the wide_type seq. for each position 
def WT_AA(num):
    return sequence[num]
##Read mutation file (from Zhang Z) & get the single/double/triple mutation data
def read_mut_data(file):
    single_data = []
    double_data = []
    #triple_data = []
    Cwt = None
    with open(file,'r+') as f:
        for line in f:
            spt = line.strip().split('\t')
            mutNum = int(spt[0])
            if mutNum==0:
                Cwt = float(spt[1])
            elif mutNum==1:
                single_data.append(spt[1:7])  ##New data 1:5
            elif mutNum==2:
                double_data.append(spt[1:11])  ##1:9
            #elif mutNum==3:
            #    triple_data.append(spt[1:12])  #1:10
    return Cwt, single_data, double_data #, triple_data   
def Read_mut_dataframe(file):
    #Cwt, single_data, double_data, triple_data = read_mut_data(file)  
    Cwt, single_data, double_data = read_mut_data(file)  
    single = pd.DataFrame(single_data, 
        columns=['Pos','Mut','str1','input','output','fitness']).drop(columns=['str1'])
    double = pd.DataFrame(double_data, 
        columns=['Pos1','Pos2','Mut1','Mut2','input','output','str1','fitness1','fitness2','fitness']).drop(columns=['str1'])
    #triple = pd.DataFrame(triple_data, 
    #    columns=['Pos1','Pos2','Pos3','Mut','input','output','fitness1','fitness2','fitness3','str1','str2','fitness']).drop(columns=['str1','str2'])
    single = single.astype({'Pos': int, 'Mut': str, 'input': float, 'output': float, 'fitness': float})
    double = double.astype({'Pos1': int, 'Pos2': int, 'Mut1': str, 'Mut2': str, 'input': float, 'output': float, 'fitness1': float, 'fitness2': float, 'fitness': float})
    #triple = triple.astype({'Pos1': int, 'Pos2': int, 'Pos3': int, 'Mut1': str, 'Mut2': str, 'Mut3': str, 'input': float, 'output': float, 'fitness1': float, 'fitness2': float, 'fitness3': float, 'fitness': float})
    double = double.sort_values(['Pos1','Pos2'],ascending = [True,True])
    if min(double['Pos1'])!=0:   ##keep the first pos is 0
        double['Pos1'] -= 1
        double['Pos2'] -= 1
    double = double.reset_index(drop=True)
    ####put the native seq into the double data
    double['wt1'] = double.apply(lambda x:WT_AA(x['Pos1']),axis=1)
    double['wt2'] = double.apply(lambda x:WT_AA(x['Pos2']),axis=1) 
    print("Mut_data read successfully:")
    print("lenght of single/double data: ",single.shape[0],double.shape[0]) 
    return Cwt,single,double #,triple

########## one hot encode for sequences ##########
def one_hot_encode(sequences):
    enc = OneHotEncoder(categories=[list('ACGT')], handle_unknown='ignore')
    seq_encode = enc.fit_transform(np.array(list(sequences)).reshape(-1, 1)).toarray()
    return seq_encode
####### Just consider base pairs: AT/TA/GC/CG/GT/TG
def bpcheck(n1,n2):
    valid_pairs = {('G', 'C'), ('C', 'G'), ('A', 'T'), ('T', 'A'), ('G', 'T'), ('T', 'G')}
    return (n1, n2) in valid_pairs
##### generate gaussian pdf
def normfun(x,mu,sigma):
    pdf = np.exp(-((x - mu)**2) / (2* sigma**2)) / (sigma * np.sqrt(2*np.pi))
    return pdf

#### Filter the double data
# delete outliers based on boxplot of input/output
def box_outlier(data,lower_output_cut=5):
    df = data.copy(deep=True)
    out_index = set()
    for col in df.columns[[4,5]]:     # input/output/fitness1/fitness2/fitness  
        Q1 = df[col].quantile(q=0.25)       # lower quartile
        Q3 = df[col].quantile(q=0.75)       # upper quartile
        IQR = Q3 - Q1
        low_whisker = max(lower_output_cut, Q1 - 1.5 * IQR)  # Lower whisker edge
        up_whisker = Q3 + 25 * IQR  # Upper whisker edge
        #print(low_whisker,up_whisker)
        rule = (df[col] > up_whisker) | (df[col] < low_whisker) 
        out_index.update(df.index[rule].tolist()) 
    df.drop(out_index, inplace=True)
    df = df[(df['output'] > lower_output_cut) & (df['input'] > lower_output_cut)]
    df = df.reset_index(drop=True)
    print("Data_shape Before/After Filter (In fitting): ",data.shape[0],' -> ',df.shape[0])
    return df
def Keep_WC_pairing(df):
    valid_wt_pairs = {('A', 'U'), ('G', 'C'), ('C', 'G'), ('U', 'A'), ('U', 'G'), ('G', 'U')}
    valid_mut_pairs = {('A', 'T'), ('G', 'C'), ('C', 'G'), ('T', 'A'), ('T', 'G'), ('G', 'T')}
    df = df[df[['wt1', 'wt2']].apply(tuple, axis=1).isin(valid_wt_pairs) & df[['Mut1', 'Mut2']].apply(tuple, axis=1).isin(valid_mut_pairs)]
    return df
def Filter_data(data):
    df = data.copy(deep=True)
    df = Keep_WC_pairing(df)
    #df.drop(df[(df['fitness']>df['fitness1']) & (df['fitness']>df['fitness2'])].index,inplace=True)
    #df.drop(df[df['fitness']<0.01].index,inplace=True)
    df.drop(df[(df['fitness']>5*(df['fitness1']+df['fitness2']))].index,inplace=True)
    print("Data Shape Before/After filter (In prediction): ",data.shape[0],' -> ',df.shape[0])
    return df.reset_index(drop=True)
###Build the X,Y dataset for SVR regression
def encode_data(data):
    tempdata=one_hot_encode(data)
    temp_dataframe=pd.DataFrame(tempdata,columns=list('acgt'))
    return (temp_dataframe)
def XY_dataset(data):
    df = data.copy(deep=True)
    ## Build training set: I-fitness~fitness1+fitness2; 
    #X = df.loc[:,['fitness1','fitness2']]
    #Y = df.loc[:,['fitness']]

    ## Build training set: II-fitness~fitness1+fitness2+seq_encode+position
    x_temp = df[['Pos1','Pos2','fitness1','fitness2']]
    #x_temp = df.loc[:,['fitness1','fitness2']]
    Y = df[['fitness']]
    #df_temp=df.loc[:,['fitness1','fitness2','fitness']]
    #df_temp=pd.DataFrame(MinMaxScaler().fit_transform(df_temp),columns=['fitness1','fitness2','fitness'])
    #print(df_temp.head())
    #x_temp = df_temp.iloc[:,[0,1]]
    #Y = df_temp.iloc[:,2]
    en_mut1=encode_data(df['Mut1'])  #one-hot for sequence
    en_mut2=encode_data(df['Mut2'])
    en_wt1=encode_data(df['wt1'])
    en_wt2=encode_data(df['wt2'])
    X_temp=pd.concat([x_temp,en_mut1,en_mut2,en_wt2,en_wt2],axis=1)
    standard_data=MinMaxScaler().fit_transform(X_temp)
    X=pd.DataFrame(standard_data)
    return X,Y

## using SVR() --> fitness~fitness1+fitness2 
def svr_model_fit(x,y):            
    svr_rbf = SVR(kernel='rbf', C=C_value, gamma=gamma_value)
    svr_rbf.fit(x,y.values.ravel())
    #y_pred = svr_rbf.predict(x)
    return svr_rbf

## Plot delt-fitness or delt-score distribution
def plot_distribution(df,Out):
    fig = plt.figure(figsize=(15, 5))
    plt.subplot(1, 3, 1)
    plt.plot(df['real'], df['pred'], 'o', markersize=1)
    abline(1, 0)
    plt.xlabel('Observed')
    plt.ylabel('Predicted')
    
    plt.subplot(1, 3, 2)
    plt.plot(df['real'], label="real")
    plt.plot(df['pred'], label="predict")
    plt.legend()
    plt.xlabel('Site')
    plt.ylabel('Fitness')
    
    plt.subplot(1, 3, 3)
    df.delt_m.plot(kind='hist', bins=50, color='blue', density=True)
    df.delt_m.plot(kind='kde', color='red')
    plt.xlabel('delt')
    plt.tight_layout()
    plt.savefig(Out + 'distribution_delt_fitness_svr.png', dpi=600)

def Plot_delt_distribution(deltList_temp,listA,listB,mean0,sd0,Out):
    meanA,sdA = np.mean(listA), np.std(listA)
    meanB,sdB = np.mean(listB), np.std(listB)
    ###ditribution of classified score
    fig = plt.figure(figsize = (10,5))
    plt.subplot(1,2,1)
    plt.hist(deltList_temp,color='blue',density=True)
    x = np.arange(min(deltList_temp),max(deltList_temp),0.1)
    y1 = normfun(x, mean0, sd0)
    plt.plot(x,y1,color='red',linewidth = 3,label="all")
    plt.xlabel('Delt_score')
    plt.legend()
    plt.subplot(1,2,2)
    plt.hist(listA,color='blue',density=True)
    y2 = normfun(x, meanA, sdA)
    plt.plot(x,y2,linewidth = 2,label="listA")
    plt.hist(listB,color='green',density=True)
    y3 = normfun(x, meanB, sdB)
    plt.plot(x,y3,linewidth = 2,label="listB")
    plt.xlabel('Delt_score')
    plt.legend()
    plt.tight_layout()
    plt.savefig(Out + 'distribution_delt_score.png',dpi=600)
##Training the regression model (SVR or MLP)
def Training_regression_model(data,Out):
    df0 = data.copy(deep=True)
    #df1 = Filter_data(data)
    df1 = box_outlier(data)    
    X1,Y1 = XY_dataset(df1)
    svr = svr_model_fit(X1,Y1)
    
    df0 = Filter_data(df0)
    X0,Y0 = XY_dataset(df0)
    y_pred = svr.predict(X0)
    #y_pred = svr.predict(X1)

    #calculate the delt between pred. & real
    df = df0.copy(deep=True)
    df['real'] = Y0.round(3)
    df['pred'] = y_pred.round(3)
    df['delt'] = (df['real'] - df['pred']).round(3)
    df['delt_m'] = (df['delt']/(df['pred']+0.2)).round(3)
    #print(df.head())
    #plot the relation or distribution of pred & real fitness
    plot_distribution(df,Out)
    ## output the file including pred_fitness & delt
    df2 = df[['Pos1','Pos2','wt1','wt2','Mut1','Mut2','delt','delt_m']]
    df2.to_csv(Out + 'pred_fitness.txt',sep = '\t', index = False)
    return df2
##Calculate P(fitness) based on two Gaussian distributions
def DeltToScore(ra, mean1, sd1, mean2, sd2,p1):
    p_unpaired = normfun(ra, mean1, sd1)
    p_paired = normfun(ra, mean2, sd2)
    ##p1 = 0.67/(L-1)
    return log(p_paired/(p_paired*p1 + p_unpaired*(1-p1)))
##output the score matrix
def Score_matrix(deltList,posList,L,listA,listB,Out):
    paired_prob = len(listB) / len(listA)
    print(f"p(paired): {paired_prob}   lenA: {len(listA)}   lenB: {len(listB)}")
    meanA,sdA = np.mean(listA), np.std(listA)
    meanB,sdB = np.mean(listB), np.std(listB)
    #print('meanA/sdA:',meanA,sdA,'\n','meanB/sdB:',meanB,sdB)
    score_matrix = np.full((L, L), np.nan) #[[np.NAN for i in range(L)] for j in range(L)]
    for k, delt in enumerate(deltList):
        final_score = np.mean([DeltToScore(d, meanA, sdA, meanB, sdB, paired_prob) for d in delt])
        i, j = posList[k]
        score_matrix[i, j] = final_score
        score_matrix[j, i] = final_score
    with open(Out + 'initial_score.txt', 'w') as f:
        for row in score_matrix:
            f.write('\t'.join(f"{score:.2f}" if not np.isnan(score) else "nan" for score in row) + '\n')
    val_max, val_min, val_med= np.nanmax(score_matrix), np.nanmin(score_matrix), np.nanmedian(score_matrix)
    plot_heatmap_score(Out + "Score_map",score_matrix,L,val_max,val_med)
    return score_matrix
## Calculate the score for each position (i,j):  (sum or max delt_m) 
def Contact_prediction(data,L,Out):
    df = data.copy(deep=True)
    df = Training_regression_model(df,Out)
    subarray = np.zeros((L,L))
    deltList_temp, deltList, posList = [], [], []
    with open(Out + 'score.txt', 'w') as Output_Score:
        for i in range(L-1):
            for j in range(i+3,L):
                 sub_df = df.loc[(df['Pos1'] == i)&(df['Pos2'] == j),['delt_m']]
                 score_temp = []
                 if not sub_df.empty:            
                    #subarray[i,j] = round(max(sub_df['delt_m']),3)  ##take max value as score
                     score_temp = sub_df['delt_m'].tolist()
                     score_sum = round(sum(score_temp),3) 
                     subarray[i,j] = score_sum  ##take sum of the score with AT/GC/GU at each pos1-pos2
                     deltList.append(score_temp)
                     posList.append((i,j))
                 for score in score_temp:
                        deltList_temp.append(score)
                        Output_Score.write(f"{i}\t{j}\t{score_sum}\t{score}\n")
    mean0 = np.mean(deltList_temp)
    sd0 = np.std(deltList_temp)
    print("deltList_l/mean/sd: ",len(deltList_temp),mean0,sd0)

    ##group the samples (ListA & ListB) according to whether score < mean0+?*sd0
    listA, listB = [], []
    threshold = mean0 + sd_cut*sd0
    listA = [score for score in deltList_temp if score < threshold]
    listB = [score for score in deltList_temp if score >= threshold]
    if len(listB)<5:
        print("Warning: sd_cut could be too large: ",sd_cut)
        exit(1)
    Plot_delt_distribution(deltList_temp,listA,listB,mean0,sd0,Out)
    #val_max = np.nanmax(subarray)
    #plot_heatmap_score(Outpath+"Score_native_map_delt",subarray,L,val_max,-2,1)
    score_matrix = Score_matrix(deltList,posList,L,listA,listB,Out)
    return score_matrix

## Min-Max normalize each row, respectively
def row_normalize(matrix_data):
    newdata = np.zeros(np.shape(matrix_data)) 
    Zmax,Zmin=matrix_data.max(axis=1),matrix_data.min(axis=1)
    row,col = matrix_data.shape[0],matrix_data.shape[1]
    for i in range(row):
        for j in range(col):
            newdata[i][j] = (matrix_data[i][j]-Zmin[i])/(Zmax[i]-Zmin[i]+0.000001)
    return newdata
def global_normalize(matrix_data):
    # Flatten the matrix to find global min and max
    global_max = np.nanmax(matrix_data)
    global_min = np.nanmin(matrix_data)
    # Apply normalization
    normalized_data = (matrix_data - global_min) / (global_max - global_min + 1e-6)  # Added small epsilon to avoid division by zero
    return normalized_data
## Processes the score matrix by replacing certain values, and normalizing.
def plot_matrix_boxplot(matrix,name):
    matrix = np.array(matrix)
    data = matrix.flatten()
    df = pd.DataFrame(data, columns=['Values'])
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='Values', color='lightblue')
    plt.xlabel('Values')
    plt.savefig(name +'.png', dpi=300) 
def Final_score(score_matrix, L, Out):  
    #plot_matrix_boxplot(score_matrix,Out+"score_distribution")
    val_max = np.nanmax(score_matrix)
    val_min = np.nanmin(score_matrix)
    #val = np.nanmean(score_matrix)
    val = np.nanmedian(score_matrix)
    # Replace values <= val or NaN with val_min
    score_matrix = np.where((score_matrix <= val) | (np.isnan(score_matrix)), val, score_matrix)
    # Normalize matrix
    #final_score = row_normalize(score_matrix)
    final_score = global_normalize(score_matrix)
    with open(Out + 'final_score.txt', 'w') as f:
        for row in final_score:
            score_list = [f"{score:.2f}" for score in row]
            f.write('\t'.join(score_list) + '\n')
    plot_heatmap_score(Out + "Score_map_normalized",final_score,L,1.0,0.0)
    return final_score
###Draw a dashed line on the current plot with the given slope and intercept###
def abline(slope, intercept):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, ls='--', c='black')
##Plot a heatmap of the score matrix with an optional mask###
def plot_heatmap_score(name,mat,L,val_max,val_min=-6.0,ismask=False):
    """
    Parameters:
    - name (str): The name of the output file.
    - mat (np.ndarray): The score matrix to plot.
    - L (int): The length of the matrix.
    - val_max (float): The maximum value for the heatmap color scale.
    - val_min (float): The minimum value for the heatmap color scale (default is -5.0).
    - ismask (bool): Whether to mask the lower triangle (default is False).
    """
    mask = np.zeros((L, L))
    if ismask:
        for i in range(1,L):
            mask[i, i] = True
            for j in range(i):
                if mat[i, j] != val_max:                
                    mask[i, j] = True
    #mask[np.tril_indices_from(mask)] = True  #set the below as 1
    plt.figure(figsize=(7, 5))
    #sns.set(font_scale=1.5)
    ax = sns.heatmap(mat,mask=mask,cmap="YlGn",vmin=val_min,vmax=val_max)   
    abline(1,0) 
    plt.xlim(0,L)
    plt.ylim(L,0)
    plt.xticks(rotation=90)
    #plt.yticks(range(0,L,5))
    # Add visible spines
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    plt.savefig(name +'.png', dpi=300)   
## Read the native base-base interaction derived from DSSR
## and evaluate the results
def native_bp(bp_file, L, val=1):
    with open(bp_file, "r") as native:
        lines = native.readlines()
    bplist = [line.split() for line in lines if int(line.split()[0]) < int(line.split()[2]) and int(line.split()[2]) > 0]
    pairs = np.zeros((L, L))
    for row in bplist:
        j, i = int(row[0]) - 1, int(row[2]) - 1
        pairs[i, j] = val
        pairs[j, i] = val
    return pairs
def Integrate_map(map1, map2, L):  ## map1: native
    final_map = np.zeros((L, L))
    # Copy the lower triangle from map1 to final_map
    final_map[np.tril_indices(L, -1)] = map1[np.tril_indices(L, -1)]
    # Copy the upper triangle from map2 to final_map
    final_map[np.triu_indices(L, 1)] = map2[np.triu_indices(L, 1)]
    return final_map
## Evaluation the results
def Evaluation_score(native_file,score_matrix,L,Out):
    ##put the native interaction into the score_matrix and output
    val_max, val_min = np.nanmax(score_matrix), np.nanmin(score_matrix)
    pair = native_bp(native_file,L,val_max)
    ## Evaluation
    rightuptri = [(i, j) for i in range(L) for j in range(i+1, L)]
    thlst = np.arange(val_min, val_max, 0.01).tolist()
    metrics = {'TPR': [], 'FPR': [], 'precision': [], 'F1': [], 'MCC': []}
    for threshold in thlst:
        TP = sum(score_matrix[i, j] >= threshold and pair[j, i] == val_max for i, j in rightuptri)
        FN = sum(score_matrix[i, j] <= threshold and pair[j, i] == val_max for i, j in rightuptri)
        FP = sum(score_matrix[i, j] >= threshold and pair[j, i] != val_max for i, j in rightuptri)
        TN = sum(score_matrix[i, j] <= threshold and pair[j, i] != val_max for i, j in rightuptri)

        TPR = TP / (TP + FN) if TP + FN > 0 else 0
        FPR = FP / (FP + TN) if FP + TN > 0 else 0
        precision = TP / (TP + FP) if TP + FP > 0 else 0
        F1 = 2 * TP / (2 * TP + FN + FP) if 2 * TP + FN + FP > 0 else 0
        MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 0
        metrics['TPR'].append(TPR)
        metrics['FPR'].append(FPR)
        metrics['precision'].append(precision)
        metrics['F1'].append(F1)
        metrics['MCC'].append(MCC)

    evadf = pd.DataFrame({'threshold': thlst, **metrics})
    #print(evadf)
    Max_F1 = max(evadf['F1'])
    Max_MCC = max(evadf['MCC'])
    print("Max_F1:",Max_F1,"MAX_MCC",Max_MCC)
    ########plot the results
    fig = plt.figure(figsize = (15,5))
    plt.subplot(131)
    plt.plot(evadf['FPR'], evadf['TPR'])
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.subplot(132)
    plt.plot(evadf['TPR'], evadf['precision'])
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.subplot(133)
    plt.plot(evadf['threshold'], evadf['F1'])
    plt.xlim(val_min, val_max)
    plt.xlabel('Threshold')
    plt.ylabel('F1 score')
    plt.tight_layout()
    plt.savefig(Out + 'Evaluate_score.jpg',dpi=600)
    #plt.show()
    score_map = Integrate_map(pair,score_matrix,L)
    plot_heatmap_score(Out+"Score_native_map",score_map,L,val_max, val_min,True)
    #plot_heatmap_score(Out+"Score_native_map",score_matrix,L,val_max, val_min)
    return Max_F1,Max_MCC 


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=str, help='The .ra data file (e.g., data/5TPY.var.ra)')
    parser.add_argument('sequence', type=str, help='The fasta sequence file (e.g., data/5TPY.fasta)')
    parser.add_argument('-o', '--out_path', type=str, default='./', help='The output path ')
    parser.add_argument('-n', '--native_secondary', type=str, default=None, 
        help='The path to the native secondary structure file (optional,e.g., 5TPY-native.txt) for result evaluation (e.g., F1 score  or score map)')
    parser.add_argument('-C', '--C_value', type=float, default=0.01, 
        help='Hyper-parameter C of svr_rbf')
    parser.add_argument('-g', '--gamma_value', type=float, default=0.1, 
        help='Hyper-parameter gamma of svr_rbf')
    parser.add_argument('-s', '--sd_cut', type=float, default=3.7, 
        help='Hyper-parameter sd_cut:delt_fitness > mean+sd_cut*sd --> possible pairing ')

    args = parser.parse_args()

    mut_file = args.data
    sequence_file = args.sequence
    Outpath = args.out_path + 'score/'
    # Define hyperparameters
    C_value, gamma_value, sd_cut = args.C_value, args.gamma_value, args.sd_cut
    #lower_output_cut = 5.0        #delete the sample with output/input<lower_cut
    # Create output directory if it doesn't exist
    os.makedirs(Outpath, exist_ok=True)

    sequence, seqLen = read_seq(sequence_file)  ## Read sequence
    #Cwt,single,double,triple = Read_mut_dataframe(mut_file)  ## Read .ra file
    Cwt,single,double = Read_mut_dataframe(mut_file)
    score_matrix = Contact_prediction(double,seqLen,Outpath) ## Calculate score
    final_score = Final_score(score_matrix,seqLen,Outpath)
    if args.native_secondary:
        native_file = args.native_secondary
        Max_F1,Max_MCC = Evaluation_score(native_file,score_matrix,seqLen,Outpath)
    run_time = time.time() - t0
    print("The run time: ",run_time,"s")
