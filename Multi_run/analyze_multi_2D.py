########################
#Author：YZ Shi
#Date:   2023.2
#Loc.:   SZBL
#Description：Analyze the results from multiple MC including cluster the contact maps
#v1.0: cluster MC results based on simility
########################
import sys, os, argparse
from collections import defaultdict  
from numpy import exp
import numpy  as np       
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
t0 = time.time()
######### Read contact maps. #########
def read_contact(file):
    contact = np.load(file)
    return contact
def read_multi_contact(filepath, filename, N):
    multi_contact, list_1, list_2 = [], [], []
    for i in range(1, N + 1):
        file = f'{filepath}2D/{i}/{filename}'     
        if os.path.isfile(file):
            contact = read_contact(file)
            multi_contact.append(contact)
            list_1.append(i)
        else:
            list_2.append(i)
    if len(list_2) > 0:
        print("Warning: Some maps missing in some MC results")
    
    return multi_contact,list_1,list_2

####### Compute similarity between two contact maps.######
#####Similarity is defined as the fraction of identical entries. #########
def contact_map_similiarity(mat1,mat2,cut):  
    identical_entries = np.sum(mat1 == mat2)
    total_entries = mat1.size
    similarity = identical_entries / total_entries
    #print(identical_entries,total_entries,similarity)
    return similarity >= cut
def jaccard_similarity(mat1, mat2, cut):
    """
    Calculate the Jaccard similarity between two binary matrices.
    """
    intersection = np.logical_and(mat1, mat2).sum()
    union = np.logical_or(mat1, mat2).sum()
    if union == 0:
        return 1.0  # Both matrices are empty
    similarity = intersection / union
    #print(intersection,union,similarity)
    return similarity >= cut
########cluster the maps & find the representative maps#######
def cluster_maps(contact, N_map, filepath, list2):
    cluster = []
    flag = [0] * N_map
    for i in range(N_map-1):
        if flag[i]==1:
            continue
        clust = [i + 1]
        flag[i] = 1
        n1 = len(contact[i])
        for j in range(i+1,N_map):        
            n2 = len(contact[j])
            if n1!=n2:
                raise ValueError("The two matrices have different dimensions")
            #if np.all(contact[i]==contact[j]):  ## equal
            #if contact_map_similiarity(contact[i],contact[j],0.999)==1:
            # Compare contact maps using Jaccard similarity
            if jaccard_similarity(contact[i],contact[j],0.9):
                clust.append(j + 1)
                flag[j] = 1
        cluster.append(clust)

    cluster.sort(key = len,reverse=True)
    N_clust = len(cluster)
    topN = min(5, N_clust)  #int(N_clust*0.1)
    top_list = [cluster[i][0] for i in range(topN)]
    print(top_list)   
    with open(filepath + 'cluster_2D.txt', 'w') as f:
        f.write('Typical:\t' + str(top_list) + '\n')
        f.write('Missing:\t' + str(list2) + '\n')
        for i, clust in enumerate(cluster):
            f.write(f'Cluster_{i + 1}:\t{clust}\n')
def average_map(contact,filepath):
    L = contact[0].shape[0] 
    mean_matrix = np.zeros((L, L))
    for matrix in contact:
        mean_matrix += matrix
    mean_matrix /= len(contact)
    return mean_matrix
###add line to one figure
def abline(slope, intercept):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, ls='--', c='black')
##Score_Heatmap plot
def plot_heatmap(filepath,name,mat,val_max=1,val_min=0,native=False,ismask=1):
    L = len(mat)
    mask = np.zeros((L, L))
    #mask[np.tril_indices_from(mask)] = True  #set the below as 1
    if native:
        pair_native = native_bp(native,L)
        mat_combine = Integrate_map(pair_native, mat, L)
    else:
        mat_combine = mat + mat.T - np.diag(np.diag(mat))
    if ismask==1:
        mask = (mat_combine == 0)
    plt.figure(figsize=(7, 5))
    #sns.set(font_scale=1.5)
    ax = sns.heatmap(mat_combine,mask=mask,cmap="YlGn",vmin=val_min,vmax=val_max,square=True,xticklabels=5,yticklabels=5)   
    abline(1,0) 
    plt.xlim(0,L)
    plt.ylim(L,0)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    plt.savefig(os.path.join(filepath, f"{name}.png"), dpi=300)   
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='The MC data path (e.g., ./2D)')
    parser.add_argument('Nfile', type=int, default=0, help='The number of MC simulation')
    parser.add_argument('-n', '--native_secondary', type=str, default=None, 
        help='The path to the native secondary structure file (optional,e.g., 5TPY-native.txt) for result evaluation (e.g., F1 score  or score map)')
    args = parser.parse_args()
    Filepath = args.data_path
    N_file = args.Nfile
    file_name = 'contact_minE.npy'

    contact, list1, list2 = read_multi_contact(Filepath, file_name, N_file)
    N_map = len(contact)
    print(N_map,len(list1),len(list2))
    cluster_maps(contact, N_map, Filepath, list2)
    avg_map = average_map(contact, Filepath)
    if args.native_secondary:
        native_file = args.native_secondary
        plot_heatmap(Filepath,"mean_2D_map_native",avg_map,0,1,native=native_file)    
    else:
        plot_heatmap(Filepath,"mean_2D_map",avg_map,0,1)    
    run_time = time.time() - t0
    print("The run time: ",run_time,"s")



