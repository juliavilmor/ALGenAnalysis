from MolecularAnalysis import mol, moldb
from MolecularAnalysis.analysis import plot
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import numpy as np
import seaborn as sns
from rdkit import Chem
from collections import Counter
from statistics import mean

def get_best_glide_docking_pose(csv_file, to_csv=True):
    """ I returns a df with the best poses for each ligand,
        based on the glide gscore """
 
    df = pd.read_csv(csv_file)
    df = df[df['SMILES'] != 'invalid_structure']
    idx = df.groupby('title')['r_i_glide_gscore'].idxmin()
    df_best = df.loc[idx]
    print('The best pose for each compound has been obtained. Total: %s compounds'%len(df_best))
    
    if to_csv:
        df_best.to_csv(csv_file.replace('.csv', '_best.csv'), index=False)
 
def plot_glide_scores(path_to_best_csv, insert_title, outdir, savefig=True):
    """Histogram of gscores."""

    df = pd.read_csv(path_to_best_csv)
    sns.histplot(data=df, x='r_i_glide_gscore').set(title=insert_title)
    filename = os.path.basename(path_to_best_csv).replace('.csv', '_gscores.png')
    if savefig:
        plt.savefig('%s/%s'%(outdir, filename))
        
def superimpose_histograms(list_of_csvs, list_of_labels, insert_title, out, savefig=True, legend_loc='upper right', xlim=None, ylim=None):
    """Superimpose histograms to compare."""

    plt.figure(figsize=(6,6))
    total = len(list_of_csvs)
    colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    colors = colors[2:total+1]
    colors = ['red'] + colors
    for i, data in enumerate(list_of_csvs):
        df = pd.read_csv(data)
        df = df[df['r_i_glide_gscore'] != 10000]
        average = mean(df['r_i_glide_gscore'].tolist())
        sns.histplot(data=df, x=df['r_i_glide_gscore'].astype(float), color=colors[i], label=list_of_labels[i], element='step', fill=False, bins=70)
        plt.axvline(average, color=colors[i], linestyle=":")
    plt.legend(loc=legend_loc, frameon=False)
    plt.title(insert_title)
    plt.xlabel("Glide docking score (kcal/mol)")
    plt.ylabel("Counts")
    if xlim:
        plt.xlim(xlim[0], xlim[1])
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    if savefig:
        plt.savefig(out)
        
def get_mean_glide_gscores(list_of_csv, out):
    "The mean of glide gscores"

    csv_list = []
    for csv_file in list_of_csv:
        df = pd.read_csv(csv_file)
        csv_list.append(df)
        
    ligs = df['title'].tolist()
    ligs = set(ligs)
    means = []
    smiles = []
    for lig in ligs:
        mean_gscore = df.loc[df['title'] == lig]['r_i_glide_gscore'].mean()
        smile = df.loc[df['title'] == lig]['SMILES'].tolist()
        smile = smile[0] # some smiles differ from quirality. I selected the first one (corresponding to SARS2)
        means.append(mean_gscore)
        smiles.append(smile)
    mean_df = pd.DataFrame(list(zip(ligs, smiles, means)), columns =['title',  'SMILES', 'r_i_glide_gscore'])
    print(mean_df)
    mean_df.to_csv('%s'%out, index=False)

def filter_by_glide_gscore(path_to_best_csv, gscore, outdir):
    """ It filters the obtained molecules by glide score.
        As a result, it gives you a csv file and a smi file."""

    df = pd.read_csv(path_to_best_csv)
    df_filt = df[df['r_i_glide_gscore'] < gscore]
    print("From %s molecules, %s were removed.\nThe new set contains %s molecules."%(len(df), (len(df)-len(df_filt)), len(df_filt)))

    # save filtered dataframe
    file_name = path_to_best_csv.replace('_best.csv', '_gscore_filtered.csv')
    df_filt.to_csv(file_name)

    # save smiles for the next round
    smiles = df_filt['SMILES'].tolist()
    smi_file = open('%s/genai_specific_set.smi'%outdir, 'w')
    for smile in smiles:
        smi_file.write('%s\n'%smile)
    smi_file.close()
    
def filter_by_glide_gscore_paninhibitors(list_of_csvs, outdir, gscore_global=-6.5, gscore_individual=-5.9):
    list_dfs = []
    for glide in list_of_csvs:
        virus = os.path.basename(glide).split('_')[0]
        receptor = os.path.basename(glide).split('_')[2]
        df = pd.read_csv(glide)
        df = df[df['r_i_glide_gscore'] != 10000]
        df = df.drop_duplicates(subset='title')
        df['virus'] = virus
        df['receptor'] = receptor
        list_dfs.append(df)
        gscores = df['r_i_glide_gscore'].tolist()
        #plt.hist(gscores, bins=15, alpha=0.5)
        #plt.axvline(x = -5.9, color = 'r')
        #plt.savefig('hist_gscores.png')
    all_df = pd.concat(list_dfs)
    ligs = all_df['title'].tolist()
    ligs = set(ligs)
    mean = {}
    for lig in ligs:
        # global threshold
        mean_gscore = all_df.loc[all_df['title'] == lig]['r_i_glide_gscore'].mean()
        if mean_gscore < gscore_global:

            # individual threshold
            ind_gscore = all_df.loc[all_df['title'] == lig]['r_i_glide_gscore'] < gscore_individual
            below = ind_gscore.tolist()
            if below == [True, True, True]:
                mean[lig] = mean_gscore
                smiles = all_df.loc[all_df['title'] == lig]['SMILES'].tolist()[0]
                smi_file = open('%s/specific_set.smi'%outdir, 'a')
                smi_file.write('%s\n'%smiles)
                smi_file.close()
    print('From %s molecules, %s were removed.\nThe new set contains %s molecules'%(len(ligs), (len(ligs) - len(mean)), len(mean)))

def cumulative_histograms(final_csvs, initial_csvs, list_of_labels, list_of_colors, insert_title, out, savefig=True, legend_loc='upper right', xlim=None, ylim=None):
    
    plt.figure(figsize=(6,6))
    fig, ax = plt.subplots()
    list_of_csvs = final_csvs + initial_csvs
    total = len(list_of_csvs)
    colors = list_of_colors
    for i, data in enumerate(list_of_csvs):
        df = pd.read_csv(data)
        df = df.rename(columns={'r_i_glide_gscore': 'gscore'})
        df = df[df['gscore'] != 10000]
        average = mean(df['gscore'].tolist())
        sns.histplot(data=df, x=df['gscore'].astype(float), color=colors[i], label=list_of_labels[i], element='step', fill=False, bins=70)
        plt.axvline(average, color=colors[i], linestyle=":")
    plt.legend(loc=legend_loc, frameon=False)
    plt.title(insert_title)
    plt.xlabel("Glide docking score (kcal/mol)")
    plt.ylabel("Counts")
    if xlim:
        plt.xlim(xlim[0], xlim[1])
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    if savefig:
        plt.savefig(out)
    

if __name__ == "__main__":
    
    # EXECUTIONS PAN-INHIBITOR
    
    # get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide1/docking/SARS2_7rnwA1_docking.csv', to_csv=True)
    # get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide1/docking/SARS_2gx4A1_docking.csv', to_csv=True)
    # get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide1/docking/MERS_7eneC1_docking.csv', to_csv=True)
    
    csvs = ['/home/cactus/julia/gensim/full/glide1/docking/SARS2_7rnwA1_docking_best.csv', '/home/cactus/julia/gensim/full/glide1/docking/SARS_2gx4A1_docking_best.csv', '/home/cactus/julia/gensim/full/glide1/docking/MERS_7eneC1_docking_best.csv']
    # labels = ['SARS2', 'SARS', 'MERS']
    # superimpose_histograms(list_of_csvs=csvs, list_of_labels=labels, insert_title='Gscores of outer 1', out='/home/cactus/julia/gensim/full/plots/hist_gscores_outer1.png', savefig=True, legend_loc='upper right')
    # superimpose_histograms(list_of_csvs=csvs, list_of_labels=labels, insert_title='Gscores of outer 1', out='/home/cactus/julia/gensim/full/plots/hist_gscores_outer1_zoom.png', savefig=True, legend_loc='upper right', xlim=[-8.5,-7])
    
    filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir='/home/cactus/julia/gensim/full/outer2', gscore_global=-7.5, gscore_individual=-7)
    