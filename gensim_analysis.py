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

def get_all_generated_molecules(results_dir, outname):
    """ It joins all generated molecules into one single file."""

    smiles_list = []
    results = glob.glob('%s/%s_*/%s_generated_smiles_threshold_*.csv'%(results_dir, outname, outname))
    df_list = []
    for result in results:
        nfile = os.path.basename(result)
        inner = nfile.replace('.csv', '').split('_')[-1]
        df = pd.read_csv(result)
        df['inner'] = int(inner)
        df_list.append(df)
    all_generated = pd.concat(df_list)
    all_generated = all_generated.sort_values('inner')
    all_generated.to_csv('%s/all_generated_molecules.csv'%results_dir, index=False)
    print('Total generate molecules: %s'%len(all_generated))
    

def create_table_gm_counts(results_dir, outname, save_df=False):
    """It creates a table of the gm results for each inner round.
       The table includes the molecule counts for each file."""

    # get all counts for each file
    counts_specific = {}
    counts_generated = {}
    counts_threshold = {}
    counts_success = {}
    results = glob.glob('%s/%s_*'%(results_dir, outname))
    for result in results:
        ninner = os.path.basename(result)
        if outname not in ninner: continue
        inner = ninner.split('_')[-1].strip()
        inner = int(inner)

        innerfiles = glob.glob('%s/*'%result)
        for tfile in innerfiles:
            nfile = os.path.basename(tfile)
            if '.png' in nfile: continue
            if outname not in nfile: continue

            if nfile == '%s_specific_smiles.csv'%outname: # the first one
                if ninner != '%s_1'%outname: continue
                initial_specific = pd.read_csv(tfile)
                num_in_specific = len(initial_specific)
                counts_specific[1] = num_in_specific
            if nfile == '%s_generated_smiles_merged.csv'%outname: # n
                if ninner == '%s_%s'%(outname, len(results)): continue
                specific = pd.read_csv(tfile)
                num_specific = len(specific)
                counts_specific[inner] = num_specific
            if nfile == '%s_generated_smiles.csv'%outname: # n
                generated = pd.read_csv(tfile)
                num_generated = len(generated)
                counts_generated[inner] = num_generated

            if nfile == '%s_generated_smiles_threshold_%s.csv'%(outname, (inner-1)): # n - 1
                threshold = pd.read_csv(tfile)
                num_threshold = len(threshold)
                counts_threshold[inner - 1] = num_threshold

            if nfile == 'report_%s.out'%outname: # n
                report = open(tfile, 'r')
                for line in report:
                    if line.startswith('Success'):
                        rate = line.split(': ')[-1].replace('%', '')
                        rate = round(float(rate), 2)
                        counts_success[inner] = rate

    # save it into a dataframe
    try:
        counts_dic = dict([(k,[counts_specific[k],counts_generated[k], counts_threshold[k], counts_success[k]]) for k in counts_specific])
    except:
        counts_specific = dict(sorted(counts_specific.items()))
        counts_specific.popitem()
        counts_dic = dict([(k,[counts_specific[k],counts_generated[k], counts_threshold[k], counts_success[k]]) for k in counts_specific])
    counts_dic = dict(sorted(counts_dic.items()))
    index_labels = ['specific', 'generated', 'thresholds', 'success_rate']
    df = pd.DataFrame(counts_dic, index=index_labels)
    print(df)

    if save_df:
        df.to_csv('%s/table_of_counts.csv'%results_dir)
        df_trans = df.transpose()
        df_trans.to_csv('%s/table_of_counts_trans.csv'%results_dir)

def plot_results_properties(results_dir, which_property, set_color='blue', title=True, save_plot=False):
    """It plots the trials vs. the selected propery
       of the GM results (From the previous dataframe!). """

    props_table = pd.read_csv('%s/table_of_counts.csv'%results_dir, index_col=0)
    inners = props_table.columns
    props_table.loc['inners'] = inners
    props_table = props_table.T
    plt.plot('inners', which_property, data=props_table, marker='o', alpha=0.4, color=set_color)
    if title:
        plt.title('Inners vs. %s'%which_property)
    plt.xlabel('Inners')
    if which_property == 'success_rate':
        plt.ylabel('%s '%which_property + '(%)')
    else:
        plt.ylabel('%s (counts)'%which_property)

    if save_plot:
        plt.savefig('%s/plot_counts_%s.png'%(results_dir, which_property))


def plot_all_props_in_one(results_dir, save_fig=True):
    """ It plots the previous 4 properties plots into one.
        It is a 2x2 plot."""

    fig, ax = plt.subplots(2,2, figsize=(12, 8))
    plt.subplot(221)
    plot_results_properties(results_dir, 'specific', set_color='blue')
    plt.subplot(222)
    plot_results_properties(results_dir, 'generated', set_color='orange')
    plt.subplot(223)
    plot_results_properties(results_dir, 'thresholds', set_color='green')
    plt.subplot(224)
    plot_results_properties(results_dir, 'success_rate', set_color='purple')
    plt.suptitle('Plots of properties for each generation inner loop')
    fig.tight_layout(pad=1.0)
    if save_fig:
        plt.savefig('%s/plot_all_counts.png'%results_dir)
        
def convert_csv_to_sdf_file(csv_to_convert, outdir, inner=False):
    """ It convert csv to smi file to be uploaded to the MolDB class.
    We convert them to sdf instead of smi to keep all the properties."""

    filename = os.path.basename(csv_to_convert).replace('.csv','')
    df = pd.read_csv(csv_to_convert)
    with Chem.SDWriter('%s/%s.sdf'%(outdir,filename)) as w:
        for index, row in df.iterrows():
            m = Chem.MolFromSmiles(row['smiles'])
            if inner:
                m.SetProp('_Name', str(inner)+'_'+str(row['Unnamed: 0']))
                m.SetProp('inner', str(inner))
            else:
                m.SetProp('_Name', str(row['inner'])+'_'+str(row['Unnamed: 0']))
                m.SetProp('inner', str(row['inner']))
                
            m.SetProp('SAacore', str(row['SAscore']))
            m.SetProp('QED', str(row['QED']))
            m.SetProp('mw', str(row['mw']))
            m.SetProp('logp', str(row['logp']))
            m.SetProp('tpsa', str(row['tpsa']))
            m.SetProp('nHD', str(row['nHD']))
            m.SetProp('nHA', str(row['nHA']))
            if 'tan_mean' in df:
                m.SetProp('tan_mean', str(row['tan_mean']))
            if 'tan_max' in df:
                m.SetProp('tan_max', str(row['tan_max']))
            w.write(m)
    
    print('csv converted to sdf successfully')

def _simplify_specific_sets(results_dir, outname):
    """ It extracts the repeated smiles in each specific set.
        As a result, it creates sdf files of all the specific sets
        containing only the new generated molecules."""
    
    def drop_duplicates_2_dfs(df1, df2):
        df1 = pd.read_csv(df1)
        df2 = pd.read_csv(df2)
        df = pd.concat([df1, df2])
        return df
    
    spec_files = glob.glob('%s/%s_*/%s_specific_smiles.csv'%(results_dir,outname, outname))
    total = len(spec_files)
    for i in range(total):
        i = i+1
        try:
            df_1 = '%s/%s_%s/%s_specific_smiles.csv'%(results_dir,outname, i, outname)
            df_2 = '%s/%s_%s/%s_specific_smiles.csv'%(results_dir,outname, i+1, outname)
            df = drop_duplicates_2_dfs(df_1, df_2)
            df.to_csv('%s/%s_%s/%s_specific_smiles_simple.csv'%(results_dir,outname, i, outname), index=False)
        except:
            pass # because the last inner does not hace specific set file.

def plot_modbs_tSNE_or_UMAP(list_of_sdfs, list_of_names, outdir, outname, sizes, alphas, markers, ptype='UMAP'):
    """ It plots the tSNE/UMAP of all the sets indicated in the list."""

    total = len(list_of_sdfs)
    moldb_list = []
    for file in list_of_sdfs:
        mols = moldb.MolDB(sdfDB=file, verbose=False)
        moldb_list.append(mols)
    
    plot_colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    plot_colors = plot_colors[1:total+1]

    min_dists = [0.2]
    neighbours = [100]

    if ptype == 'tSNE':
        plot.plotTSNE(dbs = moldb_list, names = list_of_names, output='%s/%s'%(results_dir, outname),\
                    random_max = 10000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors,\
                    sizes=sizes, alphas=alphas, linewidth=0, markers=markers, n_iter=1000, perplexity=30,\
                    early_exaggeration=12,learning_rate='auto')

    if ptype == 'UMAP':
        for i in range(len(min_dists)):
            for j in range(len(neighbours)):
                plot.plotUMAP(dbs = moldb_list, names = list_of_names, output='%s/%s_md%s_nn%s'%(outdir, outname, min_dists[i], neighbours[j]),\
                            random_max = 1000, delimiter = None, alg = 'Morgan4', colors = plot_colors, sizes = sizes,  alphas = alphas,\
                            min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 10000, markers = markers, figsize = (9,6), \
                            linewidth = 0)

        
if __name__ == "__main__":
    
    n = 'gensim_mt'
    
    #get_all_generated_molecules(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n)
    
    #create_table_gm_counts(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n, save_df=True)
    
    #plot_all_props_in_one(results_dir='/home/cactus/julia/gensim/full/outer1', save_fig=True)
    
    #convert_csv_to_sdf_file(csv_to_convert='/home/cactus/julia/gensim/full/outer1/all_generated_molecules.csv', outdir='/home/cactus/julia/gensim/full/outer1')
    
    #_simplify_specific_sets(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n)
    # spec_simples = glob.glob('/home/cactus/julia/gensim/full/outer1/gensim_mt_*/gensim_mt_specific_smiles_simple.csv')
    # spec_inner = [x.split('/')[-2].split('_')[-1] for x in spec_simples]
    # for i, spec_simple in enumerate(spec_simples):
    #     convert_csv_to_sdf_file(csv_to_convert=spec_simple, outdir='/home/cactus/julia/gensim/full/outer1/gensim_mt_%s'%spec_inner[i], inner=spec_inner[i])

    spec_simples = glob.glob('/home/cactus/julia/gensim/full/outer1/gensim_mt_*/gensim_mt_specific_smiles_simple.sdf')
    sdf_list = []
    rev_inners = list(range(1, 20))
    rev_inners.sort(reverse=True)
    for i in rev_inners:
        sdf_list.append('/home/cactus/julia/gensim/full/outer1/gensim_mt_%s/gensim_mt_specific_smiles_simple.sdf'%i)
    names = ['specific_19', 'specific_18', 'specific_17', 'specific_16', 'specific_15',\
            'specific_14', 'specific_13', 'specific_12', 'specific_11', 'specific_10',\
            'specific_9', 'specific_8', 'specific_7', 'specific_6', 'specific_5',\
            'specific_4', 'specific_3', 'specific_2', 'specific_1', 'specific_0']
    sizes = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]
    plot_modbs_tSNE_or_UMAP(list_of_sdfs=spec_simples, list_of_names=names, outdir='/home/cactus/julia/gensim/full/outer1/', outname='UMAP_test',\
                            sizes=sizes, alphas=alphas, markers=markers, ptype='UMAP')
