from MolecularAnalysis import mol, moldb
from MolecularAnalysis.analysis import plot
from MolecularAnalysis.fragments import fragmentation
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


def create_df_gscore_vs_tanimoto(files_dir, specific_set, virus='global', target='glide'):
    """It creates a dataframe of the gscore vs. tanimoto for each molecule in the specific set."""

    specific = open(specific_set, 'r')
    specific_smiles = [line for line in specific]
    specific_mols = [moldb.Mol(smile=smile) for smile in specific_smiles]

    out = open('%s/%s_df_gscore_tanimoto.csv'%(files_dir,virus), 'w')
    out.write('id,SMILE,gscore,max_tan,outer\n')

    files = glob.glob('%s/glide?/docking/%s_%s_best.csv'%(files_dir,virus,target))
    files.sort()
    files.append('%s/glide10/docking/%s_%s_best.csv'%(files_dir,virus,target))
    files = files[1:] # not glide0

    for i, csv_file in enumerate(files):
        outer_round = i+1
        df = pd.read_csv(csv_file)
        for index, row in df.iterrows():
            ids = row['title']
            ids = ids.replace(',', '-')
            smiles = row['SMILES']
            gscores = row['r_i_glide_gscore']
            mol1 = moldb.Mol(smile=smiles)
            maxtanimoto=0
            for mol2 in specific_mols:             
                try:
                    similarity = moldb.getMolSimilarity(mol1,mol2,alg='Morgan4')
                    if similarity > maxtanimoto: maxtanimoto=similarity
                except:
                    print('similarity computation failed')
            #print(ids, smiles, gscores, maxtanimoto, outer_round)
            out.write('%s,%s,%s,%.4f,outer_%s\n'%(ids,smiles,gscores,maxtanimoto,outer_round))
    out.close()
    
def filter_compounds(df, gscore_thresh, tan_thresh):
    """Helper function to filter compounds based on gscore and tanimoto thresholds."""
    filt = df[(df['gscore'] <= gscore_thresh) & (df['max_tan'] <= tan_thresh)]
    ids = filt['id'].tolist()
    rounds = [x.split('_')[-1] for x in filt['outer'].tolist()]
    return set(f'{round_}_{id_}' for round_, id_ in zip(rounds, ids))

def apply_thresholds(global_csv, individual_csvs, gscore_ind, gscore_glob, tan_ind, tan_glob):
    """Get compounds below thresholds of glide gscore and tanimoto similarity.
    
    Parameters:
         - global_csv: Path to the global CSV file.
         - individual_csvs: Dictionary of individual CSV files with keys 'SARS2', 'SARS', 'MERS'
         - gscore_ind: Individual glide gscore (kcal/mol) threshold
         - gscore_glob: Global glide_gscore (kcal/mol) threshold
         - tan_ind: Individual tanimoto similarity (against the specific set) threshold
         - tan_glob: Global tanimoto similarity (against the specific set) threshold
 
     NOTE: I calculated the intersection with the compound name, not the SMILE.
     SMILES can differ between different viruses and it gives us wrong result."""
    
    # First, get the global compounds
    glob_df = pd.read_csv(global_csv)
    glob_filt = glob_df[(glob_df['gscore'] <= gscore_glob) & (glob_df['max_tan'] <= tan_glob)]
    glob_ids = glob_filt['id'].tolist()
    glob_rounds = [x.split('_')[-1] for x in glob_filt['outer'].tolist()]
    glob_compounds = set(f'{round_}_{id_}' for round_, id_ in zip(glob_rounds, glob_ids))
    
    # Then, get the individual compounds
    all_ind_comp = []
    for csv in individual_csvs:
        ind_df = pd.read_csv(csv)

        # Filter the individual compounds based on the thresholds
        ind_compounds = filter_compounds(ind_df, gscore_ind, tan_ind)
        all_ind_comp.append(ind_compounds)
    ind_compounds = set.intersection(*all_ind_comp)
    
    # Get the intersection of the global and individual compounds
    final_compounds = glob_compounds.intersection(ind_compounds)
    print(f'Number of compounds below thresholds: {len(final_compounds)}')
    return final_compounds

def plot_gscore_vs_tanimoto(csv_file, outdir, outname, gscore_threshold=-8, tan_threshold=0.25, save_csv=True, plot_max_lim=[-6, 0.5]):
    """It plots a scatter plot of Glide gscore (x) vs. max tanimoto (y)
       against specific set. It also represent the threshold lines.
       It can also store the filtered df by the specified thresholds.
    """
    df = pd.read_csv(csv_file)

    # store filtered df
    filt_df = df[(df['gscore']<=gscore_threshold) & (df['max_tan']<=tan_threshold)]
    print('Filtered compounds:\n', filt_df)
    if save_csv:
        folder = os.path.dirname(csv_file)
        filt_df.to_csv('%s/df_gscore_%s_tan_%s.csv'%(folder, gscore_threshold, tan_threshold), header=['id','SMILES','gscore', 'max_tan','round'],index=False)

    # Apply max limits
    df = df[(df['gscore']<=plot_max_lim[0]) & (df['max_tan']<=plot_max_lim[1])]
    
    # Plot scatter plot
    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 8))
    
    df['outer'] = df['outer'].apply(lambda x: x.split('_')[-1])
    df['outer'] = df['outer'].astype(int)
    gdf = df.groupby('outer')
    dfs = [gdf.get_group(x) for x in gdf.groups]
    rounds = df['outer'].unique()
    total = len(dfs)
    palette = mcp.gen_color(cmap="YlGnBu",n=total)
    palette = list(reversed(palette))
    for i, group in enumerate(dfs):
        plt.scatter(x=group['gscore'].to_list(),y=group['max_tan'].to_list(),color=palette[i],label=rounds[i], alpha=0.8, s=14,marker='X')
    plt.axvline(x=gscore_threshold,color='red',ls='--',lw='1')
    plt.axhline(y=tan_threshold,color='red',ls='--',lw='1')
    plt.xlabel('Glide gscore',fontsize=16)
    plt.ylabel('Max tanimoto against Specific Set',fontsize=16)
    plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.savefig('%s/%s.pdf'%(outdir, outname))
    
def map_gscores_generated(csv_global, csv_virus, outdir, outname):
    
    global_df = pd.read_csv(csv_global)
    global_df = global_df.rename(columns={'gscore':'global_gscore'})
    global_df['id'] = global_df['outer'] + '_' + global_df['id']
    global_df['id'] = global_df['id'].str.replace('outer_', '')
    global_df = global_df[['id', 'max_tan', 'global_gscore']]
    
    for csv in csv_virus:
        virus = csv.split('/')[-1].split('_')[0]
        virus_df = pd.read_csv(csv)
        virus_df = virus_df.rename(columns={'gscore':f'gscore_{virus}', 'SMILE':f'SMILE_{virus}', 'max_tan':f'max_tan_{virus}'})
        virus_df['id'] = virus_df['outer'] + '_' + virus_df['id']
        virus_df['id'] = virus_df['id'].str.replace('outer_', '')
        virus_df = virus_df[['id', f'SMILE_{virus}', f'gscore_{virus}']]
        
        global_df = pd.merge(global_df, virus_df, on=['id'])
    
    global_df.to_csv('%s/%s.csv'%(outdir, outname), index=False)
        

def cluster_DBSCAN(csv_results, gscore_glob_thr, gscore_ind_thr, tanimoto_thr):
    
    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= gscore_glob_thr) & (df['gscore_SARS2'] <= gscore_ind_thr) & (df['gscore_SARS'] <= gscore_ind_thr) & (df['gscore_MERS'] <= gscore_ind_thr) & (df['max_tan'] <= tanimoto_thr)]
    df_filt = df_filt[df_filt['SMILE_SARS2'] != 'invalid_structure']
    # TODO: CHECK WHAT HAPPENED TO INVALID STRUCTURES
    
    df_filt['outer'] = df_filt['id'].apply(lambda x: x.split('_')[0])
    df_filt['outer'] = df_filt['outer'].astype(int)
    total_outers = len(df_filt['outer'].unique())
    
    similarity_thr = 0.2
    accum_scaffolds = []
    for outer in range(total_outers):
            
        df_outer = df_filt[df_filt['outer']==outer+1]
        smiles = df_outer['SMILE_SARS2'].tolist()
        ids = df_outer['id'].tolist()
        mol_list = [mol.Mol(smile=smile, name=ids[i]) for i, smile in enumerate(smiles)]
        mol_db = moldb.MolDB(molList=mol_list)
        
        scaffolds = []
        for i, molecule in enumerate(mol_db.mols):
            scaffold = fragmentation.getScaffold(molecule)
            scaffold = mol.Mol(rdkit=scaffold, name=molecule.name)
            scaffolds.append(scaffold)
        
        scaff = moldb.MolDB(molList=scaffolds)
        scaff_names = [scaffold.name for scaffold in scaff.mols]
        scaff_mols = scaff.mols
        print(scaff.size)
        
        scaff._get_similarity_matrix()
        np.fill_diagonal(scaff.simmatrix.values, 0)
        print(scaff.simmatrix)
        
        below_threshold = scaff.simmatrix.lt(similarity_thr).all(axis=1)
        print(below_threshold)
        
        row_names = scaff.simmatrix.index[below_threshold].tolist()
        print(row_names)
        print(len(row_names))
        
        # get index of the row names in the scaff_names list
        row_indexes = [scaff_names.index(name) for name in row_names]
        print(row_indexes)
        # get the scaff_mols with the indexes
        scaffs = [scaff_mols[i] for i in row_indexes]
        
        exit()
        accum_scaffolds += scaff.mols
    
    exit()
    smiles = df_filt['SMILE_SARS2'].tolist()
    mol_list = [mol.Mol(smile=smile) for smile in smiles]
    mol_db = moldb.MolDB(molList=mol_list)
    
    scaffolds = []
    for i, molecule in enumerate(mol_db.mols):
        scaffold = fragmentation.getScaffold(molecule)
        scaffold = mol.Mol(rdkit=scaffold, name=molecule.name)
        scaffolds.append(scaffold)
        
    scaff = moldb.MolDB(molList=scaffolds)
    print(scaff.size)
    
    
#cluster_DBSCAN(csv_results='/home/cactus/julia/gensim/selective/results.csv', gscore_glob_thr=-8, gscore_ind_thr=-7.5, tanimoto_thr=0.3)
    