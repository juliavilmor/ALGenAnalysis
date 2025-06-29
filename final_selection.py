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
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from collections import Counter
from statistics import mean
import itertools


def create_df_gscore_vs_tanimoto(files_dir, specific_set, virus='global', target='glide'):
    """It creates a dataframe of the gscore vs. tanimoto for each molecule in the specific set."""

    specific = open(specific_set, 'r')
    specific_smiles = [line for line in specific]
    specific_mols = [moldb.Mol(smile=smile) for smile in specific_smiles]
    
    out = open('%s/%s_df_gscore_tanimoto.csv'%(files_dir,virus), 'w')
    out.write('id,SMILES,gscore,max_tan,outer\n')

    files = glob.glob('%s/glide_?/docking/%s_%s_best.csv'%(files_dir,virus,target))
    files.sort()
    files2 = glob.glob('%s/glide_??/docking/%s_%s_best.csv'%(files_dir,virus,target))
    files2.sort()
    files = files + files2
    #files = files[1:] # not glide0

    for i, csv_file in enumerate(files):
        outer_round = i+1
        df = pd.read_csv(csv_file)
        for index, row in df.iterrows():
            ids = row['title']
            ids = ids.replace(',', '-')
            smiles = row['SMILES']
            gscores = row['r_i_docking_score']
            mol1 = moldb.Mol(smile=smiles)
            maxtanimoto=0
            for mol2 in specific_mols:
                try:
                    similarity = moldb.getMolSimilarity(mol1,mol2,alg='Morgan4')
                    if similarity > maxtanimoto: maxtanimoto=similarity
                except:
                    print('similarity computation failed')
            print(ids, smiles, gscores, maxtanimoto, outer_round)
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
        virus_df = virus_df.rename(columns={'gscore':f'gscore_{virus}', 'SMILES':f'SMILES_{virus}', 'max_tan':f'max_tan_{virus}'})
        virus_df['id'] = virus_df['outer'] + '_' + virus_df['id']
        virus_df['id'] = virus_df['id'].str.replace('outer_', '')
        virus_df = virus_df[['id', f'SMILES_{virus}', f'gscore_{virus}']]

        global_df = pd.merge(global_df, virus_df, on=['id'])

    global_df.to_csv('%s/%s.csv'%(outdir, outname), index=False)
    
def create_table_results(csv_file, outdir, outname):
    df = pd.read_csv(csv_file)
    df['id'] = df['outer'] + '_' + df['id']
    df['id'] = df['id'].str.replace('outer_', '')
    df = df[['id', 'SMILES', 'max_tan', 'gscore']]
    df.to_csv('%s/%s.csv'%(outdir, outname), index=False)

def get_scaffolds_db(molecular_db):
    scaffolds = []
    for i, molecule in enumerate(molecular_db.mols):
        scaffold = fragmentation.getScaffold(molecule)
        scaffold = mol.Mol(rdkit=scaffold, name=molecule.name)
        scaffolds.append(scaffold)
    scaffold_db = moldb.MolDB(molList=scaffolds)
    return scaffold_db

def cluster_DBSCAN(csv_results,
                   smi_specific,
                   gscore_glob_thr,
                   gscore_ind_thr,
                   tanimoto_thr,
                   similarity_thr):

    from sklearn.cluster import DBSCAN

    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= gscore_glob_thr) & (df['gscore_SARS2'] <= gscore_ind_thr) & (df['gscore_SARS'] <= gscore_ind_thr) & (df['gscore_MERS'] <= gscore_ind_thr) & (df['max_tan'] <= tanimoto_thr)]
    df_filt = df_filt[df_filt['SMILES_SARS2'] != 'invalid_structure']

    df_filt['outer'] = df_filt['id'].apply(lambda x: x.split('_')[0])
    df_filt['outer'] = df_filt['outer'].astype(int)
    total_outers = len(df_filt['outer'].unique())

    num_clusters = []
    x_values = []
    all_scaff_mols = []

    for outer in range(total_outers + 1):
        if outer == 0:
            specific_db = moldb.MolDB(smiDB=smi_specific)
            scaffold_db = get_scaffolds_db(specific_db)
            all_scaff_mols.extend(scaffold_db.mols)
            scaffold_db._get_similarity_matrix()
            simmatrix = scaffold_db.simmatrix
        else:
            df_outer = df_filt[df_filt['outer'] == outer]
            smiles = df_outer['SMILES_SARS2'].tolist()
            ids = df_outer['id'].tolist()
            mol_list = [mol.Mol(smile=smile, name=ids[i]) for i, smile in enumerate(smiles)]
            outer_db = moldb.MolDB(molList=mol_list)
            scaffold_db = get_scaffolds_db(outer_db)
            all_scaff_mols.extend(scaffold_db.mols)
            scaffold_db = moldb.MolDB(molList=all_scaff_mols)
            scaffold_db._get_similarity_matrix()
            simmatrix = scaffold_db.simmatrix
        clustering = DBSCAN(metric="precomputed",
                            eps=1-similarity_thr,
                            min_samples=1)
        clustering.fit(1-simmatrix)
        counter = Counter(clustering.labels_)
        num_clusters.append(len(counter.keys()))
        x_values.append(outer)

    return num_clusters, x_values
        #scaffolds = []
        #scaff = moldb.MolDB(molList=scaffolds)
        #scaff_names = [scaffold.name for scaffold in scaff.mols]
        #scaff_mols = scaff.mols
        #
        #scaff._get_similarity_matrix()
        #np.fill_diagonal(scaff.simmatrix.values, 0)
        #
        #below_threshold = scaff.simmatrix.lt(similarity_thr).all(axis=1)
        #
        #row_names = scaff.simmatrix.index[below_threshold].tolist()
        #
        ## get index of the row names in the scaff_names list
        #row_indexes = [scaff_names.index(name) for name in row_names]
        #print(row_indexes)
        ## get the scaff_mols with the indexes
        #scaffs = [scaff_mols[i] for i in row_indexes]


def plot_cluster_DBSCAN(csv_results,
                        smi_specific,
                        gscore_glob_thr,
                        gscore_ind_thr,
                        tanimoto_thr,
                        similarity_thrs,
                        outname):
    plt.figure(figsize=(10,6), dpi=500)
    for sim_thr in similarity_thrs:
        print(sim_thr)
        num_clusters, x_values = cluster_DBSCAN(csv_results=csv_results,
                                                smi_specific=smi_specific,
                                                gscore_glob_thr=gscore_glob_thr,
                                                gscore_ind_thr=gscore_ind_thr,
                                                tanimoto_thr=tanimoto_thr,
                                                similarity_thr=sim_thr)
        print(num_clusters, x_values)
        plt.plot(x_values, num_clusters, label='DBSCAN eps=%.2f'%sim_thr, marker='o')
    plt.xlabel('Affinity AL Cycle')
    plt.ylabel('DBSCAN scaffolds clusters')
    #plt.title('Scaffold evolution along outer loops')
    plt.legend()
    plt.ylim((0, 2000))
    plt.savefig(outname+'.pdf')

def new_scaffolds(csv_results,
                  smi_specific,
                  gscore_glob_thr,
                  gscore_ind_thr,
                  tanimoto_thr,
                  similarity_thr):

    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= gscore_glob_thr) & (df['gscore_SARS2'] <= gscore_ind_thr) & (df['gscore_SARS'] <= gscore_ind_thr) & (df['gscore_MERS'] <= gscore_ind_thr) & (df['max_tan'] <= tanimoto_thr)]
    df_filt = df_filt[df_filt['SMILES_SARS2'] != 'invalid_structure']
    print(df_filt)

    df_filt['outer'] = df_filt['id'].apply(lambda x: x.split('_')[0])
    df_filt['outer'] = df_filt['outer'].astype(int)
    total_outers = len(df_filt['outer'].unique())
    print(total_outers)
    print

    percs = []
    x_values = []
    all_scaff_fps = []
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=4)

    for outer in range(total_outers + 1):
        if outer == 0:
            specific_db = moldb.MolDB(smiDB=smi_specific)
            scaffold_db = get_scaffolds_db(specific_db)
            for scaffold in scaffold_db.mols:
                all_scaff_fps.append(fpgen.GetSparseCountFingerprint(scaffold.molrdkit))
            continue
        else:
            df_outer = df_filt[df_filt['outer'] == outer]
            smiles = df_outer['SMILES_SARS2'].tolist()
            ids = df_outer['id'].tolist()
            mol_list = [mol.Mol(smile=smile, name=ids[i]) for i, smile in enumerate(smiles)]
            outer_db = moldb.MolDB(molList=mol_list)
            scaffold_db = get_scaffolds_db(outer_db)
            #
            outer_scaff_fps = []
            for scaffold in scaffold_db.mols:
                outer_scaff_fps.append(fpgen.GetSparseCountFingerprint(scaffold.molrdkit))
            count = 0
            for fp in outer_scaff_fps:
                sims = DataStructs.BulkTanimotoSimilarity(fp, all_scaff_fps)
                if max(sims) < similarity_thr:
                    count += 1
            print(outer, scaffold_db.size, count, float(count/scaffold_db.size)*100)
            percs.append(float(count/scaffold_db.size)*100)
            x_values.append(outer)
            all_scaff_fps.extend(outer_scaff_fps)

    print(x_values)
    print(percs)
    return percs, x_values


def plot_new_scaffolds(csv_results,
                       smi_specific,
                       gscore_glob_thr,
                       gscore_ind_thr,
                       tanimoto_thr,
                       similarity_thrs,
                       outname):
    plt.figure()
    for sim_thr in similarity_thrs:
        percs, x_values = new_scaffolds(csv_results=csv_results,
                                        smi_specific=smi_specific,
                                        gscore_glob_thr=gscore_glob_thr,
                                        gscore_ind_thr=gscore_ind_thr,
                                        tanimoto_thr=tanimoto_thr,
                                        similarity_thr=sim_thr)
        plt.plot(x_values, percs, label='sim. thrs. = %.2f'%sim_thr, marker='o')
    plt.xlabel('Outer loop')
    plt.ylabel('/% of scaffolds with a similarity to all \n previous generated molecules < sim. thrs.')
    plt.title('Scaffold evolution along outer loops')
    plt.legend()
    plt.ylim((0, 110))
    plt.savefig(outname+'.pdf')
    
def map_ids_filtered_PAINS_ADMET_mols(csv_results, pains_admet_csv, outname):
    df = pd.read_csv(csv_results)
    pains_admet = pd.read_csv(pains_admet_csv)
    df = df.rename(columns={'SMILES_SARS':'canon_smiles'})
    df_map = pd.merge(pains_admet, df, on=['canon_smiles'])
    df_map = df_map[['id', 'canon_smiles', 'ad_Carcino', 'Predicted_Labels_Carcino',
                     'Predicted_Probabilities_Carcino', 'ad_cyp1a2_inhib',
                     'Predicted_Labels_cyp1a2_inhib', 'Predicted_Probabilities_cyp1a2_inhib',
                     'ad_herg', 'Predicted_Labels_herg', 'Predicted_Probabilities_herg',
                     'max_tan', 'global_gscore', 'gscore_SARS2', 'gscore_SARS', 'gscore_MERS']]

    df_map.to_csv(outname, index=False)
    
    
def get_metrics_generation(resdir, outer_name, inner_name, n, list_inners_per_outer):
    """Get the metrics of the whole generation:
        Ngen = number of molecules to generate.\
        Nval = number of valid molecules.\
        Nuniq = number of unique molecules.\
        Nunk = number of unknown molecules (not in the specific set).
        
        Validity = Nval/Ngen x 100\
        Uniqueness = Nuniq/Nval x 100\
        Unknown = Nunk/Nval x 100"""
    
    outers = len(list_inners_per_outer)
    inners = sum(list_inners_per_outer)
    print('Number of outers:', outers)
    print('Number of inners:', inners)
    
    #generated = [7500]*40 + [3500]*140 # this is case dependent!!
    generated = [5000]*inners 
    print(len(generated))
    print(generated)

    valid = []
    uniq = []
    unk = []

    for i in range(1,outers+1):
        for j in range(1,list_inners_per_outer[i-1]+1):
            print(i,j)
            
            file = glob.glob('%s/%s%s/%s%s_%s/*_generated_smiles.csv'%(resdir,outer_name,i,inner_name,i,j))[0]
            smiles = pd.read_csv(file)['smiles'].tolist()
            valid.append(len(smiles))
            
            all_mols = [mol.Mol(smile=x, allparamaters=True) for x in smiles]
            mols = []
            for mol1 in all_mols:
                print(mol1)
                try:
                    atoms = mol1.NumAtoms
                    mols.append(mol1)
                except:
                    print('Error')
                    continue
            
            mols_db = moldb.MolDB(molList=mols, verbose=False)
            mols_db.filterSimilarity(simt=1, alg='Morgan4',verbose=False)
            uniq.append(len(mols_db.dicDB))
            
            file2 = glob.glob('%s/%s%s/%s%s_%s/*_specific_smiles.csv'%(resdir,outer_name,i,inner_name,i,j))[0]
            specific = pd.read_csv(file2)['smiles'].tolist()
            all_spec_mols = [mol.Mol(smile=x, allparamaters=True) for x in specific]
            specific_mols = []
            for mol2 in all_spec_mols:
                try:
                    atoms = mol2.NumAtoms
                    specific_mols.append(mol2)
                except:
                    print('Error')
                    continue
            specific_mols_db = moldb.MolDB(molList=specific_mols, verbose=False)
            specific_mols_db.filterSimilarity(simt=1, alg='Morgan4',verbose=False)
            
            joint = moldb.joinMolDBs([mols_db,specific_mols_db], simt=1)
            unk.append(len(joint.dicDB)-len(specific_mols_db.dicDB))   

    print(len(valid))
    print(valid)
    print(len(uniq))
    print(uniq)
    print(len(unk))
    print(unk)

    df = pd.DataFrame({'Generated': generated, 'Valid': valid, 'Unique': uniq, 'Unknown': unk})
    df.to_csv('%s/metrics_generation_inners_tmp.csv'%resdir, index=False)
    
    validity = []
    for i in range(len(valid)):
        validity.append(valid[i]/generated[i]*100)

    uniqueness = []
    for i in range(len(valid)):
        if valid[i] == 0:
            uniqueness.append(0)
        else:
            uniqueness.append(uniq[i]/valid[i]*100)

    novelty = []
    for i in range(len(valid)):
        if uniq[i] == 0:
            novelty.append(0)
        else:
            novelty.append(unk[i]/uniq[i]*100)

    df2 = pd.DataFrame({'Generated': generated, 'Valid': valid, 'Unique': uniq, 'Unknown': unk,\
                        'Validity': validity, 'Uniqueness': uniqueness, 'Novelty': novelty})
    
    df2.to_csv('%s/metrics_generation_inners.csv'%resdir, index=False)
    print(df2)
    
def plot_metrics_generation(metrics_csv, output, list_inners_per_outer):
    df = pd.read_csv(metrics_csv)
    df['inner'] = df.index + 1
    print(df)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.plot(df['inner'], df['Validity'], marker='.', label='Validity')
    plt.plot(df['inner'], df['Uniqueness'], marker='.', label='Uniqueness')
    plt.plot(df['inner'], df['Novelty'], marker='.', label='Novelty')
    
    lines = list(itertools.accumulate(list_inners_per_outer))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
            
    plt.title('Metrics of Generation')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Percentage (%)')
    plt.legend()
    plt.xticks(np.array(range(0, sum(list_inners_per_outer), 10)))
    plt.savefig(output)
    
def calculate_mean_std_metrics_generation(metrics_csv):
    df = pd.read_csv(metrics_csv)
    print(df)
    
    mean_validity = df['Validity'].mean()
    std_validity = df['Validity'].std()
    mean_uniqueness = df['Uniqueness'].mean()
    std_uniqueness = df['Uniqueness'].std()
    mean_novelty = df['Novelty'].mean()
    std_novelty = df['Novelty'].std()
    
    print('Mean validity: %.2f +/- %.2f'%(mean_validity, std_validity))
    print('Mean uniqueness: %.2f +/- %.2f'%(mean_uniqueness, std_uniqueness))
    print('Mean novelty: %.2f +/- %.2f'%(mean_novelty, std_novelty))
    
def plot_QED_generation(resdir, outer_name, inner_name, n, list_inners_per_outer):
    "It plots the mean QED of the generated molecules per inner loop."
    
    outers = len(list_inners_per_outer)
    inners = sum(list_inners_per_outer)
    print('Number of outers:', outers)
    print('Number of inners:', inners)
    
    mean_qeds = []
    std_qeds = []
    for i in range(1,outers+1):
        for j in range(1,list_inners_per_outer[i-1]+1):
            print(i,j)
            file = glob.glob('%s/%s%s/%s%s_%s/*_generated_smiles.csv'%(resdir,outer_name,i,inner_name,i,j))[0]
            df = pd.read_csv(file)
            mean_qed = df['QED'].mean()
            mean_qeds.append(mean_qed)
            std_qed = df['QED'].std()
            std_qeds.append(std_qed)
            
    df = pd.DataFrame({'inner': range(1, inners + 1), 'mean_QED': mean_qeds, 'std_QED': std_qeds})
    df.to_csv('%s/QED_per_inner.csv'%(resdir), index=False)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.errorbar(df['inner'], df['mean_QED'], yerr=df['std_QED'], fmt='.', label='Mean QED', capsize=3, linestyle='-', color='blue', ecolor='lightblue', elinewidth=2, markeredgewidth=2)

    lines = list(itertools.accumulate(list_inners_per_outer))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
    plt.title('Mean QED of generated molecules')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Mean QED')
    plt.legend()
    plt.xticks(np.array(range(0, sum(list_inners_per_outer), 10)))
    plt.ylim((0, 1))
    plt.savefig('%s/QED_per_inner.pdf'%(resdir))
    
def plot_SA_generation(resdir, outer_name, inner_name, n, list_inners_per_outer):
    """It plots the mean SA of the generated molecules per inner loop."""
    
    outers = len(list_inners_per_outer)
    inners = sum(list_inners_per_outer)
    print('Number of outers:', outers)
    print('Number of inners:', inners)
    
    mean_sas = []
    std_sas = []
    for i in range(1,outers+1):
        for j in range(1,list_inners_per_outer[i-1]+1):
            print(i,j)
            file = glob.glob('%s/%s%s/%s%s_%s/*_generated_smiles.csv'%(resdir,outer_name,i,inner_name,i,j))[0]
            df = pd.read_csv(file)
            mean_sa = df['SAscore'].mean()
            mean_sas.append(mean_sa)
            std_sa = df['SAscore'].std()
            std_sas.append(std_sa)
            
    df = pd.DataFrame({'inner': range(1, inners + 1), 'mean_SA': mean_sas, 'std_SA': std_sas})
    df.to_csv('%s/SA_per_inner.csv'%(resdir), index=False)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.errorbar(df['inner'], df['mean_SA'], yerr=df['std_SA'], fmt='.', label='Mean SA', capsize=3, linestyle='-', color='blue', ecolor='lightblue', elinewidth=2, markeredgewidth=2)

    lines = list(itertools.accumulate(list_inners_per_outer))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
    plt.title('Mean SA of generated molecules')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Mean SA')
    plt.legend()
    plt.xticks(np.array(range(0, sum(list_inners_per_outer), 10)))
    plt.ylim((0, 10))
    plt.savefig('%s/SA_per_inner.pdf'%(resdir))
    
def plot_TA_similarity_generation(resdir, outer_name, inner_name, n, list_inners_per_outer):
    """It plots the mean TA similarity of the generated molecules per inner loop."""
    
    outers = len(list_inners_per_outer)
    inners = sum(list_inners_per_outer)
    print('Number of outers:', outers)
    print('Number of inners:', inners)
    
    mean_tas = []
    std_tas = []
    for i in range(1,outers+1):
        for j in range(1,list_inners_per_outer[i-1]+1):
            print(i,j)
            file = glob.glob('%s/%s%s/%s%s_%s/*_generated_smiles.csv'%(resdir,outer_name,i,inner_name,i,j))[0]
            df = pd.read_csv(file)
            mean_ta = df['tan_max'].mean()
            mean_tas.append(mean_ta)
            std_ta = df['tan_max'].std()
            std_tas.append(std_ta)
            
    df = pd.DataFrame({'inner': range(1, inners + 1), 'mean_TA': mean_tas, 'std_TA': std_tas})
    df.to_csv('%s/TA_per_inner.csv'%(resdir), index=False)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.errorbar(df['inner'], df['mean_TA'], yerr=df['std_TA'], fmt='.', label='Mean TA similarity', capsize=3, linestyle='-', color='blue', ecolor='lightblue', elinewidth=2, markeredgewidth=2)

    lines = list(itertools.accumulate(list_inners_per_outer))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
    plt.title('Max TA similarity of generated molecules')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Max TA similarity')
    plt.legend()
    plt.xticks(np.array(range(0, sum(list_inners_per_outer), 10)))
    plt.ylim((0, 1))
    plt.savefig('%s/TA_per_inner.pdf'%(resdir))
    
def plot_counts_inner_AL_filters(resdir, outer_name, inner_name, n, list_inners_per_outer):
    "It plots the counts of molecules applying the inner AL filters"
    
    outers = len(list_inners_per_outer)
    inners = sum(list_inners_per_outer)
    print('Number of outers:', outers)
    print('Number of inners:', inners)
    
    counts = []
    for i in range(1,outers+1):
        for j in range(1,list_inners_per_outer[i-1]+1):
            print(i,j)
            #print('%s/%s%s/%s%s_%s/*_generated_smiles_threshold.csv'%(resdir,outer_name,i,inner_name,i,j+1))
            file = glob.glob('%s/%s%s/%s%s_%s/*_generated_smiles_threshold.csv'%(resdir,outer_name,i,inner_name,i,j+1))[0]
            df = pd.read_csv(file)
            counts.append(len(df))
    print(counts)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.plot(range(1, inners + 1), counts, marker='.', label='Counts', color='blue')
    lines = list(itertools.accumulate(list_inners_per_outer))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
    plt.title('Counts of molecules applying inner AL filters')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Counts')
    plt.legend()
    plt.xticks(np.array(range(0, sum(list_inners_per_outer), 10)))
    plt.savefig('%s/counts_inner_AL_filters.pdf'%(resdir))
    