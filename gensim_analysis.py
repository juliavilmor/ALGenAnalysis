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
import itertools

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
    print('Total generated molecules thresholds: %s'%len(all_generated))
    

def create_table_gm_counts(results_dir, outname, save_df=False):
    """It creates a table of the gm results for each inner round.
       The table includes the molecule counts for each file."""

    # get all counts for each file
    counts_specific = {}
    counts_generated = {}
    counts_threshold = {}
    counts_success = {}
    results = glob.glob('%s/%s_*'%(results_dir, outname))
    print(results)
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
    print(counts_specific, counts_generated, counts_threshold, counts_success)

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
    
    # print the total number of molecules
    total_molecules = df.loc['generated'].sum()
    print('Total generated molecules: %s'%total_molecules)
    total_thresholds = df.loc['thresholds'].sum()
    print('Total molecules after thresholds: %s'%total_thresholds)

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
        df.drop_duplicates(subset=['smiles'], keep=False, inplace=True)
        return df
    
    spec_files = glob.glob('%s/%s_*/%s_specific_smiles.csv'%(results_dir,outname, outname))
    total = len(spec_files)
    for i in range(total):
        i = i+1
        if i == 1: # the initial specific set is the same
            os.system('cp %s/%s_%s/%s_specific_smiles.csv %s/%s_%s/%s_specific_smiles_simple.csv'%(results_dir,outname, i, outname, results_dir,outname, i, outname))
        elif i > 1 & i < (total-1): # the last one does not have specific set
            df_1 = '%s/%s_%s/%s_specific_smiles.csv'%(results_dir,outname, i-1, outname)
            df_2 = '%s/%s_%s/%s_specific_smiles.csv'%(results_dir,outname, i, outname)
            df = drop_duplicates_2_dfs(df_1, df_2)
            df.to_csv('%s/%s_%s/%s_specific_smiles_simple.csv'%(results_dir,outname, i, outname), index=False)

def _convert_smi_to_sdf(smi_file):
    """It converts a smi file into sdf file."""
    
    with Chem.SmilesMolSupplier(smi_file) as suppl:
        with Chem.SDWriter('%s'%smi_file.replace('.smi', '.sdf')) as w:
            for m in suppl:
                w.write(m)
    print('smi converted to sdf successfully')

def simplify_specific_sets_smi(list_spec_set):
        """ It extracts the repeated smiles in each specific set in smi format.
            The list must be ordered from the initial specific set to the last one."""
    
        total = len(list_spec_set)
        print(total)
        for i in range(total):
            if i == 0: # the initial specific set is the same
                pass
            elif i >= 1 & i < (total): # the last one does not have specific set
                # I want to extract the smiles that are not in the previous set
                print(i, i-1)
                set1 = set(open(list_spec_set[i], 'r').read().split('\n'))
                set2 = set(open(list_spec_set[i-1], 'r').read().split('\n'))
                print(len(set1), len(set2))
                simple = set1 - set2
                with open(list_spec_set[i].replace('.smi', '_simple.smi'), 'w') as f:
                    f.write('\n'.join(simple))
                print('Set %s simplified.'%list_spec_set[i])

def plot_modbs_tSNE_or_UMAP(list_of_sdfs, list_of_names, outdir, outname, sizes, alphas, markers, ptype='UMAP'):
    """ It plots the tSNE/UMAP of all the sets indicated in the list."""

    total = len(list_of_sdfs)
    moldb_list = []
    index_del = []
    for i, file in enumerate(list_of_sdfs):
        try:
            mols = moldb.MolDB(sdfDB=file, verbose=False)
            print(mols)
            moldb_list.append(mols)
            
        # is there are empy files (because there is no new generated molecules)    
        except: # remove them from other lists
            index_del.append(i)
    list_of_names = [i for j, i in enumerate(list_of_names) if j not in index_del]
    sizes  = [i for j, i in enumerate(sizes) if j not in index_del]
    alphas  = [i for j, i in enumerate(alphas) if j not in index_del]
    markers  = [i for j, i in enumerate(markers) if j not in index_del]
            
    plot_colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    plot_colors = plot_colors[1:total+1]
    plot_colors = [i for j, i in enumerate(plot_colors) if j not in index_del]

    min_dists = [0.2, 0.4, 0.6]
    neighbours = [50, 100, 200, 300]

    if ptype == 'tSNE':
        plot.plotTSNE(dbs = moldb_list, names = list_of_names, output='%s/%s'%(results_dir, outname),\
                    random_max = 10000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors,\
                    sizes=sizes, alphas=alphas, linewidth=0, markers=markers, n_iter=1000, perplexity=30,\
                    early_exaggeration=12,learning_rate='auto')

    if ptype == 'UMAP':
        for i in range(len(min_dists)):
            for j in range(len(neighbours)):
                plot.plotUMAP(dbs = moldb_list, names = list_of_names, output='%s/%s_md%s_nn%s'%(outdir, outname, min_dists[i], neighbours[j]),\
                            random_max = 50000, delimiter = None, alg = 'Morgan4', colors = plot_colors, sizes = sizes,  alphas = alphas,\
                            min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 10000, markers = markers, figsize = (9,6), \
                            linewidth = 0)
                print('UMAP with parameters md %s and nn %s done!'%(neighbours[j], min_dists[i]))

def get_smi_files_from_csv(csv_file, smi_column, outdir):
    """It gets the smiles from a csv file and saves them into a smi file."""
    
    df = pd.read_csv(csv_file)
    outname = os.path.basename(csv_file).replace('.csv', '.smi')
    with open('%s/%s'%(outdir, outname), 'w') as f:
        for index, row in df.iterrows():
            f.write('%s\n'%row[smi_column])
    print('smi file created successfully')
    
def remove_duplicates_from_sdf(sdf_file):
    """It removes duplicates from an sdf file."""
    
    mols = moldb.MolDB(sdfDB=sdf_file, verbose=False)
    mols.filterSimilarity(simt=1, alg='Morgan4',verbose=False)
    mols.saveToSdf(sdf_file.replace('.sdf', '_unique'))
    print('Duplicates removed successfully')
    print('The new sdf file contains %s molecules'%len(mols.mols))
    
def create_specific_set(glide_smi, previous_smi, outdir, outname):
    """It creates the specific set from the glide smiles and the previous smiles."""
    
    old_smi = moldb.MolDB(smiDB=previous_smi, verbose=False)
    new_smi = moldb.MolDB(smiDB=glide_smi, verbose=False)
    specific = moldb.joinMolDBs([old_smi, new_smi], simt=1)
    specific.saveToSmi('%s/%s.smi'%(outdir, outname))
    
def plot_UMAP(list_smis, list_names, outdir, outname, sizes, alphas, markers, colors=None):
    """It plots the tSNE of all the sets indicated in the list."""
    total = len(list_smis)
    mol_list = []
    for smi in list_smis:
        try:
            mol = moldb.MolDB(smiDB=smi, verbose=False)
            mol_list.append(mol)
        except Exception as e:
            print(f"Error loading SMILES file '{smi}': {e}")
    #mol_list = [moldb.MolDB(smiDB=smi, verbose=False) for smi in list_smis]
    if colors:
        colors = colors
    else:
        colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
        colors = colors[2:total+1]
        colors = colors + ['red']
    
    # If there are empty files (because there is no new generated molecules)
    # remove them from the lists
    len_moldb = [len(mol.mols) for mol in mol_list]
    index_del = [i for i, x in enumerate(len_moldb) if x == 0]
    list_names = [i for j, i in enumerate(list_names) if j not in index_del]
    sizes  = [i for j, i in enumerate(sizes) if j not in index_del]
    alphas  = [i for j, i in enumerate(alphas) if j not in index_del]
    markers  = [i for j, i in enumerate(markers) if j not in index_del]
    colors = [i for j, i in enumerate(colors) if j not in index_del]
    mol_list = [i for j, i in enumerate(mol_list) if j not in index_del]
    
    min_dists = [0.2, 0.4, 0.6]
    neighbours = [50, 100, 200]
    for i in range(len(min_dists)):
        for j in range(len(neighbours)):
            plot.plotUMAP(dbs = mol_list, names = list_names, output='%s/%s_md%s_nn%s'%(outdir, outname, min_dists[i], neighbours[j]),\
                        random_max = 50000, delimiter = None, alg = 'Morgan4', colors = colors, sizes = sizes,  alphas = alphas,\
                        min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 10000, markers = markers, figsize = (9,6), \
                        linewidth = 0)

def plot_tSNE(list_smis, list_names, outdir, outname, sizes, alphas, markers, colors=None):
    """It plots the tSNE of all the sets indicated in the list."""
    total = len(list_smis)
    mol_list = [moldb.MolDB(smiDB=smi, verbose=False) for smi in list_smis]
    if colors:
        colors = colors
    else:
        colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
        colors = colors[2:total+1]
        colors = colors + ['red']
    
    # If there are empty files (because there is no new generated molecules)
    # remove them from the lists
    len_moldb = [len(mol.mols) for mol in mol_list]
    index_del = [i for i, x in enumerate(len_moldb) if x == 0]
    list_names = [i for j, i in enumerate(list_names) if j not in index_del]
    sizes  = [i for j, i in enumerate(sizes) if j not in index_del]
    alphas  = [i for j, i in enumerate(alphas) if j not in index_del]
    markers  = [i for j, i in enumerate(markers) if j not in index_del]
    colors = [i for j, i in enumerate(colors) if j not in index_del]
    mol_list = [i for j, i in enumerate(mol_list) if j not in index_del]
    
    perplexities = [10, 20, 30, 40]
    
    for i in range(len(perplexities)):
        plot.plotTSNE(dbs = mol_list, names = list_names, output='%s/%s_pp%s'%(outdir, outname, perplexities[i]),\
                        random_max = 10000, delimiter = None, alg = 'Morgan4', colors = colors,\
                        sizes=sizes, alphas=alphas, linewidth=0, n_iter=1000, perplexity=perplexities[i],\
                        early_exaggeration=12,learning_rate='auto', markers=markers, figsize=(9,6))       
        
def plot_specific_set_evolution(results_dir, outdir, outname):
    """It plots the evolution of the specific set."""
    plt.figure()
    fig, ax = plt.subplots(figsize=(8,6))
    outers = glob.glob('%s/outer_?'%results_dir)
    outers.sort()
    outers2 = glob.glob('%s/outer_??'%results_dir)
    outers2.sort()
    outers2 = outers2[:-1]
    outers = outers + outers2

    sizes_specific = []
    inner_sizes = []
    for i, outer in enumerate(outers):
        table = pd.read_csv('%s/outer_%s/table_of_counts_trans.csv'%(results_dir, i+1))
        size_specific = table['specific'].tolist()
        sizes_specific.extend(size_specific)
        inner_size = len(table)
        inner_sizes.append(inner_size)
    x = list(range(1, len(sizes_specific)+1))
    plt.plot(x, sizes_specific, marker='.', color='blue', label='specific set')
    
    lines = list(itertools.accumulate(inner_sizes))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
        
    ax.set_yscale('log', base=10)
    ax.set_ylabel('Specific set size (log scale)')
    ax.set_xlabel('Outer AL cycle')
    
    plt.savefig('%s/%s.pdf'%(outdir, outname))

def superpose_specific_set_evolution(results_dir_1, results_dir_2, gscore_values_1, gscore_values_2, outdir, outname):
    """It superposes the evolution of the specific set."""
    plt.figure()
    fig, ax = plt.subplots(figsize=(8,6))
    outers = glob.glob('%s/outer_?'%results_dir_1)
    outers.sort()
    outers_ = glob.glob('%s/outer_??'%results_dir_1)
    outers_.sort()
    outers = outers + outers_
    outers = outers[:-1]

    sizes_specific = []
    inner_sizes = []
    for i, outer in enumerate(outers):
        table = pd.read_csv('%s/outer_%s/table_of_counts_trans.csv'%(results_dir_1, i+1))
        size_specific = table['specific'].tolist()
        sizes_specific.extend(size_specific)
        inner_size = len(table)
        inner_sizes.append(inner_size)
    x = list(range(1, len(sizes_specific)+1))
    plt.plot(x, sizes_specific, marker='.', color='seagreen', label='Regular Configuration')
    
    outers2 = glob.glob('%s/outer_?'%results_dir_2)
    outers2.sort()
    outers_ = glob.glob('%s/outer_??'%results_dir_2)
    outers_.sort()
    outers2 = outers2 + outers_
    outers2 = outers2[:-1]
    
    sizes_specific2 = []
    inner_sizes2 = []
    for i, outer in enumerate(outers2):
        table = pd.read_csv('%s/outer_%s/table_of_counts_trans.csv'%(results_dir_2, i+1))
        size_specific = table['specific'].tolist()
        sizes_specific2.extend(size_specific)
        inner_size = len(table)
        inner_sizes2.append(inner_size)
    x2 = list(range(1, len(sizes_specific2)+1))
    plt.plot(x2, sizes_specific2, marker='.', color='royalblue', label='Ablated Configuration')
    plt.legend(loc='upper right')
    
    lines1 = list(itertools.accumulate(inner_sizes))
    lines2 = list(itertools.accumulate(inner_sizes2))
    
    for line in lines2:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)

    ax.set_yscale('log', base=10)
    ax.set_ylabel('Specific set size (log scale)')
    ax.set_xlabel('Chemical AL cycle')

    # Create a secondary y-axis with the gscore values
    ax2 = ax.twinx()
    ax2.set_ylabel('Global Docking Score Threshold')
    ax2.plot(lines1, gscore_values_1, marker='.', color='lightgreen', label='Regular Configuration', linewidth=5, alpha=0.5)
    ax2.plot(lines2, gscore_values_2, marker='.', color='cornflowerblue', label='Ablated Configuration', linewidth=5, alpha=0.5)
    ax2.set_ylim(-9, -6)
    
    plt.legend(loc='lower right')
    plt.savefig('%s/%s.pdf'%(outdir, outname))

def plot_specific_set_evolution_dockingscores(results_dir, gscore_values, outdir, outname):
        plt.figure()
        fig, ax = plt.subplots(figsize=(10,6), dpi=300)
        outers = glob.glob('%s/outer_?'%results_dir_1)
        outers.sort()
        outers_ = glob.glob('%s/outer_??'%results_dir_1)
        outers_.sort()
        outers = outers + outers_
        outers = outers[:-1]

        sizes_specific = []
        inner_sizes = []
        for i, outer in enumerate(outers):
            table = pd.read_csv('%s/outer_%s/table_of_counts_trans.csv'%(results_dir, i+1))
            size_specific = table['specific'].tolist()
            sizes_specific.extend(size_specific)
            inner_size = len(table)
            inner_sizes.append(inner_size)
        x = list(range(1, len(sizes_specific)+1))
        plt.plot(x, sizes_specific, marker='.', color='royalblue', label='specific set IL1beta')

        lines = list(itertools.accumulate(inner_sizes))

        for line in lines:
            plt.axvline(line, color='black', linestyle=':', alpha=0.5)

        ax.set_yscale('log', base=10)
        ax.set_ylabel('Specific set size (log scale)')
        ax.set_xlabel('Inner round')

        # Create a secondary y-axis with the gscore values
        ax2 = ax.twinx()
        ax2.set_ylabel('Average Docking score threshold')
        ax2.plot(lines, gscore_values, marker='.', color='cornflowerblue', label='Docking score threshold', linewidth=5, alpha=0.5)
        ax2.set_ylim(-9, -6)

        plt.legend(loc='lower right')
        plt.savefig('%s/%s.png'%(outdir, outname))

def tensordti_csv_format(csv, outer):
    """
    The output should be in the following format:
    id,smiles

    Args:
        csv (str): all_generated_molecules.csv
        outer (int): outer round number
    """
    df = pd.read_csv(csv)
    df = df.rename(columns={'Unnamed: 0': 'id'})
    df['id'] = str(outer) + '_' + df['inner'].astype(str) + '_' + df['id'].astype(str)
    df = df[['id', 'smiles']]
    
    df = df.drop_duplicates(subset=['smiles'])
    
    df.to_csv(csv.replace('.csv', '_tensordti.csv'), index=False)
    
def remove_duplicates_from_csv(csv_file):
    """It removes duplicates from a csv file."""
    
    df = pd.read_csv(csv_file)
    mol_list = df['smiles'].tolist()
    mol_list = [mol.Mol(smile=smile) for smile in mol_list]
    mols = moldb.MolDB(molList=mol_list, verbose=False)
    mols.filterSimilarity(simt=1, alg='Morgan4', verbose=False)
    unique_smiles = mols.smiles
    
    unique_df = df[df['smiles'].isin(unique_smiles)]
    unique_df.to_csv(csv_file.replace('.csv', '_unique.csv'), index=False)
    
    print('From %s molecules, %s unique molecules were found.'%(len(df), len(unique_df)))
    
    
if __name__ == "__main__":
    
    n = 'gensim_mt'
    
    #get_all_generated_molecules(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n)
    
    #create_table_gm_counts(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n, save_df=True)
    
    #plot_all_props_in_one(results_dir='/home/cactus/julia/gensim/full/outer1', save_fig=True)
    
    #convert_csv_to_sdf_file(csv_to_convert='/home/cactus/julia/gensim/full/outer1/all_generated_molecules.csv', outdir='/home/cactus/julia/gensim/full/outer1')
    
    # _simplify_specific_sets(results_dir='/home/cactus/julia/gensim/full/outer1', outname=n)
    # spec_simples = glob.glob('/home/cactus/julia/gensim/full/outer1/gensim_mt_*/gensim_mt_specific_smiles_simple.csv')
    # spec_inner = [x.split('/')[-2].split('_')[-1] for x in spec_simples]
    # for i, spec_simple in enumerate(spec_simples):
    #     convert_csv_to_sdf_file(csv_to_convert=spec_simple, outdir='/home/cactus/julia/gensim/full/outer1/gensim_mt_%s'%spec_inner[i], inner=spec_inner[i])

    # _convert_smi_to_sdf(smi_file='/home/cactus/julia/gensim/full/outer1/full_init_spec_set.smi')

    spec_simples = glob.glob('/home/cactus/julia/gensim/full/outer1/gensim_mt_*/gensim_mt_specific_smiles_simple.sdf')
    sdf_list = []
    rev_inners = list(range(1, len(spec_simples)+1))
    rev_inners.sort(reverse=True)
    for i in rev_inners:
        sdf_list.append('/home/cactus/julia/gensim/full/outer1/gensim_mt_%s/gensim_mt_specific_smiles_simple.sdf'%i)
    names = ['specific_19', 'specific_18', 'specific_17', 'specific_16', 'specific_15',\
            'specific_14', 'specific_13', 'specific_12', 'specific_11', 'specific_10',\
            'specific_9', 'specific_8', 'specific_7', 'specific_6', 'specific_5',\
            'specific_4', 'specific_3', 'specific_2', 'specific_1', 'specific_0']
    sizes = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]
    plot_modbs_tSNE_or_UMAP(list_of_sdfs=sdf_list, list_of_names=names, outdir='/home/cactus/julia/gensim/full/outer1/plots', outname='UMAP_outer1',\
                            sizes=sizes, alphas=alphas, markers=markers, ptype='UMAP')
