from MolecularAnalysis import mollib, plotlib
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


########## GENAI BLOCK ###########

def join_all_molecule_files(results_dir, file_names_to_join):
    """ It joins all files with the same name.
        Ex: Join all threshold files from trial_x folders"""

    # get all molecule smiles
    smiles_list = []
    results = glob.glob('%s/*'%results_dir)
    for result in results:
       trial = glob.glob('%s/*'%result)
       for tfile in trial:
           nfile = os.path.basename(tfile)
           if file_names_to_join not in nfile: continue
           file_to_read = pd.read_csv(tfile)
           smiles = file_to_read['smiles'].values.tolist()
           smiles_list.extend(smiles)
    # save into new file
    new_file = open('%s/all_%s.smi'%(results_dir, file_names_to_join), 'w')
    for smile in smiles_list:
        new_file.write('%s\n'%smile)
    new_file.close()

def create_table_gmn_counts(results_dir, save_df=False):
    """It creates a table of the gmn results for each round.
       The table includes the molecule counts for each file."""

    # get all counts for each file
    counts_specific = {}
    counts_generated = {}
    counts_threshold = {}
    counts_success = {}
    results = glob.glob('%s/*'%results_dir)
    for result in results:
        nametrial = os.path.basename(result)
        if 'trial' not in nametrial: continue
        ntrial = nametrial.split('_')[-1].strip()
        try:
            ntrial = int(ntrial)
        except:
            pass # just in case there are more non-desirable files
        #ntrial = int(float(ntrial))
        trial = glob.glob('%s/*'%result)
        for tfile in trial:
            nfile = os.path.basename(tfile)
            if 'trial' not in nfile: continue
            #print(nfile)
            if nfile == 'trial_specific_smiles.csv': # the first one
                if nametrial != 'trial_1': continue
                initial_specific = pd.read_csv(tfile)
                num_in_specific = len(initial_specific)
                counts_specific[1] = num_in_specific
            if nfile == 'trial_generated_smiles_merged.csv': # n
                if nametrial == 'trial_11': continue
                specific = pd.read_csv(tfile)
                num_specific = len(specific)
                counts_specific[ntrial] = num_specific
            if nfile == 'trial_generated_smiles.csv': # n
                generated = pd.read_csv(tfile)
                num_generated = len(generated)
                counts_generated[ntrial] = num_generated

            if nfile == 'trial_generated_smiles_threshold.csv': # n - 1
                threshold = pd.read_csv(tfile)
                num_threshold = len(threshold)
                counts_threshold[ntrial - 1] = num_threshold

            if nfile == 'report_trial.out': # n
                report = open(tfile, 'r')
                for line in report:
                    if line.startswith('Success'):
                        rate = line.split(': ')[-1].replace('%', '')
                        rate = round(float(rate), 2)
                        counts_success[ntrial] = rate

    # put into a dataframe
    try:
        counts_dic = dict([(k,[counts_specific[k],counts_generated[k], counts_threshold[k], counts_success[k]]) for k in counts_specific])
    except:
        counts_specific = dict(sorted(counts_specific.items()))
        counts_specific.popitem()
        counts_dic = dict([(k,[counts_specific[k],counts_generated[k], counts_threshold[k], counts_success[k]]) for k in counts_specific])
    counts_dic = dict(sorted(counts_dic.items()))
    index_labels = ['specific_set', 'generated_set', 'thresholds', 'success_rate(%)']
    df = pd.DataFrame(counts_dic, index=index_labels)
    print(df)

    if save_df:
        df.to_csv('%s/table_of_results.csv'%results_dir)

def plot_results_properties(results_dir, which_property, set_color='blue', title=True, save_plot=False):
    """It plots the trials vs. the selected propery
       of the GMN results (From the previous dataframe!!). """

    props_table = pd.read_csv('%s/table_of_results.csv'%results_dir, index_col=0)
    trials = props_table.columns
    props_table.loc['trials'] = trials
    props_table = props_table.T
    plt.plot('trials', which_property, data=props_table, marker='o', alpha=0.4, color=set_color)
    if title:
        plt.title('Trials vs. %s'%which_property)
    plt.xlabel('Trials')
    if which_property == 'success_rate(%)':
        plt.ylabel(which_property)
    else:
        plt.ylabel('%s (counts)'%which_property)

    if save_plot:
        plt.savefig('%s/plot_trial_%s.png'%(results_dir, which_property))


def plot_all_props_in_one(results_dir, save_fig=True):
    """ It plots the previous 4 properties plots into one.
        It is a 2x2 plot."""

    fig, ax = plt.subplots(2,2)
    plt.subplot(221)
    plot_results_properties(results_dir, 'specific_set', set_color='blue')
    plt.subplot(222)
    plot_results_properties(results_dir, 'generated_set', set_color='orange')
    plt.subplot(223)
    plot_results_properties(results_dir, 'thresholds', set_color='green')
    plt.subplot(224)
    plot_results_properties(results_dir, 'success_rate(%)', set_color='purple')
    plt.suptitle('Plots of properties for each generation epoch')
    fig.tight_layout(pad=1.0)
    if save_fig:
        plt.savefig('%s/plot_all_properties.png'%results_dir)

def _create_smiles_txt_file(results_dir, file_to_convert):
    """ It create smiles txt files to be uploaded to the molDB class.
    -------------NOOOOT ANYMOREEEE, UPLOAD MOLLIB-----------------"""

    file_to_convert = '%s/%s'%(results_dir, file_to_convert)

    smi_dir = '%s/smi_files'%results_dir
    if not os.path.isdir(smi_dir):
         os.system('mkdir %s'%smi_dir)
    if '.smi' in file_to_convert:
        filename = os.path.basename(file_to_convert).replace('.smi', '')
        smi_file = open(file_to_convert, 'r')
        new_file = open('%s/%s.txt'%(smi_dir, filename), 'w')
        for i, line in enumerate(smi_file):
            smile = line.strip()
            new_file.write('%s %s %s\n'%(smile, smile, i))
        new_file.close()

    if '.csv' in file_to_convert:
        filename = os.path.basename(file_to_convert).replace('.csv','')
        csv_file = pd.read_csv(file_to_convert)
        smiles = initial_specific.iloc[:,0].values.tolist()
        new_file = open('%s/%s.txt'%(smi_dir, filename), 'w')
        for i, smile in enumerate(smiles):
            new_file.write('%s %s %s\n'%(smile, smile, i))
        new_file.close()

def _convert_csv_to_smi_file(csv_to_convert, outdir, smi_column=1):
    """ It convert csv to smi file to be uploaded to the molDB class."""

    filename = os.path.basename(csv_to_convert).replace('.csv','')
    csv_file = pd.read_csv(csv_to_convert)
    smiles = csv_file.iloc[:,smi_column].values.tolist()
    new_file = open('%s/%s.smi'%(outdir, filename), 'w')
    for i, smile in enumerate(smiles):
        new_file.write('%s\n'%smile)
    new_file.close()


def plot_tSNE_before_and_after(results_dir):
    """It plots the tSNE of the results of a GMN round,
       before and after."""

    specific_1 = mollib.MolDB(smiDB='%s/smi_files/trial_specific_smiles.smi'%results_dir, verbose=False)
    generated = mollib.MolDB(smiDB='%s/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    glide_filtered = mollib.MolDB(smiDB='%s/glide_filtered.smi'%results_dir, verbose=False)


    plotlib.plot_tSNE(dbs=[specific_1, generated, glide_filtered], names=['specific_1','generated', 'glide_filtered'],output='%s/tSNE_before_and_after'%results_dir, random_max = 10000, delimiter = None, fpsalg = 'Morgan4',colors = ['lightcoral', 'lightseagreen', 'teal'], sizes=[0.7, 0.7,  0.7], alphas = [0.7, 0.2, 0.9], linewidths = [0, 0, 0], n_iter=10000, perplexity=20, early_exaggeration=12, learning_rate='auto')


def _create_specific_set_smi_files(results_dir):
    """ It creates smi files of all specific sets
        to upload then to molDB class."""

    smi_dir = '%s/smi_files'%results_dir
    if not os.path.isdir(smi_dir):
         os.system('mkdir %s'%smi_dir)

    results = glob.glob('%s/*'%results_dir)
    for result in results:
        nametrial = os.path.basename(result)
        if 'trial' not in nametrial: continue
        ntrial = nametrial.split('_')[-1].strip()
        try:
            ntrial = int(ntrial)
        except:
            pass # just in case there are more non-desirable files
        trial = glob.glob('%s/*'%result)
        for tfile in trial:
            nfile = os.path.basename(tfile)
            if 'trial' not in nfile: continue
            #print(nfile)
            if nfile == 'trial_specific_smiles.csv': # the first one
                if nametrial != 'trial_1': continue
                initial_specific = pd.read_csv(tfile)
                smiles = initial_specific.iloc[:,1].values.tolist()
                new_file = open('%s/specific_1.smi'%(smi_dir), 'w')
                for i, smile in enumerate(smiles):
                    new_file.write('%s\n'%(smile))
                new_file.close()
            if nfile == 'trial_generated_smiles_merged.csv': # n
                if nametrial == 'trial_11': continue
                specific = pd.read_csv(tfile)
                smiles = specific.iloc[:,0].values.tolist()
                new_file = open('%s/specific_%s.smi'%(smi_dir, ntrial), 'w')
                for i, smile in enumerate(smiles):
                    new_file.write('%s\n'%(smile))
                new_file.close()

def plot_specific_sets_tSNE(results_dir):
    """From the previous files...
       It plots the tSNE of all the specific sets.

       NOTE: at each trial the molecules are superimposed. To avoid it,
       use the following function (simplified.....)"""

    specific_1 = mollib.MolDB(txtDB='%s/smi_files/specific_1.txt'%results_dir, verbose=False)
    specific_2 = mollib.MolDB(txtDB='%s/smi_files/specific_2.txt'%results_dir, verbose=False)
    specific_3 = mollib.MolDB(txtDB='%s/smi_files/specific_3.txt'%results_dir, verbose=False)
    specific_4 = mollib.MolDB(txtDB='%s/smi_files/specific_4.txt'%results_dir, verbose=False)
    specific_5 = mollib.MolDB(txtDB='%s/smi_files/specific_5.txt'%results_dir, verbose=False)
    specific_6 = mollib.MolDB(txtDB='%s/smi_files/specific_6.txt'%results_dir, verbose=False)
    specific_7 = mollib.MolDB(txtDB='%s/smi_files/specific_7.txt'%results_dir, verbose=False)
    specific_8 = mollib.MolDB(txtDB='%s/smi_files/specific_8.txt'%results_dir, verbose=False)
    specific_9 = mollib.MolDB(txtDB='%s/smi_files/specific_9.txt'%results_dir, verbose=False)
    specific_10 = mollib.MolDB(txtDB='%s/smi_files/specific_10.txt'%results_dir, verbose=False)

    plot_colors = mcp.gen_color(cmap="BuPu",n=10)
    plotlib.plot_tSNE(dbs=[specific_10, specific_9, specific_8, specific_7, specific_6, specific_5, specific_4, specific_3, specific_2, specific_1], names=['specific_10','specific_9', 'specific_8', 'specific_7', 'specific_6', 'specific_5', 'specific_4', 'specific_3', 'specific_2', 'specific_1'],output='%s/tSNE_specific_sets'%results_dir, random_max = 10000, delimiter = None, fpsalg = 'Morgan4',colors = plot_colors, sizes=[0.5,0.5, 0.6, 0.6, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7], alphas = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], n_iter=1000, perplexity=30, early_exaggeration=12, learning_rate='auto')

def _simplify_specific_sets(results_dir):
    """ It extracts the repeated smiles in each specific set.
        As a result, it creates smi files of all the specific sets
        containing only the new generated molecules."""

    if not os.path.isdir('%s/smi_files'%results_dir):
        os.system('mkdir %s/smi_files'%results_dir)

    smi_files = glob.glob('%s/trial_*/trial_specific_smiles.csv'%results_dir)
    smiles = {}
    for i, specific in enumerate(smi_files):
        num = specific.split('/')[2].split('_')[-1]
        smi_file = pd.read_csv(specific)
        smile = smi_file.iloc[:,1].values.tolist()
        smiles[num] = set(smile)

    smiles = dict(sorted(smiles.items()))

    specifics = []
    specific_0 = smiles['1']
    specifics.append(specific_0)
    specific_1 = smiles['2'] - smiles['1']
    specifics.append(specific_1)
    specific_2 = smiles['3'] - smiles['2']
    specifics.append(specific_2)
    specific_3 = smiles['4'] - smiles['3']
    specifics.append(specific_3)
    specific_4 = smiles['5'] - smiles['4']
    specifics.append(specific_4)
    specific_5 = smiles['6'] - smiles['5']
    specifics.append(specific_5)
    specific_6 = smiles['7'] - smiles['6']
    specifics.append(specific_6)
    specific_7 = smiles['8'] - smiles['7']
    specifics.append(specific_7)
    #specific_8 = smiles['9'] - smiles['8']
    #specifics.append(specific_8)
    #specific_9 = smiles['10'] - smiles['9']
    #specifics.append(specific_9)

    for i, specific in enumerate(specifics):
        new_file = open('%s/smi_files/simp_specific_%d.smi'%(results_dir, i), 'w')
        for j, smile in enumerate(specific):
            new_file.write('%s\n'%(smile))
        new_file.close()

def plot_specific_simplified_sets_tSNE(results_dir):
    """From the previous files...
       It plots the tSNE of all the simplified specific sets."""

    specific_0 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_0.smi'%results_dir, verbose=False)
    specific_1 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_1.smi'%results_dir, verbose=False)
    specific_2 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_2.smi'%results_dir, verbose=False)
    specific_3 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_3.smi'%results_dir, verbose=False)
    specific_4 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_4.smi'%results_dir, verbose=False)
    specific_5 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_5.smi'%results_dir, verbose=False)
    specific_6 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_6.smi'%results_dir, verbose=False)
    specific_7 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_7.smi'%results_dir, verbose=False)
    #specific_8 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_8.smi'%results_dir, verbose=False)
    #specific_9 = mollib.MolDB(smiDB='%s/smi_files/simp_specific_9.smi'%results_dir, verbose=False)

    plot_colors = mcp.gen_color(cmap="YlGnBu", n=9)
    plot_colors = plot_colors[1:9]
    sizes = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    linewidths = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    markers = ["o", "o", "o", "o", "o", "o", "o", "*"]

    #plotlib.plot_tSNE(dbs=[specific_9, specific_8, specific_7, specific_6, specific_5, specific_4, specific_3, specific_2, specific_1, specific_0], names=['specific_9', 'specific_8', 'specific_7', 'specific_6', 'specific_5', 'specific_4', 'specific_3', 'specific_2', 'specific_1', 'specific_0'],output='%s/tSNE_specific_simp_sets'%results_dir, random_max = 10000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes=sizes, alphas=alphas, linewidths=linewidths, markers=markers, n_iter=1000, perplexity=30, early_exaggeration=12,learning_rate='auto')

    min_dists = [0.2]
    neighbours = [100, 150, 200]

    for i in range(len(min_dists)):
        for j in range(len(neighbours)):
            plotlib.plot_UMAP(dbs=[specific_7, specific_6, specific_5, specific_4, specific_3, specific_2, specific_1, specific_0], names=['specific_7', 'specific_6', 'specific_5', 'specific_4', 'specific_3', 'specific_2', 'specific_1', 'specific_0'], output='%s/UMAP_inner_specifics_md%s_nn%s'%(results_dir, min_dists[i], neighbours[j]), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes,  alphas = alphas, min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 10000, markers = markers, figsize = (9,6), linewidths = linewidths)

def plot_tSNE_for_generation_round(results_dir, biased=False):
    """It plots the tSNE of the gernated compounds at each round (outer generation)."""

    if biased:
        initial_specific = mollib.MolDB(smiDB='%s/round_1/smi_files/simp_specific_0.smi'%results_dir, verbose=False)
        generated_1 = mollib.MolDB(smiDB='%s/round_1/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_2 = mollib.MolDB(smiDB='%s/round_2_biased/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_3 = mollib.MolDB(smiDB='%s/round_3_biased/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_4 = mollib.MolDB(smiDB='%s/round_4_biased/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        outname = 'tSNE_generated_compounds_biased'
    else:
        initial_specific = mollib.MolDB(smiDB='%s/round_1/Mpro_initial_specificset.smi'%results_dir, verbose=False)
        generated_1 = mollib.MolDB(smiDB='%s/round_1/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_2 = mollib.MolDB(smiDB='%s/round_2/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_3 = mollib.MolDB(smiDB='%s/round_3_reducedta/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        generated_4 = mollib.MolDB(smiDB='%s/round_4/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
        outname = 'tSNE_generated_compounds_unbiased'

    plot_colors = mcp.gen_color(cmap="YlGnBu", n=6)
    plot_colors = plot_colors[1:6]
    sizes = [0.4, 0.4, 0.4, 0.4, 0.6]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9]
    linewidths = [0.0, 0.0, 0.0, 0.0, 0.0]
    markers = ["o", "o", "o", "o", "*"]

    #plotlib.plot_tSNE(dbs=[generated_4, generated_3, generated_2, generated_1, initial_specific], names=['generated_4', 'generated_3', 'generated_2', 'generated_1', 'initial_specific'],output='%s/plots/%s'%(results_dir, outname), random_max = 10000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes, alphas = alphas, linewidths = linewidths, markers = markers, n_iter=10000, perplexity=30, early_exaggeration=12,learning_rate='auto')

    #plotlib.plot_UMAP(dbs=[generated_4, generated_3, generated_2, generated_1, initial_specific], names=['generated_4', 'generated_3', 'generated_2', 'generated_1', 'initial_specific'], output='%s/plots/UMAP_1'%(results_dir), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes, alphas = alphas, min_dist = 0.1, n_neighbors = 100, n_epochs = 1000, markers = markers, figsize = (8,8), linewidths = linewidths)

    min_dists = [0, 0.2, 0.4]
    neighbours = [100, 150, 200, 250]

    for i in range(len(min_dists)):
        for j in range(len(neighbours)):
            plotlib.plot_UMAP(dbs=[generated_4, generated_3, generated_2, generated_1, initial_specific], names=['outer4', 'outer3', 'outer2', 'outer1', 'specific'], output='%s/plots/UMAP_final_md%s_nn%s_full'%(results_dir, min_dists[i], neighbours[j]), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes,  alphas = alphas, min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 10000, markers = markers, figsize = (8,6), linewidths = linewidths)


def plot_specific_set_evolution(results_dir, biased=True):
    """It plots the evolution of the specific set over
       the outer generations."""

    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 6))
    sizes_specific = []
    if biased:
        batches = glob.glob('%s/round_*_biased_new'%(results_dir))
        batch1 = '%s/round_1'%results_dir
        batches.append(batch1)
    else:
        batches = glob.glob('%s/round_?'%(results_dir))
    batches.sort()
    for j, batch in enumerate(batches):
        table = pd.read_csv('%s/table_of_results.csv'%batch, index_col=0)
        size_specific = table.loc[['specific_set']].values.tolist()[0]
        sizes_specific.extend(size_specific)
    x = list(range(1,len(sizes_specific)+1))
    plt.plot(x,sizes_specific,label='specific set', marker='.')
    for i in range(0,50,10):
        plt.axvline(x = i, color = 'black',linestyle='--')
    ax.legend(prop={'size': 12})
    plt.xlabel('Inner generation loops',fontsize=18)
    plt.ylabel('Specific set size',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig('%s/plots/specific_set_evolution.pdf'%results_dir)



####### GLIDE BLOCK ##########

def get_best_glide_docking_pose_for_ligand(path_to_csv):
    """ I returns only the best poses for each ligand,
        based on the glide gscore property."""

    results = pd.read_csv(path_to_csv)
    ids = results['title'].values.tolist()
    ids = set(ids) # unique
    best = {}
    for id in ids:
        # min glide score
        glide_score = results.loc[results['title'] == id]['r_i_glide_gscore'].min()
        best[id] = glide_score
    filts = []
    for k, v in best.items():
        filtered = results.loc[(results['title'] == k)&(results['r_i_glide_gscore'] == v)]
        filts.append(filtered)
    best_df = pd.concat(filts)

    file_name = path_to_csv.replace('.csv', '_best.csv')
    best_df.to_csv(file_name)

def plot_glide_scores(path_to_best_csv, insert_title, savefig=True):
    """Histogram of gscores."""

    df = pd.read_csv(path_to_best_csv)
    sns.histplot(data=df, x='r_i_glide_gscore').set(title=insert_title)
    filename = path_to_best_csv.replace('.csv', '_gscores.png')
    if savefig:
        plt.savefig(filename)

def superimpose_histograms(list_of_csvs, list_of_labels, insert_title, out, savefig=True, legend_loc='upper right'):
    """Superimpose histograms to compare.
       Ex: gscores comparison inital vs. generated."""

    plt.figure(figsize=(6,6))
    colors = mcp.gen_color(cmap="YlGnBu", n=6)
    colors = colors[1:6]
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
    if savefig:
        plt.savefig(out)

def zoom_in_histogram(threshold, out, savefig=True):
    """ If zooms in the previous plot by indicating the x limit."""

    #csvs = ['smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_4.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_3.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_2.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_1.csv','initial_set/Mpro_6xbgA1_receptor_SARS2_Mpro_inhibitors_DB_pv_best.csv']
    #csvs = ['smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_totalgen.csv', 'initial_set/Mpro_6xbgA1_receptor_SARS2_Mpro_inhibitors_DB_pv_best.csv']
    csvs = ['paninhibitor/plots/glides/global_sel_glide_gscores_4.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_3.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_2.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_1.csv', 'paninhibitor/plots/glides/global_initial_gscores.csv']
    superimpose_histograms(list_of_csvs=csvs, list_of_labels=['outer4', 'outer3', 'outer2', 'outer1', 'specific'], insert_title='Gscore of initial and resulting compounds', out='paninhibitor/plots/gscore_initial_and_resulting.png', savefig=False, legend_loc='upper left')
    #superimpose_histograms(list_of_csvs=csvs, list_of_labels=['generated', 'specific'], insert_title='Gscore of initial and resulting compounds', out='smallmol_specific_set/plots/     gscore_initial_and_resulting_totalgen.png', savefig=False, legend_loc='upper left')
    plt.xlim(-10, threshold)
    #plt.ylim(0, 290)
    if savefig:
        plt.savefig(out)

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
        plt.hist(gscores, bins=15, alpha=0.5)
        plt.axvline(x = -5.9, color = 'r')
        plt.savefig('hist_gscores.png')
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


def _get_global_glide_gscores_paninhibitor(SARS2_csv, SARS_csv, MERS_csv, out):
    "The mean of glide gscores"

    SARS2 = pd.read_csv(SARS2_csv)
    SARS = pd.read_csv(SARS_csv)
    MERS = pd.read_csv(MERS_csv)
    df = pd.concat([SARS2, SARS, MERS])
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
    mean_df.to_csv('%s'%out)



########## FINAL VALIDATION BLOCK ##########

def create_df_gscore_vs_tanimoto(files_dir, gscore_limit=-7, tan_limit=0.3, biased=False):
    """It creates a dataframe to use it for plotting in the next function:
       A scatter plot of Glide gscore (x) vs. max tanimoto (y)

       files_dir: dir where all glide csv results are located.
       limits: to store data until the specified limits"""

    specific = open('%s/specific_set.smi'%files_dir)
    specific_smiles = [line for line in specific]
    specific_mols = [mollib.Mol(smile=smile) for smile in specific_smiles]

    out = open('%s/global_sel_df_gscore_%s_tan_%s.csv'%(files_dir, gscore_limit, tan_limit), 'w')
    out.write('id,SMILE,gscore,max_tan,round\n')
    if biased:
        files = glob.glob('%s/glide_docking_biased*'%files_dir)
    else:
        files = glob.glob('%s/global_sel_glide_gscores_*'%files_dir)
    files.sort()
    for csv_file in files:
        gen_round = os.path.basename(csv_file).replace('.csv', '').split('_')[-1]
        df = pd.read_csv(csv_file)
        for index, row in df.iterrows():
            ids = row['title']
            ids = ids.split('.')[-1]
            smiles = row['SMILES']
            gscores = row['r_i_glide_gscore']
            if gscores > gscore_limit: continue
            mol1 = mollib.Mol(smile=smiles)
            maxtanimoto=0
            for mol2 in specific_mols:
                try:
                    similarity = mollib.get_MolSimilarity(mol1,mol2,fingerprint='Morgan4')
                    if similarity > maxtanimoto: maxtanimoto=similarity
                except:
                    print('similarity computation failed')
            #print(maxtanimoto)
            if maxtanimoto <= tan_limit:
                out.write('%s,%s,%s,%.4f,round_%s\n'%(ids,smiles,gscores,maxtanimoto,gen_round))
    out.close()

def plot_gscore_vs_tanimoto(csv_file, outdir, gscore_threshold=-8, tan_threshold=0.25):
    """It plots a scatter plot of Glide gscore (x) vs. max tanimoto (y)
       against specific set. It also represent the threshold lines.
       It also stores the filtered df by the specified thresholds.

       csv_file: csv obtained in the previous function.
       Thresholds: indicate where to draw the threshold lines."""

    df = pd.read_csv(csv_file)

    # store filtered df
    folder = os.path.dirname(csv_file)
    filt_df = df[(df['gscore']<=gscore_threshold) & (df['max_tan']<=tan_threshold)]
    filt_df.to_csv('%s/df_gscore_%s_tan_%s.csv'%(folder, gscore_threshold, tan_threshold), header=['id','SMILES','gscore', 'max_tan','round'],index=False)
    print('Filtered compounds:\n', filt_df)

    # plot scatter plot
    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 8))
    palette = mcp.gen_color(cmap="YlGnBu",n=7)[1:]
    palette = list(reversed(palette[:-1]))+palette[-1:]
    gdf = df.groupby('round')
    dfs = [gdf.get_group(x) for x in gdf.groups]
    rounds = df['round'].unique()
    for i, group in enumerate(dfs):
        plt.scatter(x=group['gscore'].to_list(),y=group['max_tan'].to_list(),color=palette[i],label=rounds[i], alpha=0.8, s=14,marker='X')
    plt.axvline(x=gscore_threshold,color='red',ls='--',lw='1')
    plt.axhline(y=tan_threshold,color='red',ls='--',lw='1')
    plt.xlabel('Glide gscore',fontsize=16)
    plt.ylabel('Max tanimoto against Specific Set',fontsize=16)
    plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.savefig('%s/scatter_gscore_vs_tanimoto.pdf'%outdir)


def get_mol_formal_charges(filtered_df):
    """It filter by uniqueness and calculates the formal charge
       of the molecules."""

    mols_list = []
    df = pd.read_csv(filtered_df)
    for index, row in df.iterrows():
        smiles = row['SMILES']
        mols = mollib.Mol(smile=smiles)
        mols_list.append(mols)

    # filter uniqueness
    moldb = mollib.MolDB(molList=mols_list, verbose=False)
    moldb.filter_similarity(simthreshold=1,fingerprint='Morgan4',verbose=False)
    moldb.save_tosmi('%s/unique_molecules.smi'%os.path.dirname(filtered_df))

    smi_file = open('%s/unique_molecules.smi'%os.path.dirname(filtered_df), 'r')
    dfs = []
    for line in smi_file:
        line = line.rstrip()
        df_smi = df[df['SMILES'] == line]
        dfs.append(df_smi)
    unique = pd.concat(dfs)
    unique.drop_duplicates(subset=['SMILES'])
    outname = filtered_df.replace('.csv', '_unique.csv')
    unique.to_csv(outname)
    print(unique)
    exit()

    # calculate formal charge
    smi_file = open('%s/unique_molecules.smi'%os.path.dirname(filtered_df), 'r')
    mol_charge = {}
    charges = []
    for line in smi_file:
        m = Chem.MolFromSmiles(line)
        charge = Chem.GetFormalCharge(m)
        charges.append(charge)
        mol_charge[line] = charge
    count_dict = dict(Counter(charges).items())
    print(count_dict)
    print(len(mol_charge))


def _map_PELE_glide(glide_csv, PELE_csv):
    """Create a df with all PELE and Glide results."""

    glide = pd.read_csv(glide_csv)
    PELE = pd.read_csv(PELE_csv)
    molecules = []
    for index, row in glide.iterrows():
        ids = row['id']
        rounds = row['round']
        molecule = '%s_lig_%s'%(rounds, ids)
        molecules.append(molecule)
    glide['molecule'] = molecules
    print(glide)
    print(PELE)
    df = glide.merge(PELE, on='molecule')
    print(df)
    outdir = os.path.dirname(PELE_csv)
    df.to_csv('%s/PELE_Glide_analysis.csv'%outdir)

def arreglo(glide_csv):
    df = pd.read_csv(glide_csv)
    names = []
    for index, row in df.iterrows():
         molecules = row['molecule']
         mae_names = 'conn%s'%molecules
         names.append(mae_names)
    df['mae_name'] = names
    df.to_csv('PELE_Glide_mae.csv')

def plot_gscore_vs_PELE_BFE(csv_file, outdir, BFE_threshold=-61.98, csv_ref='smallmol_specific_set/PELE_unbiased/validation_PELE_glide.csv'):
    """Scatter plot Glide gscore vs. PELE BFE."""

    #df_ref = pd.read_csv(csv_ref)
    #print(df_ref)
    df = pd.read_csv(csv_file)
    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 8))
    #plt.scatter(x=df_ref.gscore,y=df_ref.all_BFE, color='red', alpha=0.8, s=16, marker='*')
    palette = mcp.gen_color(cmap="YlGnBu",n=6)[1:]
    palette = list(reversed(palette[:-1]))+palette[-1:]
    gdf = df.groupby('round')
    dfs = [gdf.get_group(x) for x in gdf.groups]
    rounds = df['round'].unique()
    for i, group in enumerate(dfs):
        plt.scatter(x=group['gscore'].to_list(),y=group['all_BFE'].to_list(),color=palette[i],label=rounds[i], alpha=0.8, s=14,marker='X')
    plt.axhline(y=BFE_threshold,color='red',ls='--',lw='1')
    plt.xlabel('Glide docking score (kcal/mol)',fontsize=16)
    plt.ylabel('PELE BFE (kcal/mol)',fontsize=16)
    plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.savefig('%s/PELE_Glide_scatter.pdf'%outdir)

    # store filtered df
    filt_df = df[(df['all_BFE']<=BFE_threshold)]
    outname = csv_file.replace('.csv', '_filt_final.csv')
    filt_df.to_csv('%s'%(outname))
    print('FINAL COMPOUNDS:\n', filt_df)


# COMPARISON UNBIASED VS. BIASED

def plot_UMAP_biased_unbiased(results_dir):
    """It plots the UMAP of the gernated compounds at each outer round,
       for biased and unbiased generations."""

    specific = mollib.MolDB(smiDB='%s/round_1/smi_files/simp_specific_0.smi'%results_dir, verbose=False)
    biased_1 = mollib.MolDB(smiDB='%s/round_1/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    biased_2 = mollib.MolDB(smiDB='%s/round_2_biased_new/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    biased_3 = mollib.MolDB(smiDB='%s/round_3_biased_new/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    biased_4 = mollib.MolDB(smiDB='%s/round_4_biased_new/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    unbiased_1 = mollib.MolDB(smiDB='%s/round_1/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    unbiased_2 = mollib.MolDB(smiDB='%s/round_2/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    unbiased_3 = mollib.MolDB(smiDB='%s/round_3/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)
    unbiased_4 = mollib.MolDB(smiDB='%s/round_4/all_trial_generated_smiles_threshold.smi'%results_dir, verbose=False)

    cool_colors = mcp.gen_color(cmap="YlGnBu", n=6)
    cool_colors = cool_colors[1:6]
    warm_colors = mcp.gen_color(cmap="YlOrRd", n=5)
    warm_colors = warm_colors[1:5]
    plot_colors = cool_colors[0], warm_colors[0], cool_colors[1], warm_colors[1], cool_colors[2], warm_colors[2], cool_colors[3], warm_colors[3], cool_colors[4]
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    linewidths = [0.0, 0.0, 0.0, 0.0, 0.0]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o"]

    plotlib.plot_UMAP(dbs=[unbiased_4, biased_4, unbiased_3, biased_3, unbiased_2, biased_2, unbiased_1, biased_1, specific], names=['outer4',  'biased_outer4', 'outer3', 'biased_outer3', 'outer2',  'biased_outer2', 'outer1', 'biased_outer1', 'specific'], output='%s/plots/UMAP_biased_unbiased'%(results_dir), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes, alphas = alphas, min_dist = 0.1, n_neighbors = 100, n_epochs = 1000, markers =  markers, figsize = (8,8), linewidths = linewidths)


# COMPARISON OF DIFFERENT GENERATION TYPES
def plot_UMAP_comparison_2dbs(db1, db2, out, titles, palette = "YnGnBu"):
    """It plots the UMAP of the generated compounds for paninhibitor design.
    It compares the ones generated with full specific set vs. the ones
    generated with selective specific set."""


    smidb1 = mollib.MolDB(smiDB=db1, verbose=False)
    smidb2 = mollib.MolDB(smiDB=db2, verbose=False)

    plot_colors = mcp.gen_color(cmap=palette, n=2)

    min_dists = [0, 0.2, 0.4]
    neighbours = [100, 150, 200, 250]

    for i in range(len(min_dists)):
        for j in range(len(neighbours)):
            plotlib.plot_UMAP(dbs=[smidb1, smidb2], names=titles, output='%s_md%s_nn%s'%(out, min_dists[i], neighbours[j]), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = [0.4, 0.4], alphas = [0.9, 0.9], min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 1000, markers = ["o", "o"], figsize = (8,6), linewidths = [0,0])

def plot_UMAP_full_selective(pan_dir, out):
    """Comparison UMAP full vs. selective specific set.
    Everything in one plot."""

    full_specific = mollib.MolDB(smiDB='%s/round_1/Mpro_initial_specificset.smi'%pan_dir, verbose=False)
    full_1 = mollib.MolDB(smiDB='%s/round_1/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    full_2 = mollib.MolDB(smiDB='%s/round_2/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    full_3 = mollib.MolDB(smiDB='%s/round_3_reducedta/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    full_4 = mollib.MolDB(smiDB='%s/round_4/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    sel_specific = mollib.MolDB(smiDB='%s/round_1_selective/initial_specific_set.smi'%pan_dir, verbose=False)
    sel_1 = mollib.MolDB(smiDB='%s/round_1_selective/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    sel_2 = mollib.MolDB(smiDB='%s/round_2_selective/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    sel_3 = mollib.MolDB(smiDB='%s/round_3_selective_reducedta/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)
    sel_4 = mollib.MolDB(smiDB='%s/round_4_selective/all_trial_generated_smiles_threshold.smi'%pan_dir, verbose=False)

    cool_colors = mcp.gen_color(cmap="YlGnBu", n=6)
    cool_colors = cool_colors[1:6]
    warm_colors = mcp.gen_color(cmap="YlOrRd", n=6)
    warm_colors = warm_colors[1:6]
    plot_colors = cool_colors[0], warm_colors[0], cool_colors[1], warm_colors[1], cool_colors[2], warm_colors[2], cool_colors[3], warm_colors[3], cool_colors[4], warm_colors[4]
    sizes = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.3]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    linewidths = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "*", "*"]

    min_dists = [0.8, 1]
    neighbours = [25, 50]

    for i in range(len(min_dists)):
         for j in range(len(neighbours)):
             plotlib.plot_UMAP(dbs=[full_4, sel_4, full_3, sel_3, full_2, sel_2, full_1, sel_1, full_specific, sel_specific], names=['full_outer4',  'sel_outer4', 'full_outer3', 'sel_outer3', 'full_outer2',  'sel_outer2', 'full_outer1', 'sel_outer1', 'full_specific', 'sel_specific'], output='%s/UMAP_full_selective_md%s_nn%s'%(out, min_dists[i], neighbours[j]), random_max = 1000, delimiter = None, fpsalg = 'Morgan4', colors = plot_colors, sizes = sizes, alphas = alphas, min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 1000, markers = markers, figsize = (8,6), linewidths = linewidths)



if __name__ == "__main__":

    #join_all_molecule_files('paninhibitor/round_4_selective', 'trial_generated_smiles_threshold')

    #create_table_gmn_counts('paninhibitor/round_4_selective', save_df=True)

    #plot_all_props_in_one('paninhibitor/round_4_selective')

    #plot_specific_sets_tSNE('full_specific_set/round_1')

    #_simplify_specific_sets('paninhibitor/round_4_selective')
    #plot_specific_simplified_sets_tSNE('paninhibitor/round_4_selective')

    #_convert_csv_to_smi_file('smallmol_specific_set/round_1/trial_1/trial_specific_smiles.csv', outdir='smallmol_specific_set/round_1/smi_files')

    #plot_tSNE_before_and_after('smallmol_specific_set/round_1')

    #plot_tSNE_for_generation_round('paninhibitor', biased=False)

    #get_best_glide_docking_pose_for_ligand('paninhibitor/glide_docking_4/docking/output_models/MERS_7eneC1_prep/MERS_7eneC1_prep_all_trial_generated_smiles_threshold.csv')

    #filter_by_glide_gscore(path_to_best_csv='smallmol_specific_set/biased_glide_5/docking/output_models/Mpro_6xbgA1_receptor/Mpro_6xbgA1_receptor_Mpro_smallmol_generated_threshold_best.csv', gscore=-6.5, outdir='smallmol_specific_set/biased_glide_5/')

    #csvs = glob.glob('paninhibitor/glide_docking_4/docking/output_models/*/*_best.csv')
    #filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir='paninhibitor/glide_docking_4/')

    #csvs = ['smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_4.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_3.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_2.csv', 'smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_1.csv','initial_set/Mpro_6xbgA1_receptor_SARS2_Mpro_inhibitors_DB_pv_best.csv']
    #csvs = ['smallmol_specific_set/plots/gscore_tanimoto_unbiased/glide_docking_totalgen.csv', 'initial_set/Mpro_6xbgA1_receptor_SARS2_Mpro_inhibitors_DB_pv_best.csv']
    csvs = ['paninhibitor/plots/glides/global_sel_glide_gscores_4.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_3.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_2.csv', 'paninhibitor/plots/glides/global_sel_glide_gscores_1.csv', 'paninhibitor/plots/glides/global_initial_gscores.csv']
    #superimpose_histograms(list_of_csvs=csvs, list_of_labels=['outer4', 'outer3', 'outer2', 'outer1', 'specific'], insert_title='Gscore of initial and resulting compounds', out='paninhibitor/plots/global_sel_hist_gscore_bins.pdf')

    #zoom_in_histogram(threshold=-6.5, out='paninhibitor/plots/global_sel_hist_gscore_zoom_bins.pdf', savefig=True)

    #plot_specific_set_evolution('smallmol_specific_set', biased=False)

    #create_df_gscore_vs_tanimoto('paninhibitor/plots/glides', gscore_limit=-5.9, tan_limit=0.5)

    #plot_gscore_vs_tanimoto('paninhibitor/plots/gscore_tan/global_sel_df_gscore_-5.9_tan_0.5.csv', outdir='paninhibitor/plots/gscore_tan', gscore_threshold=-7, tan_threshold=0.3)

    #get_mol_formal_charges('smallmol_specific_set/plots/gscore_tanimoto_biased/df_gscore_-7_tan_0.3.csv')

    #_map_PELE_glide(glide_csv='smallmol_specific_set/plots/gscore_tanimoto_unbiased/df_gscore_-7_tan_0.3_unique.csv', PELE_csv='smallmol_specific_set/candidates/PELE_unmodified/candidates_unmod.csv')
    #plot_gscore_vs_PELE_BFE('smallmol_specific_set/PELE_unbiased/PELE_Glide_analysis.csv', outdir='smallmol_specific_set/plots/')

    #plot_UMAP_biased_unbiased('smallmol_specific_set')

    #arreglo('smallmol_specific_set/candidates/PELE_unmodified/PELE_Glide_analysis.csv')

    #_get_global_glide_gscores_paninhibitor(SARS2_csv = 'paninhibitor/plots/glides/SARS2_selective_glide_docking_4.csv', SARS_csv = 'paninhibitor/plots/glides/SARS_selective_glide_docking_4.csv', MERS_csv = 'paninhibitor/plots/glides/MERS_selective_glide_docking_4.csv', out = 'paninhibitor/plots/glides/global_sel_glide_gscores_4.csv')

    plot_UMAP_comparison_2dbs(db1='paninhibitor/all_generated_selective.smi', db2='paninhibitor/all_generated_1target.smi', out='paninhibitor/plots/comp_1targ_3targ/UMAP_generated_sel_1targ_3targ', titles = ["1target", "3targets"], palette = "winter")

    #plot_UMAP_full_selective(pan_dir='paninhibitor', out='paninhibitor/plots')


