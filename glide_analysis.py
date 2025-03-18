from MolecularAnalysis import mol, moldb
from MolecularAnalysis.analysis import plot
import glob
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import numpy as np
import seaborn as sns
from rdkit import Chem
from collections import Counter
from statistics import mean

def create_glide_docking_folder(destination_path, template_path="templates/glide_template", ligands_file_path=None):
    """
    Creates a new directory with the structure and files copied from a template directory.
    
    destination_path: Path where the new glide docking directory should be created.
    template_path: Path to the template directory containing the required structure and files.
    ligands_file_path: Path to the 'all_generated_molecules_unique.sdf' file to be copied into the ligands folder.
    """
    
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template directory '{template_path}' does not exist.")
    
    if os.path.exists(destination_path):
        raise FileExistsError(f"Destination directory '{destination_path}' already exists.")
    
    shutil.copytree(template_path, destination_path)
    print(f"Glide docking folder created at: {destination_path}")
    
    # Copy the ligands file if provided
    if ligands_file_path:
        if not os.path.exists(ligands_file_path):
            raise FileNotFoundError(f"Ligands file '{ligands_file_path}' does not exist.")
        
        ligands_dest = os.path.join(destination_path, "ligands", "all_generated_molecules_unique.sdf")
        os.makedirs(os.path.dirname(ligands_dest), exist_ok=True)
        shutil.copy(ligands_file_path, ligands_dest)
        print(f"Ligands file copied to: {ligands_dest}")
    
def create_glide_run_script(destination_path, glide_files_path):
    """
    Creates a SLURM run script in the specified destination directory.
    
    destination_path: Path where the run script should be created.
    glide_files_path: Path where the Glide docking files are located.
    """
    script_content = f"""#!/bin/bash
#SBATCH --job-name=glide_docking
#SBATCH --output=logs/glide_%j.out
#SBATCH --error=logs/glide_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --qos=bsc_ls
#SBATCH --time=48:00:00
#SBATCH --constraint=schrodinger
##SBATCH -x amd02

export LD_PRELOAD=/lib64/libcrypto.so.1.1.1k
export PATH=$PATH:/gpfs/projects/bsc72/Programs/schrodinger2024-1
export SCHRODINGER=/gpfs/projects/bsc72/Programs/schrodinger2024-1

module load gcc openmpi greasy/2.2.3

echo "Loaded modules: $(module list)"
echo "Using SCHRODINGER directory: $SCHRODINGER"

cd {glide_files_path}/ligands/
$SCHRODINGER/ligprep -retain_i -ph 7.4 -pht 0.5 -bff 16 -g -s 4 -epik -isd all_generated_molecules_unique.sdf -osd all_generated_molecules_prep.mae -HOST localhost:32
cd ..

wait

cd docking/
$SCHRODINGER/glide MERS_7eneC1.in -adjust -HOST localhost:32
$SCHRODINGER/glide SARS_2gx4A1.in -adjust -HOST localhost:32
$SCHRODINGER/glide SARS2_7rnwA1.in -adjust -HOST localhost:32
"""
    
    script_path = os.path.join(destination_path, "run_glide.sh")
    with open(script_path, "w") as script_file:
        script_file.write(script_content)
    
    print(f"Run script created at: {script_path}")
    
def create_glide_run_script_from_template(template_path, destination_path, glide_files_path):
    """
    Creates a SLURM run script based on a template file, replacing placeholders with actual values.
    
    template_path: Path to the template run script.
    destination_path: Path where the modified run script should be created.
    glide_files_path: Path where the Glide docking files are located.
    """
    
    # Read the template file
    with open(template_path, "r") as template_file:
        script_content = template_file.read()
    
    # Replace placeholders with actual values
    script_content = script_content.replace("{GLIDE_FILES_PATH}", glide_files_path)
    
    # Define the output script path
    script_path = os.path.join(destination_path, "run_glide.sh")
    
    # Write the modified script to the destination path
    with open(script_path, "w") as script_file:
        script_file.write(script_content)
    
    print(f"Run script created at: {script_path}")

def get_best_glide_docking_pose(csv_file, to_csv=True):
    """ I returns a df with the best poses for each ligand,
        based on the glide gscore """
 
    df = pd.read_csv(csv_file, on_bad_lines='skip')
    df = df[df['SMILES'] != 'invalid_structure']
    idx = df.groupby('title')['r_i_docking_score'].idxmin()
    df_best = df.loc[idx]
    print('The best pose for each compound has been obtained. Total: %s compounds'%len(df_best))
    
    if to_csv:
        df_best.to_csv(csv_file.replace('.csv', '_best.csv'), index=False)
 
def plot_glide_scores(path_to_best_csv, insert_title, outdir, savefig=True):
    """Histogram of gscores."""

    df = pd.read_csv(path_to_best_csv)
    plt.figure(figsize=(10, 8), dpi=200)
    sns.histplot(data=df, x='r_i_docking_score').set(title=insert_title)
    filename = os.path.basename(path_to_best_csv).replace('.csv', '_gscores.png')
    if savefig:
        plt.savefig('%s/%s'%(outdir, filename))
        
def superimpose_histograms(list_of_csvs, list_of_labels, insert_title, out, savefig=True, legend_loc='upper right', xlim=None, ylim=None):
    """Superimpose histograms to compare."""

    #plt.figure(figsize=(8, 6), dpi=200)
    plt.figure(figsize=(6, 6), dpi=200)
    total = len(list_of_csvs)
    colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    colors = colors[2:total+1]
    colors = ['red'] + colors
    for i, data in enumerate(list_of_csvs):
        df = pd.read_csv(data)
        df = df[df['r_i_docking_score'] != 10000]
        average = mean(df['r_i_docking_score'].tolist())
        sns.histplot(data=df, x=df['r_i_docking_score'].astype(float), color=colors[i], label=list_of_labels[i], element='step', fill=False, bins=70)
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
        mean_gscore = df.loc[df['title'] == lig]['r_i_docking_score'].mean()
        smile = df.loc[df['title'] == lig]['SMILES'].tolist()
        smile = smile[0] # some smiles differ from quirality. I selected the first one (corresponding to SARS2)
        means.append(mean_gscore)
        smiles.append(smile)
    mean_df = pd.DataFrame(list(zip(ligs, smiles, means)), columns =['title',  'SMILES', 'r_i_docking_score'])
    print(mean_df)
    mean_df.to_csv('%s'%out, index=False)

def filter_by_glide_docking_score(path_to_best_csv, gscore, outdir):
    """ It filters the obtained molecules by glide score.
        As a result, it gives you a csv file and a smi file."""

    df = pd.read_csv(path_to_best_csv)
    df_filt = df[df['r_i_docking_score'] < gscore]
    print("From %s molecules, %s were removed.\nThe new set contains %s molecules."%(len(df), (len(df)-len(df_filt)), len(df_filt)))

    # save filtered dataframe
    file_name = path_to_best_csv.replace('_best.csv', '_gscore_filtered.csv')
    df_filt.to_csv(file_name)

    # save smiles for the next round
    smiles = df_filt['SMILES'].tolist()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
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
        df = df[df['r_i_docking_score'] != 10000]
        df = df.drop_duplicates(subset='title')
        df['virus'] = virus
        df['receptor'] = receptor
        list_dfs.append(df)
        gscores = df['r_i_docking_score'].tolist()
        #plt.hist(gscores, bins=15, alpha=0.5)
        #plt.axvline(x = -5.9, color = 'r')
        #plt.savefig('hist_gscores.png')
    all_df = pd.concat(list_dfs)
    ligs = all_df['title'].tolist()
    ligs = set(ligs)
    mean = {}
    for lig in ligs:
        # global threshold
        mean_gscore = all_df.loc[all_df['title'] == lig]['r_i_docking_score'].mean()
        if mean_gscore < gscore_global:

            # individual threshold
            ind_gscore = all_df.loc[all_df['title'] == lig]['r_i_docking_score'] < gscore_individual
            below = ind_gscore.tolist()
            if below == [True, True, True]:
                mean[lig] = mean_gscore
                smiles = all_df.loc[all_df['title'] == lig]['SMILES'].tolist()[0]
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                smi_file = open('%s/specific_set.smi'%outdir, 'a')
                smi_file.write('%s\n'%smiles)
                smi_file.close()
    print('From %s molecules, %s were removed.\nThe new set contains %s molecules'%(len(ligs), (len(ligs) - len(mean)), len(mean)))

def cumulative_histograms(final_csvs, initial_csvs, list_of_labels, list_of_colors, insert_title, out, savefig=True, legend_loc='upper right', xlim=None, ylim=None):
    
    fig, ax = plt.subplots(figsize=(6,6))
    list_of_csvs = final_csvs + initial_csvs
    total = len(list_of_csvs)
    colors = list_of_colors
    for i, data in enumerate(list_of_csvs):
        df = pd.read_csv(data)
        df = df.rename(columns={'r_i_docking_score': 'gscore'})
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
    