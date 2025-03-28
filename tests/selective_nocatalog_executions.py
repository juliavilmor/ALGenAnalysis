import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from gensim_analysis import *
from glide_analysis import *
from final_selection import *
from mycolorpy import colorlist as mcp

if __name__ == "__main__":

    n = 'gensim_nocatalog'
    resdir = '/home/cactus/julia/gensim/selective_nocatalog_pretrained'

    # GENSIM ANALYSIS
    """
    outer_round = 10
    resdir = resdir + '/outer_%d'%outer_round
    
    get_all_generated_molecules(results_dir=resdir, outname=n)

    create_table_gm_counts(results_dir=resdir, outname=n, save_df=True)
    plot_all_props_in_one(results_dir=resdir, save_fig=True)

    convert_csv_to_sdf_file(csv_to_convert='%s/all_generated_molecules.csv'%resdir, outdir=resdir)
    remove_duplicates_from_sdf(sdf_file='%s/all_generated_molecules.sdf'%resdir)
    """
    
    # RUN GLIDE DOCKING
    """
    glide_round = 10
    create_glide_docking_folder(destination_path='%s/glide_%d'%(resdir, glide_round),\
                                template_path='/home/cactus/julia/gensim/ALGenAnalysis/templates/glide_Mpro_multitarget',\
                                ligands_file_path='%s/outer_%d/all_generated_molecules_unique.sdf'%(resdir, glide_round))
    create_glide_run_script(destination_path='.',\
                            glide_files_path='%s/glide_%d'%(resdir, glide_round))
    """
    
    # GLIDE ANALYSIS
    """
    glide_round = 10
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS2_7rnwA1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS_2gx4A1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/MERS_7eneC1.csv'%(resdir, glide_round))
    
    csvs = ['%s/glide_%d/docking/SARS2_7rnwA1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/SARS_2gx4A1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/MERS_7eneC1_best.csv'%(resdir, glide_round)]
    
    resdir = resdir + '/outer_%d'%glide_round
    filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir=resdir.replace(str(glide_round), str(glide_round+1)), gscore_global=-8.3, gscore_individual=-7.8)
    """
    
    
    # TEST: Apply catalogs to this no catalog restuls
    """
    from molecular_filters.molecular_filters import filter_PAINS, filter_Brenk, filter_NIH, filter_REOS
    generated_molecules = glob.glob('/home/cactus/julia/gensim/selective_nocatalog_pretrained/outer_*/specific_set.smi')
    print(len(generated_molecules))
    molecules = [x for sublist in [open(x).read().splitlines() for x in generated_molecules] for x in sublist]
    print(len(molecules))
    
    valid_smiles = []
    for smile in molecules:
        smile = filter_PAINS(smile)
        if smile is None: continue
        smile = filter_Brenk(smile)
        if smile is None: continue
        smile = filter_NIH(smile)
        if smile is None: continue
        smile = filter_REOS(smile)
        if smile is None: continue
        valid_smiles.append(smile)
    print('Total valid smiles: ', len(valid_smiles))
    print(valid_smiles)
    """
    
    # get global gscores
    """
    for i in range(1,11):
        csvs = ['%s/glide_%s/docking/SARS2_7rnwA1_best.csv'%(resdir,i),\
                '%s/glide_%s/docking/SARS_2gx4A1_best.csv'%(resdir,i),\
                '%s/glide_%s/docking/MERS_7eneC1_best.csv'%(resdir,i)]
        get_mean_glide_gscores(list_of_csv=csvs, out='%s/glide_%s/docking/global_glide_best.csv'%(resdir,i))
        print('Glide %s done.'%i)
    """
    
    # PLOT HISTOGRAMS OF GLIDE DOCKING SCORES
    """
    virus = 'global'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = 'glide'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    glide_csvs = ['/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target),\
                    '%s/glide_1/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_2/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_3/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_4/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_5/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_6/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_7/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_8/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_9/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_10/docking/%s_%s_best.csv'%(resdir,virus,target)]
    labels = ['initial specific', 'outer1', 'outer2', 'outer3', 'outer4', 'outer5', 'outer6', 'outer7', 'outer8', 'outer9', 'outer10']
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores.png'%(resdir,virus), savefig=True, legend_loc='upper right')
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores_zoom.png'%(resdir,virus), savefig=True, legend_loc='upper left', xlim=[-10.5,-7.5], ylim=[0, 500])
    """
    
    # TABLE OF GSCORES vs. TANIMOTO + THRESHOLD COUNTS
    """
    virus = 'MERS'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = '7eneC1'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    create_df_gscore_vs_tanimoto(files_dir=resdir, specific_set='%s/sel_init_spec_set.smi'%resdir, virus=virus, target=target)
    """
    
    """
    tables = ['%s/SARS2_df_gscore_tanimoto.csv'%resdir,\
                '%s/SARS_df_gscore_tanimoto.csv'%resdir,\
                '%s/MERS_df_gscore_tanimoto.csv'%resdir]
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    for i in range(7):
        apply_thresholds(global_csv='%s/global_df_gscore_tanimoto.csv'%resdir, individual_csvs=tables, 
                         gscore_ind=ind_gscores[i], gscore_glob=glob_gscores[i], tan_ind=0.3, tan_glob=0.3)
    """
    
    # PLOT SPECIFIC SET EVOLUTION
    """
    plot_specific_set_evolution(results_dir=resdir,\
                                outdir='%s/plots'%resdir, outname='specific_set_evolution')
    """
    
    # PLOT TSNEs
    # by outer loop
    """
    sdf_list = glob.glob('%s/outer_?/specific_outer_?.smi'%resdir)
    sdf_list.sort()
    sdf_list.append('%s/outer_10/specific_outer_10.smi'%resdir)
    sdf_list.insert(0, '%s/sel_init_spec_set.smi'%resdir)   
    simplify_specific_sets_smi(list_spec_set=sdf_list)
    """
    """
    sdf_list = glob.glob('%s/outer_?/specific_outer_?_simple.smi'%resdir)
    sdf_list.sort(reverse=True)
    sdf_list.insert(0, '%s/outer_10/specific_outer_10_simple.smi'%resdir)
    sdf_list.insert(0, '%s/outer_11/specific_set.smi'%resdir)
    sdf_list.append('%s/sel_init_spec_set.smi'%resdir)

    names = ['outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]
    
    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='%s/plots'%resdir, outname='UMAP_spec_sets',\
              sizes=sizes, alphas=alphas, markers=markers)
    """
    
    # PLOT GSCORE VS. TANIMOTO
    """
    virus = 'global'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    plot_gscore_vs_tanimoto(csv_file='/home/cactus/julia/gensim/selective/%s_df_gscore_tanimoto.csv'%virus,\
                            outdir='/home/cactus/julia/gensim/selective/plots', outname='%s_scatter_gscore_tanimoto'%virus,\
                            gscore_threshold=-8, tan_threshold=0.3,\
                            save_csv=False, plot_max_lim=[-7, 0.5])
    """
    
    # PLOT CUMULATIVE HISTOGRAMS OF GSCORES
    """
    virus = 'global'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = 'glide'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    cumulative_histograms(final_csvs=['/home/cactus/julia/gensim/full/%s_df_gscore_tanimoto.csv'%virus, \
                          '/home/cactus/julia/gensim/selective/%s_df_gscore_tanimoto.csv'%virus],\
                          initial_csvs=['/home/cactus/julia/gensim/full/glide0/docking/%s_%s_docking_best.csv'%(virus,target),\
                          '/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['generated FULL', 'generated SELECTIVE', 'initial FULL', 'initial SELECTIVE'],\
                          list_of_colors = ['lightgreen', 'cornflowerblue', 'red', 'darkred'], insert_title='Glide docking score',\
                          out='/home/cactus/julia/gensim/Mpro_GMN/plots/%s_cum_hist_gscores.png'%virus, savefig=True, legend_loc='upper right', xlim=None, ylim=None)
    cumulative_histograms(final_csvs=['/home/cactus/julia/gensim/full/%s_df_gscore_tanimoto.csv'%virus, \
                          '/home/cactus/julia/gensim/selective/%s_df_gscore_tanimoto.csv'%virus],\
                          initial_csvs=['/home/cactus/julia/gensim/full/glide0/docking/%s_%s_docking_best.csv'%(virus,target),\
                          '/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['generated FULL', 'generated SELECTIVE', 'initial FULL', 'initial SELECTIVE'],\
                          list_of_colors = ['lightgreen', 'cornflowerblue', 'red', 'darkred'], insert_title='Glide docking score',\
                          out='/home/cactus/julia/gensim/Mpro_GMN/plots/%s_cum_hist_gscores_zoom.png'%virus, savefig=True, legend_loc='upper left', xlim=[-10.5, -7.5], ylim=[0, 2500])
    """
    
    # PLOT SPECIFIC SET EVOLUTION (JUST FOR SELECTIVE)
    """
    def superpose_specific_set_evolution(results_dir, gscore_values, outdir, outname):        
        plt.figure()
        fig, ax = plt.subplots(figsize=(10,6), dpi=300)
        outers = glob.glob('%s/outer_?'%results_dir)
        outers.sort()
        outers = outers + ['%s/outer_10'%results_dir]

        sizes_specific = []
        inner_sizes = []
        for i, outer in enumerate(outers):
            table = pd.read_csv('%s/outer_%s/table_of_counts_trans.csv'%(results_dir, i+1))
            size_specific = table['specific'].tolist()
            sizes_specific.extend(size_specific)
            inner_size = len(table)
            inner_sizes.append(inner_size)
        x = list(range(1, len(sizes_specific)+1))
        plt.plot(x, sizes_specific, marker='.', color='royalblue', label='specific set FULL')
        
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
    
    superpose_specific_set_evolution(results_dir=resdir,\
                                     gscore_values=[-7.5, -7.6, -7.7, -7.8, -7.9, -8.0, -8.1, -8.2, -8.3, -8.3],\
                                     outdir='%s/plots'%resdir, outname='specific_set_evolution_dscores')
    """
    
    # Count compounds under glide threshold INITIAL specific set
    """
    def apply_gscore_thresholds(global_csv, individual_csvs, gscore_ind, gscore_glob):
        
        # First, get the global compounds
        glob_df = pd.read_csv(global_csv)
        glob_filt = glob_df[(glob_df['r_i_docking_score'] <= gscore_glob)]
        glob_ids = glob_filt['title'].tolist()
        glob_compounds = set(glob_ids)
        
        # Then, get the individual compounds
        all_ind_comp = []
        for csv in individual_csvs:
            ind_df = pd.read_csv(csv)

            # Filter the individual compounds based on the thresholds
            ind_filt = ind_df[(ind_df['r_i_docking_score'] <= gscore_ind)]
            ind_ids = ind_filt['title'].tolist()
            ind_compounds = set(ind_ids)
            all_ind_comp.append(ind_compounds)
        ind_compounds = set.intersection(*all_ind_comp)
        
        # Get the intersection of the global and individual compounds
        final_compounds = glob_compounds.intersection(ind_compounds)
        print(f'Number of compounds below thresholds: {len(final_compounds)}')
        return final_compounds
    
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    for i in range(7):
        apply_gscore_thresholds(global_csv='/home/cactus/julia/gensim/selective/glide0/docking/global_glide_best.csv',\
                                individual_csvs=['/home/cactus/julia/gensim/selective/glide0/docking/SARS2_7rnwA1_best.csv',\
                                '/home/cactus/julia/gensim/selective/glide0/docking/SARS_2gx4A1_best.csv',\
                                '/home/cactus/julia/gensim/selective/glide0/docking/MERS_7eneC1_best.csv'],\
                                gscore_ind=ind_gscores[i], gscore_glob=glob_gscores[i])
    """
    
    # GET TABLE OF RESULTS
    """
    map_gscores_generated(csv_global='%s/global_df_gscore_tanimoto.csv'%resdir,\
                            csv_virus=['%s/SARS2_df_gscore_tanimoto.csv'%resdir,\
                            '%s/SARS_df_gscore_tanimoto.csv'%resdir,\
                            '%s/MERS_df_gscore_tanimoto.csv'%resdir],\
                            outdir=resdir, outname='results')
    
    """
    
    # CLUSTER DBSCANS
    """
    plot_cluster_DBSCAN(csv_results='%s/results.csv'%resdir,
                        smi_specific='%s/sel_init_spec_set.smi'%resdir,
                        gscore_glob_thr=-8,
                        gscore_ind_thr=-7.5,
                        tanimoto_thr=0.3,
                        similarity_thrs=[0.7, 0.6, 0.5, 0.4, 0.3],
                        outname='%s/plots/cluster_dbscans_selective'%resdir)

    plot_new_scaffolds(csv_results='%s/results.csv'%resdir,
                       smi_specific='%s/sel_init_spec_set.smi'%resdir,
                       gscore_glob_thr=-8,
                       gscore_ind_thr=-7.5,
                       tanimoto_thr=0.3,
                       similarity_thrs=[0.7, 0.6, 0.5, 0.4, 0.3],
                       outname='%s/plots/perc_newscaffolds_selective'%resdir)
    """
    
    # APPLY THRESHOLDS
    """
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    df = pd.read_csv('/home/cactus/julia/gensim/selective/results.csv')
    for i in range(7):
        df_filt = df[(df['global_gscore'] <= glob_gscores[i]) & (df['gscore_SARS2'] <= ind_gscores[i]) & (df['gscore_SARS'] <= ind_gscores[i]) & (df['gscore_MERS'] <= ind_gscores[i])  & (df['max_tan'] <= 0.3)]
        print(len(df_filt))
    """       
  
    # FILTER CSV RESULTS FOR ALEXIS
    """
    indv = -7.5
    glo = -8
    csv_results = '%s/results.csv'%resdir
    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= glo) & (df['gscore_SARS2'] <= indv) & (df['gscore_SARS'] <= indv) & (df['gscore_MERS'] <= indv)  & (df['max_tan'] <= 0.3)]
    print(df_filt)
    df_filt.to_csv('%s/results_filt.csv'%resdir, index=False)
    """
    
    # GET SMILES FOR FILTERED PAINS ADMET MOLECULES
    """
    get_smi_files_from_csv(csv_file='/home/cactus/julia/gensim/selective/final_output_selective_highPAINS_ADMET.csv',
                           smi_column='canon_smiles', outdir='/home/cactus/julia/gensim/selective')
    """
    
    # PLOT THE FILTERED PAINS ADMET MOLECULES TO THE UMAP
    """
    sdf_list = glob.glob('/home/cactus/julia/gensim/selective/outer?/sel_spec_set_outer?_simple.smi')
    sdf_list.sort(reverse=True)
    sdf_list.insert(0, '/home/cactus/julia/gensim/selective/outer10/sel_spec_set_outer10_simple.smi')
    sdf_list.insert(0, '/home/cactus/julia/gensim/selective/outer11/specific_set.smi')
    sdf_list.append('/home/cactus/julia/gensim/selective/sel_init_spec_set.smi')
    sdf_list.append('/home/cactus/julia/gensim/selective/PAINS_ADMET/final_output_selective_highPAINS_ADMET.smi')

    names = ['outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific', 'filtered']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o","o", "o", "o", "o", "o", "o", "o", "o", "*", "X"]
    total = len(sdf_list)
    colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    colors = colors[3:total+1]
    colors = colors + ['red'] + ['fuchsia']

    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/selective/plots', outname='UMAP_spec_sets_filtered_poster',\
              sizes=sizes, alphas=alphas, markers=markers, colors=colors)
    plot_tSNE(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/selective/plots', outname='tSNE_spec_sets_filtered',\
                sizes=sizes, alphas=alphas, markers=markers, colors=colors)
    """
    
    # MAP IDS TO THE FINAL OUTPUT PAINS ADMET CSV FILE
    """
    map_ids_filtered_PAINS_ADMET_mols('/home/cactus/julia/gensim/selective/results_filt_selective.csv',
                                      '/home/cactus/julia/gensim/selective/final_output_selective_highPAINS_ADMET.csv',
                                      '/home/cactus/julia/gensim/selective/final_output_selective_highPAINS_ADMET_mapped.csv')
    """
    
    # PLOT SUMMARY OF GENERATION
    
    df = pd.read_csv('%s/summary_selective_nocatalog.csv'%resdir)
    # Transpose the dataframe
    df = df.T
    # add index as a column
    df.reset_index(inplace=True)
    # add the first row as the header and remove it
    df.columns = df.iloc[0]
    df = df.drop(0)
    df = df.rename(columns={"Unnamed: 0": "outer_round"})
    print(df)

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=(8, 10), dpi=200)

    # Stacked bar plot
    bars1 = ax1.bar(df["outer_round"], df["generated"], label="Generated", color="#A3BBAD")
    bars2 = ax1.bar(df["outer_round"], df["inner_al"], label="Inner Active Learning", color="#357266")
    bars3 = ax1.bar(df["outer_round"], df["outer_al"], label="Outer Active Learning", color="#0E3B43")

    # Annotate bars with percentages
    for i in range(1, len(df)+1):
        # Inner Active Learning percentage annotation
        ax1.text(df["outer_round"][i], df["outer_al"][i] + (df["inner_al"][i]), 
                f"{df['perc_inner_al'][i]:.1f}%", ha='center', va='baseline', color="white", fontsize=8)

        # Outer Active Learning percentage annotation
        ax1.text(df["outer_round"][i], df["outer_al"][i] / 2, 
                f"{df['perc_outer_al'][i]:.1f}%", ha='center', va='baseline', color="white", fontsize=8)

    # Labels for left y-axis
    ax1.set_xlabel("Outer Round")
    ax1.set_ylabel("Number of Molecules")
    ax1.set_title("Summary of Generation")

    # Second y-axis for gscore thresholds
    ax2 = ax1.twinx()
    ax2.plot(df["outer_round"], df["glob_gsscore_thr"], label="Global gscore thr.", color="darkred", marker="s", linestyle="dotted")
    ax2.plot(df["outer_round"], df["ind_gscore_thr"], label="Indiv. gscore thr.", color="red", marker="s", linestyle="dotted")
    ax2.set_ylim(-10, -6)
    
    # Labels for right y-axis
    ax2.set_ylabel("Docking Score Threshold")

    # Legends
    ax1.legend(loc="upper center")
    ax2.legend(loc="upper right")

    # Show plot
    plt.savefig('%s/plots/summary_generation.png'%resdir)
