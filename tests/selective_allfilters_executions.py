import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from gensim_analysis import *
from glide_analysis import *
from final_selection import *
from mycolorpy import colorlist as mcp

if __name__ == "__main__":
    
    n = 'gensim_mt_sel'
    resdir = '/home/cactus/julia/gensim/selective_allfilters_pretrained'
    
    # GENSIM ANALYSIS
    """
    outer_round = 13
    resdir = resdir + '/outer_%d'%outer_round
    
    get_all_generated_molecules(results_dir=resdir, outname=n)
    
    create_table_gm_counts(results_dir=resdir, outname=n, save_df=True)
    plot_all_props_in_one(results_dir=resdir, save_fig=True)
    
    convert_csv_to_sdf_file(csv_to_convert='%s/all_generated_molecules.csv'%resdir, outdir=resdir)
    remove_duplicates_from_sdf(sdf_file='%s/all_generated_molecules.sdf'%resdir)
    """
    
    # RUN GLIDE DOCKING
    """
    glide_round = 13
    create_glide_docking_folder(destination_path='%s/glide_%d'%(resdir, glide_round),\
                                template_path='/home/cactus/julia/gensim/ALGenAnalysis/templates/glide_Mpro_multitarget',\
                                ligands_file_path='%s/outer_%d/all_generated_molecules_unique.sdf'%(resdir, glide_round))
    create_glide_run_script(destination_path='.',\
                            glide_files_path='%s/glide_%d'%(resdir,glide_round))
    """
    
    # GLIDE ANALYSIS
    """
    glide_round = 13
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS2_7rnwA1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS_2gx4A1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/MERS_7eneC1.csv'%(resdir, glide_round))
    
    csvs = ['%s/glide_%d/docking/SARS2_7rnwA1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/SARS_2gx4A1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/MERS_7eneC1_best.csv'%(resdir, glide_round)]
    
    resdir = resdir + '/outer_%d'%glide_round
    filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir=resdir.replace(str(glide_round), str(glide_round+1)), gscore_global=-8.3, gscore_individual=-7.8)
    """
    
    # COMPARE ALLFILTERS VS. NOCATALOG GENERATIONS
    """
    sdf_list = ['/home/cactus/julia/gensim/selective_allfilters/outer_1/all_generated_molecules_unique.sdf',\
                '/home/cactus/julia/gensim/selective_nocatalog/outer_1/all_generated_molecules_unique.sdf']
    names = ['generated_allfilters', 'generated_nocatalog']
    sizes = [0.3, 0.3]
    alphas = [0.6, 0.6]
    markers = ["o", "o"]
    
    plot_modbs_tSNE_or_UMAP(list_of_sdfs=sdf_list, list_of_names=names, outdir='/home/cactus/julia/gensim/selective_allfilters', outname='UMAP_allfilters_nocatalog',\
                            sizes=sizes, alphas=alphas, markers=markers, ptype='UMAP')
    """
    
    # COMPARE ALLFILTERS VS. NOCATALOG AFTER GLIDE THRESHOLD
    """
    smi_list = ['/home/cactus/julia/gensim/selective_allfilters/outer_2/specific_set.smi',\
                '/home/cactus/julia/gensim/selective_nocatalog/outer_2/specific_set.smi']
    names = ['allfilters', 'nocatalog']
    sizes = [0.3, 0.3]
    alphas = [0.6, 0.6]
    markers = ["o", "o"]
    plot_UMAP(list_smis=smi_list, list_names=names, outdir='/home/cactus/julia/gensim/selective_allfilters', outname='UMAP_alfilters_nocatalog_glide',\
              sizes=sizes, alphas=alphas, markers=markers)
    """
    
    # get global gscores
    """
    for i in range(11,14):
        csvs = ['%s/glide_%s/docking/SARS2_7rnwA1_best.csv'%(resdir, i),\
                '%s/glide_%s/docking/SARS_2gx4A1_best.csv'%(resdir, i),\
                '%s/glide_%s/docking/MERS_7eneC1_best.csv'%(resdir, i)]
        get_mean_glide_gscores(list_of_csv=csvs, out='%s/glide_%s/docking/global_glide_best.csv'%(resdir,i))
        print('Glide %s done.'%i)
    """
    
    # PLOT HISTOGRAMS OF GLIDE DOCKING SCORES
    """
    virus = 'SARS2'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = '7rnwA1'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
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
                    '%s/glide_10/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_11/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_12/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_13/docking/%s_%s_best.csv'%(resdir,virus,target),]
    labels = ['initial', 'outer1', 'outer2', 'outer3', 'outer4', 'outer5', 'outer6', 'outer7', 'outer8', 'outer9', 'outer10',\
            'outer11', 'outer12', 'outer13']
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores_extended.png'%(resdir,virus), savefig=True, legend_loc='upper right')
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores_zoom_extended.png'%(resdir,virus), savefig=True, legend_loc='upper left', xlim=[-10.5,-7.5], ylim=[0, 500])
    """
    
    # TABLE OF GSCORES vs. TANIMOTO + THRESHOLD COUNTS
    """
    virus = 'SARS2'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = '7rnwA1'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
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
    
    # PLOT SPECIFIC SET EVOLUTION SUPERPOSED FOR ALLFILTERS VS. NOCATALOG
    """
    superpose_specific_set_evolution(results_dir_1='/home/cactus/julia/gensim/selective_allfilters_pretrained',\
                                     results_dir_2='/home/cactus/julia/gensim/selective_nocatalog_pretrained',\
                                     gscore_values_1=[-7.5, -7.6, -7.7, -7.8, -7.9, -8.0, -8.1, -8.2, -8.2, -8.3, -8.3, -8.3, -8.3],\
                                     gscore_values_2=[-7.5, -7.6, -7.7, -7.8, -7.9, -8.0, -8.1, -8.2, -8.3, -8.3, -8.3, -8.4, -8.4, -8.4, -8.4],\
                                     outdir='%s/plots'%resdir, outname='specific_set_evolution_allfilters_nocatalog_extended')
    """
    
    # PLOT TSNEs
    # by outer loop
    """
    sdf_list = glob.glob('%s/outer_?/specific_outer_?.smi'%resdir)
    sdf_list.sort()
    sdf_list2 = glob.glob('%s/outer_??/specific_outer_??.smi'%resdir)
    sdf_list2.sort()
    sdf_list += sdf_list2
    sdf_list.insert(0, '%s/sel_init_spec_set.smi'%resdir)   
    simplify_specific_sets_smi(list_spec_set=sdf_list)
    """
    
    """
    sdf_list = glob.glob('%s/outer_?/specific_outer_?_simple.smi'%resdir)
    sdf_list.sort(reverse=True)
    sdf_list2 = glob.glob('%s/outer_??/specific_outer_??_simple.smi'%resdir)
    sdf_list2.sort(reverse=True)
    sdf_list = sdf_list2 + sdf_list
    sdf_list.insert(0, '%s/outer_14/specific_set.smi'%resdir)
    sdf_list.append('%s/sel_init_spec_set.smi'%resdir)

    names = ['outer13', 'outer12', 'outer11', 'outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]
    
    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='%s/plots'%resdir, outname='UMAP_spec_sets_extended',\
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
    
    cumulative_histograms(final_csvs=['/home/cactus/julia/gensim/selective_allfilters_pretrained/%s_df_gscore_tanimoto.csv'%virus, \
                          '/home/cactus/julia/gensim/selective_nocatalog_pretrained/%s_df_gscore_tanimoto.csv'%virus],\
                          initial_csvs=['/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['Generated Regular', 'Generated Ablated', 'Initial Specific'],\
                          list_of_colors = ['lightgreen', 'cornflowerblue', 'red'], insert_title='Glide docking score',\
                          out='/home/cactus/julia/gensim/selective_allfilters_pretrained/plots/%s_cum_hist_gscores_extended.pdf'%virus,\
                          savefig=True, legend_loc='upper right', xlim=None, ylim=None)
    
    cumulative_histograms(final_csvs=['/home/cactus/julia/gensim/selective_allfilters_pretrained/%s_df_gscore_tanimoto.csv'%virus, \
                          '/home/cactus/julia/gensim/selective_nocatalog_pretrained/%s_df_gscore_tanimoto.csv'%virus],\
                          initial_csvs=['/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['Generated Regular', 'Generated Ablated', 'Initial specific'],\
                          list_of_colors = ['lightgreen', 'cornflowerblue', 'red'], insert_title='Glide docking score',\
                          out='/home/cactus/julia/gensim/selective_allfilters_pretrained/plots/%s_cum_hist_gscores_zoom_extended.pdf'%virus,\
                          savefig=True, legend_loc='upper left', xlim=[-10.5, -8], ylim=[0, 2000])
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
                        outname='%s/plots/cluster_dbscans_selective_extended'%resdir)

    plot_new_scaffolds(csv_results='%s/results.csv'%resdir,
                       smi_specific='%s/sel_init_spec_set.smi'%resdir,
                       gscore_glob_thr=-8,
                       gscore_ind_thr=-7.5,
                       tanimoto_thr=0.3,
                       similarity_thrs=[0.7, 0.6, 0.5, 0.4, 0.3],
                       outname='%s/plots/perc_newscaffolds_selective_extended'%resdir)
    """
    
    # APPLY THRESHOLDS
    """
    glob_gscores = [-7, -7.5, -8, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8, -8, -8, -8]
    df = pd.read_csv('%s/results.csv'%resdir)
    for i in range(8):
        df_filt = df[(df['global_gscore'] <= glob_gscores[i]) & (df['gscore_SARS2'] <= ind_gscores[i]) & (df['gscore_SARS'] <= ind_gscores[i]) & (df['gscore_MERS'] <= ind_gscores[i])]
        print(len(df_filt))
    """       
    
    # FILTER CSV RESULTS FOR ALEXIS
    """
    indv = -8
    glo = -8
    csv_results = '%s/results.csv'%resdir
    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= glo) & (df['gscore_SARS2'] <= indv) & (df['gscore_SARS'] <= indv) & (df['gscore_MERS'] <= indv)]
    print(df_filt)
    df_filt.to_csv('%s/results_filt_8.csv'%resdir, index=False)
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

    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/selective/plots', outname='UMAP_spec_sets_filtered',\
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
    
    df = pd.read_csv('%s/summary_selective_allfilters.csv'%resdir)
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
