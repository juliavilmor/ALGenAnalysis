from gensim_analysis import *
from glide_analysis import *
from final_selection import *

if __name__ == "__main__":
    
    n = 'gensim_mt'
    resdir = '/home/cactus/julia/gensim/full/outer10'
    
    # GENSIM ANALYSIS
    """
    get_all_generated_molecules(results_dir=resdir, outname=n)
    
    create_table_gm_counts(results_dir=resdir, outname=n, save_df=True)
    
    plot_all_props_in_one(results_dir=resdir, save_fig=True)
    
    convert_csv_to_sdf_file(csv_to_convert='%s/all_generated_molecules.csv'%resdir, outdir=resdir)
    remove_duplicates_from_sdf(sdf_file='%s/all_generated_molecules.sdf'%resdir)
    """
    
    # GLIDE ANALYSIS
    """
    get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide10/docking/SARS2_7rnwA1_docking.csv')
    get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide10/docking/SARS_2gx4A1_docking.csv')
    get_best_glide_docking_pose(csv_file='/home/cactus/julia/gensim/full/glide10/docking/MERS_7eneC1_docking.csv')
    
    csvs = ['/home/cactus/julia/gensim/full/glide10/docking/SARS2_7rnwA1_docking_best.csv',\
            '/home/cactus/julia/gensim/full/glide10/docking/SARS_2gx4A1_docking_best.csv',\
            '/home/cactus/julia/gensim/full/glide10/docking/MERS_7eneC1_docking_best.csv']
    filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir=resdir.replace('10', '11'), gscore_global=-8, gscore_individual=-7.5)
    """
    
    # get global gscores
    """
    for i in range(11):
        csvs = ['/home/cactus/julia/gensim/full/glide%s/docking/SARS2_7rnwA1_docking_best.csv'%i,\
                '/home/cactus/julia/gensim/full/glide%s/docking/SARS_2gx4A1_docking_best.csv'%i,\
                '/home/cactus/julia/gensim/full/glide%s/docking/MERS_7eneC1_docking_best.csv'%i]
        get_mean_glide_gscores(list_of_csv=csvs, out='/home/cactus/julia/gensim/full/glide%s/docking/global_glide_docking_best.csv'%i)
        print('Glide %s done.'%i)
    """
    
    # PLOT HISTOGRAMS OF GLIDE DOCKING SCORES
    """
    virus = 'SARS2'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = '7rnwA1'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    glide_csvs = ['/home/cactus/julia/gensim/full/glide0/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide1/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide2/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide3/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide4/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide5/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide6/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide7/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide8/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide9/docking/%s_%s_docking_best.csv'%(virus,target),\
                    '/home/cactus/julia/gensim/full/glide10/docking/%s_%s_docking_best.csv'%(virus,target)]
    labels = ['initial', 'outer1', 'outer2', 'outer3', 'outer4', 'outer5', 'outer6', 'outer7', 'outer8', 'outer9', 'outer10']
    superimpose_histograms(list_of_csvs=glide_csvs,  list_of_labels=labels, insert_title='Glide docking score', out='/home/cactus/julia/gensim/full/plots/%s_hist_gscores.png'%virus, savefig=True, legend_loc='upper right')
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='/home/cactus/julia/gensim/full/plots/%s_hist_gscores_zoom.png'%virus, savefig=True, legend_loc='upper left', xlim=[-10,-7.5], ylim=[0, 250])
    """
    
    # TABLE OF THRESHOLD COUNTS
    """
    virus = 'global'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = 'glide'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    create_df_gscore_vs_tanimoto(files_dir='/home/cactus/julia/gensim/full/', specific_set='/home/cactus/julia/gensim/full/full_init_spec_set.smi', virus=virus, target='%s_docking'%target)
    """
    
    """
    tables = ['/home/cactus/julia/gensim/full/SARS2_df_gscore_tanimoto.csv',\
                '/home/cactus/julia/gensim/full/SARS_df_gscore_tanimoto.csv',\
                '/home/cactus/julia/gensim/full/MERS_df_gscore_tanimoto.csv']
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    for i in range(7):
        apply_thresholds(global_csv='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.csv', individual_csvs=tables, 
                         gscore_ind=ind_gscores[i], gscore_glob=glob_gscores[i], tan_ind=1, tan_glob=1)
    """
    
    # PLOT SPECIFIC SET EVOLUTION
    """
    plot_specific_set_evolution(results_dir='/home/cactus/julia/gensim/full',\
                                outdir='/home/cactus/julia/gensim/full/plots', outname='specific_set_evolution')
    superpose_specific_set_evolution(results_dir_1='/home/cactus/julia/gensim/full',\
                                     results_dir_2='/home/cactus/julia/gensim/selective',\
                                     gscore_values_1=[-7.5, -7.6, -7.7, -7.7, -7.8, -7.9, -7.9, -7.9, -8.0, -8.0],\
                                     gscore_values_2=[-7.5, -7.6, -7.7, -7.8, -7.9, -7.9, -8.0, -8.1, -8.2, -8.2],\
                                     outdir='/home/cactus/julia/gensim/Mpro_GMN/plots', outname='specific_set_evolution_full_sel')
    """
    
    # PLOT UMAP OF SPECIFIC SETS
    """
    sdf_list = glob.glob('/home/cactus/julia/gensim/full/outer?/full_spec_set_outer?.smi')
    sdf_list.sort()
    sdf_list.append('/home/cactus/julia/gensim/full/outer10/full_spec_set_outer10.smi')
    sdf_list.insert(0, '/home/cactus/julia/gensim/full/full_init_spec_set.smi')   
    simplify_specific_sets_smi(list_spec_set=sdf_list)
    """
    
    """
    sdf_list = glob.glob('/home/cactus/julia/gensim/full/outer?/full_spec_set_outer?_simple.smi')
    sdf_list.sort(reverse=True)
    sdf_list.insert(0, '/home/cactus/julia/gensim/full/outer10/full_spec_set_outer10_simple.smi')
    sdf_list.insert(0, '/home/cactus/julia/gensim/full/outer11/specific_set.smi')
    sdf_list.append('/home/cactus/julia/gensim/full/full_init_spec_set.smi')

    names = ['outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]
    
    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/full/plots', outname='UMAP_spec_sets',\
              sizes=sizes, alphas=alphas, markers=markers)
    """
    
    # PLOT GSCORE VS. TANIMOTO
    """
    virus = 'MERS'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    plot_gscore_vs_tanimoto(csv_file='/home/cactus/julia/gensim/full/%s_df_gscore_tanimoto.csv'%virus,\
                            outdir='/home/cactus/julia/gensim/full/plots', outname='%s_scatter_gscore_tanimoto'%virus,\
                            gscore_threshold=-8, tan_threshold=0.3,\
                            save_csv=False, plot_max_lim=[-7, 0.5])
    """
    
    # FINAL SELECTION FOR PELE
    """
    tables = ['/home/cactus/julia/gensim/full/SARS2_df_gscore_tanimoto.csv',\
                '/home/cactus/julia/gensim/full/SARS_df_gscore_tanimoto.csv',\
                '/home/cactus/julia/gensim/full/MERS_df_gscore_tanimoto.csv']
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    for i in range(7):
        apply_thresholds(global_csv='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.csv', individual_csvs=tables, 
                         gscore_ind=ind_gscores[i], gscore_glob=glob_gscores[i], tan_ind=0.3, tan_glob=0.3)
    """
    
    # Count compounds under glide threshold INITIAL specific set
    """
    def apply_gscore_thresholds(global_csv, individual_csvs, gscore_ind, gscore_glob):
        
        # First, get the global compounds
        glob_df = pd.read_csv(global_csv)
        glob_filt = glob_df[(glob_df['r_i_glide_gscore'] <= gscore_glob)]
        glob_ids = glob_filt['title'].tolist()
        glob_compounds = set(glob_ids)
        
        # Then, get the individual compounds
        all_ind_comp = []
        for csv in individual_csvs:
            ind_df = pd.read_csv(csv)

            # Filter the individual compounds based on the thresholds
            ind_filt = ind_df[(ind_df['r_i_glide_gscore'] <= gscore_ind)]
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
        apply_gscore_thresholds(global_csv='/home/cactus/julia/gensim/full/glide0/docking/global_glide_docking_best.csv',\
                                individual_csvs=['/home/cactus/julia/gensim/full/glide0/docking/SARS2_7rnwA1_docking_best.csv',\
                                '/home/cactus/julia/gensim/full/glide0/docking/SARS_2gx4A1_docking_best.csv',\
                                '/home/cactus/julia/gensim/full/glide0/docking/MERS_7eneC1_docking_best.csv'],\
                                gscore_ind=ind_gscores[i], gscore_glob=glob_gscores[i])
    """
    
    # GET TABLE OF RESULTS
    """
    map_gscores_generated(csv_global='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.csv',\
                      csv_virus=['/home/cactus/julia/gensim/full/SARS2_df_gscore_tanimoto.csv',\
                      '/home/cactus/julia/gensim/full/SARS_df_gscore_tanimoto.csv',\
                      '/home/cactus/julia/gensim/full/MERS_df_gscore_tanimoto.csv'],\
                      outdir='/home/cactus/julia/gensim/full', outname='results')
    """
    
    # CLUSTER DBSCAN
    """
    plot_cluster_DBSCAN(csv_results='/home/cactus/julia/gensim/full/results.csv',
                        smi_specific='/home/cactus/julia/gensim/full/full_init_spec_set.smi',
                        gscore_glob_thr=-8,
                        gscore_ind_thr=-7.5,
                        tanimoto_thr=0.3,
                        similarity_thrs=[0.7, 0.6, 0.5, 0.4, 0.3],
                        outname='plots/cluster_dbscans_full') 

    plot_new_scaffolds(csv_results='/home/cactus/julia/gensim/full/results.csv',
                       smi_specific='/home/cactus/julia/gensim/full/full_init_spec_set.smi',
                       gscore_glob_thr=-8,
                       gscore_ind_thr=-7.5,
                       tanimoto_thr=0.3,
                       similarity_thrs=[0.7, 0.6, 0.5, 0.4, 0.3],
                       outname='plots/perc_newscaffolds_full')    
    """
    
    # APPLY THRESHOLDS
    """
    glob_gscores = [-7, -7.5, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8.5, -9, -9.5]
    df = pd.read_csv('/home/cactus/julia/gensim/full/results.csv')
    for i in range(7):
        df_filt = df[(df['global_gscore'] <= glob_gscores[i]) & (df['gscore_SARS2'] <= ind_gscores[i]) & (df['gscore_SARS'] <= ind_gscores[i]) & (df['gscore_MERS'] <= ind_gscores[i])  & (df['max_tan'] <= 0.3)]
        print(len(df_filt))
    """
    
    # FILTER CSV RESULTS FOR ALEXIS
    """
    indv = -7
    glo = -7.5
    csv_results = '/home/cactus/julia/gensim/full/results.csv'
    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= glo) & (df['gscore_SARS2'] <= indv) & (df['gscore_SARS'] <= indv) & (df['gscore_MERS'] <= indv)  & (df['max_tan'] <= 0.3)]
    print(df_filt)
    df_filt.to_csv('/home/cactus/julia/gensim/full/results_filt_full.csv', index=False)
    """
    
    # GET SMILES FOR FILTERED PAINS ADMET MOLECULES
    """
    get_smi_files_from_csv(csv_file='/home/cactus/julia/gensim/full/final_output_full_highPAINS_ADMET.csv',
                           smi_column='canon_smiles', outdir='/home/cactus/julia/gensim/full')
    """
    
    # PLOT THE FILTERED PAINS ADMET MOLECULES TO THE UMAP
    """
    sdf_list = glob.glob('/home/cactus/julia/gensim/full/outer?/full_spec_set_outer?_simple.smi')
    sdf_list.sort(reverse=True)
    sdf_list.insert(0, '/home/cactus/julia/gensim/full/outer10/full_spec_set_outer10_simple.smi')
    sdf_list.insert(0, '/home/cactus/julia/gensim/full/outer11/specific_set.smi')
    sdf_list.append('/home/cactus/julia/gensim/full/full_init_spec_set.smi')
    sdf_list.append('/home/cactus/julia/gensim/full/PAINS_ADMET/final_output_full_highPAINS_ADMET.smi')

    names = ['outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific', 'filtered']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*", "X"]
    total = len(sdf_list)
    colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    colors = colors[3:total+1]
    colors = colors + ['red'] + ['fuchsia']

    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/full/plots', outname='UMAP_spec_sets_filtered',\
              sizes=sizes, alphas=alphas, markers=markers, colors=colors)
    plot_tSNE(list_smis=sdf_list, list_names=names, outdir='/home/cactus/julia/gensim/full/plots', outname='tSNE_spec_sets_filtered',\
              sizes=sizes, alphas=alphas, markers=markers, colors=colors)
    """
    
    # MAP IDS TO THE FINAL OUTPUT PAINS ADMET CSV FILE
    """
    map_ids_filtered_PAINS_ADMET_mols('/home/cactus/julia/gensim/full/results_filt_full.csv',
                                      '/home/cactus/julia/gensim/full/final_output_full_highPAINS_ADMET.csv',
                                      '/home/cactus/julia/gensim/full/final_output_full_highPAINS_ADMET_mapped.csv')
    """
    
    # GET INTERACTION FINGERPRINTS
    from MolecularAnalysis.utils.interactions import InteractionFingerprints
    SARS2_pdbfile = '/home/cactus/julia/gensim/selective/glide1/structures/SARS2_7rnwA1_prep.pdb'
    SARS2_sdffile = '/home/cactus/julia/gensim/full/maestro_molecules/SARS2_ligands_filtered_full.sdf'
    SARS_pdbfile = '/home/cactus/julia/gensim/selective/glide1/structures/SARS_2gx4A1_prep.pdb'
    SARS_sdffile = '/home/cactus/julia/gensim/full/maestro_molecules/SARS_ligands_filtered_full.sdf'
    MERS_pdbfile = '/home/cactus/julia/gensim/selective/glide1/structures/MERS_7eneC1_prep.pdb'
    MERS_sdffile = '/home/cactus/julia/gensim/full/maestro_molecules/MERS_ligands_filtered_full.sdf'
    
    int_fps_SARS2 = InteractionFingerprints(SARS2_pdbfile, SARS2_sdffile)
    int_fps_SARS2.to_pickle('/home/cactus/julia/gensim/full/int_fingerprints/SARS2_int_fps.pkl')
    int_fps_SARS2.plot_barcode('/home/cactus/julia/gensim/full/plots/SARS2_int_fps_barcode.png')
    int_fps_SARS2.plot_clustermap('/home/cactus/julia/gensim/full/plots/SARS2_int_fps_clustermap.png')
    int_fps_SARS2.get_int_perc(plot='/home/cactus/julia/gensim/full/plots/SARS2_int_fps_perc_heatmap.png')
    int_fps_SARS2.get_int_res_perc(plot='/home/cactus/julia/gensim/full/plots/SARS2_int_fps_res_perc_heatmap.png')
    int_fps_SARS2.plot_pair_similarity_dist('/home/cactus/julia/gensim/full/plots/SARS2_int_fps_pair_similarity_dist.png')
    int_fps_SARS2.get_numdiff_interactions()
    
    int_fps_SARS = InteractionFingerprints(SARS_pdbfile, SARS_sdffile)
    int_fps_SARS.to_pickle('/home/cactus/julia/gensim/full/int_fingerprints/SARS_int_fps.pkl')
    int_fps_SARS.plot_barcode('/home/cactus/julia/gensim/full/plots/SARS_int_fps_barcode.png')
    int_fps_SARS.plot_clustermap('/home/cactus/julia/gensim/full/plots/SARS_int_fps_clustermap.png')
    int_fps_SARS.get_int_perc(plot='/home/cactus/julia/gensim/full/plots/SARS_int_fps_perc_heatmap.png')
    int_fps_SARS.get_int_res_perc(plot='/home/cactus/julia/gensim/full/plots/SARS_int_fps_res_perc_heatmap.png')
    int_fps_SARS.plot_pair_similarity_dist('/home/cactus/julia/gensim/full/plots/SARS_int_fps_pair_similarity_dist.png')
    int_fps_SARS.get_numdiff_interactions()
    
    int_fps_MERS = InteractionFingerprints(MERS_pdbfile, MERS_sdffile)
    int_fps_MERS.to_pickle('/home/cactus/julia/gensim/full/int_fingerprints/MERS_int_fps.pkl')
    int_fps_MERS.plot_barcode('/home/cactus/julia/gensim/full/plots/MERS_int_fps_barcode.png')
    int_fps_MERS.plot_clustermap('/home/cactus/julia/gensim/full/plots/MERS_int_fps_clustermap.png')
    int_fps_MERS.get_int_perc(plot='/home/cactus/julia/gensim/full/plots/MERS_int_fps_perc_heatmap.png')
    int_fps_MERS.get_int_res_perc(plot='/home/cactus/julia/gensim/full/plots/MERS_int_fps_res_perc_heatmap.png')
    int_fps_MERS.plot_pair_similarity_dist('/home/cactus/julia/gensim/full/plots/MERS_int_fps_pair_similarity_dist.png')
    int_fps_MERS.get_numdiff_interactions()