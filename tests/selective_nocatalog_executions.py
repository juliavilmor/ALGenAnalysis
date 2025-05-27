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
    outer_round = 15
    resdir = resdir + '/outer_%d'%outer_round

    get_all_generated_molecules(results_dir=resdir, outname=n)

    create_table_gm_counts(results_dir=resdir, outname=n, save_df=True)
    plot_all_props_in_one(results_dir=resdir, save_fig=True)

    convert_csv_to_sdf_file(csv_to_convert='%s/all_generated_molecules.csv'%resdir, outdir=resdir)
    remove_duplicates_from_sdf(sdf_file='%s/all_generated_molecules.sdf'%resdir)
    """

    # RUN GLIDE DOCKING
    """
    glide_round = 15
    create_glide_docking_folder(destination_path='%s/glide_%d'%(resdir, glide_round),\
                                template_path='/home/cactus/julia/gensim/ALGenAnalysis/templates/glide_Mpro_multitarget',\
                                ligands_file_path='%s/outer_%d/all_generated_molecules_unique.sdf'%(resdir, glide_round))
    create_glide_run_script(destination_path='.',\
                            glide_files_path='%s/glide_%d'%(resdir, glide_round))
    """

    # GLIDE ANALYSIS
    """
    glide_round = 15
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS2_7rnwA1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/SARS_2gx4A1.csv'%(resdir, glide_round))
    get_best_glide_docking_pose(csv_file='%s/glide_%d/docking/MERS_7eneC1.csv'%(resdir, glide_round))

    csvs = ['%s/glide_%d/docking/SARS2_7rnwA1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/SARS_2gx4A1_best.csv'%(resdir, glide_round),\
            '%s/glide_%d/docking/MERS_7eneC1_best.csv'%(resdir, glide_round)]

    resdir = resdir + '/outer_%d'%glide_round
    filter_by_glide_gscore_paninhibitors(list_of_csvs=csvs, outdir=resdir.replace(str(glide_round), str(glide_round+1)), gscore_global=-8.4, gscore_individual=-7.9)
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
    for i in range(1,16):
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
                    '%s/glide_10/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_11/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_12/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_13/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_14/docking/%s_%s_best.csv'%(resdir,virus,target),\
                    '%s/glide_15/docking/%s_%s_best.csv'%(resdir,virus,target)]
    # labels = ['initial specific', 'outer1', 'outer2', 'outer3', 'outer4', 'outer5', 'outer6', 'outer7', 'outer8', 'outer9', 'outer10',\
    #         'outer11', 'outer12', 'outer13', 'outer14', 'outer15']
    labels = ['Fixed set', 'Generated 1', 'Generated 2', 'Generated 3', 'Generated 4', 'Generated 5', 'Generated 6',\
            'Generated 7', 'Generated 8', 'Generated 9', 'Generated 10', 'Generated 11', 'Generated 12',\
            'Generated 13', 'Generated 14', 'Generated 15']
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores_extended.pdf'%(resdir,virus), savefig=True, legend_loc='upper right')
    superimpose_histograms(list_of_csvs=glide_csvs, list_of_labels=labels, insert_title='Glide docking score', out='%s/plots/%s_hist_gscores_zoom_extended.pdf'%(resdir,virus), savefig=True, legend_loc='upper left', xlim=[-10,-7.5], ylim=[0, 700])
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
    sdf_list.insert(0, '%s/outer_16/specific_set.smi'%resdir)
    sdf_list.append('%s/sel_init_spec_set.smi'%resdir)

    names = ['outer15', 'outer14', 'outer13', 'outer12', 'outer11', 'outer10', 'outer9', 'outer8', 'outer7', 'outer6', 'outer5', 'outer4', 'outer3', 'outer2', 'outer1', 'initial specific']
    sizes = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*"]

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

    # PLOT CUMULATIVE HISTOGRAMS OF GSCORES (JUST NOCATALOG INITIAL VS. GENERATED)
    """
    virus = 'global'      # Select: 'SARS2', 'SARS', 'MERS', 'global'
    target = 'glide'   # Select: '7rnwA1', '2gx4A1', '7eneC1', 'glide'
    cumulative_histograms(final_csvs=['%s/%s_df_gscore_tanimoto.csv'%(resdir,virus)],\
                          initial_csvs=['/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['generated', 'initial'],\
                          list_of_colors = ['cornflowerblue', 'red'], insert_title='Glide docking score',\
                          out='%s/plots/%s_cum_hist_gscores_extended.png'%(resdir,virus), savefig=True, legend_loc='upper right', xlim=None, ylim=None)
    cumulative_histograms(final_csvs=['%s/%s_df_gscore_tanimoto.csv'%(resdir,virus)],\
                          initial_csvs=['/home/cactus/julia/gensim/selective/glide0/docking/%s_%s_best.csv'%(virus,target)],\
                          list_of_labels = ['generated', 'initial'],\
                          list_of_colors = ['cornflowerblue', 'red'], insert_title='Glide docking score',\
                          out='%s/plots/%s_cum_hist_gscores_zoom_extended.png'%(resdir,virus), savefig=True, legend_loc='upper right', xlim=[-10,-7.5], ylim=[0,1200])
    """

    # PLOT SPECIFIC SET EVOLUTION (JUST FOR THIS ONE)
    """
    def superpose_specific_set_evolution(results_dir, gscore_values, outdir, outname):
        plt.figure()
        fig, ax = plt.subplots(figsize=(10,6), dpi=300)
        outers = glob.glob('%s/outer_?'%results_dir)
        outers.sort()
        outers2 = glob.glob('%s/outer_??'%results_dir)
        outers2.sort()
        outers = outers + outers2
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
                                     gscore_values=[-7.5, -7.6, -7.7, -7.8, -7.9, -8.0, -8.1, -8.2, -8.3, -8.3, -8.3, -8.4, -8.4, -8.4, -8.4],\
                                     outdir='%s/plots'%resdir, outname='specific_set_evolution_dscores_extended')
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

    glob_gscores = [-7, -7.5, -8, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8, -8, -8, -8]
    for i in range(8):
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
                        gscore_ind_thr=-8,
                        tanimoto_thr=1,
                        similarity_thrs=[0.7, 0.6, 0.5],
                        outname='%s/plots/cluster_dbscans_selective_extended2'%resdir)
    
    plot_new_scaffolds(csv_results='%s/results.csv'%resdir,
                       smi_specific='%s/sel_init_spec_set.smi'%resdir,
                       gscore_glob_thr=-8,
                       gscore_ind_thr=-7.5,
                       tanimoto_thr=1,
                       similarity_thrs=[0.7, 0.6, 0.5],
                       outname='%s/plots/perc_newscaffolds_selective_extended'%resdir)
    """
    
    # APPLY THRESHOLDS
    """
    from molecular_filters.molecular_filters import filter_PAINS, filter_Brenk, filter_NIH, filter_REOS
    glob_gscores = [-7, -7.5, -8, -8, -8.5, -9, -9.5, -10]
    ind_gscores = [-6.5, -7, -7.5, -8, -8, -8, -8, -8]
    df = pd.read_csv('%s/results.csv'%resdir)
    for i in range(8):
        df_filt = df[(df['global_gscore'] <= glob_gscores[i]) & (df['gscore_SARS2'] <= ind_gscores[i]) & (df['gscore_SARS'] <= ind_gscores[i]) & (df['gscore_MERS'] <= ind_gscores[i])]
        print(len(df_filt))

        # # APPLY CATALOG FILTERS
        # df = df_filt
        # df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_PAINS)
        # df.dropna(subset=['SMILES_SARS2'], inplace=True)
        # print(df.shape)
        # df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_Brenk)
        # df.dropna(subset=['SMILES_SARS2'], inplace=True)
        # print(df.shape)
        # df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_NIH)
        # df.dropna(subset=['SMILES_SARS2'], inplace=True)
        # print(df.shape)
        # df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_REOS)
        # df.dropna(subset=['SMILES_SARS2'], inplace=True)
        # print(df.shape)
        # print(df)
        # print('--------------------')
    """

    # FILTER CSV RESULTS
    """
    indv = -8
    glo = -8
    csv_results = '%s/results.csv'%resdir
    df = pd.read_csv(csv_results)
    df_filt = df[(df['global_gscore'] <= glo) & (df['gscore_SARS2'] <= indv) & (df['gscore_SARS'] <= indv) & (df['gscore_MERS'] <= indv)]
    print(df_filt)
    df_filt.to_csv('%s/results_filt_8.csv'%resdir, index=False)
    """

    # FILTER THE CSV RESULTS BY CATALOGUES
    """
    from molecular_filters.molecular_filters import filter_PAINS, filter_Brenk, filter_NIH, filter_REOS
    df = pd.read_csv('%s/results_filt_8.csv'%resdir)
    df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_PAINS)
    df.dropna(subset=['SMILES_SARS2'], inplace=True)
    print(df.shape)
    df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_Brenk)
    df.dropna(subset=['SMILES_SARS2'], inplace=True)
    print(df.shape)
    df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_NIH)
    df.dropna(subset=['SMILES_SARS2'], inplace=True)
    print(df.shape)
    df['SMILES_SARS2'] = df['SMILES_SARS2'].apply(filter_REOS)
    df.dropna(subset=['SMILES_SARS2'], inplace=True)
    print(df.shape)
    df.to_csv('%s/results_filt_8_catalogues.csv'%resdir, index=False)
    """

    # PLOT UMAP WITH FINAL SELECTED MOLECULES
    """
    sdf_list = glob.glob('%s/outer_?/specific_outer_?_simple.smi'%resdir)
    sdf_list.sort(reverse=True)
    sdf_list2 = glob.glob('%s/outer_??/specific_outer_??_simple.smi'%resdir)
    sdf_list2.sort(reverse=True)
    sdf_list = sdf_list2 + sdf_list
    sdf_list.insert(0, '%s/outer_16/specific_set.smi'%resdir)
    sdf_list.append('%s/sel_init_spec_set.smi'%resdir)

    sdf_list.append('%s/results_filt_8_catalogues.smi'%resdir)

    names = ['Generated 15', 'Generated 14', 'Generated 13', 'Generated 12', 'Generated 11', 'Generated 10', 'Generated 9',\
            'Generated 8', 'Generated 7', 'Generated 6', 'Generated 5', 'Generated 4', 'Generated 3', 'Generated 2',\
            'Generated 1', 'Initial specific', 'Candidates']
    sizes = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    alphas = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    markers = ["o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "o", "*", "X"]
    total = len(sdf_list)
    colors = mcp.gen_color(cmap="YlGnBu", n=total+1)
    colors = colors[3:total+1]
    colors = colors + ['red'] + ['fuchsia']

    plot_UMAP(list_smis=sdf_list, list_names=names, outdir='%s/plots'%resdir, outname='UMAP_finalsel',\
              sizes=sizes, alphas=alphas, markers=markers, colors=colors)
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
    """
    df = pd.read_csv('%s/summary_selective_nocatalog_extended.csv'%resdir)
    # Transpose the dataframe
    df = df.T
    # add index as a column
    df.reset_index(inplace=True)
    # add the first row as the header and remove it
    df.columns = df.iloc[0]
    df = df.drop(0)
    df = df.rename(columns={"Unnamed: 0": "outer_round"})
    outers = df['outer_round'].tolist()
    df['outer_round'] = [x.split('_')[1] for x in outers]
    print(df)

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=(12, 8), dpi=500)

    # Stacked bar plot
    bars1 = ax1.bar(df["outer_round"], df["generated"], label="Generated", color="#A3BBAD")
    bars2 = ax1.bar(df["outer_round"], df["inner_al"], label="Chemical AL", color="#357266")
    bars3 = ax1.bar(df["outer_round"], df["outer_al"], label="Affinity AL", color="#0E3B43")

    # Annotate bars with percentages
    for i in range(1, len(df)+1):
        # Inner Active Learning percentage annotation
        ax1.text(df["outer_round"][i], df["outer_al"][i] + (df["inner_al"][i]),
                f"{df['perc_inner_al'][i]:.1f}%", ha='center', va='baseline', color="white", fontsize=10)

        # Outer Active Learning percentage annotation
        ax1.text(df["outer_round"][i], df["outer_al"][i] / 2,
                f"{df['perc_outer_al'][i]:.1f}%", ha='center', va='baseline', color="white", fontsize=10)

    # Labels for left y-axis
    ax1.set_xlabel("Affinity AL Cycle")
    ax1.set_ylabel("Number of Molecules")
    ax1.set_title("Summary of Generation")

    # Second y-axis for gscore thresholds
    ax2 = ax1.twinx()
    ax2.plot(df["outer_round"], df["glob_gsscore_thr"], label="Global", color="darkred", marker="s", linestyle="dotted")
    ax2.plot(df["outer_round"], df["ind_gscore_thr"], label="Individual", color="red", marker="s", linestyle="dotted")
    ax2.set_ylim(-10, -6)

    # Labels for right y-axis
    ax2.set_ylabel("Docking Score Threshold")

    # Legends
    ax1.legend(loc="upper center")
    ax2.legend(loc="upper right")

    # Show plot
    plt.savefig('%s/plots/summary_generation_extended.pdf'%resdir)
    """

    # PLOT THE METRICS OF GENERATION --> NO (this is wrong)
    """
    df = pd.read_csv('%s/summary_selective_nocatalog_extended_metrics.csv'%resdir, index_col=0)
    df = df.T
    #df.reset_index(inplace=True)
    df['cycle'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    df.index = df['cycle']
    df = df[['validity', 'perc_inner_al', 'uniqueness', 'perc_outer_al', 'novelty']]
    df = df.rename(columns={"validity": "Validity", "perc_inner_al": "Chemical AL",\
                        "uniqueness": "Uniqueness", "perc_outer_al": "Affinity AL",\
                        "novelty": "Novelty"})
    print(df)

    plt.figure(figsize=(10, 6), dpi=500)
    for column in df.columns:
        plt.plot(df.index, df[column], marker='o', label=column)

    plt.title('Metrics of Generation')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Percentage (%)')
    plt.legend()
    plt.savefig('%s/plots/metrics_generation.png'%resdir)
    """

    # PLOT METRICS VALIDITY, UNIQUENESS AND NOVELTY FOR ALL INNER ROUNDS
    from MolecularAnalysis import moldb, mol

    """
    generated = [7500]*40 + [3500]*140
    print(len(generated))
    print(generated)

    valid = []
    uniq = []
    unk = []
    for i in range(1,41):
        print(i)
        smiles = pd.read_csv('%s/outer_1/%s_%s/gensim_nocatalog_generated_smiles.csv'%(resdir,n,i))['smiles'].tolist()
        valid.append(len(smiles))
        all_mols = [mol.Mol(smile=x, allparamaters=True) for x in smiles]
        mols = []
        for mol1 in all_mols:
            try:
                atoms = mol1.NumAtoms
                mols.append(mol1)
            except:
                print('Error')
                continue

        mols_db = moldb.MolDB(molList=mols, verbose=False)
        mols_db.filterSimilarity(simt=1, alg='Morgan4',verbose=False)
        uniq.append(len(mols_db.dicDB))
        
        specific = pd.read_csv('%s/outer_1/%s_%s/gensim_nocatalog_specific_smiles.csv'%(resdir,n,i))['smiles'].tolist()
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
        print(len(joint.dicDB))
        print(len(specific_mols_db.dicDB))
        print(len(mols_db.dicDB))
        print(len(specific_mols_db.dicDB) + len(mols_db.dicDB))
        print(len(joint.dicDB) - len(specific_mols_db.dicDB))

        unk.append(len(joint.dicDB)-len(specific_mols_db.dicDB))

    for i in range(2,16):
        for j in range(1,11):
            print(i,j)
            
            smiles = pd.read_csv('%s/outer_%s/%s_%s/gensim_nocatalog_generated_smiles.csv'%(resdir,i,n,j))['smiles'].tolist()
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
            
            specific = pd.read_csv('%s/outer_%s/%s_%s/gensim_nocatalog_specific_smiles.csv'%(resdir,i,n,j))['smiles'].tolist()
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
    
    df = pd.DataFrame({'Generated': generated, 'Valid': valid, 'Unique': uniq})
    df.to_csv('%s/metrics_generation_inners_tmp1.csv'%resdir, index=False)
    
    df = pd.DataFrame({'Generated': generated, 'Valid': valid, 'Unique': uniq, 'Unknown': unk})
    df.to_csv('%s/metrics_generation_inners_tmp2.csv'%resdir, index=False)
    
    """
    

    # validity is valid/generated * 100
    # uniqueness is uniq/valid * 100
    # novelty is unk/uniq * 100
    """
    generated = pd.read_csv('%s/metrics_generation_inners_tmp2.csv'%resdir)['Generated'].tolist()
    valid = pd.read_csv('%s/metrics_generation_inners_tmp2.csv'%resdir)['Valid'].tolist()
    uniq = pd.read_csv('%s/metrics_generation_inners_tmp2.csv'%resdir)['Unique'].tolist()
    unk = pd.read_csv('%s/metrics_generation_inners_tmp2.csv'%resdir)['Unknown'].tolist()
    
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
    """
    
    # PLOT THE METRICS VALIDITY, UNIQUENESS AND NOVELTY FOR ALL INNER ROUNDS
    """
    df = pd.read_csv('%s/metrics_generation_inners.csv'%resdir)
    print(df)
    
    plt.figure(figsize=(15, 5), dpi=500)
    plt.plot(df.index, df['Validity'], marker='.', label='Validity')
    plt.plot(df.index, df['Uniqueness'], marker='.', label='Uniqueness')
    plt.plot(df.index, df['Novelty'], marker='.', label='Novelty')
    
    outer_sizes = [40, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
    lines = list(itertools.accumulate(outer_sizes))
    for line in lines:
        plt.axvline(line, color='black', linestyle=':', alpha=0.5)
            
    plt.title('Metrics of Generation')
    plt.xlabel('Affinity AL cycle')
    plt.ylabel('Percentage (%)')
    plt.legend()
    plt.xticks(np.array(range(0, 190, 10)))
    plt.savefig('%s/plots/metrics_generation_inners.pdf'%resdir)
    """
    
    # CALCULATE MEANS AND STD FOR VALIDITY, UNIQUENESS AND NOVELTY
    """
    df = pd.read_csv('%s/metrics_generation_inners.csv'%resdir)
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
    
    """
    
    # CLUSTER SCAFFOLDS DBSCAN BUT WITH THRESHOLDS APPLIED
    #"""
    def cluster_DBSCAN(csv_results,
                   smi_specific,
                   gscore_glob_thr_list,
                   gscore_ind_thr_list,
                   similarity_thr):

        from sklearn.cluster import DBSCAN

        df = pd.read_csv(csv_results)

        df['outer'] = df['id'].apply(lambda x: x.split('_')[0])
        df['outer'] = df['outer'].astype(int)
        total_outers = len(df['outer'].unique())

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
                df_filt = df[(df['global_gscore'] <= gscore_glob_thr_list[outer-1]) & (df['gscore_SARS2'] <= gscore_ind_thr_list[outer-1]) & (df['gscore_SARS'] <= gscore_ind_thr_list[outer-1]) & (df['gscore_MERS'] <= gscore_ind_thr_list[outer-1])]                
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
    
    def plot_cluster_DBSCAN(csv_results,
                        smi_specific,
                        gscore_glob_thr_list,
                        gscore_ind_thr_list,
                        similarity_thrs,
                        outname):
        plt.figure(figsize=(10,6), dpi=500)
        for sim_thr in similarity_thrs:
            print(sim_thr)
            num_clusters, x_values = cluster_DBSCAN(csv_results=csv_results,
                                                    smi_specific=smi_specific,
                                                    gscore_glob_thr_list=gscore_glob_thr_list,
                                                    gscore_ind_thr_list=gscore_ind_thr_list,
                                                    similarity_thr=sim_thr)
            print(num_clusters, x_values)
            plt.plot(x_values, num_clusters, label='DBSCAN eps=%.2f'%sim_thr, marker='o')
        plt.xlabel('Affinity AL Cycle')
        plt.ylabel('DBSCAN scaffolds clusters')
        plt.legend()
        plt.savefig(outname+'.pdf')
        
    plot_cluster_DBSCAN(csv_results='%s/results.csv'%resdir,
                        smi_specific='%s/sel_init_spec_set.smi'%resdir,
                        gscore_glob_thr_list=[-7.5, -7.6, -7.7, -7.8, -7.9, -8.0, -8.1, -8.2, -8.3, -8.3, -8.3, -8.4, -8.4, -8.4, -8.4],
                        gscore_ind_thr_list=[-7, -7.1, -7.2, -7.3, -7.4, -7.5, -7.6, -7.7, -7.8, -7.8, -7.8, -7.9, -7.9, -7.9, -7.9],
                        similarity_thrs=[0.7, 0.6, 0.5],
                        outname='%s/plots/cluster_dbscans_selective_extended_thrs'%resdir)
    #"""