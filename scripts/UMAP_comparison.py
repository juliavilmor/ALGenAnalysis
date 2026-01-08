from MolecularAnalysis import moldb
from MolecularAnalysis.analysis import plot
from mycolorpy import colorlist as mcp
from gensim_analysis import get_smi_files_from_csv


def plot_UMAP_comparison_2dbs(db1, db2, out, titles):

     smidb1 = moldb.MolDB(smiDB=db1, verbose=False)
     smidb2 = moldb.MolDB(smiDB=db2, verbose=False)

     plot_colors = ['lightgrey', 'cornflowerblue']

     min_dists = [0.2, 0.4, 0.6, 0.8]
     neighbours = [100, 200, 300]

     for i in range(len(min_dists)):
         for j in range(len(neighbours)):
             plot.plotUMAP(dbs=[smidb1, smidb2], names=titles, output='%s_md%s_nn%s'%(out, min_dists[i], neighbours[j]),\
                           random_max = 1000, delimiter = None, colors = plot_colors, sizes = [0.4, 0.4], alphas = [0.9, 0.9],\
                           min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 1000, markers = ["o", "o"],\
                           figsize = (8,6), linewidth = 0)
             
def plot_UMAP_comparison_4dbs(db1, db2, db3, db4, out, titles):
    
    smidb1 = moldb.MolDB(smiDB=db1, verbose=False)
    smidb2 = moldb.MolDB(smiDB=db2, verbose=False)
    smidb3 = moldb.MolDB(smiDB=db3, verbose=False)
    smidb4 = moldb.MolDB(smiDB=db4, verbose=False)

    plot_colors = ['lightgreen', 'cornflowerblue', 'seagreen', 'royalblue']

    min_dists = [0.2, 0.4, 0.6]
    neighbours = [100, 200]

    for i in range(len(min_dists)):
        for j in range(len(neighbours)):
            plot.plotUMAP(dbs=[smidb1, smidb2, smidb3, smidb4], names=titles, output='%s_md%s_nn%s'%(out, min_dists[i], neighbours[j]),\
                        random_max = 1000, delimiter = None, colors = plot_colors, sizes = [0.2, 0.2, 0.3, 0.3], alphas = [0.9, 0.9, 0.9, 0.9],\
                        min_dist = min_dists[i], n_neighbors = neighbours[j], n_epochs = 1000, markers = ["o", "o", "*", "*"],\
                        figsize = (8,6), linewidth = 0)

# FULL VS. SELECTIVE
# Comparison SPECIFIC SETS ['seagreen', 'royalblue']
"""
plot_UMAP_comparison_2dbs(db1='/home/cactus/julia/gensim/full/full_init_spec_set.smi',\
                          db2='/home/cactus/julia/gensim/selective/sel_init_spec_set.smi',\
                          out='plots/UMAP_spec_set_full_selective', titles = ["FULL", "SELECTIVE"])
"""

# Comparison GENERATED MOLECULES ['lightgreen', 'cornflowerblue']
# get_smi_files_from_csv(csv_file='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.csv', outdir='/home/cactus/julia/gensim/full/')
# get_smi_files_from_csv(csv_file='/home/cactus/julia/gensim/selective/global_df_gscore_tanimoto.csv', outdir='/home/cactus/julia/gensim/selective/')
"""
plot_UMAP_comparison_2dbs(db1='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.smi',\
                          db2='/home/cactus/julia/gensim/selective/global_df_gscore_tanimoto.smi',\
                          out='plots/UMAP_generated_full_selective', titles = ["FULL", "SELECTIVE"])
"""
# Comparison all together
"""
plot_UMAP_comparison_4dbs(db1='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.smi',\
                          db2='/home/cactus/julia/gensim/selective/global_df_gscore_tanimoto.smi',\
                          db3='/home/cactus/julia/gensim/full/full_init_spec_set.smi',\
                          db4='/home/cactus/julia/gensim/selective/sel_init_spec_set.smi',\
                          out='plots/UMAP_all_full_sel', titles = ["FULL generated", "SELECTIVE generated", "FULL specific set", "SELECTIVE specific set"])
"""

# 1 TARGET VS. 3 TARGETS
# comparison FULL ['lightgrey', 'lightgreen']
# plot_UMAP_comparison_2dbs(db1='/home/cactus/julia/Mpro_GMN/paninhibitor/all_generated_1target.smi',\
#                           db2='/home/cactus/julia/gensim/full/global_df_gscore_tanimoto.smi',\
#                           out='plots/UMAP_1tar_3tarFULL', titles = ["1target", "3targets FULL"])

# comparison SELECTIVE ['lightgrey', 'cornflowerblue']
plot_UMAP_comparison_2dbs(db1='/home/cactus/julia/Mpro_GMN/paninhibitor/all_generated_1target.smi',\
                          db2='/home/cactus/julia/gensim/selective/global_df_gscore_tanimoto.smi',\
                          out='plots/UMAP_1tar_3tarSELECTIVE', titles = ["1target", "3targets SELECTIVE"])