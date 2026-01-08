import schrodinger.structure as struct
import pandas as pd
import glob
import os

def rename_molecules_maegz(maegz_file, output_name, virus):
    
    writer = struct.StructureWriter(output_name)
    
    m = struct.StructureReader(maegz_file)
    for i, s in enumerate(m):
        if i == 0: # because the first structure is the receptor
            writer.append(s)
        else:
            title = s.title
            s.title = title + '_' + virus
            writer.append(s)
        
    m.close()
    writer.close()
    
if __name__ == '__main__':

    """ Example usage:
        $SCHRODINGER/run python3 rename_maegz.py
    """
    resdir = '/home/cactus/julia/gensim/selective_allfilters_pretrained'
    virus = 'MERS' # Choose from: 'SARS2', 'SARS', 'MERS'
    maegz_file = '%s/%s_ligands_filtered_regular.maegz'%(resdir,virus)
    outname = '%s/%s_ligands_filtered_regular_virus.maegz'%(resdir,virus)

    rename_molecules_maegz(maegz_file, outname, virus)
    print('Done!')