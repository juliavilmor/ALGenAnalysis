import schrodinger.structure as struct
import pandas as pd
import glob
import os


def extract_ligand_best_poses(csv_file, maegz_files, output_dir, outname, virus):
    """This function extracts the glide ligand poses from maegz files based on the gscore in the csv file."""

    df = pd.read_csv(csv_file)
    print(df)
    writer = struct.StructureWriter('%s/%s'%(output_dir,outname))

    for index, row in df.iterrows():
        print('Processing row %d/%d: %s'%(index+1, len(df), row['id']))
        
        # Extract the outer round and the ids from the row
        outer = row['id'].split('_')[0]
        ids = row['id'].split('_')[1:]
        ids = '_'.join(ids)              # Title in csv file
        ids = ids.split('-')[0]          # just in case there are duplicates      
        gscores = row['gscore_%s'%virus] # Gscore in csv file
        
        # Find the corresponding maegz file for outer round
        maegz_file = maegz_files[int(outer)-1]
        print('\t', maegz_file)

        m = struct.StructureReader(maegz_file)
        
        match_found = False  # Track if match is found

        for i, s in enumerate(m):
            if i == 0: continue # because the first structure is the receptor
            props = s.property.keys()
            title = s.title
            title = title.split('.')[-1]                 # Title in maegz file
            title = title.split(',')[0]                  # just in case there are duplicates
            st_gscore = s.property['r_i_docking_score']  # Gscore in maegz file
            
            # If the title matches the id and the gscore matches
            #if title == ids and round(st_gscore, 5) == round(gscores, 5):
            if title == ids:   
                print('\t', title, ids)
                print('\t', round(gscores, 5), round(st_gscore, 5))
            
                # change maegz title
                s.title = row['id'].split('-')[0]

                # write the new maegz file
                writer.append(s)
                
                match_found = True
                break
            
        if not match_found:
            print('WARNING: No match found for %s (%s) in %s'%(ids, gscores, maegz_file))
            
        m.close()
    writer.close()
    print('Structures saved to: %s/%s'%(output_dir,outname))


if __name__ == '__main__':

    """ Example usage:
        $SCHRODINGER/run python3 filter_maegz.py
    """

    resdir = '/home/cactus/julia/gensim/selective_nocatalog_pretrained'
    virus = 'MERS' # Choose from: 'SARS2', 'SARS', 'MERS'
    target = '7eneC1' # Choose from: '7rnwA1', '2gx4A1', '7eneC1'
    csv_file = '%s/results_filt_8_catalogues.csv'%resdir
    maegz_files = glob.glob('%s/glide_?/docking/%s_%s_pv_best.maegz'%(resdir,virus, target))
    maegz_files.sort()
    maegz_files2 = glob.glob('%s/glide_??/docking/%s_%s_pv_best.maegz'%(resdir, virus, target))
    maegz_files2.sort()
    maegz_files = maegz_files + maegz_files2
    outdir = resdir
    outname = '%s_ligands_filtered.maegz'%virus

    extract_ligand_best_poses(csv_file, maegz_files, outdir, outname, virus)
    print('Done!')


