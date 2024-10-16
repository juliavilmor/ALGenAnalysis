import schrodinger.structure as struct
import pandas as pd
import glob
import os


def extract_ligand_best_poses(csv_file, maegz_files, output_dir, outname):
    """This function extracts the glide ligand poses from maegz files based on the gscore in the csv file."""

    df = pd.read_csv(csv_file)
    writer = struct.StructureWriter('%s/%s'%(output_dir,outname))

    for index, row in df.iterrows():
        outer = row['id'].split('_')[0]
        print(outer)
        # Find the corresponding maegz file for outer round
        maegz_file = maegz_files[int(outer)-1]
        print(maegz_file)

        m = struct.StructureReader(maegz_file)

        for s in m:
            props = s.property.keys()
            title = s.title
            title = title.split('.')[-1]    # Title in maegz file

            ids = row['id'].split('_')[1:]
            ids = '_'.join(ids)             # Title in csv file
            gscores = row['gscore_SARS2']

            if title != ids: continue
            print(title, ids)
            for prop in props:

                if prop != 'r_i_glide_gscore': continue
                st_gscore = s.property[prop]
                st_gscore = round(st_gscore, 2)
                #print(title, st_gscore)
                if st_gscore != round(gscores, 2): continue
                print(round(gscores,2), st_gscore)

                # change maegz title
                s.title = row['id']

                writer.append(s)
        m.close()
    writer.close()


if __name__ == '__main__':

    csv_file = '/home/cactus/julia/gensim/full/final_output_full_highPAINS_ADMET_mapped.csv'
    maegz_files = glob.glob('/home/cactus/julia/gensim/full/glide?/docking/SARS2_*.maegz')
    maegz_files.sort()
    maegz_files = maegz_files[1:]
    maegz_files.append('/home/cactus/julia/gensim/full/glide10/docking/SARS2_7rnwA1_docking_pv.maegz')
    outdir = '/home/cactus/julia/gensim/full/'
    outname = 'SARS2_ligands_filtered.maegz'

    extract_ligand_best_poses(csv_file, maegz_files, outdir, outname)


