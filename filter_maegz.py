import schrodinger.structure as struct
import pandas as pd
import glob
import os


def extract_ligand_best_poses(csv_file, maegz_files, output_dir, outname, virus):
    """This function extracts the glide ligand poses from maegz files based on the gscore in the csv file."""

    df = pd.read_csv(csv_file)
    writer = struct.StructureWriter('%s/%s'%(output_dir,outname))

    for index, row in df.iterrows():
        outer = row['id'].split('_')[0]
        print('Outer: ', outer)
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
            gscores = row['gscore_%s'%virus]

            if title != ids: continue
            print(title, ids)
            for prop in props:

                if prop != 'r_i_glide_gscore': continue
                st_gscore = s.property[prop]
                if round(st_gscore, 3) != round(gscores, 3): continue
                print(round(gscores, 3), round(st_gscore, 3))

                # change maegz title
                s.title = row['id']
                
                writer.append(s)
        m.close()
    writer.close()


if __name__ == '__main__':

    virus = 'MERS' # Choose from: 'SARS2', 'SARS', 'MERS'
    target = '7eneC1' # Choose from: '7rnwA1', '2gx4A1', '7eneC1'
    csv_file = '/home/cactus/julia/gensim/full/final_output_full_highPAINS_ADMET_mapped.csv'
    maegz_files = glob.glob('/home/cactus/julia/gensim/full/glide?/docking/%s_*.maegz'%virus)
    maegz_files.sort()
    maegz_files = maegz_files[1:]
    maegz_files.append('/home/cactus/julia/gensim/full/glide10/docking/%s_%s_docking_pv.maegz'%(virus, target))
    outdir = '/home/cactus/julia/gensim/full/'
    outname = '%s_ligands_filtered.maegz'%virus

    extract_ligand_best_poses(csv_file, maegz_files, outdir, outname, virus)


