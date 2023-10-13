import pandas as pd
import os
import glob

outer_dfs = glob.glob('*_df_gscore_-5.9_tan_0.5.csv')
global_dfs = glob.glob('global_*-6.5*')
dfs = outer_dfs + global_dfs

def apply_thresholds(gscore_ind, gscore_glob, tan_ind, tan_glob, sel=False):
    """Get compounds below thresholds of glide gscore and tanimoto similarity.

    Parameters:
        - gscore_ind: Individual glide gscore (kcal/mol) threshold
        - gscore_glob: Global glide_gscore (kcal/mol) threshold
        - tan_ind: Individual tanimoto similarity (against the specific set) threshold
        - tan_glob: Global tanimoto similarity (against the specific set) threshold

    NOTE: I calculated the intersection with the compound name, not the SMILE.
    SMILES can differ between different viruses and it gives us wrong result."""

    all_comp = []

    for df in dfs:
        if sel == False:
            if "sel" in df: continue
        else:
            if "sel" not in df: continue

        if "global" in df:
            glob_df = pd.read_csv(df)
            glob_filt = glob_df[(glob_df['gscore']<=gscore_glob) & (glob_df['max_tan']<=tan_glob)]
            glob_id = glob_filt['id'].tolist()
            glob_id = [x.split(':')[-1] for x in glob_id]
            glob_round = glob_filt['round'].tolist()
            glob_round = [x.split('_')[-1] for x in glob_round]
            glob_comp = []
            for i in range(len(glob_id)):
                comp = '%s_%s'%(glob_round[i], glob_id[i])
                glob_comp.append(comp)
            #glob_comp = glob_filt['SMILE'].tolist()
            glob_comp = set(glob_comp)
            print('Globals: ', len(glob_comp))
            all_comp.append(glob_comp)
        else:
            ind_df = pd.read_csv(df)
            ind_filt = ind_df[(ind_df['gscore']<=gscore_ind) & (ind_df['max_tan']<=tan_ind)]
            ind_id = ind_filt['id'].tolist()
            ind_id = [x.split(':')[-1] for x in ind_id]
            ind_round = ind_filt['round'].tolist()
            ind_round = [x.split('_')[-1] for x in ind_round]
            ind_comp = []
            for i in range(len(ind_id)):
                comp = '%s_%s'%(ind_round[i], ind_id[i])
                ind_comp.append(comp)
            #ind_comp = filt_df['SMILE'].tolist()
            ind_comp = set(ind_comp)
            if "SARS2" in df:
                SARS2_comp = ind_comp
                print('SARS2: ', len(SARS2_comp))
                all_comp.append(SARS2_comp)
            elif "SARS" in df:
                SARS_comp = ind_comp
                print('SARS: ', len(SARS_comp))
                all_comp.append(SARS_comp)
            elif "MERS" in df:
                MERS_comp = ind_comp
                print('MERS: ', len(MERS_comp))
                all_comp.append(MERS_comp)

    inter = set.intersection(*map(set,all_comp))
    print(len(inter))
    print(sorted(inter))


# TRY DIFFERENT THRESHOLDS
apply_thresholds(gscore_ind=-6.5, gscore_glob=-7, tan_ind=0.3, tan_glob=0.3, sel=False)




