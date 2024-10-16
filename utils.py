import pandas as pd

def set_unique(set1):
    df1 = pd.read_csv(set1)
    print('Total: ', len(df1))
    smi1 = set(df1['smiles'].tolist())
    print('Unique: ', len(smi1))
    
    # f = open('unique_smiles.smi', 'w')
    # for smi in smi1:
    #     f.write('%s\n'%smi)
    # f.close()
        
def gensim_set_comparison(set1_csv, set2_csv):
    
    df1 = pd.read_csv(set1_csv)
    df2 = pd.read_csv(set2_csv)
    smi1 = set(df1['smiles'].tolist())
    smi2 = set(df2['smiles'].tolist())
    diff = smi1.difference(smi2)
    inter = smi1.intersection(smi2)
    
    print('Total set 1: %s'%len(df1))
    print('Total set 2: %s'%len(df2))
    print('Total difference: %s'%len(diff))
    print('Total intersection: %s\n'%len(inter))
    
    print('The difference is: %s'%diff)
    print('The intersection is: %s'%inter)
    
if __name__ == "__main__":
    
    specific_1 = '/home/cactus/julia/gensim/selective/outer1/gensim_mt_sel_1/gensim_mt_sel_specific_smiles.csv'
    gen_thres_1 = '/home/cactus/julia/gensim/selective/outer1/gensim_mt_sel_2/gensim_mt_sel_generated_smiles_threshold_1.csv'
    specific_2 = '/home/cactus/julia/gensim/selective/outer1/gensim_mt_sel_2/gensim_mt_sel_specific_smiles.csv'
    
    set_unique(set1=specific_1)
    gensim_set_comparison(set1_csv=specific_2, set2_csv=specific_1)
    
    
    