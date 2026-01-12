from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from chembl_structure_pipeline import standardizer
import useful_rdkit_utils as uru
from rdkit.Chem.Draw import MolsToGridImage
import pandas as pd
from collections import Counter
import glob

def smiles_cleaning(smiles):
    try:
        remover = SaltRemover()
        mol = Chem.MolFromSmiles(smiles)
        mol = remover.StripMol(mol, dontRemoveEverything=True)
        mol = rdMolStandardize.Cleanup(mol)
        # mol = standardizer.standardize_mol(mol) # Finally, do NOT apply standardizer!!
        return Chem.MolToSmiles(mol)
    except:
        print('Error on molecule with SMILES:', smiles)
        return None

def filter_PAINS(smiles):
    mol = Chem.MolFromSmiles(smiles)
    params_pains = FilterCatalogParams()
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog_pains = FilterCatalog(params_pains)
    flag = catalog_pains.HasMatch(mol) # Check if there is a matching PAINS
    if flag:
        description = [entry.GetDescription() for entry in catalog_pains.GetMatches(mol)]
        print('PAINS found on molecule with SMILES:', smiles, ' --> ', description)
        return None
    else:
        return Chem.MolToSmiles(mol)

def filter_Brenk(smiles):
    mol = Chem.MolFromSmiles(smiles)
    params_unwanted = FilterCatalogParams()
    params_unwanted.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog_unwanted = FilterCatalog(params_unwanted)
    flag = catalog_unwanted.HasMatch(mol) # Check if there is a matching unwanted substructure
    if flag:
        description = [entry.GetDescription() for entry in catalog_unwanted.GetMatches(mol)]
        print('Brenk found on molecule with SMILES:', smiles, ' --> ', description)
        return None
    else:
        return Chem.MolToSmiles(mol)
    
def filter_NIH(smiles):
    mol = Chem.MolFromSmiles(smiles)
    params_nih = FilterCatalogParams()
    params_nih.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    catalog_nih = FilterCatalog(params_nih)
    flag = catalog_nih.HasMatch(mol) # Check if there is a matching NIH
    if flag:
        description = [entry.GetDescription() for entry in catalog_nih.GetMatches(mol)]
        print('NIH found on molecule with SMILES:', smiles, ' --> ', description)
        return None
    else:
        return Chem.MolToSmiles(mol)
    
def filter_REOS(smiles):
    """Adds Pat Walter's Chembl Filters.
        The filters are based on the REOS (Rapid Elimination Of Swill).
        The Catalog CHEMBL contains all these filters:
        CHEMBL = CHEMBL_Glaxo | CHEMBL_Dundee | CHEMBL_BMS | CHEMBL_SureChEMBL |
        CHEMBL_MLSMR | CHEMBL_Inpharmatica | CHEMBL_LINT"""
        
    mol = Chem.MolFromSmiles(smiles)
    params_reos = FilterCatalogParams()
    params_reos.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL)
    catalog_reos = FilterCatalog(params_reos)
    flag = catalog_reos.HasMatch(mol) # Check if there is a matching REOS
    if flag:
        description = [entry.GetDescription() for entry in catalog_reos.GetMatches(mol)]
        print('REOS found on molecule with SMILES:', smiles, ' --> ', description)
        return None
    else:
        return Chem.MolToSmiles(mol)
    
def filter_rings_8_atoms(smiles):
    """Filter the molecules that have rings with 8 or more atoms"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:        
        if len(ring) >= 8:
            print('Ring with 8 or more atoms found on molecule with SMILES:', smiles)
            return None
    return Chem.MolToSmiles(mol)
    
def filter_3_more_fused_rings(smiles):
    """Filter the molecules that have 3 or more fused rings"""
    mol = Chem.MolFromSmiles(smiles)
    rings = mol.GetRingInfo().AtomRings()
    # Flatten the list and count occurrences, then get items with count > 1
    fused_rings_idxs =  [idx for idx, count in Counter([i for t in rings for i in t]).items() if count > 1]
    if len(fused_rings_idxs) >= 4:
        # Check if these idxs are in the same ring to be considered all fused together
        intersected_rings = []
        for i, ring1 in enumerate(rings):
            for j, ring2 in enumerate(rings):
                if i < j:
                    intersection =  set(ring1).intersection(set(ring2))
                    if len(intersection) > 1:
                        intersected_rings.append([i, j])
        intersected_rings_idxs = [i for t in intersected_rings for i in t]
        unique = set(intersected_rings_idxs)
        if len(unique) < len(intersected_rings_idxs):
            print('3 or more fused rings found on molecule with SMILES:', smiles)
            return None
        else:
            return Chem.MolToSmiles(mol)
    else:
        return Chem.MolToSmiles(mol)  

def filter_3_4_rings_fused(smiles):
    """Filter the molecules that have 3 or 4-atom rings fused with other rings"""
    mol = Chem.MolFromSmiles(smiles)
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 3 or len(ring) == 4:
            other_rings = [r for r in rings if r != ring]
            other_rings = set([i for t in other_rings for i in t])
            fused_rings_idxs = set(ring).intersection(other_rings)
            if len(fused_rings_idxs) > 1:
                print('Fused 3-4-atom rings found on molecule with SMILES:', smiles)
                return None
    return Chem.MolToSmiles(mol)

def filter_4_ring_with_substitutions(smiles):
    """Filter the molecules that have 4-atom rings with more than 2 substitutions"""
    mol = Chem.MolFromSmiles(smiles)
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 4:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            bonds = [len(atom.GetNeighbors()) for atom in atoms]
            if sum([1 for b in bonds if b > 2]) > 2:
                print('4-atom ring with more than 2 substitutions found on molecule with SMILES:', smiles)
                return None
    return Chem.MolToSmiles(mol)
            
def filter_more_than_2_7rings(smiles):
    """Filter the molecules that have 2 or more 7-membered rings"""
    mol = Chem.MolFromSmiles(smiles)
    rings = mol.GetRingInfo().AtomRings()
    count = 0
    for ring in rings:
        if len(ring) >= 7:
            count += 1
    if count >= 2:
        return None
    return Chem.MolToSmiles(mol)

def filter_consecutive_nitrogens_ring(smiles):
    """Filter the molecules that have consecutive nitrogens in a ring"""
    mol = Chem.MolFromSmiles(smiles)
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        atoms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if 'N' not in atoms: continue
        for i in range(len(ring)-1):
            if mol.GetAtomWithIdx(ring[i]).GetAtomicNum() == 7 and mol.GetAtomWithIdx(ring[i+1]).GetAtomicNum() == 7:
                print('Consecutive nitrogens in a ring found on molecule with SMILES:', smiles)
                return None
        if mol.GetAtomWithIdx(ring[0]).GetAtomicNum() == 7 and mol.GetAtomWithIdx(ring[-1]).GetAtomicNum() == 7:
            print('Consecutive oxygens in a ring found on molecule with SMILES:', smiles)
            return None
    return Chem.MolToSmiles(mol)    
    
def filter_more_than_4_chiral_centers(smiles):
    """Filter the molecules that have more than 4 chiral centers"""
    mol = Chem.MolFromSmiles(smiles)
    chiral = Chem.FindMolChiralCenters(mol)
    if len(chiral) > 4:
        print('More than 4 chiral centers found on molecule with SMILES:', smiles)
        return None
    return Chem.MolToSmiles(mol)

def filter_manual_SMARTS_pattens(smiles):
    """Filter the molecules that have manual SMARTS patterns"""
    patterns = ['[C-]']
    mol = Chem.MolFromSmiles(smiles)
    for pattern in patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            print('Manual SMARTS pattern found on molecule with SMILES:', smiles, ' --> ', pattern)
            return None
    return Chem.MolToSmiles(mol)

def is_trifurcated(
    mol,
    min_branch_size=3,
    exclude_rings=True
    ):
    """
    Returns the mol object if it contains a valid trifurcation center,
    otherwise returns None.

    Parameters
    ----------
    mol : RDKit Mol
        Molecule to analyze.
    min_branch_size : int
        Minimum branch size (including branch-start).
    exclude_rings : bool
        Whether to ignore trifurcations where the center is in a ring.
    """
    # SMARTS pattern for trifurcation (heavy atom connected to 3 branching atoms)
    TRIFURCATION_SMARTS = "[*;!H0]([$([*])])([$([*])])([$([*])])"
    TRIFURCATION_PATTERN = Chem.MolFromSmarts(TRIFURCATION_SMARTS)

    if mol is None:
        return None

    matches = mol.GetSubstructMatches(TRIFURCATION_PATTERN)
    if not matches:
        return None

    # DFS helper function to count branch size
    def count_branch_atoms(mol, start_idx, forbidden):
        visited = set(forbidden)
        stack = [start_idx]
        visited.add(start_idx)
        count = 0

        while stack:
            idx = stack.pop()
            count += 1
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in visited:
                    visited.add(nidx)
                    stack.append(nidx)

        return count

    # Evaluate each match
    for match in matches:
        center_idx = match[0]
        branch_idxs = match[1:]
        center_atom = mol.GetAtomWithIdx(center_idx)

        # Skip rings if required
        if exclude_rings and center_atom.IsInRing():
            continue

        # Center must have at least 3 heavy neighbors
        heavy_neighbors = [
            a for a in center_atom.GetNeighbors()
            if a.GetAtomicNum() != 1
        ]
        if len(heavy_neighbors) < 3:
            continue

        # Branches must be distinct
        if len(set(branch_idxs)) != 3:
            continue

        # Check branches independently
        valid = True
        for i, bidx in enumerate(branch_idxs):
            forbidden = {center_idx} | set(branch_idxs[:i] + branch_idxs[i+1:])
            size = count_branch_atoms(mol, bidx, forbidden)

            if size < min_branch_size:
                valid = False
                break

        if valid:
            return mol  # It IS trifurcated

    return None  # No valid trifurcation found

if __name__ == '__main__':
    
    molecule_file = glob.glob('/home/cactus/julia/gensim/selective_new/outer_1_test2/gensim_mt_sel_*/gensim_mt_sel_generated_smiles.csv')
    smiles_df = pd.concat([pd.read_csv(file) for file in molecule_file])

    #  COUNT NUMBER OF MOLECULES BELOW EACH FILTER (CONSECUTIVELY)
    def count_molecules_filters_consecutively(df):
        counts = {}
        
        df['smiles'] = df['smiles'].apply(filter_rings_8_atoms)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_rings_8_atoms'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_3_more_fused_rings)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_3_more_fused_rings'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_3_4_rings_fused)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_3_4_rings_fused'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_4_ring_with_substitutions)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_4_ring_with_substitutions'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_more_than_2_7rings)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_more_than_2_7rings'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_consecutive_nitrogens_ring)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_consecutive_nitrogens_ring'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_more_than_4_chiral_centers)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_more_than_4_chiral_centers'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_manual_SMARTS_pattens)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_manual_SMARTS_pattens'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_PAINS)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_PAINS'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_Brenk)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_Brenk'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_NIH)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_NIH'] = len(df)
        df['smiles'] = df['smiles'].apply(filter_REOS)
        df.dropna(subset=['smiles'], inplace=True)
        counts['filter_REOS'] = len(df)
        df = df[df['QED'] > 0.8]
        counts['QED'] = len(df)
        df = df[df['SAscore'] < 3]
        counts['SAscore'] = len(df)
        df = df[df['tan_max'] < 0.4]
        counts['tan_max'] = len(df)
        
        return counts
        
    count_molecules = count_molecules_filters_consecutively(smiles_df)
    print(count_molecules)
    count_molecules.to_csv('count_filters_consecutive.csv', index=False)    
    exit()    
        
    # plot QED vs MW
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    corr, _ = stats.pearsonr(smiles_df['QED'], smiles_df['mw'])
    sns.lmplot(x='QED', y='mw', data=smiles_df, scatter_kws={'alpha':0.4}, line_kws={'color': 'red'})
    plt.annotate(f'correlation: {corr:.2f}', xy=(0.5, 0.5), xycoords='axes fraction')
    plt.savefig('QED_vs_MW.png')
    """

    # COUNT NUMBER OF MOLECULES BELOW EACH FILTER (INDEPENDENTLY)
    """
    def count_below_threshold(df, qed_threshold, sascore_threshold, tan_max_threshold):

        qed_count = df[df['QED'] > qed_threshold].shape[0]
        sascore_count = df[df['SAscore'] < sascore_threshold].shape[0]
        tan_max_count = df[df['tan_max'] < tan_max_threshold].shape[0]

        return {
            'QED_below_threshold': qed_count,
            'SAscore_below_threshold': sascore_count,
            'tan_max_below_threshold': tan_max_count
        }
        
    count_df = count_below_threshold(smiles_df, 0.6, 3, 0.4)
    print(count_df)
    exit()
    filters = {'filter_rings_8_atoms': filter_rings_8_atoms,
            'filter_3_more_fused_rings': filter_3_more_fused_rings,
            'filter_3_4_rings_fused': filter_3_4_rings_fused,
            'filter_4_ring_with_substitutions': filter_4_ring_with_substitutions,
            'filter_more_than_2_7rings': filter_more_than_2_7rings,
            'filter_consecutive_nitrogens_ring': filter_consecutive_nitrogens_ring,
            'filter_more_than_4_chiral_centers': filter_more_than_4_chiral_centers,
            'filter_manual_SMARTS_pattens': filter_manual_SMARTS_pattens,
            'filter_PAINS': filter_PAINS,
            'filter_Brenk': filter_Brenk,
            'filter_NIH': filter_NIH,
            'filter_REOS': filter_REOS}

    def apply_functions_and_count(functions_dict):
        counts = []
        for name, function in functions_dict.items():
            smiles_df['smiles'] = smiles_df['smiles'].apply(function)
            smiles_df.dropna(subset=['smiles'], inplace=True)
            counts.append({"Function": name, "Molecules": len(smiles_df)})
        df = pd.DataFrame(counts)
        return df
        
    count_df = apply_functions_and_count(filters)
    count_df.to_csv('count_filters.csv', index=False)
    """

    # TEST HERE
    #df = pd.read_csv('../full/PAINS_ADMET/final_output_full_highPAINS_ADMET_mapped.csv') # Result: 38 --> 9
    #df = pd.read_csv('../selective/PAINS_ADMET/final_output_selective_highPAINS_ADMET_mapped.csv') # Result: 27 --> 8
    df = pd.read_csv('../full/results_filt_full.csv') # Result: 1750 --> 480 --> 170
    #df = pd.read_csv('../selective/results_filt_selective.csv') # Result: 3865 --> 1972 --> 901
    print(df)

    # Apply these filters to the dataset
    df['clean_smiles'] = df['SMILES_SARS'].apply(smiles_cleaning)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_rings_8_atoms)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_3_more_fused_rings)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_3_4_rings_fused)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_4_ring_with_substitutions)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_more_than_2_7rings)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_consecutive_nitrogens_ring)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_more_than_4_chiral_centers)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_manual_SMARTS_pattens)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_PAINS)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_Brenk)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_NIH)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    df['clean_smiles'] = df['clean_smiles'].apply(filter_REOS)
    df.dropna(subset=['clean_smiles'], inplace=True)
    print(df.shape)
    print(df)

    exit()

    # PLOT MOLECULES
    clean_smiles = df['clean_smiles'].tolist()
    mol_list = [Chem.MolFromSmiles(x) for x in clean_smiles]
    ids = df['id'].tolist()

    img = MolsToGridImage(
            mol_list, molsPerRow=5, subImgSize=(300,300), legends=ids
    )
    img.save('test_mol_props.png')
    exit()

    # Apply these filters to the specific dataset
    specific_full = open('/home/cactus/julia/gensim/selective/sel_init_spec_set.smi', 'r').read().splitlines()
    valid_smiles = []
    for smile in specific_full:
        print(smile)
        smile = smiles_cleaning(smile)
        if smile is None: continue
        smile = filter_rings_8_atoms(smile)
        if smile is None: continue
        smile = filter_3_more_fused_rings(smile)
        if smile is None: continue
        smile = filter_3_4_rings_fused(smile)
        if smile is None: continue
        smile = filter_4_ring_with_substitutions(smile)
        if smile is None: continue
        smile = filter_more_than_2_7rings(smile)
        if smile is None: continue
        smile = filter_consecutive_nitrogens_ring(smile)
        if smile is None: continue
        smile = filter_more_than_4_chiral_centers(smile)
        if smile is None: continue
        smile = filter_manual_SMARTS_pattens(smile)
        if smile is None: continue
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