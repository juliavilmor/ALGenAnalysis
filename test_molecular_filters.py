from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from chembl_structure_pipeline import standardizer
import useful_rdkit_utils as uru
from rdkit.Chem.Draw import MolsToGridImage
import pandas as pd
from collections import Counter

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
            print('Manual SMARTS pattern found on molecule with SMILES:', smiles)
            return None
    return Chem.MolToSmiles(mol)


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
clean_smiles = df['clean_smiles'].tolist()
mol_list = [Chem.MolFromSmiles(x) for x in clean_smiles]
ids = df['id'].tolist()

img = MolsToGridImage(
        mol_list, molsPerRow=5, subImgSize=(300,300), legends=ids
)
img.save('test_mol_props.png')