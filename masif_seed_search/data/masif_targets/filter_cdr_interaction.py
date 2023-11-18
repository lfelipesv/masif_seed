# usage: python filter_cdr_interaction.py <target_chain>

import os
import sys

import Bio.PDB
import numpy as np
import pandas as pd

def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two):
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

# Get target name and pdb file
target_name = sys.argv[1]
target_pdb_path = os.path.join('targets', target_name, 'out_peptides', target_name)
target_pdb = os.path.join(target_pdb_path, target_name + '.pdb')

# Path to the designs and original synthetic structures
design_path = os.path.join('targets', target_name, 'out_peptides', target_name, 'designs')
design_pdbs_path = os.path.join(design_path, 'pdbs')
design_scores_path = os.path.join(design_path, 'scores')
synthetic_path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/external/synthetic_data/predictions_flat' # hard coded!

# Complexes directory
complexes_path = os.path.join(design_path, 'complexes')

# Create output directory if does not exists
new_path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/target-conditioned-design/search-results/gfd2_C'
outdir = os.path.join(new_path, 'filtered_complexes')
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# Get total number of pdb files in design_pdbs_path
total_pdbs = len([name for name in os.listdir(complexes_path) if os.path.isfile(os.path.join(complexes_path, name))])

# split by underscore, antigen chain name is second element of list
target_chain = target_name.split('_')[1]
heavy_chain = 'H'
light_chain = 'L'

# CDR region start and end positions for each chain (according to Chothia numbering scheme)
cdr_regions = {
    'heavy': [('CDR1', 26, 35), ('CDR2', 50, 65), ('CDR3', 95, 102)],
    'light': [('CDR1', 24, 34), ('CDR2', 50, 56), ('CDR3', 89, 97)]
}

# Loop through each pdb in design_pdb_path
i = 1
for pdb_file in os.listdir(complexes_path):
    print('Processing pdb {}/{}'.format(i, total_pdbs))
    if pdb_file.endswith('.pdb'):
        pdb_filepath = os.path.join(complexes_path, pdb_file)
        structure = Bio.PDB.PDBParser().get_structure(pdb_file, pdb_filepath)
        model = structure[0]
        heavy_dist_matrix = calc_dist_matrix(model[heavy_chain], model[target_chain])
        heavy_contact_map = heavy_dist_matrix < 12.0
        light_dist_matrix = calc_dist_matrix(model[light_chain], model[target_chain])
        light_contact_map = light_dist_matrix < 12.0
        # find CDR regions in heavy and light chains
        heavy_cdr_regions = []
        light_cdr_regions = []
        # renumber from IMGT to Chothia scheme and save in tmp.pdb
        out_data_path = 'tmp.pdb'
        cmd = f'python ImmunoPDB.py -i {pdb_filepath} -o {out_data_path}'
        exit_code = os.system(cmd)
        if exit_code != 0:
            print(f'renumbering failed for {pdb_filepath}. scheme Chothia')
            break
        # extract CDR regions from tmp.pdb
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("tmp", out_data_path)
        heavy_chain_cdrs = {name: (int(start), int(end)) for name, start, end in cdr_regions['heavy']}
        light_chain_cdrs = {name: (int(start), int(end)) for name, start, end in cdr_regions['light']}
        heavy_cdrs_mask = np.zeros(len(structure[0][heavy_chain]), dtype=bool)
        light_cdrs_mask = np.zeros(len(structure[0][light_chain]), dtype=bool)
        # fill mask with True values for CDR regions
        # notice that there is a mismatch between Chothia numbers and id's in the mask
        for chain in structure[0]:
            id_ = 0
            for residue in chain:
                if chain.id == heavy_chain:
                    for name, (start, end) in heavy_chain_cdrs.items():
                        if residue.id[1] in range(start, end + 1):
                            heavy_cdrs_mask[id_] = True
                elif chain.id == light_chain:
                    for name, (start, end) in light_chain_cdrs.items():
                        if residue.id[1] in range(start, end + 1):
                            light_cdrs_mask[id_] = True
                id_ += 1
        # check if there is any contact between CDR regions and target chain
        if np.any(heavy_contact_map[heavy_cdrs_mask, :]) and np.any(light_contact_map[light_cdrs_mask, :]):
            print(f'CDR regions of {pdb_file} are in contact with target chain')
            out_path = os.path.join(outdir, pdb_file)
            os.system(f'cp {pdb_filepath} {out_path}')
    i += 1

# In part 2 loop through filtered designs and combine their csv files
csv_files = []
pdb_files = []
for pdb_file in os.listdir(outdir):
    if pdb_file.endswith('.pdb'):
        synthetic_pdb_name_list = pdb_file.split('_')
        synthetic_pdb_name_list = synthetic_pdb_name_list[:-2]
        synthetic_pdb_name = '_'.join(synthetic_pdb_name_list)
        for csv_file in os.listdir(design_scores_path):
            if csv_file.startswith(synthetic_pdb_name):
                csv_files.append(csv_file)
                pdb_files.append(pdb_file)

# combine csv files into one csv file
with open(os.path.join(outdir, 'filtered_scores.csv'), 'w') as outfile:
    for fname in csv_files:
        with open(os.path.join(design_scores_path, fname)) as infile:
            # skip header unless it is the first file
            if fname == csv_files[0]:
                for line in infile:
                    outfile.write(line)
            else:
                # skip header
                next(infile)
                for line in infile:
                    outfile.write(line)

# Load csv file and add pdb_file column using pdb_files list
df = pd.read_csv(os.path.join(outdir, 'filtered_scores.csv'))
df['pdb_filename'] = pdb_files
df.to_csv(os.path.join(outdir, 'filtered_scores.csv'), index=False)