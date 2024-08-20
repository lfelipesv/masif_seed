# usage: python substitute_and_create_complexes.py <target_chain>

import os
import sys

import Bio.PDB

# function to align two pdbs and save in a tmp.pdb file
def align_to_target(reference_pdb, sample_pdb):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    reference_structure = pdb_parser.get_structure("reference", reference_pdb)
    sample_structure = pdb_parser.get_structure("sample", sample_pdb)
    reference_model = reference_structure[0]
    sample_model = sample_structure[0]
    reference_atoms = []
    sample_atoms = []
    for reference_chain in reference_model:
        for reference_residue in reference_chain:
            reference_atoms.append(reference_residue['CA'])
    for sample_chain in sample_model:
        for sample_residue in sample_chain:
            sample_atoms.append(sample_residue['CA'])
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(reference_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    io.save("tmp.pdb")

# function to combine two pdbs and save in outdir
def combine_pdbs(target_pdb, synthetic_pdb, synthetic_pdb_name, target_name, outdir):
    out_name = synthetic_pdb_name.split('.')[0] + '_' + target_name + '.pdb'
    out_path = os.path.join(outdir, out_name)
    # read target pdb and synthetic pdb using biopython
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    target_structure = pdb_parser.get_structure("target", target_pdb)
    synthetic_structure = pdb_parser.get_structure("synthetic", synthetic_pdb)
    # get target model and synthetic model
    target_model = target_structure[0]
    synthetic_model = synthetic_structure[0]
    # detach chain from target_name and add to synthetic model
    # Get a list of the chains in a structure
    chains = list(target_model.get_chains())
    # Detach this chain from structure
    chains[0].detach_parent()
    # Add it onto structure1
    synthetic_model.add(chains[0])
    # save combined pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(synthetic_structure)
    io.save(out_path)

def change_chain_name(pdb_name, heavy_chain, light_chain):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("tmp", pdb_name)
    for chain in structure.get_chains():
        if chain.id == heavy_chain:
            chain.id = 'H'
        elif chain.id == light_chain:
            chain.id = 'L'
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_name)

# Get target name and pdb file
target_name = sys.argv[1]
target_pdb_path = os.path.join('targets', target_name, 'out_peptides', target_name)
target_pdb = os.path.join(target_pdb_path, target_name + '.pdb')

# Path to the designs and original synthetic structures
design_path = os.path.join('targets', target_name, 'out_peptides', target_name, 'designs')
design_pdbs_path = os.path.join(design_path, 'pdbs')
sabdab_path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/preprocessing_sabdab/output/ab_HL_sabdab' # hard coded!

# Create output directory if does not exists
outdir = os.path.join(design_path, 'complexes')
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# Get total number of pdb files in design_pdbs_path
total_pdbs = len([name for name in os.listdir(design_pdbs_path) if os.path.isfile(os.path.join(design_pdbs_path, name))])

# Loop through each pdb in design_pdb_path
i = 1
for pdb_file in os.listdir(design_pdbs_path):
    print('Processing pdb {}/{}'.format(i, total_pdbs))
    if pdb_file.endswith('.pdb'):
        # get name of the original sabdab pdb file: first part of the name before the underscore
        aux_name = pdb_file.split('_')[0]
        # find sabdab pdb file in sabdab_path: starts with sabdab_pdb_name
        for file in os.listdir(sabdab_path):
            if file.startswith(aux_name):
                sabdab_pdb_name = file
        # sabdab_pdb_name = pdb_file.split('_')[0] + '.pdb'
        pdb_filename = pdb_file.split('.')[0]
        # get the path to the sabdab pdb file
        sabdab_pdb_path = os.path.join(sabdab_path, sabdab_pdb_name)
        # align synthetic pdb to target pdb
        reference_pdb = os.path.join(design_pdbs_path, pdb_file)
        try:
            align_to_target(reference_pdb, sabdab_pdb_path)
            # get chain names from sabdab_pdb_name 
            # light chain is the last character of the string
            # heavy chain is the second to last character of the string
            heavy_chain = sabdab_pdb_name.split('.')[0][-2]
            light_chain = sabdab_pdb_name.split('.')[0][-1]
            # change chain names in tmp.pdb and save again to tmp.pdb using biopython
            change_chain_name('tmp.pdb', heavy_chain, light_chain)
            # combine synthetic and target pdbs
            # sabdab aligned pdb is in tmp.pdb
            combined_pdb = combine_pdbs(target_pdb, 'tmp.pdb', pdb_filename, target_name, outdir)
        except:
            print('Error in pdb {}'.format(pdb_file))
    i += 1

