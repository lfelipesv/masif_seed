#!/usr/bin/python
import Bio
from Bio.PDB import * 
import sys
import importlib
import os
import shutil

from default_config.masif_opts import masif_opts
# Local includes
from input_output.protonate import protonate

# if len(sys.argv) <= 2: 
#     print("Usage: " + sys.argv[0] + " PDBID_A_B" + " DATA_DIR")
#     print("A or B are the chains to include in this pdb.")
#     sys.exit(1)

if len(sys.argv) <= 1: 
    print("Usage: " + sys.argv[0] + " PDBID_A_B")
    print("A or B are the chains to include in this pdb.")
    sys.exit(1)

if not os.path.exists(masif_opts['raw_pdb_dir']):
    os.makedirs(masif_opts['raw_pdb_dir'])

if not os.path.exists(masif_opts['tmp_dir']):
    os.mkdir(masif_opts['tmp_dir'])

in_fields = sys.argv[1].split('_')
pdb_id = in_fields[0]
filename = sys.argv[1]
# path = sys.argv[2]
# path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/target-conditioned-design/masif_seed/synthetic_cdrs_renamed/'
# path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/target-conditioned-design/masif_seed/synthetic_renamed/'
path = '/home/lfelipesv/Desktop/Felipe/Projects/protein/target-conditioned-design/masif_seed/masif/data/masif_peptides/ab_A_sabdab/'

print(path + filename + '.pdb')

# copy sys.argv[1] from sys.argv[2] to masif_opts['raw_pdb_dir'] using shutil.copyfile
# shutil.copy(path + '/' + filename + '.pdb', masif_opts['raw_pdb_dir']) 

##### Protonate with reduce, if hydrogens included.
# - Always protonate as this is useful for charges. If necessary ignore hydrogens later.
protonated_file = masif_opts['raw_pdb_dir'] + pdb_id + ".pdb"
protonate(path + filename + '.pdb', protonated_file)
pdb_filename = protonated_file