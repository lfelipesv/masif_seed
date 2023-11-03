masif_root=$(git rev-parse --show-toplevel)/masif
masif_source=$masif_root/source
export PYTHONPATH=$PYTHONPATH:$masif_source
export masif_seed_root
export masif_source
PDB_ID=$(echo $1| cut -d"_" -f1)
CHAIN1=$(echo $1| cut -d"_" -f2)
CHAIN2=$(echo $1| cut -d"_" -f3)
python -W ignore $masif_source/data_preparation/create_file_and_copy_data.py $1
python -W ignore $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $1