import os
import argparse

if __name__ == "__main__":
    '''
    Create txt file with pdb ids from data folder
    '''

    # Parse Input Parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, help='path to data folder')
    parser.add_argument('--outfile', type=str, help='output file name')
    args = parser.parse_args()

    # Loop through each pdb files in data path folder
    # and create txt file with filename
    # one file per line
    i = 1
    length = len(os.listdir(args.path))
    with open(args.outfile, "w") as f:
        for filename in os.listdir(args.path):
            if filename.endswith(".pdb"):
                # split filename to get pdb id
                file_id = filename.split(".")[0]
                f.write(file_id + "\n")

            # print how many files renamed
            print(f'{i}/{length} files added to file list')
            i += 1