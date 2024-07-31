# This script takes in a list of 1-deme treeSeq files (input), and
# outputs an excel file with their pi values

# sys arg 1 = current_dir
# sys arg 2 = mu
# sys arg 3 = N
# sys arg 4 = r
# sys arg 5 = output_csv file

# Import modules
import msprime
import pyslim
import tskit
import numpy as np
import os
import sys
import subprocess
import csv


# Define function to return list of pi values given list of panmictic treeSeq files
def calc_pi(files):
    # files: TreeSeq files

    # List to store pi values from each TreeSeq file
    pi = [None] * len(files)

    # Iterate over files to calculate FST
    for i in range(0, len(files)):

        # Load in current file
        current_file = tskit.load(files[i])

       # Simulate neutral mutations for current file
        current_file = msprime.sim_mutations(current_file,
                                             rate=float(sys.argv[2]),
                                             model=msprime.SLiMMutationModel(
                                                 type=0),
                                             keep=False)

        #### Calculate pi ####

        # Filter test tree and calculate FST for it
        # filtered_tree = current_file.delete_sites(sites_to_rm)

        # Calculate average pairwise diversity
        pi[i] = current_file.diversity(windows=[0., 350.,  10351., 10700.])[1]

        # Return iteration completed
        print(i)

    return (pi)

# Filter entries so it only contains TreeSeq files of particular N, r
def filter_files_r(file):
    if file[-5:] == 'trees' and (N in file) and (r in file):
        return True
    return False


def get_t(files):
    # Get selection coefficients
    t = [None]*len(treeSeqFiles)

    for i in range(0, len(treeSeqFiles)):
        t[i] = treeSeqFiles[i].split('_')[3].split('=')[1]
    return (t)


if __name__ == '__main__':

    # Define current directory
    current_dir = sys.argv[1]
    os.chdir(current_dir)

    # Define N and r to filter files
    N = 'N=' + str(sys.argv[3]) + '_'
    r = 'r=' + str(sys.argv[4]) + '_'

    # Get files and sort
    entries = os.listdir(current_dir)
    
    treeSeqFiles = sorted(list(filter(filter_files_r, entries)))

    # Calculate pi_s and pi_t for TreeSeq files and store
    pi = calc_pi(treeSeqFiles)

    # Get selection coefficient and migration rate of file
    t = get_t(treeSeqFiles)

    # Write into csv w/ name
    with open(str(sys.argv[5]), 'w') as f:
        writer = csv.writer(f)

        for i in range(0, len(s)):
            writer.writerow([t[i], pi[i]])

