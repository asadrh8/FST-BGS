# This script takes in a list of treeSeq files (input), and
# outputs an excel file with all their pi_s, pi_t, and FST values

# sys arg 1 = directory – containing treeSeq Files
# sys arg 2 = mu – neutral
# sys arg 3 = r – between deleterious and neutral site
# sys arg 4 = output_csv file

# Import modules
import msprime
import pyslim
import tskit
import numpy as np
import os
import sys
import subprocess
import csv

# Define function to calculate pi_s, pi_t, and FST for TreeSeq of given configuration
# (10^4 neutral sites flanked by 350 selected sites on each side, 10 demes)


def calc_pi_FST(files):
    # files: TreeSeq files

    # List to store pi_s, pi_t, and FST values from each TreeSeq file
    pi_s = [None] * len(files)
    pi_t = [None] * len(files)
    FST = [None] * len(files)

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

        # Get pi_s and pi_t using calc_pi_s_pi_t
        pi_s[i], pi_t[i] = calc_pi_s_pi_t(current_file)

        #### Apply MAF Filter of MAF < 0.05 to calculate FST ####
        filtered_tree = MAF_filter(current_file)

        # Get pi_s and pi_t with MAF < 0.05 filter to calculate FST
        pi_s_MAF, pi_t_MAF = calc_pi_s_pi_t(filtered_tree)

        # Calculate the average pairwise diversity between subpopulations
        num_subpops = 10
        pi_b_rep = (num_subpops / (num_subpops - 1)) * \
            pi_t_MAF - (1 / (num_subpops - 1)) * pi_s_MAF

        # Calculate FST
        FST[i] = 1 - pi_s_MAF / pi_b_rep

        # Return iteration completed
        print(i)

    return ([pi_s, pi_t, FST])


def calc_pi_s_pi_t(file):
    # Input: TreeSeq File with neutral mutations added through msprime

    # Calculate average pairwise diversity within subpopulations
    num_subpops = 10
    pi_subpop = [None] * num_subpops

    # Get pi for each subpopulation
    for subpop in range(1, num_subpops + 1):

        pi_subpop[subpop - 1] = file.diversity(sample_sets=file.samples(population=subpop),
                                               windows=[0., 350., 10351., 10700.])[1]

    # Average subpopulation pairwise diversity
    pi_s = sum(pi_subpop)/num_subpops

    # Calculate total average pairwise diversity
    pi_t = file.diversity(
        windows=[0., 350., 10351., 10700.])[1]

    return ([pi_s, pi_t])


def MAF_filter(current_file):
    # Input TreeSeq File with neutral mutations added through msprime

    #### Apply MAF Filter of MAF < 0.05 to calculate FST ####

    # Get genotype matrix
    genotype_matrix = current_file.genotype_matrix()

    # Num sites
    num_sites = len(genotype_matrix[:, 0])

    # Get list of sites with MAF < 0.05 to remove
    sites_to_rm = []

    # Total number of individuals
    N_total = 10**4

    for s in range(0, num_sites):

        # Get frequency of most common allele = highest frequency allele / num total individuals
        major_allele_freq = max(np.bincount(
            genotype_matrix[s, :])) / N_total

        if (1 - major_allele_freq) < 0.05:
            sites_to_rm.append(s)

    # Filter test tree and calculate FST for it
    filtered_tree = current_file.delete_sites(sites_to_rm)

    return (filtered_tree)

# Filter entries so it only contains TreeSeq files of particular N, r


def filter_files(file):
    if file[-5:] == 'trees' and (r in file):
        return True
    return False

# Get selection coefficient of each TreeSeq


def get_t_and_m(files):
    # Get selection coefficients
    t = [None]*len(files)
    m = [None]*len(files)

    for i in range(0, len(files)):
        t[i] = float(files[i].split('_')[3].split('=')[1])
        m[i] = float(files[i].split('_')[4].split('=')[1])

    return ([s, m])


if __name__ == '__main__':

    # Define current directory
    current_dir = sys.argv[1]
    os.chdir(current_dir)

    # Define r to filter files
    r = 'r=' + str(sys.argv[3]) + '_'

    # Get files and sort
    entries = os.listdir(current_dir)
    treeSeqFiles = sorted(list(filter(filter_files, entries)))

    # Calculate pi_s, pi_t, and FST for TreeSeq files and store
    pi_s, pi_t, FST = calc_pi_FST(treeSeqFiles)

    # Get selection coefficient and migration rate of file
    t, m = get_t_and_m(treeSeqFiles)

    # Write into csv w/ name
    with open(str(sys.argv[4]), 'w') as f:
        writer = csv.writer(f)

        for i in range(0, len(t)):
            writer.writerow([t[i], m[i], pi_s[i], pi_t[i], FST[i]])
