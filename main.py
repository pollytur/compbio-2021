#!/usr/bin/python3

from occurence import *
from processing import *
from comparison import *
import matplotlib.pyplot as plt
import math as m
import sys


# filepath for the fasta files
filename_initial = sys.argv[1]
filename_variant = sys.argv[2]

# filename_initial = "./../files_for_project/salmonella-enterica.reads.fna"
# filename_variant = "./../files_for_project/salmonella-enterica-variant.reads.fna"

# N is the window size and T is the threshold for the error filtering
N = 17
T = 10

# sampling rate for testing purposes

rate = 1

if __name__ == "__main__":

    ##############################################################################
    # Processing FASTA files and filtering out the errors

    # storing the k-mers built from the reads of the non mutant salmonella in a dictionnary
    print("Calculating occurences for the 'original' sequence...")
    initial_dic = occurence_counter(filename_initial, N, rate)
    print("Original sequence processed.")
    print("Filtering out errors for 'original' sequence dictionary..")
    initial_dic = filter_out(initial_dic, threshold = T)
    print("Filtering done for original.")


    # storing the k-mers built from the reads of the mutant salmonella in a dictionnary
    print("Calculating occurences for the variant sequence...")
    variant_dic = occurence_counter(filename_variant, N, rate)
    print("Variant sequence processed.")
    print("Filtering out errors for variant sequence dictionary..")
    variant_dic = filter_out(variant_dic, threshold = T)
    print("Filtering done for variant.")


    ##############################################################################
    # Finding sequences in the variant which are candidates for having mutations

    print("Extracting all sequences in the variant which cannot be found in the original sequence and writing them in a file...")
    unique_sequences_dict = comparison_dictionaries(initial_dic, variant_dic)
    print("There are {:d} {:d}-mers that are in the variant's reads which are not in the initial one.".format(len(unique_sequences_dict), N))


    ##############################################################################
    # Clustering candidate sequences 

    print("Clustering sequences...")
    Clusters = naive_clustering(list(unique_sequences_dict.keys()))
    
    # Removing small clusters
    temp = [cluster for cluster in Clusters if len(cluster) >= N]
    Clusters = temp

    # Sorting
    Clusters_sorted = [order_cluster(cluster) for cluster in Clusters]
    Clusters = Clusters_sorted

    print("\n\n")
    print("#####################################################")
    print("\n\n")


    # Reconstructing sequences
    reconstructed_sequences = []
    mutated_sequences = []
    for cluster in Clusters:
        reconstructed_sequence = reconstruct(cluster)
        mutated_sequence = fast_hamming(cluster,initial_dic, N)
        if mutated_sequence == cluster[int(m.floor(len(cluster)/2))]:
            continue
        reconstructed_sequences.append(reconstructed_sequence)
        mutated_sequences.append(mutated_sequence)

        # we retrieve the sequence which contains all the k-mers containing the mutated nucleotides
        #this original sequence is of the order of 2*k
        original_sequence = gene_substitution(reconstructed_sequence, cluster[int(m.floor(len(cluster)/2))], mutated_sequence)

        # the cluster may be in the wrong order if this happen we flip flip it and we redo the process
        if original_sequence == reconstructed_sequence:
            reconstructed_sequences.remove(reconstructed_sequence)
            cluster = cluster[::-1]
            reconstructed_sequence = reconstruct(cluster)
            reconstructed_sequences.append(reconstructed_sequence)
            original_sequence = gene_substitution(reconstructed_sequence, cluster[int(m.floor(len(cluster)/2))], mutated_sequence)
            
    
        # initial_read = find_sequence_in_reads(filename_initial, original_sequence)
        # variant_read = find_sequence_in_reads(filename_variant, reconstructed_sequence)
        print("The sequence:                     ", mutated_sequence)
        print("has been mutated to the sequence: ", cluster[int(m.floor(len(cluster)/2))])
        print("The mutations which have taken place are: ")
        print("\n\n")
        snps = sequence_comparison(mutated_sequence, cluster[int(m.floor(len(cluster)/2))])
        for snp in snps:
            print(snp[0], " to ", snp[1], " at position ", snp[2])
        print("\n\n")

   
