# from occurence import file_extractor, dictionary_writer
from occurence import *
from processing import *
from collections import defaultdict
import numpy as np
import math as m
import itertools as it

file_path_1 = "outputs/dic_initial_filtered.txt"
file_path_2 = "outputs/dic_variant_filtered.txt"



def comparison_dictionaries(dic_initial, dic_variant):
    diff_dict = defaultdict(lambda: 0)
    for key in dic_variant.keys():
        if key not in dic_initial.keys():
            diff_dict[key] = dic_variant[key]
    return diff_dict


def Levenshtein(sequence1, sequence2):
    # Initializing distance matrix
    D = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    for i in range(len(sequence1) + 1):
        D[i][0] = i
    for i in range(len(sequence2) + 1):
        D[0][i] = i

    # Calculating distances between the prefixes
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            if sequence1[i-1] == sequence2[j-1]:
                D[i][j] = D[i-1][j-1]
            else:
                D[i][j] = min(D[i-1][j], D[i][j-1], D[i-1][j-1]) + 1
    return D[len(sequence1)][len(sequence2)]


def naive_clustering(unique_sequences):
    assert len(unique_sequences) != 0
    unique_sequences_copy = unique_sequences
    # Vector of clusters
    Clusters = []
    Empty = False
    while not Empty:
        # Make a new cluster
        cluster = [unique_sequences_copy[0]]
        unique_sequences_copy.remove(unique_sequences_copy[0])
        Cluster_done = False
        # Fill the cluster with elements which link together(have a Levenshtein distance of 2)
        while not Cluster_done:
            j = 0
            Found = False
            while j < len(unique_sequences_copy):
                indices = []
                # Find indices of simillar sequences
                for seq in cluster:
                    if Levenshtein(seq, unique_sequences_copy[j]) == 2.:
                        indices.append(j)
                        Found = True
                # Add sequences to cluster
                for k in indices:
                    cluster.append(unique_sequences_copy[k])
                    unique_sequences_copy.remove(unique_sequences_copy[k])  
                j += 1
            # Check if cluster is done
            if not Found:
                Cluster_done = True
        Clusters.append(cluster)
        # Checking if list is empty
        if not unique_sequences:
            Empty = True
    return Clusters


def order_cluster(cluster):
    # Find some starting sequence
    close_seqs = np.zeros(len(cluster))
    for i in range(len(cluster)):
        seq = cluster[i]
        for seq_compare in set(cluster) - set([seq]):
            if Levenshtein(seq, seq_compare) == 2.:
                close_seqs[i] += 1
    idcs = list(np.where(close_seqs == 1))
    start_seq = cluster[int(idcs[0][0])]

    # Order sequences
    ordered_cluster = [start_seq]
    for i in range(len(cluster)):
        for seq in set(cluster) - set(ordered_cluster):
            if Levenshtein(seq, ordered_cluster[i]) == 2.:
                ordered_cluster.append(seq)
    return ordered_cluster[::-1]


def reconstruct(cluster):
    # Assumed to be already ordered
    reconstructed_sequence = cluster[0]
    for i in range(1,len(cluster)):
        reconstructed_sequence += cluster[i][-1]
    return reconstructed_sequence


def create_alphabet(word_length):
    letters = []
    for i in range(word_length):
        letters.append('A')
        letters.append('T')
        letters.append('G')
        letters.append('C')
    dic = list(it.combinations(letters, word_length))
    return dic


def fast_hamming(cluster, initial_dic, N):
    seq = cluster[int(m.floor(len(cluster)/2))]
    word_length = len(cluster) - N + 1
    letters_alphabet = create_alphabet(word_length)
    index = int(m.floor(len(seq)/2))
    print("Initial sequence: ", seq)
    for word in letters_alphabet:
        seq = list(seq)
        for i in range(word_length):
            seq[index - int(m.floor(word_length/2)) + i] = word[i]
        seq = ''.join(seq)
        if seq in initial_dic:
            print("We got it: ", seq)
            return seq
    print('No mutation found for that cluster')
    return seq


def gene_substitution(sequence, mutated_sequence, original_sequence):
    """
    replaces the original_sequence chunk of sequence by mutated_sequence
    """
    return sequence.replace(mutated_sequence, original_sequence)


def sequence_comparison(seq_1, seq_2):
    """
    compares 2 iterable and stores their differences in a list of tuples
    """
    assert len(seq_1) == len(seq_2)
    res = []
    for k, s_1, s_2 in enumerate(zip(seq_1, seq_2)):
        if s_1 != s_2:
            res.append((s_1, s_2, k))
    return res