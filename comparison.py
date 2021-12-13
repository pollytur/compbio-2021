# from occurence import file_extractor, dictionary_writer
from occurence import *
from processing import *
from collections import defaultdict
import numpy as np
import math

file_path_1 = "outputs/dic_initial_filtered.txt"
file_path_2 = "outputs/dic_variant_filtered.txt"

def create_alphabet():
    dic = defaultdict(int)
    seq = {'A', 'B', 'C'}
    letters = ['A', 'T', 'G', 'C']
    for letter_zero in letters:
        for letter_one in letters:
            for letter_two in letters:
                print(seq)
                seq = list(seq)
                seq[0] = letter_zero
                seq[1] = letter_one
                seq[2] = letter_two
                seq = ''.join(seq)
                print(seq)
                dic[seq] = 1
    print(dic)
    return dic

def fast_hamming(sequence, initial_dic):
    seq = sequence
    three_letters_alphabet = create_alphabet()
    index = int(math.floor(len(sequence)/2))
    print("Initial sequence: ", sequence)
    for word in three_letters_alphabet:
        seq = list(sequence)
        seq[index - 1] = word[0]
        seq[index] = word[1]
        seq[index + 1] = word[2]
        seq = ''.join(seq)
        print("")
        print("Word: ", word)
        print("Modified sequence: ", sequence)
        if seq in initial_dic:
            print("We got it: ", seq)
    return seq

def comparison(file_path_1, file_path_2):
    initial_dic = file_extractor(file_path_1)
    variant_dic = file_extractor(file_path_2)
    return set(initial_dic) - set(variant_dic)

def dictionary_writer(dic, filename):
    """
    writes the dictionnary obtained with the previous function in a file
    named filename. The write is made in a sorted order sarting with the less occurente sequences
    """
    with open(filename, "w") as f:
        for w in reversed(sorted(dic, key=dic.get, reverse=True)):
            f.write(w + " " + str(dic[w]) + "\n")

def comparison_dictionaries(dic_initial, dic_variant):
    # res = []
    diff_dict = defaultdict(lambda: 0)
    for key in dic_variant.keys():
        if key not in dic_initial.keys():
            diff_dict[key] = dic_variant[key]
            # res.append(key)
    return diff_dict

def write_unique_sequences(unique_sequences):
    with open("outputs/unique_sequences.txt","w") as f:
        for seq in unique_sequences:
            f.write(seq + "\n")

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

if __name__ == "__main__":
    fast_hamming('AACCCGTCCCCTCGCCC', filter_out(file_extractor("./outputs/initial_dic_17.txt"), 17))


#################################NOT USED###############################################3

# BUILD part of PAM algorithm which returns a list of medoids around which the data will be clustered
# def BUILD(unique_sequences_dict_enum,D,clustersN,N):
#     # Obtaining first medoid
#     D_sums = np.zeros(N)
#     for i in range(N):
#         D_sums[i] = np.sum(D[i])
#     m1_index = np.argmin(D_sums)

#     # List of medoids
#     medoids = [(m1_index,unique_sequences_dict_enum[m1_index][1])]

#     # Finding other medoids
#     for i in range(clustersN - 1):
#         # Contrivution matrix
#         C = np.zeros((N,N))
#         C_sums = np.zeros(N)
#         for i, key1 in set(unique_sequences_dict_enum) - set(medoids):
#             for j, key2 in set(unique_sequences_dict_enum) - set(medoids):
#                 # Finding distance to closest medoid
#                 med_dist = []
#                 for k in range(len(medoids)):
#                     med_dist.append(D[j][medoids[k][0]])
#                 Dj = min(med_dist)
    
#                 # Calculating contribution
#                 C[j][i] = max(Dj - Levenshtein(key1,key2),0)
#             C_sums[i] = np.sum(C[:][i])
#         med_idx = np.argmax(C_sums)
#         medoids.append((med_idx,unique_sequences_dict_enum[med_idx][1]))
#     return medoids


# def SWAP(unique_sequences_dict_enum,D,medoids,N):
#     not_medoids = list(set(unique_sequences_dict_enum) - set(medoids))
#     while True:
#         T = np.full((N,N),np.inf)
#         for not_medoid in not_medoids: #h
#             for medoid in medoids: #i
#                 # Indices
#                 i = medoid[0]
#                 h = not_medoid[0]
#                 Tih = 0
#                 for j, elem_j in set(unique_sequences_dict_enum) - set(not_medoid) - set(medoid):
#                     # Distances
#                     Dh = D[h][j]
#                     Di = D[i][j]
#                     D_other_medoids = [D[elem[0]][j] for elem in set(medoids)-set([medoid])]
#                     Dmj = min(D_other_medoids)
#                     # Different cases which will determine contribution of voter
#                     if Dh > Dmj and Di > Dmj:
#                         Cjih = 0
#                     elif Di <= Dmj:
#                         if Dh < Dmj:
#                             Cjih = Dh - Di
#                         else:
#                             Cjih = Dmj - Di
#                     elif Dh <= Dmj and Dmj < Di:
#                         Cjih = Dh - Dmj 
#                     # Updating sum
#                     Tih += Cjih
#                 T[i][h] = Tih
#         i,h = np.unravel_index(T.argmin(),T.shape)
#         if T[i][h] >= 0:
#             break
#         else:
#             elem_h = unique_sequences_dict_enum[h] #not_medoid
#             elem_i = unique_sequences_dict_enum[i] #medoid
#             # Move old not_medoid to medoids
#             not_medoids.remove(elem_h)
#             medoids.append(elem_h)
#             # Move old medoid to not_meoids
#             medoids.remove(elem_i)
#             not_medoids.append(elem_i)
#     return medoids    


# def PAM(unique_sequences_dict, clustersN):
#     N = len(unique_sequences_dict)
#     # Initializing dissimilarity matrix
#     D = np.zeros((N,N))
#     unique_sequences_dict_enum = list(enumerate(unique_sequences_dict))

#     for i, key1 in unique_sequences_dict_enum:
#         for j, key2 in unique_sequences_dict_enum:
#             if i == j:
#                 D[i][j] = 0
#             else:
#                 D[i][j] = Levenshtein(key1,key2)
    
#     # BUILD
#     medoids = BUILD(unique_sequences_dict_enum, D, clustersN, N)

#     # SWAP
#     medoids = SWAP(unique_sequences_dict_enum, D, medoids, N)

#     return medoids

# def cluster_sequence(medoids, non_medoid_sequences):
#     clusterN = len(medoids)
#     non_medoid_sequences_copy = non_medoid_sequences
#     Clusters = []
#     for i in range(clusterN):
#         cluster = [medoids[i]]
#         Done = False
#         k = 0
#         while not Done:
#             j = 0
#             Found = False
#             while j < len(non_medoid_sequences_copy):
#                 if Levenshtein(cluster[k], non_medoid_sequences_copy[j]) == 2.:
#                     cluster.append(non_medoid_sequences_copy[j])
#                     non_medoid_sequences_copy.remove(non_medoid_sequences_copy[j])
#                     Found = True
#                 else:
#                     j += 1
#             if Found:
#                 k += 1
#             else:
#                 Done = True
#         Clusters.append(cluster)
#     return Clusters
###########################################################################################


if __name__ == "__maian__":
    unique_sequences = file_extractor("outputs/orig_main_unique_sequences_dict10.txt")
    unique_sequences_list = list(unique_sequences.keys())

    # print("Calculating medoids for clustering...")
    # medoid_sequences = PAM(unique_sequences,4)
    # medoid_sequences = [elem[1] for elem in medoid_sequences]
    # non_medoid_sequences = list(set(unique_sequences) - set(medoid_sequences))
    # Clusters = cluster_sequence(medoid_sequences, non_medoid_sequences)

    Clusters = naive_clustering(unique_sequences_list)
    for cluster in Clusters:
        cluster = sorted(cluster, reverse=True)

    print("Clustering done. Writing results to the file 'outputs/clustered_unique_sequences.txt'")
    # Writing clustered sequences to file
    with open("outputs/clustered_unique_sequences.txt","w") as f:
        for i in range(len(Clusters)):
            for sequence in Clusters[i]:
                f.write(sequence + " ===> " + str(i+1) + "\n")
            f.write("--------------------------------------------------\n")



    # # Attaching an index for each cluster(medoid)
    # medoid_sequences_enumed = list(enumerate(medoid_sequences))

    # print("Assigning each sequence to the cluster with the lowest Levenshtein distance...")
    # # Assigning each sequence to the cluster of the medoid which is 'closest' 
    # non_medoid_sequences_enumed = []
    # for sequence in non_medoid_sequences:
    #     distances = [Levenshtein(sequence,medoid) for medoid in medoid_sequences]
    #     non_medoid_sequences_enumed.append((np.argmin(distances),sequence))

    # # Making final list of clustering and sorting it
    # clustered_sequences = medoid_sequences_enumed + non_medoid_sequences_enumed
    # clustered_sequences.sort(key=lambda tup: tup[0])

    # print("Clustering done. Writing results to the file 'outputs/clustered_unique_sequences.txt'")
    # # Writing clustered sequences to file
    # with open("outputs/clustered_unique_sequences.txt","w") as f:
    #     for cluster,sequence in clustered_sequences:
    #         f.write(sequence + " ===> " + str(cluster) + "\n")


    # print("Extracting all sequences in the variant which cannot be found in the original and writing them in a file...")
    # unique_sequences,unique_sequences_dict = comparison_dictionaries(initial_dic_filtered, variant_dic_filtered)
    # write_unique_sequences(unique_sequences)
    # dictionary_writer(unique_sequences_dict,"outputs/unique_sequences_dict.txt")
    # print("Done. Written in file ./outputs/unique_sequences.txt")
    # # print(len(unique_sequences))

    # smaller_kmers = align_kmers(unique_sequences)
    # print(smaller_kmers)