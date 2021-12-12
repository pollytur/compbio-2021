from occurence import occurence_counter
from processing import filter_out, file_extractor
from comparison import *
from stats import draw_hist
import matplotlib.pyplot as plt

filename_initial = "./../files_for_project/salmonella-enterica.reads.fna"
filename_variant = "./../files_for_project/salmonella-enterica-variant.reads.fna"
N = 17
T = 10
if __name__ == "__main__":

    ##############################################################################
    # Processing FASTA files and filtering out the errors
    print("Calculating occurences for the 'original' sequence...")
    initial_dic = occurence_counter(filename_initial, N, 1.)
    dictionary_writer(initial_dic, "./outputs/initial_dic_{:d}.txt".format(N))
    # initial_dic = file_extractor("./outputs/initial_dic_{:d}.txt".format(N))
    print("Original sequence processed.")
    print("Filtering out errors for both 'original' sequence dictionary..")
    initial_dic = filter_out(initial_dic, threshold = T)
    print("Filtering done for original.")

    print("Calculating occurences for the variant sequence...")
    variant_dic = occurence_counter(filename_variant, N, 1.)
    dictionary_writer(variant_dic, "./outputs/variant_dic_{:d}.txt".format(N))
    # variant_dic = file_extractor("./outputs/variant_dic_{:d}.txt".format(N))
    print("Variant sequence processed.")
    print("Filtering out errors for both variant sequence dictionary..")
    variant_dic = filter_out(variant_dic, threshold = T)
    print("Filtering done for variant.")

    ##############################################################################
    # Finding sequences in the variant which are candidates for having mutations
    print("Extracting all sequences in the variant which cannot be found in the original sequence and writing them in a file...")
    unique_sequences_dict = comparison_dictionaries(initial_dic, variant_dic)
    dictionary_writer(unique_sequences_dict,"outputs/main_unique_sequences_dict.txt")
    print("Done. Written in file ./outputs/main_unique_sequences.txt")
    print("There are {:d} unique sequences.".format(len(unique_sequences_dict)))


    ##############################################################################
    # Clustering
    print("Clustering sequences...")
    Clusters = naive_clustering(list(unique_sequences_dict.keys()))
    for cluster in Clusters:
        cluster = sorted(cluster)

    print("Clustering done. Writing results to the file 'outputs/clustered_unique_sequences.txt'")
    # Writing clustered sequences to file
    with open("outputs/clustered_unique_sequences_N{:d}_T{:d}.txt".format(N,T),"w") as f:
        for i in range(len(Clusters)):
            for sequence in Clusters[i]:
                f.write(sequence + " ===> " + str(i+1) + "\n")
            f.write("--------------------------------------------------\n")



    # =========================================================================================
    # print("Making histograms...")
    # plt.figure()
    # plt.hist(initial_dic.values(), range(1,max(initial_dic.values())))
    # plt.title("Initial sequence")
    
    # plt.figure()
    # plt.hist(variant_dic.values(), range(1,max(variant_dic.values())))
    # plt.title("Variant sequence")
    # print("Done.")
    # plt.show()
