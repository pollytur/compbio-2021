import Bio.SeqIO
import numpy as np
from collections import defaultdict

#N is the size of the window
#N must be a divisor of the length of the chunk(150bp)
N = 100

filename_initial = "./../files_for_project/salmonella-enterica.reads.fna"
filename_variant = "./../files_for_project/salmonella-enterica-variant.reads.fna"

def occurence_counter(file_path, N, percentage = 1.0):
    """
    returns a dictionnary where the keys are the N-sized sequences of nucleotides
    and the values are the number of occurences of each sequence in the batch of samples
    """
    dic = defaultdict(int)
    #we parse the .fna file
    records = list(Bio.SeqIO.parse(file_path,"fasta"))
    for record in records[:int(len(records)*percentage)]:
        #for each sample we extract the nucleotide sequence in the SEQ variable
        SEQ = str(record.seq)
        #we go through the sample with our N-sized window, the stride for this parsing is 1
        for i in range(len(SEQ) - N):
            seq = SEQ[i:i + N]
            #if the key was not present in the dictionnary we add it and we initialize its value to 1
            #if the key was present we just add 1 to its value
            dic[seq] += 1
    return dic

#we write the result of the previpous function in a file because we want to separate the extraction from
#the processing as the extraction is pretty long (arround 10sec for our example)
def dictionary_writer(dic, filename):
    """
    writes the dictionnary obtained with the previous function in a file
    named filename. The write is made in a sorted order sarting with the less occurente sequences
    """
    with open(filename, "w") as f:
        for w in reversed(sorted(dic, key=dic.get, reverse=True)):
            f.write(w + " " + str(dic[w]) + "\n")

file_path_1 = "outputs/dic_initial.txt"
file_path_2 = "outputs/dic_variant.txt"

if __name__ == "__main__":
    print("Calculating occurences for the 'original' sequence...")
    D_1 = occurence_counter(filename_initial, N, 0.25)
    print("Original sequence processed.")
    # print("Calculating occurences for the variant sequence...")
    # D_2 = occurence_counter(filename_variant, N, 0.25)
    # print("Variant sequence processed.")
    
    print("Writing dictionaries to files...")
    dictionary_writer(D_1, file_path_1)
    # dictionary_writer(D_2, file_path_2)
    print("Dictionaries written to files. Pre-processing is complete.")




        
        
