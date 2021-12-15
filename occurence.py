import Bio.SeqIO
import numpy as np
from collections import defaultdict


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






        
        
