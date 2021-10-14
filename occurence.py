import Bio.SeqIO

#N is the size of the window
#N must divide the length of the samples
N = 10

filename = "salmonella-enterica.reads.fna"


def occurence_counter(file_path, N):
    """
    returns a dictionnary where the keys are the N-sized sequences of nucleotides
    and the values are the number of occurences of each sequence in the batch of samples
    """
    dic = {}
    #we parse the .fna file
    for i,record in enumerate(Bio.SeqIO.parse(file_path,"fasta")):
        #for each sample we extract the nucleotide sequence in the SEQ variable
        SEQ = str(record.seq)
        #we go through the sample with our N-sized window, the stride for this parsing is 1
        counter = 0
        while(counter + N <= len(SEQ)):
            seq = SEQ[counter:counter + N]
            #if the key was not present in the dictionnary we add it and we initialize its value to 1
            #if the key was present we just add 1 to its value
            if seq in dic.keys():
                dic[seq] += 1
            else:
                dic[seq] = 1
            counter += 1
    return dic

file_path = "outputs\dic.txt"

#we write the result of the previpous function in a file because we want to separate the extraction from
#the processing as the extraction is pretty long (arround 10sec for our example)
def dictionnary_writer(dic, filename):
    """
    writes the dictionnary obtained with the previous function in a file
    named filename. The write is made in a sorted order sarting with the less occurente sequences
    """
    with open(filename, "w") as f:
        for w in reversed(sorted(dic, key=dic.get, reverse=True)):
            f.write(w + " " + str(dic[w]) + "\n")

if __name__ == "__main__":
    D = occurence_counter(filename, N)
    dictionnary_writer(D, file_path)




        
        
