import Bio.SeqIO

#N must divide the length of the samples
N = 10
dic = {}

for i,record in enumerate(Bio.SeqIO.parse("salmonella-enterica.reads.fna","fasta")):
    ID = str(record.id)
    SEQ = str(record.seq)
    counter = 0
    while(counter + N <= len(SEQ)):
        seq = SEQ[counter:counter + N]
        if seq in dic.keys():
            dic[seq] += 1
        else:
            dic[seq] = 1
        counter += 1

with open("dic.txt", "w") as f:
    for w in reversed(sorted(dic, key=dic.get, reverse=True)):
        f.write(w + " " + str(dic[w]) + "\n")






        
        
