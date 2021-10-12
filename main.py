from Bio import SeqIO
from collections import Counter
# todo add click here

if __name__ == '__main__':
    path = "salmonella-enterica.reads.fna"
    fasta_sequences = SeqIO.parse(path, "fasta")
    window = 5
    small_seqs = []
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
        for i in range(window, len(sequence)):
            small_seqs.append(sequence[i-window: window])

    fasta_counter = Counter(small_seqs)
    n = 50
    m = 500
    n_least_common = fasta_counter.most_common()[:-n-1:-1]
    m_most_common = fasta_counter.most_common(m)

    for rare in n_least_common:
        for popular in m_most_common:
            similarity_counter = 0
            if popular == "":
                continue
            print("Rare = ", rare)
            print("Popular = ", popular)
            for i in range(window):
                if(rare[i] == popular[i]):
                    similarity_counter += 1

            if similarity_counter == window - 1:
                print("Rare = ", rare)
                print("Popular = ", popular)


