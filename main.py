from Bio import SeqIO
from collections import Counter
# todo add click here hello
a = 9090
if __name__ == '__main__':
    path = "salmonella-enterica.reads.fna"
    fasta_sequences = SeqIO.parse(path, "fasta")
    window = 10
    small_seqs = []
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
        for i in range(window, len(sequence)):
            small_seqs.append(sequence[i-window: window])

    fasta_counter = Counter(small_seqs)
    n = 50
    m = 50000
    sim_counter = 0
    n_least_common = fasta_counter.most_common()[:-n-1:-1]
    m_most_common = fasta_counter.most_common(m)[1:m]
    print(m_most_common)
    for rare in n_least_common:
        for popular in m_most_common:
            similarity_counter = 0

            if len(popular[0]) < window or len(rare[0]) < window:
                continue

            for i in range(window):
                if rare[0][i] == popular[0][i]:
                    similarity_counter += 1
                #print("similarity_counter = ", similarity_counter)

            if similarity_counter == window - 1:
                print("Rare = ", rare)
                print("Popular = ", popular)
                sim_counter += 1

    print("sim_counter = ", sim_counter)


