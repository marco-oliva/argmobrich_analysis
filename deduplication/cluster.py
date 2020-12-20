from Bio import SeqIO

from sklearn.cluster import KMeans

import numpy as np

import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    fasta_file = sys.argv[1]

    #Iterate through every read. Accumulate number of reads while recording read length
    read_lengths = ([], [])
    outlier_read_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        read_len = len(record.seq)
        read_lengths[0].append(read_len)
        read_lengths[1].append(record)


    X = np.array(read_lengths[0])

    num_clusters = 200
    kmeans = KMeans(n_clusters=num_clusters).fit(X.reshape((X.shape[0],1)))

    for i in range(0, num_clusters):
        sequences = []
        for j in range(0, kmeans.labels_.shape[0]):
            if kmeans.labels_[j] == i:
                sequences.append(read_lengths[1][j])

        output_filename = os.path.splitext(os.path.basename(fasta_file))[0] + '_' + str(i) + ".fasta"
        with open(output_filename, 'w') as output_fasta:
            SeqIO.write(sequences, output_fasta, "fasta")
