import numpy as np
import os.path, gzip, errno
from Bio import SeqIO
from sklearn.cluster import KMeans
import argparse


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python 2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-r', help='Reads File (FASTA/FASTQ)', type=str, dest='reads_file', required=True)
    parser.add_argument('-o', help='Output directory', type=str, dest='out_dir', required=True)
    parser.add_argument('-n', help='Number of clusters to generate', type=int, dest='num_clusters', default=200)
    args = parser.parse_args()

    reads_file = args.reads_file
    output_dir = args.out_dir
    num_clusters = args.num_clusters

    fastq_extensions = ['fastq', 'fq', 'FASTQ', 'FQ', 'fq.gz', 'fastq.gz', 'FQ.gz', 'FASTQ.gz']
    fasta_extensions = ['fasta', 'fa', 'FASTA', 'FA', 'fa.gz', 'fasta.gz', 'FA.gz', 'FASTA.gz']

    file_type = ''
    for p_ext in fastq_extensions:
        if (reads_file.endswith(p_ext)):
            file_type = 'fastq'
    for p_ext in fasta_extensions:
        if (reads_file.endswith(p_ext)):
            file_type = 'fasta'

    if file_type == '':
        print('File has to be either fastq or fasta')
        exit()

    # Iterate through every read. Accumulate number of reads while recording read length
    read_lengths = ([], [])

    if (reads_file.endswith('.gz')):
        print('Opening gzipped file')
        file_handler = gzip.open(reads_file, 'rt')
    else:
        print('Opening uncompressed file')
        file_handler = open(reads_file, 'rt')

    for record in SeqIO.parse(file_handler, file_type):
        read_len = len(record.seq)
        read_lengths[0].append(read_len)
        read_lengths[1].append(record)

    X = np.array(read_lengths[0])

    if (len(read_lengths[0]) < num_clusters):
        num_clusters = int(len(read_lengths[0]) / 10)
    print("Num of cluster used: {}".format(num_clusters))

    kmeans = KMeans(n_clusters=num_clusters).fit(X.reshape((X.shape[0], 1)))

    fasta_files = list()
    for i in range(0, num_clusters):
        sequences = []
        for j in range(0, kmeans.labels_.shape[0]):
            if kmeans.labels_[j] == i:
                sequences.append(read_lengths[1][j])

        output_filename = output_dir + '/' + os.path.splitext(os.path.basename(reads_file))[0] + '_' + str(i) + '.fasta'
        fasta_files.append(output_filename)
        with open(output_filename, 'w') as output_fasta:
            SeqIO.write(sequences, output_fasta, 'fasta')

    for file_name in fasta_files:
        print(file_name)

if __name__ == "__main__":
    main()


