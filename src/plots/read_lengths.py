import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import gzip
import numpy as np


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def reject_outliers(data, m = 2.):
    data = np.array(data)
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]


def main():
    parser = argparse.ArgumentParser(description='Plot read lengths')
    parser.add_argument('-s', help='Standard Deviations', default=2, dest='std_dev', type=int)
    parser.add_argument('-b', help='Number of bins', default=100, dest='bins', type=int)
    parser.add_argument('input', help='Input Fastq file', dest='reads_file', required=True)
    args = parser.parse_args()

    file_path = args.input
    if is_gz_file(file_path):
        handle = gzip.open(file_path)
    else:
        handle = open(file_path)

    read_lengths = dict()
    for record in SeqIO.parse(handle, "fastq"):
        read_lengths[record.name] = len(record.seq)

    read_lengths_list = reject_outliers([float(v) for v in read_lengths.values()], args.std_dev)

    # Plot
    plt.hist(read_lengths_list, bins=args.bins)
    plt.ylabel('Count')
    plt.xlabel('Read Length')


if __name__ == '__main__':
    main()
