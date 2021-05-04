from Bio import SeqIO

import matplotlib.pyplot as plt
from scipy.stats import kurtosis
from scipy.stats import skew

import os.path
import statistics
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    fastq_file = sys.argv[1]

    #Iterate through every read. Accumulate number of reads while recording read length
    num_reads = 0
    read_lengths = []
    #outlier_read_lengths = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        num_reads += 1

        read_len = len(record.seq)
        #if read_len < 10000:
        read_lengths.append(read_len)
        #else:
            #outlier_read_lengths.append(read_len)
        if read_len == 46438:
            print(record.id)
            sys.exit(0)


    print(os.path.basename(fastq_file) + " number of reads: " + str(num_reads))
    print(os.path.basename(fastq_file) + " mean: " + str(statistics.mean(read_lengths)))
    print(os.path.basename(fastq_file) + " range: " + str((min(read_lengths), max(read_lengths))))
    print(os.path.basename(fastq_file) + " max elem: " + str(read_lengths.index(max(read_lengths))))
    print(os.path.basename(fastq_file) + " std dev: " + str(statistics.stdev(read_lengths)))
    print(os.path.basename(fastq_file) + " variance: " + str(statistics.variance(read_lengths)))
    print(os.path.basename(fastq_file) + " skew: " + str(skew(read_lengths)))
    print(os.path.basename(fastq_file) + " kurtosis: " + str(kurtosis(read_lengths)))
    print("=================================")

    #Plot histogram
    fig, ax = plt.subplots(tight_layout=True)
    ax.set_title(os.path.basename(fastq_file))
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Frequency")
    bin_tuple = ax.hist(read_lengths, bins=500)

    plt.savefig(os.path.basename(fastq_file) + ".png")

    #Plot upper histogram
    #fig, ax = plt.subplots(tight_layout=True)
    #ax.set_title(os.path.basename(fastq_file))
    #ax.set_xlabel("Read length (bp)")
    #ax.set_ylabel("Number of reads")
    #ax.hist(outlier_read_lengths, bins="fd", range=(min(outlier_read_lengths), max(outlier_read_lengths)))

    #plt.savefig(os.path.basename(fastq_file) + "_upper.png")

    sys.exit(0)
