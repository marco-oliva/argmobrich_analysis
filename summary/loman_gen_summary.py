from Bio import SeqIO

import matplotlib.pyplot as plt
import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 3):
        sys.exit(1)

    if sys.argv[2] == "1":
        fastq_file = sys.argv[1] + ".fastq"

        #Iterate through every read. Accumulate number of reads while recording read length
        num_reads = 0
        read_lengths = []
        outlier_read_lengths = []
        for record in SeqIO.parse(fastq_file, "fastq"):
            num_reads += 1

            read_len = len(record.seq)
            if read_len < 10000:
                read_lengths.append(read_len)
            else:
                outlier_read_lengths.append(read_len)


        print(os.path.basename(fastq_file) + ": " + str(num_reads))

    else:
        num_reads = 0
        read_lengths = []
        outlier_read_lengths = []
        for i in range(1, int(sys.argv[2])+1):
            if i < 10:
                fastq_file = sys.argv[1] + ".part-0" + str(i) + ".fastq"
            else:
                fastq_file = sys.argv[1] + ".part-" + str(i) + ".fastq"

            #Iterate through every read. Accumulate number of reads while recording read length
            for record in SeqIO.parse(fastq_file, "fastq"):
                num_reads += 1

                read_len = len(record.seq)
                if read_len < 10000:
                    read_lengths.append(read_len)
                else:
                    outlier_read_lengths.append(read_len)

            print(os.path.basename(fastq_file) + ": " + str(num_reads))

    #Plot histogram
    fig, ax = plt.subplots(tight_layout=True)
    ax.set_title(os.path.basename(fastq_file))
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.hist(read_lengths, bins="fd", range=(min(read_lengths), max(read_lengths)))

    plt.savefig(os.path.basename(fastq_file) + ".png")

    #Plot upper histogram
    fig, ax = plt.subplots(tight_layout=True)
    ax.set_title(os.path.basename(fastq_file))
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.hist(outlier_read_lengths, bins="fd", range=(min(outlier_read_lengths), max(outlier_read_lengths)))

    plt.savefig(os.path.basename(fastq_file) + "_upper.png")

    sys.exit(0)
