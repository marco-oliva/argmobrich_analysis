from Bio import SeqIO

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    params = []
    with open(sys.argv[1], 'r') as input_handle:
        for line in input_handle:
            line_params = line.split(',')
            plot_title = line_params[0]
            fastq_file = line_params[1][:-1]
            params.append((plot_title, fastq_file))

    for param_pair in params:
        plot_title = param_pair[0]
        fastq_file = param_pair[1]

        #Iterate through every read. Accumulate number of reads while recording read length
        num_reads = 0
        read_lengths = []
        inset_read_lengths = []
        for record in SeqIO.parse(fastq_file, "fastq"):
            num_reads += 1

            read_len = len(record.seq)
            if read_len < 8000:
                read_lengths.append(read_len)
            else:
                read_lengths.append(read_len)
                inset_read_lengths.append(read_len)


        #Plot histogram
        fig, ax = plt.subplots()
        ax.set_title(plot_title)
        ax.set_xlabel("Read length (nt)")
        ax.set_ylabel("Frequency (Number of Reads)")
        ax.hist(read_lengths, bins=500)


        #Setup inset histogram
        #Make rect
        inset_ax = plt.axes([0,0,1,1])

        #Set pos and relative size
        inset_pos = InsetPosition(ax, [0.3,0.5,0.6,0.4])
        inset_ax.set_axes_locator(inset_pos)
        #Indicate what portion of axes is covered by inset
        mark_inset(ax, inset_ax, loc1=3, loc2=4, fc="none", ec="0.5")
        inset_ax.hist(inset_read_lengths, bins=100)

        plt.savefig(plot_title + ".png")

    sys.exit(0)
