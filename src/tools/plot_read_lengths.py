#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import gzip
import numpy as np
from common import *


def main():
    parser = argparse.ArgumentParser(description='Plot read lengths')
    parser.add_argument('-s', help='Reject outliers. Standard Deviations', default=0, dest='std_dev', type=int)
    parser.add_argument('-b', help='Number of bins', default=100, dest='bins', type=int)
    parser.add_argument('-o', help="Output file", default='', dest='out_file', type=str)
    parser.add_argument('--title', help='Plot title', default='', dest='title', type=str)
    parser.add_argument('--inset', help="Add inset. Left endpoint of the interval", dest='inset', default=0, type=int)
    parser.add_argument('--inset-position', help='Inset position', default='0.3,0.5,0.6,0.4', dest='inset_pos', type=str)
    parser.add_argument('input', help='Input Fastq file')
    args = parser.parse_args()

    file_path = args.input
    if is_gz_file(file_path):
        handle = gzip.open(file_path, 'rt')
    else:
        handle = open(file_path, 'rt')

    read_lengths = list()
    inset_read_lengths = list()
    for record in SeqIO.parse(handle, "fastq"):
        read_length = len(record.seq)
        read_lengths.append(read_length)
        if args.inset != 0 and read_length > args.inset:
            inset_read_lengths.append(read_length)

    if args.std_dev != 0:
        read_lengths = reject_outliers(read_lengths, args.std_dev)
        max_read_length = max(read_lengths)
        inset_read_lengths = [i_read_length for i_read_length in inset_read_lengths if i_read_length <= max_read_length]

    # Plot histogram
    fig, ax = plt.subplots()
    ax.set_title(args.title)
    ax.set_xlabel("Read length (nt)")
    ax.set_ylabel("Frequency (Number of Reads)")
    ax.hist(read_lengths, bins=args.bins)

    if args.inset != 0 and len(inset_read_lengths) > 0:
        # Setup inset histogram
        # Make rect
        inset_ax = plt.axes([0, 0, 1, 1])

        # Set pos and relative size
        inset_pos_list_str = args.inset_pos.split(',')
        inset_pos_list_float = [float(s) for s in inset_pos_list_str]
        inset_pos = InsetPosition(ax, inset_pos_list_float)
        inset_ax.set_axes_locator(inset_pos)

        # Indicate what portion of axes is covered by inset
        mark_inset(ax, inset_ax, loc1=3, loc2=4, fc="none", ec="0.5")
        inset_ax.hist(inset_read_lengths, bins=args.bins)

    # Save
    if args.out_file != '':
        plt.savefig(args.out_file)
    else:
        plt.savefig(file_path + '.hist.pdf')


if __name__ == '__main__':
    main()
