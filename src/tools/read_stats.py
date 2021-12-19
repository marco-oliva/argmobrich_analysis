#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import statistics
import json
from scipy.stats import kurtosis
from scipy.stats import skew
import gzip
from common import *

def main():
    parser = argparse.ArgumentParser(description='Get reads statistics from reads file')
    parser.add_argument('input', help='Input Fastq file')
    args = parser.parse_args()

    TELS_statistics = dict()

    file_path = args.input
    if is_gz_file(file_path):
        handle = gzip.open(file_path, 'rt')
    else:
        handle = open(file_path, 'rt')

    read_lengths = list()
    for record in SeqIO.parse(handle, "fastq"):
        read_length = len(record.seq)
        read_lengths.append(read_length)

    TELS_statistics['NUM_OF_READS'] = str(len(read_lengths))
    TELS_statistics['READ_LENGTH_MEAN'] = str(statistics.mean(read_lengths))
    TELS_statistics['READ_LENGTH_MEDIAN'] = str(statistics.median(read_lengths))
    TELS_statistics['READ_LENGTH_RANGE'] = str((min(read_lengths), max(read_lengths)))
    if len(read_lengths) >= 2:
        TELS_statistics['READ_LENGTH_STD_DEV'] = str(statistics.stdev(read_lengths))
        TELS_statistics['READ_LENGTH_VARIANCE'] = str(statistics.variance(read_lengths))
        TELS_statistics['READ_LENGTH_SKEW'] = str(skew(read_lengths))
        TELS_statistics['READ_LENGTH_KURTOSIS'] = str(kurtosis(read_lengths))
    else:
        TELS_statistics['READ_LENGTH_STD_DEV'] = 0
        TELS_statistics['READ_LENGTH_VARIANCE'] = 0
        TELS_statistics['READ_LENGTH_SKEW'] = 0
        TELS_statistics['READ_LENGTH_KURTOSIS'] = 0

    print(json.dumps(TELS_statistics))


if __name__ == '__main__':
    main()
