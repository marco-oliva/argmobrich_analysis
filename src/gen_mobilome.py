from Bio import SeqIO
import pysam
import argparse
import csv
import configparser
import numpy as np
import json
import os
from common import *


def long_reads_strategy(config):
    mges_gene_lengths = dict()
    mge_combined_reference_fasta_filename = config['DATABASE']['MGES']

    for rec in SeqIO.parse(mge_combined_reference_fasta_filename, "fasta"):
        mges_gene_lengths[rec.name] = len(rec.seq)

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    reads_aligned_to_args = dict()
    if config['INPUT']['ARGS_SAM_FILE'] != '':
        args_gene_lengths = dict()
        args_reference_fasta_filename = config['DATABASE']['MEGARES']

        for rec in SeqIO.parse(args_reference_fasta_filename, "fasta"):
            args_gene_lengths[rec.name] = len(rec.seq)

        args_sam_file = pysam.AlignmentFile(config['INPUT']['ARGS_SAM_FILE'], 'r')

    mges_sam_file = pysam.AlignmentFile(config['INPUT']['MGES_SAM_FILE'], 'r')

    gene_dict = dict()
    reads_aligned = set()
    # Iterate through every aligned segment
    for read in mges_sam_file.fetch():
        if read.is_unmapped:
            continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary:
            continue

        # check coverage
        if (float(read.reference_length) / mges_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_MGE_THRESHOLD']):
            if not read.reference_name in gene_dict:
                gene_dict[read.reference_name] = 1
            else:
                gene_dict[read.reference_name] += 1
            reads_aligned.add(read.query_name)

    # Prepare rows of tsv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['MGES_' + stat_name, stat_value])

    csv_rows.append(['Mobimole'])

    gene_riches = [(header, gene_dict[header]) for header in gene_dict]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(gene_riches)])

    # Column headers
    csv_rows.append(["MGE Header", "Num Reads"])

    # Output individual MGEs and their conts
    for gene_count_tuple in sorted(gene_riches, key=lambda gene_count_tuple: gene_count_tuple[1], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

def short_reads_stratedy(config):
    mge_genes = dict()
    mge_combined_reference_fasta_filename = config['DATABASE']['MGES']
    for rec in SeqIO.parse(mge_combined_reference_fasta_filename, "fasta"):
        mge_genes[rec.name] = np.zeros(len(rec.seq))

    mges_sam_file = pysam.AlignmentFile(config['INPUT']['MGES_SAM_FILE'], 'r')

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    gene_hits = dict()
    reads_aligned_per_gene = dict()
    # Iterate through every read. Accumulate number of reads aligned and number of alignments per aclame mge
    for read in mges_sam_file.fetch():
        if read.is_unmapped:
            continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary:
            continue

        for i in range(read.reference_start, read.reference_end):
            mge_genes[read.reference_name][i] = 1

        if read.reference_name not in gene_hits:
            gene_hits[read.reference_name] = 1
            reads_aligned_per_gene[read.reference_name] = set()
            reads_aligned_per_gene[read.reference_name].add(read.query_name)
        else:
            gene_hits[read.reference_name] += 1
            reads_aligned_per_gene[read.reference_name].add(read.query_name)

    # check coverage
    covered_genes = set()
    reads_aligned = set()
    for mge_gene, coverage_vector in mge_genes.items():
        if float(sum(coverage_vector) / len(coverage_vector)) > float(config['MISC']['GLOBAL_MGE_THRESHOLD']):
            covered_genes.add(mge_gene)
            reads_aligned.update(reads_aligned_per_gene[mge_gene])

    # Prepare rows of tsv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['MGES_' + stat_name, stat_value])

    csv_rows.append(['Mobimole'])
    covered_gene_richness = [(header, gene_hits[header]) for header in covered_genes]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(covered_gene_richness)])

    # Column headers
    csv_rows.append(["MGE Header", "Num Reads"])

    # Output individual MGEs and their counts
    for gene_count_tuple in sorted(covered_gene_richness, key=lambda gene_count_tuple: gene_count_tuple[1], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

def main():
    parser = argparse.ArgumentParser(description='Colocalizations Finder.')
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-m', help='MGEs alignment file', dest='mges_sam', required=True)
    parser.add_argument('-a', help='ARGs alignment file', dest='args_sam', required=True)
    parser.add_argument('-o', help='Output file prefix', dest='out_prefix', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    logger = init_logger()

    config['INPUT'] = dict()
    config['INPUT']['MGES_SAM_FILE'] = args.mges_sam
    config['INPUT']['ARGS_SAM_FILE'] = args.args_sam
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.reads_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.reads_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    config['OUTPUT'] = dict()
    config['OUTPUT']['OUTPUT_PREFIX'] = args.out_prefix
    config['OUTPUT']['OUT_DIR'] = os.path.dirname(os.path.abspath(config['OUTPUT']['OUTPUT_PREFIX']))

    strategy = config['MISC']['MOBILOME_STRATEGY']
    logger.info('Strategy: {}'.format(strategy))
    if strategy == 'SHORT':
        short_reads_stratedy(config)
    elif strategy == 'LONG':
        long_reads_strategy(config)
    else:
        logger.error("MOBILOME_STRATEGY must be in [LONG, SHORT], {} not supported".format(strategy))


if __name__ == "__main__":
    main()