from Bio import SeqIO
import pysam
import argparse
import csv
import configparser
import numpy as np
from common import *


def long_reads_strategy(config):
    # Get aclame, iceberg, and plasmid finder lengths for coverage
    aclame_gene_lengths = {}
    iceberg_gene_lengths = {}
    pf_gene_lengths = {}
    aclame_reference_fasta_filename = config['DATABASE']['ACLAME']
    iceberg_reference_fasta_filename = config['DATABASE']['ICEBERG']
    pf_reference_fasta_filename = config['DATABASE']['PLASMIDS']

    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        aclame_gene_lengths[rec.name] = len(rec.seq)

    for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        iceberg_gene_lengths[rec.name] = len(rec.seq)

    for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
        pf_gene_lengths[rec.name] = len(rec.seq)

    aclame_sam_file = pysam.AlignmentFile(config['INPUT']['ACLAME_SAM_FILE'], 'r')
    iceberg_sam_file = pysam.AlignmentFile(config['INPUT']['ICEBERG_SAM_FILE'], 'r')
    pf_sam_file = pysam.AlignmentFile(config['INPUT']['PLASMIDS_SAM_FILE'], 'r')

    reads_aligned = set()
    gene_dict = {}
    # Iterate through every read. Accumulate number of reads aligned and number of alignments per aclame mge
    for read in aclame_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        # check coverage
        if (float(read.reference_length) / aclame_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_MGE_THRESHOLD']):
            if (not read.reference_name in gene_dict):
                gene_dict[read.reference_name] = ["aclame", 1]
            else:
                gene_dict[read.reference_name][1] += 1

            if (not read.query_name in reads_aligned):
                reads_aligned[read.query_name] = True

    # Iterate through every read. Accumulate number of reads aligned and number of alignments per iceberg mge
    for read in iceberg_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        # check coverage
        if (float(read.reference_length) / iceberg_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_MGE_THRESHOLD']):
            if (not read.reference_name in gene_dict):
                gene_dict[read.reference_name] = ["iceberg", 1]
            else:
                gene_dict[read.reference_name][1] += 1

            if (not read.query_name in reads_aligned):
                reads_aligned[read.query_name] = True

    # Iterate through every read. Accumulate number of reads aligned and number of alignments per plasmidfinder mge
    for read in pf_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        # check coverage
        if (float(read.reference_length) / pf_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_MGE_THRESHOLD']):
            if (not read.reference_name in gene_dict):
                gene_dict[read.reference_name] = ["pf", 1]
            else:
                gene_dict[read.reference_name][1] += 1

            if (not read.query_name in reads_aligned):
                reads_aligned[read.query_name] = True

    # Prepare rows of tsv
    # csv_rows = [[os.path.basename(sys.argv[1])[0:6]]] #First row is library label
    csv_rows = list()
    gene_riches = [(header, gene_dict[header][0], gene_dict[header][1]) for header in gene_dict]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(gene_riches)])

    # Output how many reads aligned to an MGE
    csv_rows.append(["Number of Reads Aligned to MGE Databases:", len(reads_aligned)])

    # Column headers
    csv_rows.append(["MGE Database", "MGE Header", "Num Reads"])

    # Output individual MGEs and their conts
    for gene_count_tuple in sorted(gene_riches, key=lambda gene_count_tuple: gene_count_tuple[2], reverse=True):
        csv_rows.append([gene_count_tuple[1], gene_count_tuple[0], gene_count_tuple[2]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

def short_reads_stratedy(config):
    mge_genes = {}

    aclame_reference_fasta_filename = config['DATABASE']['ACLAME']
    iceberg_reference_fasta_filename = config['DATABASE']['ICEBERG']
    pf_reference_fasta_filename = config['DATABASE']['PLASMIDS']

    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        mge_genes[('ACLAME', rec.name)] = np.zeros(len(rec.seq))

    for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        mge_genes[('ICEBERG', rec.name)] = np.zeros(len(rec.seq))

    for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
        mge_genes[('PLASMID_FINDER', rec.name)] = np.zeros(len(rec.seq))

    aclame_sam_file = pysam.AlignmentFile(config['INPUT']['ACLAME_SAM_FILE'], 'r')
    iceberg_sam_file = pysam.AlignmentFile(config['INPUT']['ICEBERG_SAM_FILE'], 'r')
    pf_sam_file = pysam.AlignmentFile(config['INPUT']['PLASMIDS_SAM_FILE'], 'r')

    gene_hits = dict()
    reads_aligned = set()
    # Iterate through every read. Accumulate number of reads aligned and number of alignments per aclame mge
    for read in aclame_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        for i in range(read.reference_start, read.reference_end):
            mge_genes[('ACLAME', read.reference_name)][i] = 1

        if (not read.reference_name in gene_hits):
            gene_hits[read.reference_name] = 1
        else:
            gene_hits[read.reference_name] += 1

    # Iterate through every read. Accumulate number of reads aligned and number of alignments per iceberg mge
    for read in iceberg_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        for i in range(read.reference_start, read.reference_end):
            mge_genes[('ICEBERG', read.reference_name)][i] = 1

        if (not read.reference_name in gene_hits):
            gene_hits[read.reference_name] = 1
        else:
            gene_hits[read.reference_name] += 1

    # Iterate through every read. Accumulate number of reads aligned and number of alignments per plasmidfinder mge
    for read in pf_sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name in reads_aligned:
            continue
        else:
            reads_aligned.add(read.query_name)

        for i in range(read.reference_start, read.reference_end):
            mge_genes[('PLASMID_FINDER', read.reference_name)][i] = 1

        if (not read.reference_name in gene_hits):
            gene_hits[read.reference_name] = 1
        else:
            gene_hits[read.reference_name] += 1

    # check coverage
    covered_genes = set()
    for mge_gene, coverage_vector in mge_genes.items():
        if (float(sum(coverage_vector) / len(coverage_vector)) > float(config['MISC']['GLOBAL_MGE_THRESHOLD'])):
            covered_genes.add(mge_gene)

    # Prepare rows of tsv
    csv_rows = list()
    covered_gene_richness = [(db_name, header, gene_hits[header]) for db_name, header in covered_genes]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(covered_gene_richness)])

    # Column headers
    csv_rows.append(["MGE Database", "MGE Header", "Num Reads"])

    # Output individual MGEs and their counts
    for gene_count_tuple in sorted(covered_gene_richness, key=lambda gene_count_tuple: gene_count_tuple[2], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1], gene_count_tuple[2]])

    # Write csv
    filename_prefix = config['OUTPUT']['OUTPUT_PREFIX']
    with open(filename_prefix + '_' + config['MISC']['MOBILOME_STRATEGY'] + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

def main():
    parser = argparse.ArgumentParser(description='Colocalizations Finder.')
    parser.add_argument('-p', help='Plasmids Finder alignment file', dest='plasmids_sam', required=True)
    parser.add_argument('-i', help='Iceberg alignment file', dest='iceberg_sam', required=True)
    parser.add_argument('-a', help='ACLAME alignment file', dest='aclame_sam', required=True)
    parser.add_argument('-o', help='Output file prefix', dest='out_prefix', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    logger = init_logger()

    config['INPUT'] = dict()
    config['INPUT']['PLASMIDS_SAM_FILE'] = args.plasmids_sam
    config['INPUT']['ICEBERG_SAM_FILE'] = args.iceberg_sam
    config['INPUT']['ACLAME_SAM_FILE'] = args.aclame_sam
    config['OUTPUT'] = dict()
    config['OUTPUT']['OUTPUT_PREFIX'] = args.out_prefix

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