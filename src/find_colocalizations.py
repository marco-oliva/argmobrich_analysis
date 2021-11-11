#!/usr/bin/env python3
import configparser
import sys
from Bio import SeqIO
import pysam
import itertools
import csv
import argparse
import gzip
import logging
import itertools

def populate_megares_ontology(config):
    hierarchy_dict = {}
    class_samples = {}
    mech_samples = {}
    group_samples = {}

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    with open(config['DATABASE']['MEGARES_ONTOLOGY'], 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            cl = row[1]
            mech = row[2]
            group = row[3]

            # Set up hiearachy dict. This will be our tree structure
            if not cl in hierarchy_dict:
                hierarchy_dict[cl] = {}

            if not mech in hierarchy_dict[cl]:
                hierarchy_dict[cl][mech] = []

            if not group in hierarchy_dict[cl][mech]:
                hierarchy_dict[cl][mech].append(group)

            # Make sure each label of each category is represented so we can color nodes appropriately
            class_samples[cl] = set()
            mech_samples[mech] = set()
            group_samples[group] = set()

            # FIll in our dict
            megares_ontology[row[0]] = {"class": cl,
                                        "mech": mech,
                                        "group": group
                                        }
    return megares_ontology


# Check if valid alignements
def not_overlapping(intervals_list, left_boundary=0, right_boundary=0):
    sorted_intervals = sorted(intervals_list)

    if sorted_intervals[0][0] < left_boundary:
        return False
    if sorted_intervals[-1][1] > right_boundary:
        return False

    curr_right_lim = 0
    for interval in sorted_intervals:
        if interval[0] < curr_right_lim:
            return False
        curr_right_lim = interval[1]
    return True


def get_colocalizations(config, reads_file_path, to_megares_path, to_aclme_path, to_iceberg_path, to_plasmids_path,
                        to_kegg_path, skip_begin, skip_end):
    logger = logging.getLogger()

    # Get read lengths
    logger.info("Reading FASTQ")
    reads_length = dict()
    if reads_file_path.endswith('gz'):
        reads_file_handle = gzip.open(reads_file_path, 'rt')
    else:
        reads_file_handle = open(reads_file_path, 'rt')

    fasta_exts = [".fasta", ".fa", ".fna", ".FASTA", ".FA", ".FNA"]
    fastq_exts = [".fastq", ".fq", ".FASTQ", ".FQ"]

    is_fastq = False
    for ext_fq in fastq_exts:
        if ext_fq in reads_file_path:
            is_fastq = True

    is_fasta = False
    if not is_fastq:
        for ext_fa in fasta_exts:
            if ext_fa in reads_file_path:
                is_fasta = True

    if is_fasta:
        file_format = "fasta"
    elif is_fastq:
        file_format = "fastq"
    else:
        print("Couldn't recognize reads file format")
        sys.exit(1)

    for record in SeqIO.parse(reads_file_handle, file_format):
        reads_length[record.name] = len(record.seq)

    reads_file_handle.close()

    # Get AMR genes lengths for coverage
    logger.info("Reading ARGS DB")
    megares_gene_lengths = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    # Get aclame, iceberg, and plasmid finder lengths for coverage
    logger.info("Reading MGES DB")
    mge_gene_lengths = dict()
    mge_gene_lengths['aclme'] = dict()
    aclame_reference_fasta_filename = config['DATABASE']['ACLAME']
    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        mge_gene_lengths['aclme'][rec.name] = len(rec.seq)

    mge_gene_lengths['iceberg'] = dict()
    iceberg_reference_fasta_filename = config['DATABASE']['ICEBERG']
    for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        mge_gene_lengths['iceberg'][rec.name] = len(rec.seq)

    mge_gene_lengths['plasmids'] = dict()
    plasmids_reference_fasta_filename = config['DATABASE']['PLASMIDS']
    for rec in SeqIO.parse(plasmids_reference_fasta_filename, "fasta"):
        mge_gene_lengths['plasmids'][rec.name] = len(rec.seq)

    kegg_gene_lengths = dict()
    kegg_reference_fasta_filename = config['DATABASE']['KEGG']
    for rec in SeqIO.parse(kegg_reference_fasta_filename, "fasta"):
        kegg_gene_lengths[rec.name] = len(rec.seq)

    # Open aligned to megares sam
    logger.info("Reading ARGS alignment file")
    amr_positions = dict()
    read_to_amr = dict()
    amr_to_generated_bases = dict()
    amr_length = dict()
    with pysam.AlignmentFile(to_megares_path, "rb") as to_megares_samfile:
        for read in to_megares_samfile:
            if (not read.is_unmapped and ((not read.is_secondary) and (not read.is_supplementary))):
                # Check coverage
                if ((read.reference_length / (megares_gene_lengths[read.reference_name])) > float(config['MISC']['GLOBAL_AMR_THRESHOLD'])):
                    if (read.query_name not in read_to_amr):
                        read_to_amr[read.query_name] = list()
                        amr_positions[read.query_name] = list()
                    amr_positions[read.query_name].append([read.query_alignment_start, read.query_alignment_end])
                    read_to_amr[read.query_name].append(read.reference_name)

                    amr_length[read.reference_name] = megares_gene_lengths[read.reference_name]

                    if (read.reference_name not in amr_to_generated_bases):
                        amr_to_generated_bases[read.reference_name] = read.reference_length
                    else:
                        amr_to_generated_bases[read.reference_name] += read.reference_length

    # Open aligned to Kegg
    logger.info("Reading MGES alignment files")
    kegg_positions = dict()
    read_to_kegg = dict()
    with pysam.AlignmentFile(to_kegg_path, "rb") as to_kegg_samfile:
        for read in to_kegg_samfile:
            if (not read.is_unmapped):
                # Check coverage
                if ((read.reference_length / (kegg_gene_lengths[read.reference_name])) > float(config['MISC']['GLOBAL_KEGG_THRESHOLD'])):
                    if (read.query_name not in read_to_kegg):
                        read_to_kegg[read.query_name] = list()
                        kegg_positions[read.query_name] = list()
                    kegg_positions[read.query_name].append([read.query_alignment_start, read.query_alignment_end])
                    read_to_kegg[read.query_name].append(read.reference_name)

    # Mobile elements dict
    mge_positions = dict()
    read_to_mges = dict()
    mge_alignment_files = {'aclme': to_aclme_path, 'iceberg': to_iceberg_path, 'plasmids': to_plasmids_path}
    for mge_db_name, mge_alignment_file_path in mge_alignment_files.items():
        mge_positions[mge_db_name] = dict()
        read_to_mges[mge_db_name] = dict()

        with pysam.AlignmentFile(mge_alignment_file_path, "rb") as mge_alignment_file:
            for read in mge_alignment_file:
                if (not read.is_unmapped):

                    # Check coverage
                    gene_length = mge_gene_lengths[mge_db_name][read.reference_name]

                    if ((read.reference_length / gene_length) > float(config['MISC']['GLOBAL_MGE_THRESHOLD'])):
                        if (read.query_name not in read_to_mges[mge_db_name]):
                            read_to_mges[mge_db_name][read.query_name] = list()
                            mge_positions[mge_db_name][read.query_name] = list()
                        mge_positions[mge_db_name][read.query_name].append(
                            [read.query_alignment_start, read.query_alignment_end])
                        read_to_mges[mge_db_name][read.query_name].append(read.reference_name)

    genes_lists = dict()
    genes_list_csvs = reads_file_path + config['EXTENSION']['GENES_LIST']
    with open(genes_list_csvs, 'w') as genes_list_handle:
        writer = csv.writer(genes_list_handle)
        header = ['Read Name', 'AMR Genes', 'MGE Genes', 'KEGG genes']
        writer.writerow(header)

        for read, amr_list in read_to_amr.items():
            # Amr genes list
            amr_genes_list = list()
            for idx, amr_name in enumerate(amr_list):
                amr_genes_list.append([amr_name, amr_positions[read][idx], 'amr'])

            # Kegg genes list
            kegg_genes_list = list()
            if (read in read_to_kegg):
                for idx, kegg_name in enumerate(read_to_kegg[read]):
                    kegg_genes_list.append([kegg_name, kegg_positions[read][idx], 'kegg'])

            # MGE genes list
            mge_genes_list = list()
            for db_name, read_to_mge_db in read_to_mges.items():
                if (read in read_to_mge_db):
                    for idx, mge_name in enumerate(read_to_mge_db[read]):
                        mge_genes_list.append([mge_name, mge_positions[db_name][read][idx], 'mge'])

            genes_lists[read] = list()
            genes_lists[read].extend(amr_genes_list)
            genes_lists[read].extend(kegg_genes_list)
            genes_lists[read].extend(mge_genes_list)

            amr_genes_concatenated = ''
            if len(amr_genes_list) > 0:
                amr_genes_concatenated = ';'.join(list(itertools.chain(*amr_genes_list)))

            mge_genes_concatenated = ''
            if len(mge_genes_list) > 0:
                mge_genes_concatenated = ';'.join(list(itertools.chain(*mge_genes_list)))

            kegg_genes_concatenated = ''
            if len(kegg_genes_list) > 0:
                kegg_genes_concatenated = ';'.join(list(itertools.chain(*kegg_genes_list)))

            row = [read, amr_genes_concatenated, mge_genes_concatenated, kegg_genes_concatenated]
            writer.writerow(row)

    # Candidate colocalizations
    candidate_colocalizations = {k: v for k, v in genes_lists.items() if len(v) >= 2}

    colocalizations = dict()
    for read, candidate_colocalization_list in candidate_colocalizations.items():
        if (read not in reads_length):
            logger.error("{} not in read lengths.".format(read))
            continue
        for i in range(1, len(candidate_colocalization_list)):
            for coloc in itertools.combinations(candidate_colocalization_list[1:], i):
                intervals = list()
                intervals.append(candidate_colocalization_list[0][1])
                for gene in coloc:
                    intervals.append(gene[1])

                if not_overlapping(intervals, skip_begin, reads_length[read] - (skip_end + 1)):
                    if read not in colocalizations:
                        colocalizations[read] = list()
                    colocalization = list()
                    colocalization.append(candidate_colocalization_list[0])
                    colocalization.extend(coloc)

                    colocalization.sort(key=lambda x: x[1][0])
                    colocalizations[read].append(colocalization)

    arg_set = set()
    for read, args_list in read_to_amr.items():
        for arg in args_list:
            arg_set.add(arg)

    mge_set = set()
    for read, mges_list in read_to_mges.items():
        for mge in mges_list:
            mge_set.add(mge)

    logger.info("MGES: {}\tARGS: {}".format(len(mge_set), len(arg_set)))
    return arg_set, mge_set, colocalizations

############################################################
# Compute values

def main():

    parser = argparse.ArgumentParser(description='Colocalizations Finder.')
    parser.add_argument('-k', help='Kegg alignment file', dest='kegg_sam', required=True)
    parser.add_argument('-p', help='Plasmids Finder alignment file', dest='plasmids_sam', required=True)
    parser.add_argument('-i', help='Iceberg alignment file', dest='iceberg_sam', required=True)
    parser.add_argument('-m', help='Megares alignment file', dest='megares_sam', required=True)
    parser.add_argument('-a', help='ACLAME alignment file', dest='aclame_sam', required=True)
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-e', help='Skip n bases at the end of the read', dest='skip_end', default=0, type=int)
    parser.add_argument('-b', help='Skip n bases at the beginning of the read', dest='skip_begin', default=0, type=int)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    amr_set, mge_set, sample_colocalizations = get_colocalizations(config,
                                                                   args.reads_file,
                                                                   args.megares_sam,
                                                                   args.aclame_sam,
                                                                   args.iceberg_sam,
                                                                   args.plasmids_sam,
                                                                   args.kegg_sam,
                                                                   args.skip_begin,
                                                                   args.skip_end)

    logger = logging.getLogger()
    logger.info("Found {} reads with colocalizations".format(len(sample_colocalizations)))

    # Output colocalizations to stdout
    csvfile = sys.stdout
    writer = csv.writer(csvfile)
    header = ['read', 'ARG', 'ARG positions', 'MGE(s)', 'MGE(s) positions', 'KEGG(s)', 'KEGG(s) positions']
    writer.writerow(header)
    for read, possible_genes_list in sample_colocalizations.items():
        for genes_list in possible_genes_list:
            row = list()
            row.append(read)
            # ARG
            arg = ""
            arg_position = ""
            for gene in genes_list:
                if gene[2] == 'amr':
                    arg = gene[0]
                    arg_position = "{}:{}".format(gene[1][0], gene[1][1])

            # Mges
            mges = ""
            mges_positions = ""
            for gene in genes_list:
                if gene[2] == 'mge':
                    if mges != "":
                        mges += ';'
                        mges_positions += ';'
                    mges += gene[0]
                    mges_positions += "{}:{}".format(gene[1][0], gene[1][1])

            # Kegg
            keggs = ""
            keggs_positions = ""
            for gene in genes_list:
                if gene[2] == 'kegg':
                    if keggs != "":
                        keggs += ';'
                        keggs_positions += ';'
                    keggs += gene[0]
                    keggs_positions += "{}:{}".format(gene[1][0], gene[1][1])

            row.extend([arg, arg_position, mges, mges_positions, keggs, keggs_positions])
            writer.writerow(row)


if __name__ == "__main__":
    main()