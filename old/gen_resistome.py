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
    # Get megares lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    sam_file = pysam.AlignmentFile(config['INPUT']['SAM_FILE'], 'r')

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    ontology_filename = config['DATABASE']['MEGARES_ONTOLOGY']
    with open(ontology_filename, 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            # Fill in our dict
            megares_ontology[row[0]] = {"class": row[1],
                                        "mechanism": row[2],
                                        "group": row[3]
                                        }

    # Get reads lengths
    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    # Iterate through every read. Accumulate number of reads while recording read length
    reads_aligned = set()
    for read in sam_file.fetch():
        if read.is_unmapped:
            continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary:
            continue

        # check coverage
        if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > float(
                config['MISC']['GLOBAL_AMR_THRESHOLD']):
            classname = megares_ontology[read.reference_name]["class"]
            mech = megares_ontology[read.reference_name]["mechanism"]
            group = megares_ontology[read.reference_name]["group"]

            # update gene dict
            if read.reference_name not in gene_dict:
                gene_dict[read.reference_name] = 1
            else:
                gene_dict[read.reference_name] += 1

            # update class dict
            if classname not in class_dict:
                class_dict[classname] = 1
            else:
                class_dict[classname] += 1

            # update mechanism dict
            if mech not in mech_dict:
                mech_dict[mech] = 1
            else:
                mech_dict[mech] += 1

            # update group dict
            if group not in group_dict:
                group_dict[group] = 1
            else:
                group_dict[group] += 1

            reads_aligned.add(read.query_name)

    # Prepare rows of diversity csv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['ARG_' + stat_name, stat_value])

    csv_rows.append(['Resistome'])
    csv_rows.append(
        ["MEGARes Gene Header", "Num Reads", "Group", "Num Reads", "Mechanism", "Num Reads", "Class",
         "Num Reads"])
    start = len(csv_rows)
    for gene in sorted(gene_dict, key=lambda gene: gene_dict[gene], reverse=True):
        csv_rows.append([gene, gene_dict[gene]])

    i = start
    for group in sorted(group_dict, key=lambda group: group_dict[group], reverse=True):
        csv_rows[i].extend([group, group_dict[group]])
        i += 1

    i = start
    for mech in sorted(mech_dict, key=lambda mech: mech_dict[mech], reverse=True):
        csv_rows[i].extend([mech, mech_dict[mech]])
        i += 1

    i = start
    for class_name in sorted(class_dict, key=lambda class_name: class_dict[class_name], reverse=True):
        csv_rows[i].extend([class_name, class_dict[class_name]])
        i += 1

    # Write diversity tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_dict), len(class_dict), len(mech_dict), len(group_dict)])

    # Write richness tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_richness.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

def short_reads_stratedy(config):
    # Get megares gene for coverage
    megares_genes = {}
    megares_reference_fasta_filename = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_genes[rec.name] = np.zeros(len(rec.seq))

    sam_file = pysam.AlignmentFile(config['INPUT']['SAM_FILE'], 'r')

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    ontology_filename = config['DATABASE']['MEGARES_ONTOLOGY']
    with open(ontology_filename, 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            # FIll in our dict
            megares_ontology[row[0]] = {"class": row[1],
                                        "mechanism": row[2],
                                        "group": row[3]
                                        }

    reads_lengths = dict()
    with open(config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['READS_LENGTH'], 'rt') as reads_lengths_json_fp:
        reads_lengths = json.load(reads_lengths_json_fp)

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    reads_aligned_per_gene = dict()
    # Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file.fetch():
        if read.is_unmapped:
            continue

        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary:
            continue

        for i in range(read.reference_start, read.reference_end):
            megares_genes[read.reference_name][i] = 1

        classname = megares_ontology[read.reference_name]["class"]
        mech = megares_ontology[read.reference_name]["mechanism"]
        group = megares_ontology[read.reference_name]["group"]

        # update gene dict
        if (not read.reference_name in gene_dict):
            gene_dict[read.reference_name] = 1
            reads_aligned_per_gene[read.reference_name] = set()
            reads_aligned_per_gene[read.reference_name].add(read.query_name)
        else:
            gene_dict[read.reference_name] += 1
            reads_aligned_per_gene[read.reference_name].add(read.query_name)

        # update class dict
        if classname not in class_dict:
            class_dict[classname] = 1
        else:
            class_dict[classname] += 1

        # update mechanism dict
        if mech not in mech_dict:
            mech_dict[mech] = 1
        else:
            mech_dict[mech] += 1

        # update group dict
        if group not in group_dict:
            group_dict[group] = 1
        else:
            group_dict[group] += 1

    # check coverage
    covered_genes = set()
    reads_aligned = set()
    for megares_gene, coverage_vector in megares_genes.items():
        if (float(sum(coverage_vector) / len(coverage_vector)) > float(config['MISC']['GLOBAL_AMR_THRESHOLD'])):
            covered_genes.add(megares_gene)
            reads_aligned.update(reads_aligned_per_gene[megares_gene])


    # Get only covered genes
    gene_covered_dict = dict()
    class_covered_dict = dict()
    mech_covered_dict = dict()
    group_covered_dict = dict()

    for gene in covered_genes:
        classname = megares_ontology[gene]["class"]
        mech = megares_ontology[gene]["mechanism"]
        group = megares_ontology[gene]["group"]

        gene_covered_dict[gene]         = gene_dict[gene]
        class_covered_dict[classname]   = class_dict[classname]
        mech_covered_dict[mech]         = mech_dict[mech]
        group_covered_dict[group]       = group_dict[group]


    # Prepare rows of diversity csv
    csv_rows = list()
    csv_rows.append(['Statistics'])
    arg_containing_reads_stats = reads_statistics(reads_aligned, reads_lengths)
    for stat_name, stat_value in arg_containing_reads_stats.items():
        csv_rows.append(['ARG_' + stat_name, stat_value])

    csv_rows.append(['Resistome'])
    csv_rows.append(
        ["MEGARes Gene Header", "Num Reads", "Group", "Num Reads", "Mechanism", "Num Reads", "Class",
         "Num Reads"])
    start = len(csv_rows)
    for gene in sorted(gene_covered_dict, key=lambda gene: gene_covered_dict[gene], reverse=True):
        csv_rows.append([gene, gene_covered_dict[gene]])

    i = start
    for group in sorted(group_covered_dict, key=lambda group: group_covered_dict[group], reverse=True):
        csv_rows[i].extend([group, group_covered_dict[group]])
        i += 1

    i = start
    for mech in sorted(mech_covered_dict, key=lambda mech: mech_covered_dict[mech], reverse=True):
        csv_rows[i].extend([mech, mech_covered_dict[mech]])
        i += 1

    i = start
    for class_name in sorted(class_covered_dict, key=lambda class_name: class_covered_dict[class_name], reverse=True):
        csv_rows[i].extend([class_name, class_covered_dict[class_name]])
        i += 1

    # Write diversity tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_covered_dict), len(class_covered_dict), len(mech_covered_dict), len(group_covered_dict)])

    # Write richness tsv
    with open(config['OUTPUT']['OUTPUT_PREFIX'] + '_' + config['MISC']['RESISTOME_STRATEGY'] + "_amr_richness.csv",'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)


def main():
    parser = argparse.ArgumentParser(description='Compute resistome')
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-s', help='Alignment file', dest='sam_file', required=True)
    parser.add_argument('-o', help='Output Prefix', dest='out_prefix', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    logger = init_logger()

    config['INPUT'] = dict()
    config['INPUT']['SAM_FILE'] = args.sam_file
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.reads_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.reads_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    config['OUTPUT'] = dict()
    config['OUTPUT']['OUTPUT_PREFIX'] = args.out_prefix
    config['OUTPUT']['OUT_DIR'] = os.path.dirname(os.path.abspath(config['OUTPUT']['OUTPUT_PREFIX']))

    strategy = config['MISC']['RESISTOME_STRATEGY']
    logger.info('Strategy: {}'.format(strategy))
    if strategy == 'SHORT':
        short_reads_stratedy(config)
    elif strategy == 'LONG':
        long_reads_strategy(config)
    else:
        logger.error("RESISTOME_STRATEGY must be in [LONG, SHORT], {} not supported".format(strategy))


if __name__ == "__main__":
    main()

