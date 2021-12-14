#!/usr/bin/env python
# coding: utf-8

import pysam
from Bio import SeqIO
import csv
import glob
import re
import os
import pandas as pd
import argparse
import configparser

from common import *


# Hardcoded stuff
base_path = '/home/marco/Downloads/Noyes_TELS_I/ZymoBIOMICS.STD.refseq.v2/Genomes'
alignment_dir = base_path + '/' + 'alignment_files'

MEGARES_ext = '_ato_MEGARES.sam'
MGES_ext = '_ato_MGES.sam'
KEGG_ext = '_ato_KEGG.sam'

reference_files_list = ['Bacillus_subtilis_complete_genome.fasta',
                        'Cryptococcus_neoformans_draft_genome.fasta',
                        'Enterococcus_faecalis_complete_genome.fasta',
                        'Escherichia_coli_complete_genome.fasta',
                        'Lactobacillus_fermentum_complete_genome.fasta',
                        'Listeria_monocytogenes_complete_genome.fasta',
                        'Pseudomonas_aeruginosa_complete_genome.fasta',
                        'Saccharomyces_cerevisiae_draft_genome.fasta',
                        'Salmonella_enterica_complete_genome.fasta',
                        'Staphylococcus_aureus_complete_genome.fasta']

alignment_files_base_path = '/home/marco/Downloads/Noyes_TELS_I/MOCK_ato_WGS'

sample_files_base_path = '/home/marco/Downloads/Noyes_TELS_I/'
sample_name = 'deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq'

COLOC_DISTANCE = 500


def read_db_lengths(db_path):
    out_dict = dict()
    with open(db_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            out_dict[record.name] = len(record.seq)
            
    return out_dict

class Gene:

    MGE_TYPE_STR = 'MGE'
    ARG_TYPE_STR = 'ARG'
    KEGG_TYPE_STR = 'KEGG'

    def __init__(self, gene_name, gene_type, start_pos, end_pos):
        self.name = gene_name
        self.type = gene_type
        self.start_pos = start_pos
        self.end_pos = end_pos

    def __str__(self):
        return self.name + ',' + self.type + ',' + str(self.start_pos) + ',' + str(self.end_pos)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.start_pos == other.start_pos and self.end_pos == other.end_pos

    def __lt__(self, other):
        if self.start_pos < other.start_pos:
            return True
        if self.end_pos < other.end_pos:
            return True
        return False

    def __le__(self, other):
        if self.start_pos <= other.start_pos:
            return True
        if self.end_pos <= other.end_pos:
            return True
        return False

    def __gt__(self, other):
        return not (self.__le__(other))

    def __ge__(self, other):
        return not (self.__lt__(other))

    def length(self):
        return self.end_pos - self.start_pos

    def is_arg(self):
        return self.type == Gene.ARG_TYPE_STR

    def is_mge(self):
        return self.type == Gene.MGE_TYPE_STR

    def is_kegg(self):
        return self.type == Gene.KEGG_TYPE_STR


def read_sam(file_name, genes_lengths, threshold, gene_type='unknown'):
    alignments_list = list()
    with pysam.AlignmentFile(file_name) as sam_file:
        for read in sam_file:
            if read.is_unmapped or read.is_secondary:
                continue
            if float(read.query_alignment_length) / float(genes_lengths[read.query_name]) > threshold:
                alignments_list.append(Gene(read.query_name, gene_type, read.reference_start, read.reference_end))

    return alignments_list


def gen_resistome(config, organism_name, genes_list):
    megares_ontology, _ = read_megares_ontology(config)

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    for gene in genes_list:
        if gene.is_arg():
            classname = megares_ontology[gene.name]["class"]
            mech = megares_ontology[gene.name]["mechanism"]
            group = megares_ontology[gene.name]["group"]

            # update gene dict
            if gene.name not in gene_dict:
                gene_dict[gene.name] = 1
            else:
                gene_dict[gene.name] += 1

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

    # Prepare rows of diversity csv
    csv_rows = list()
    csv_rows.append(
        ["MEGARes Gene Header", "Num Matches", "Group", "Num Matches", "Mechanism", "Num Matches", "Class",
         "Num Matches"])
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
    with open(config['OUTPUT']['OUT_FILE_PREFIX'] + organism_name + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_dict), len(class_dict), len(mech_dict), len(group_dict)])

    # Write richness tsv
    with open(config['OUTPUT']['OUT_FILE_PREFIX'] + organism_name + "_amr_richness.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)
        
def gen_mobilome(config, organism_name, genes_list):
    gene_dict = dict()
    
    for gene in genes_list:
        if gene.is_mge():
            if gene.name not in gene_dict:
                gene_dict[gene.name] = 1
            else:
                gene_dict[gene.name] += 1

    # Prepare rows of tsv
    csv_rows = list()
    gene_riches = [(header, gene_dict[header]) for header in gene_dict]

    # Output how many different MGEs are in the data, this is actually diversity!
    csv_rows.append(["MGE Richness:", len(gene_riches)])

    # Column headers
    csv_rows.append(["MGE Header", "Num Matches"])

    # Output individual MGEs and their conts
    for gene_count_tuple in sorted(gene_riches, key=lambda gene_count_tuple: gene_count_tuple[1], reverse=True):
        csv_rows.append([gene_count_tuple[0], gene_count_tuple[1]])

    # Write csv
    with open(config['OUTPUT']['OUT_FILE_PREFIX'] + organism_name + "_mobilome.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)


def gen_colocalizations(config, organism_name, genes_list):
    gene_it = 1
    colocalizations = list()
    if len(genes_list) == 0:
        return 
    curr_colocalization = [genes_list[0]]
    while (gene_it < len(genes_list)):
        if genes_list[gene_it].start_pos - curr_colocalization[-1].end_pos < COLOC_DISTANCE:
            curr_colocalization.append(genes_list[gene_it])
        else:
            if len(curr_colocalization) > 1:
                colocalizations.append(curr_colocalization)
                curr_colocalization = [genes_list[gene_it]]
            else:
                curr_colocalization = [genes_list[gene_it]]
        gene_it += 1
    
    
    # Output colocalizations
    with open(config['OUTPUT']['OUT_FILE_PREFIX'] + organism_name + '_colocalizations.csv', 'wt') as csvfile:
        writer = csv.writer(csvfile)
        header = ['ARG', 'ARG positions', 'MGE(s)', 'MGE(s) positions']
        writer.writerow(header)
        
        for coloc in colocalizations:
            row = list()
            mges_list = list()
            args_list = list()
            for gene in coloc:
                if gene.is_arg():
                    args_list.append(gene)
                elif gene.is_mge():
                    mges_list.append(gene)
                if len(args_list) != 0:
                    args = ""
                    args_positions = ""
                    for gene in args_list:
                        if args != "":
                            args += ';'
                            args_positions += ';'
                        args += gene.name
                        args_positions += "{}:{}".format(gene.start_pos, gene.end_pos)

                if len(mges_list) != 0:
                    mges = ""
                    mges_positions = ""
                    for gene in mges_list:
                        if mges != "":
                            mges += ';'
                            mges_positions += ';'
                        mges += gene.name
                        mges_positions += "{}:{}".format(gene.start_pos, gene.end_pos)
            if len(args_list) != 0 and len(mges_list) != 0:
                row.extend([args, args_positions, mges, mges_positions])
                writer.writerow(row)
                

def merge_annotations(A, B):
    out = list()
    A_it = 0
    B_it = 0
    while A_it < len(A) and B_it < len(B):
        if A[A_it] <= B[B_it]:
            out.append(A[A_it])
            A_it += 1
        else:
            out.append(B[B_it])
            B_it += 1

    while A_it < len(A):
        out.append(A[A_it])
        A_it += 1
    while B_it < len(B):
        out.append(B[B_it])
        B_it += 1
        
    return out


# Get per read position
def get_positions(file_name, read_set):
    reads_positions = dict()
    with pysam.AlignmentFile(file_name) as sam_file:
        for read in sam_file:
            if read.is_unmapped or read.is_secondary:
                continue
            if read.query_name in read_set:
                reads_positions[read.query_name] = (read.reference_start, read.reference_end)
    return reads_positions


# It is a linear scan at the moment, make it at least log. It's a sorted array of ranges.
def search_annotation(annotation, start, end):
    ann_it = 0
    genes_in_range = list()
    while ann_it < len(annotation):
        current_gene = annotation[ann_it]
        if current_gene.start_pos > start:
            genes_in_range.append(current_gene)
        elif current_gene.end_pos > start:
            genes_in_range.append(current_gene)
        if current_gene.end_pos > end:
            break
        ann_it += 1
    return genes_in_range

def main():
    parser = argparse.ArgumentParser(description='WGS Analysis, many hardcoded things')
    parser.add_argument('-c', help='Config file', dest='config_path', default='./config_test.ini')
    parser.add_argument('-o', help="Output file", default='', dest='out_file', type=str)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    config['OUTPUT'] = dict()
    config['OUTPUT']['OUT_FILE_PREFIX'] = os.path.abspath(args.out_file)

    mges_lengths = read_db_lengths(config['DATABASE']['MGES'])
    args_lengths = read_db_lengths(config['DATABASE']['MEGARES'])

    complete_genomes_annotations = dict()
    for file_name in reference_files_list:
        print('Working on {}'.format(file_name))
        mges_list = read_sam(alignment_dir + '/' + file_name + MGES_ext, mges_lengths, 0.5, Gene.MGE_TYPE_STR)
        args_list = read_sam(alignment_dir + '/' + file_name + MEGARES_ext, args_lengths, 0.8, Gene.ARG_TYPE_STR)

        genome_annotation = list()

        # Args
        for arg in sorted(args_list):
            if len(genome_annotation) == 0:
                genome_annotation.append(arg)
                continue

            if not arg.start_pos < genome_annotation[-1].end_pos:
                genome_annotation.append(arg)

        # MGEs
        mge_annotation = list()
        mges_it = 0
        sorted_mges_list = sorted(mges_list)
        for arg in genome_annotation:
            while (mges_it < len(sorted_mges_list)) and (sorted_mges_list[mges_it].end_pos < arg.start_pos):
                if len(mge_annotation) == 0:
                    mge_annotation.append(sorted_mges_list[mges_it])
                    continue

                if not sorted_mges_list[mges_it].start_pos < mge_annotation[-1].end_pos:
                    mge_annotation.append(sorted_mges_list[mges_it])

                mges_it += 1

        # Merge ARGs and MGEs
        complete_genome_annotation = merge_annotations(genome_annotation, mge_annotation)
        complete_genomes_annotations[file_name] = complete_genome_annotation

        # Extract resistome, mobilome and colocalizations and output to file
#        gen_resistome(config, file_name, complete_genome_annotation)
#        gen_mobilome(config, file_name, complete_genome_annotation)
#        gen_colocalizations(config, file_name, complete_genome_annotation)


    # Get sample colocalizations
    S_colocalizations_list = list()
    S_reads_set = set()
    S_colocalizations_path = '/home/marco/Downloads/Noyes_TELS_I/deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq_colocalizations.csv'
    with open(S_colocalizations_path, 'rt') as coloc_csv_hanlde:
        csv_reader = csv.reader(coloc_csv_hanlde)
        header = next(csv_reader)
        for row in csv_reader:
            read = row[0]
            args_list = row[1].split(';')
            args_pos = row[2].split(';')
            mges_list = row[3].split(';')
            mges_pos = row[4].split(';')

            S_colocalizations_genes_list = list()
            for idx, arg in enumerate(args_list):
                S_colocalizations_genes_list.append(
                    Gene(arg, Gene.ARG_TYPE_STR, args_pos[idx].split(':')[0], args_pos[idx].split(':')[1]))

            for idx, mge in enumerate(mges_list):
                S_colocalizations_genes_list.append(
                    Gene(mge, Gene.MGE_TYPE_STR, mges_pos[idx].split(':')[0], mges_pos[idx].split(':')[1]))

            S_colocalizations_list.append((read, sorted(S_colocalizations_genes_list)))
            S_reads_set.add(read)

    # Get per read positions
    reads_positions = dict()
    for reference_file in reference_files_list:
        sam_file_name = alignment_files_base_path + '/' + sample_name + '_ato_' + reference_file + '.sam'
        reads_positions[reference_file] = dict()
        reads_positions[reference_file].update(get_positions(sam_file_name, S_reads_set))


    # Search colocalizations in reference genomes
    for read, colocalization in S_colocalizations_list:
        # Go trough all the references and look for the colocalization
        for reference_file in reference_files_list:
            ref_annotation = complete_genomes_annotations[reference_file]
            if read in reads_positions[reference_file]:
                read_pos = reads_positions[reference_file][read]
                genes = search_annotation(ref_annotation, read_pos[0], read_pos[1])

                adjusted_colocalization = list()
                for gene in colocalization:
                    adjusted_colocalization.append(Gene(gene.name,
                                                        gene.type,
                                                        int(gene.start_pos) + int(read_pos[0]),
                                                        int(gene.end_pos) + int(read_pos[0])))

                # Check if genes and colocalization align
                print(reference_file[:-6] + '\n{' + str(genes) + '}\n{' + str(adjusted_colocalization) + '}\n====================')


if __name__ == '__main__':
    main()
