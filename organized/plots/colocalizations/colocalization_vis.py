#!/usr/bin/env python
# coding: utf-8

import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import seaborn as sns
import csv
import requests

# Read megares ontology
hierarchy_dict = {}
class_samples = {}
mech_samples = {}
group_samples = {}

#Create ontology dictionary from MEGARes ontology file
megares_ontology = {}
ontology_filename = "/home/noyes046/jsettle/argmobrich/MEGARESONTOLOGY.tsv"
with open(ontology_filename, 'r') as ontology_tsv:
    ontology_reader = csv.reader(ontology_tsv, delimiter='\t')
    for row in ontology_reader:
        #Skip column names
        if row[0] == "header":
            continue


        cl =  row[1]
        mech = row[2]
        group = row[3]

        #Set up hiearachy dict. This will be our tree structure
        if not cl in hierarchy_dict:
            hierarchy_dict[cl] = {}

        if not mech in hierarchy_dict[cl]:
            hierarchy_dict[cl][mech] = []

        if not group in hierarchy_dict[cl][mech]:
            hierarchy_dict[cl][mech].append(group)

        #Make sure each label of each category is represented so we can color nodes appropriately
        class_samples[cl] = set()
        mech_samples[mech] = set()
        group_samples[group] = set()

        #FIll in our dict
        megares_ontology[row[0]] = { "class"        : cl,
                                     "mech"         : mech,
                                     "group"        : group
                                   }

# Check if valid alignements
def not_overlapping(intervals_list):
    sorted_intervals = sorted(intervals_list)
    curr_right_lim = 0
    for interval in sorted_intervals:
        if (interval[0] < curr_right_lim):
            return False
        curr_right_lim = interval[1]
    return True

def get_colocalizations(reads_file_path, to_megares_path, to_aclme_path, to_iceberg_path, to_plasmids_path, to_kegg_path):
    # Get read lengths
    reads_length = dict()
    with open(reads_file_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            reads_length[record.id] = len(record.seq)

    tot_bases_in_sample = sum(reads_length.values())

    # Get AMR genes lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    # Get aclame, iceberg, and plasmid finder lengths for coverage
    mge_gene_lengths = dict()
    mge_gene_lengths['aclme'] = dict()
    aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        mge_gene_lengths['aclme'][rec.name] = len(rec.seq)

    mge_gene_lengths['iceberg'] = dict()
    iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
    for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        mge_gene_lengths['iceberg'][rec.name] = len(rec.seq)

    mge_gene_lengths['plasmids'] = dict()
    plasmids_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"
    for rec in SeqIO.parse(plasmids_reference_fasta_filename, "fasta"):
        mge_gene_lengths['plasmids'][rec.name] = len(rec.seq)

    # Get Kegg
    kegg_gene_lengths = dict()
    kegg_reference_fasta_filename = "/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/cargo_genes/kegg_prokaryotes.fasta"
    for rec in SeqIO.parse(kegg_reference_fasta_filename, "fasta"):
        kegg_gene_lengths[rec.name] = len(rec.seq)

    # Thresholds
    amr_threshold  = 0.8
    mge_threshold  = 0.8
    kegg_threshold = 0.8

    # Open aligned to megares sam
    amr_positions = dict()
    read_to_amr = dict()
    amr_to_generated_bases = dict()
    amr_length = dict()
    with pysam.AlignmentFile(to_megares_path, "rb") as to_megares_samfile:
        for read in to_megares_samfile: 
            if (not read.is_unmapped and ((not read.is_secondary) and (not read.is_supplementary))):
            #if (not read.is_unmapped):
                # Check coverage
                if ((read.reference_length / (megares_gene_lengths[read.reference_name])) > amr_threshold):
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

    # Open aligned to Kegg Prokaryotes
    kegg_positions = dict()
    read_to_kegg = dict()
    with pysam.AlignmentFile(to_kegg_path, "rb") as to_kegg_samfile:
        for read in to_kegg_samfile: 
            if (not read.is_unmapped):
                # Check coverage
                if ((read.reference_length / (kegg_gene_lengths[read.reference_name])) > kegg_threshold):
                    if (read.query_name not in read_to_kegg):
                        read_to_kegg[read.query_name] = list()
                        kegg_positions[read.query_name] = list()
                    kegg_positions[read.query_name].append([read.query_alignment_start, read.query_alignment_end])
                    read_to_kegg[read.query_name].append(read.reference_name) 

    # Mobile elements dict
    mge_positions = dict()
    read_to_mges = dict()
    mge_alignment_files = {'aclme' : to_aclme_path, 'iceberg' : to_iceberg_path, 'plasmids' : to_plasmids_path}
    for mge_db_name, mge_alignment_file_path in mge_alignment_files.items():
        mge_positions[mge_db_name] = dict()
        read_to_mges[mge_db_name] = dict()

        with pysam.AlignmentFile(mge_alignment_file_path, "rb") as mge_alignment_file:
            for read in mge_alignment_file:
                if (not read.is_unmapped):

                    # Check coverage
                    gene_length = mge_gene_lengths[mge_db_name][read.reference_name]

                    if ((read.reference_length / gene_length) > mge_threshold):
                        if (read.query_name not in read_to_mges[mge_db_name]):
                            read_to_mges[mge_db_name][read.query_name] = list()
                            mge_positions[mge_db_name][read.query_name] = list()
                        mge_positions[mge_db_name][read.query_name].append([read.query_alignment_start, read.query_alignment_end])
                        read_to_mges[mge_db_name][read.query_name].append(read.reference_name) 
                        
                        
    # Get all possible valid colocalizations
    genes_lists = dict()
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

    # Candidate colocalizations
    candidate_colocalizations = {k: v for k, v in genes_lists.items() if len(v) >= 2}

    colocalizations = dict()
    for read, candidate_colocalization_list in candidate_colocalizations.items():
        for i in range(1, len(candidate_colocalization_list)):
            for coloc in itertools.combinations(candidate_colocalization_list[1:], i):
                intervals = list()
                intervals.append(candidate_colocalization_list[0][1])
                for gene in coloc:
                    intervals.append(gene[1])

                if (not_overlapping(intervals)):
                    if (read not in colocalizations):
                        colocalizations[read] = list()
                    colocalization = list()
                    colocalization.append(candidate_colocalization_list[0])
                    colocalization.extend(coloc)

                    colocalization.sort(key=lambda x: x[1][0])
                    colocalizations[read].append(colocalization)   
    
    # Extract closest ARG-MGE or MGE-ARG colocalizations, closest pair
    arg_mge_colocalizations = dict()
    for read, colocalizations_list in colocalizations.items():
        arg = None
        for gene in colocalizations_list[0]:
            if (gene[2] == 'amr'):
                arg = gene
                break

        closest_mge = None
        for gene in colocalizations_list[0]:
            if (gene[2] == 'mge'):
                if (closest_mge is None):
                    closest_mge = gene
                else:
                    if (min(abs(gene[1][0] - arg[1][1]), abs(gene[1][1] - arg[1][0])) < 
                        min(abs(closest_mge[1][0] - arg[1][1]), abs(closest_mge[1][1] - arg[1][0]))):
                        closest_mge = gene

        if ((arg is None) or (closest_mge is None)):
            continue
        else:
            if (arg[1][0] < closest_mge[1][0]):
                arg_mge_colocalizations[read] = [arg, closest_mge]
            else: 
                arg_mge_colocalizations[read] = [closest_mge, arg]
            continue
            
    # Colocalization incidence function values

    # Compute the values for the X axis
    X_axis_per_read = dict()
    for read, colocalization in arg_mge_colocalizations.items():
        distance = colocalization[1][1][0] - colocalization[0][1][1]
        X_axis_per_read[read] = distance / reads_length[read]

    # Compute ARGs frequencies
    amr_frequencies = dict()
    for read, amr_list in read_to_amr.items():
        if (amr_list[0] in amr_frequencies.keys()):
            amr_frequencies[amr_list[0]] = amr_frequencies[amr_list[0]] + 1
        else:
            amr_frequencies[amr_list[0]] = 1

    # Compute colocalized ARGs frequencies
    amr_col_frequencies = dict()
    for read, _ in arg_mge_colocalizations.items():
        amr = read_to_amr[read][0]
        if (amr in amr_col_frequencies.keys()):
            amr_col_frequencies[amr] = amr_col_frequencies[amr] + 1
        else:
            amr_col_frequencies[amr] = 1

    # Compute Y axis values per read
    Y_axis_per_read = dict()
    for read, _ in arg_mge_colocalizations.items():
        amr  = read_to_amr[read][0]
        Fac  = amr_col_frequencies[amr]
        Fanc = amr_frequencies[amr]
        SeqA = amr_to_generated_bases[amr]
        SeqT = tot_bases_in_sample
        Al   = amr_length[amr]

        Y_axis_per_read[read] = (Fac/Fanc) * (SeqA/SeqT) * (1/Al)

    x = list()
    y = list()
    for read, _ in arg_mge_colocalizations.items():
        x.append(X_axis_per_read[read])
        y.append(Y_axis_per_read[read])
     
    arg_set = set()
    for read, args_list in read_to_amr.items():
        for arg in args_list:
            arg_set.add(arg)
        
    mge_set = set()
    for read, mges_list in read_to_mges.items():
        for mge in mges_list:
            mge_set.add(mge)
    
    return reads_length, arg_set, mge_set, colocalizations, x, y


def get_paths(sample_letter):
    base_reads       = '/panfs/roc/groups/11/noyes046/jsettle/argmobrich/analysis/deduplication/deduplicated_sequel-demultiplex.1896_{}01.ccs.fastq'
    base_to_megares  = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/{}01_ato_megares.sam'
    base_to_aclame   = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/{}01_ato_aclame.sam'
    base_to_iceberg  = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/{}01_ato_iceberg.sam'
    base_to_plasmids = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/{}01_ato_plasmids.sam'
    base_to_kegg     = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/{}01_ato_kegg_prokaryotes.sam'
    reads_file_path  = base_reads.format(sample_letter)
    to_megares_path  = base_to_megares.format(sample_letter)
    to_aclme_path    = base_to_aclame.format(sample_letter)
    to_iceberg_path  = base_to_iceberg.format(sample_letter)
    to_plasmids_path = base_to_plasmids.format(sample_letter)
    to_kegg_path     = base_to_kegg.format(sample_letter)
    return reads_file_path, to_megares_path, to_aclme_path, to_iceberg_path, to_plasmids_path, to_kegg_path


# Compute values 
sample_name = 'B'
reads_length, amr_set, mge_set, sample_colocalizations, x, y = get_colocalizations(*get_paths(sample_name))

# Output colocalizations csv
with open('{}_colocalizations.csv'.format(sample_name), 'w') as csvfile:
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
                if (gene[2] == 'amr'):
                    arg = gene[0]
                    arg_position="{}:{}".format(gene[1][0], gene[1][1])

            # Mges
            mges = ""
            mges_positions = ""
            for gene in genes_list:
                if (gene[2] == 'mge'):
                    if (mges != ""):
                        mges += ';'
                        mges_positions += ';'
                    mges += gene[0]
                    mges_positions += "{}:{}".format(gene[1][0], gene[1][1])

            # Kegg
            keggs = ""
            keggs_positions = ""
            for gene in genes_list:
                if (gene[2] == 'kegg'):
                    if (keggs != ""):
                        keggs += ';'
                        keggs_positions += ';'
                    keggs += gene[0]
                    keggs_positions += "{}:{}".format(gene[1][0], gene[1][1])


            row.extend([arg, arg_position, mges, mges_positions, keggs, keggs_positions])
            writer.writerow(row)

# Print colocalizations that have at least one kegg, one amr, one mge
candidate_colocalizations = dict()
no_kegg = 0
no_kegg_lim = 4
for read, colocalizations in sample_colocalizations.items():
    for genes_list in colocalizations:
        contains_kegg = False
        contains_mge = False
        contains_arg = False
        for gene in genes_list:
            if (gene[2] == 'amr'):
                contains_arg = True
            if (gene[2] == 'mge'):
                contains_mge = True
            if (gene[2] == 'kegg'):
                contains_kegg = True
        if ((contains_arg and contains_mge and contains_kegg)):
            if (reads_length[read] < 10000):
                candidate_colocalizations[read] = genes_list


# Read annotations
annotations_df = pd.read_excel("../heatmaps/TELS_MGE_ACCESSIONS.xlsx")

mge_genes_names = dict()
mge_genes_types = dict()
mge_accessions = set()

for index, row in annotations_df.iterrows():
    mge_genes_names[str(row['ACCESSIONS'])] = str(row['MGE NAME'])
    mge_genes_types[str(row['MGE NAME'])] = str(row['MGE TYPE'])
    mge_accessions.add(str(row['ACCESSIONS']))


# Plot candidate visualizations
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

padding = 2
height = 2
curr_y = 1+ padding
plt.figure(figsize=(16,9))
currentAxis = plt.gca()
use_reads_lengths = list()
for read, genes_list in candidate_colocalizations.items():
    read_length = reads_length[read]
    use_reads_lengths.append(read_length)
    read_r = Rectangle((0, curr_y), read_length, height, facecolor='blue', alpha=0.8, edgecolor='black', label='Read')
    currentAxis.add_patch(read_r)
    
    lastx = 0
    for gene in genes_list:
        if (gene[2] == 'amr'):
            color = 'red'
            label = 'ARG'
            element_name = megares_ontology[gene[0]]['group']
        elif (gene[2] == 'mge'):
            color = 'green'
            label = 'MGE'
            if gene[0] in mge_genes_names:
                element_name = mge_genes_names[gene[0]]
            else: 
                print('Da annotare:', gene[0])
                element_name = 'ANNOTAME QUESTO'
        else:
            color = 'orange'
            label = 'Kegg'
            gene_annotation = requests.get("http://rest.kegg.jp/find/genes/{}".format(gene[0]))
            element_name = str(gene_annotation.text).split(';')[0].split('\t')[1]
            
        length = gene[1][1] - gene[1][0]
        start = gene[1][0]
        element_r = Rectangle((gene[1][0], curr_y), length, height, facecolor=color, alpha=0.8, edgecolor='black', label=label)
        currentAxis.add_patch(element_r)
        rx, ry = element_r.get_xy()
        cx = rx + element_r.get_width()/2.0
        cy = ry + element_r.get_height()/2.0
        currentAxis.annotate(element_name, (cx, cy - (height/2 + 1)), color='black', fontsize=8, ha='center', va='center')
    
    curr_y += height + padding
    
sns.despine(top=True, right=True, left=True, bottom=False)
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
plt.xlim([0, max(use_reads_lengths)])
plt.ylim([0, curr_y])
plt.xlabel('Read Length')

#remove duplicates from legends
handles, labels = currentAxis.get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
    if label not in newLabels:
        newLabels.append(label)
        newHandles.append(handle)

plt.legend(newHandles, newLabels)
plt.savefig('colocalizations_{}01.svg'.format(sample_name))



