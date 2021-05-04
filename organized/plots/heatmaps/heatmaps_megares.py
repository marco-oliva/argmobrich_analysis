#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt
import itertools
import csv
import os


# In[2]:


#Get megares lengths for coverage
megares_gene_lengths = {}
megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
    megares_gene_lengths[rec.name] = len(rec.seq)


# In[3]:


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


# In[4]:


# Samples
base_path = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/heatmaps/megares'
samples_names = { 
    'A01' : 'fecal 2kb (TELS)', 
    'B01' : 'fecal 5kb (TELS)', 
    'C01' : 'fecal 8kb (TELS)', 
    'MOCK_D01' : 'mock 2kb (TELS)',
    'MOCK_E01' : 'mock 5kb (TELS)',
    'MOCK_F01' : 'mock 8kb (TELS)',
    'ERR3152366' : 'mock (GridION)',
    'ERR3152367' : 'mock (PromethION)'
}

sam_files = []
for sample_name, _ in samples_names.items():
    sam_files.append(pysam.AlignmentFile('{}/{}_ato_megares.sam'.format(base_path, sample_name), 'r'))


# In[5]:


for sam_file in sam_files:
    #Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file:
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if "RequiresSNPConfirmation" in read.reference_name:
            continue
        

        #check coverage
        if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:
            cl = megares_ontology[read.reference_name]["class"]
            mech = megares_ontology[read.reference_name]["mech"]
            group = megares_ontology[read.reference_name]["group"]

            filename = str(sam_file.filename, 'utf-8')
            sample_name = os.path.splitext(filename)[0].split('/')[-1].split('_ato_')[0]

            class_samples[cl].add(samples_names[sample_name])
            mech_samples[mech].add(samples_names[sample_name])
            group_samples[group].add(samples_names[sample_name])

    sam_file.close()


# In[6]:


represented_hierarchy = {}
for cl in class_samples:
    if len(class_samples[cl]) > 0:
        represented_hierarchy[cl] = {}
        for mech in hierarchy_dict[cl]:
            if len(mech_samples[mech]) > 0:
                represented_hierarchy[cl][mech] = []


# In[7]:


bovine_sample_names = [samples_names['A01'], samples_names['B01'], samples_names['C01']]
mock_sample_names = [samples_names['MOCK_D01'], samples_names['MOCK_E01'], samples_names['MOCK_F01']]
loman_sample_names = [samples_names['ERR3152366'], samples_names['ERR3152367']]

#Categories for y axis are mechanisms 
mechtick_vals = []
mechtick_text = []
classtick_vals = []
classtick_text = []
left = 0
right = 0
mech_names = []
for cl in represented_hierarchy:
    left = len(mech_names)
    right = left + len(represented_hierarchy[cl])
    classtick_pos = (right+left-1)/2.0
    classtick_vals.append(classtick_pos)
    classtick_text.append(cl)
    i = left+1
    for mech in represented_hierarchy[cl]:
        current_tick_pos = (left + i - 1)/2.0
        mechtick_vals.append(current_tick_pos)
        mechtick_text.append(mech)

        mech_names.append(mech)
        i += 1
        left += 1

num_mechs = len(mech_names)

#Just want 2 colors for presence/absence
color_scheme = [[0.0, 'rgb(198, 198, 198)'], 
                [0.5, 'rgb(198, 198, 198)'],
                [0.5, 'rgb(194, 59, 34)'],
                [1.0, 'rgb(194, 59, 34)']]

label_heatmap_matrix = []
i = 1
for cl in represented_hierarchy:
    val = i / float(len(mech_names) + 1.0)
    for mech in represented_hierarchy[cl]:
        label_heatmap_matrix.append([val])
    i+=1

j = 1
bovine_presence_matrix = []
for cl in represented_hierarchy:
    val = j / float(len(mech_names)+1.0)
    for mech in represented_hierarchy[cl]:
        bovine_presence_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
        for i in range(0, len(bovine_presence_matrix_row)):
            if bovine_sample_names[i] in mech_samples[mech]:
                bovine_presence_matrix_row[i] = val
        bovine_presence_matrix.append(bovine_presence_matrix_row)
    j+=1

j = 1
mock_presence_matrix = []
for cl in represented_hierarchy:
    val = j / float(len(mech_names)+1.0)
    for mech in represented_hierarchy[cl]:
        mock_presence_matrix_row = [0 for x in range(0, len(mock_sample_names))]
        for i in range(0, len(mock_presence_matrix_row)):
            if mock_sample_names[i] in mech_samples[mech]:
                mock_presence_matrix_row[i] = val
        mock_presence_matrix.append(mock_presence_matrix_row)
    j+=1

j = 1
loman_presence_matrix = []
for cl in represented_hierarchy:
    val = j / float(len(mech_names)+1.0)
    for mech in represented_hierarchy[cl]:
        loman_presence_matrix_row = [0 for x in range(0, len(loman_sample_names))]
        for i in range(0, len(loman_presence_matrix_row)):
            if loman_sample_names[i] in mech_samples[mech]:
                loman_presence_matrix_row[i] = val
        loman_presence_matrix.append(loman_presence_matrix_row)
    j+=1


# In[20]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

custom_color_list = ["grey", "lightpink", "palegreen", "peachpuff", "powderblue", "plum", "lightsalmon", "burlywood", "tomato", "chocolate", "gold", "salmon", "orchid"]
n_bins = len(mech)
custom_color = LinearSegmentedColormap.from_list('custom', custom_color_list, N=n_bins)

color_max = np.max(label_heatmap_matrix)
color_min = 0


fig, axs = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1, 3, 3, 2]}, figsize=(32, 18))
axs[0].imshow(label_heatmap_matrix, aspect='auto', vmin=color_min, vmax=color_max,  cmap=custom_color)
axs[1].imshow(bovine_presence_matrix, aspect='auto', vmin=color_min, vmax=color_max, cmap=custom_color)
axs[2].imshow(mock_presence_matrix, aspect='auto', vmin=color_min, vmax=color_max, cmap=custom_color)
axs[3].imshow(loman_presence_matrix, aspect='auto', vmin=color_min, vmax=color_max, cmap=custom_color)

# Use the pyplot interface to change just one subplot...
plt.sca(axs[0])
plt.yticks(classtick_vals, classtick_text, color='black', fontsize=20)
plt.xticks([], [], color='black')

plt.sca(axs[1])
plt.xticks([0,1,2], ["fecal 2kb (TELS)", "fecal 5kb (TELS)", "fecal 8kb (TELS)"], color='black', rotation=-45, ha="left", fontsize=20)
plt.yticks([], [], color='black')

plt.sca(axs[2])
plt.xticks([0,1,2], ["mock 2kb (TELS)", "mock 5kb (TELS)", "mock 8kb (TELS)"], color='black', rotation=-45, ha="left", fontsize=20)
plt.yticks([], [], color='black')

plt.sca(axs[3])
plt.xticks([0,1], ["mock (GridION)", "mock (PromethION)"], color='black', rotation=-45, ha="left", fontsize=20)
plt.yticks(mechtick_vals, mechtick_text, color='black')
plt.tick_params(axis='y', labelleft=False, labelright=True, left=False, right=True, labelsize=15)
plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.75)
plt.savefig('heatmap_megares.svg')


# In[ ]:




