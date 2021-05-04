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
import pandas as pd
import math


# In[2]:


#Get aclame, iceberg, and plasmid finder lengths for coverage
aclame_gene_lengths = {}
#iceberg_gene_lengths = {}
pf_gene_lengths = {}
aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
#iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
pf_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
    aclame_gene_lengths[rec.name] = len(rec.seq)

#for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
    #iceberg_gene_lengths[rec.name] = len(rec.seq)

for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
    pf_gene_lengths[rec.name] = len(rec.seq)

#Get informative names for aclame from excel
plasmid_ref_df = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=0, header=1, usecols="A,F", index_col=0)
plasmid_ref_df_2 = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=1, header=0, usecols="A,F", index_col=0)
#print(plasmid_ref_df)
plasmid_ref_df = plasmid_ref_df.append(plasmid_ref_df_2)
#print(plasmid_ref_df)
prophage_ref_df = pd.read_excel("aclame_genes_prophages_0.4.xls", header=1, usecols="A,F", index_col=0)


# In[3]:


# Samples
base_path = '/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/tmp/heatmaps/mge'
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
    sam_files.append(pysam.AlignmentFile('{}/{}_mobilome.sam'.format(base_path, sample_name), 'r'))


# In[4]:


# Read annotations
annotations_df = pd.read_excel("TELS_MGE_ACCESSIONS.xlsx")

mge_genes_names = dict()
mge_genes_types = dict()
mge_accessions = set()

for index, row in annotations_df.iterrows():
    mge_genes_names[str(row['ACCESSIONS'])] = str(row['MGE NAME'])
    mge_genes_types[str(row['MGE NAME'])] = str(row['MGE TYPE'])
    mge_accessions.add(str(row['ACCESSIONS']))


# In[5]:


samples = {}

for sam_file in sam_files:
    filename = str(sam_file.filename, 'utf-8')
    sample_name = os.path.splitext(filename)[0].split('/')[-1].split('_mobilome')[0]
    
    #Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
            
        gene_db_name = ''
        translated_name = ''
            
        #check coverage
        if (read.reference_name in aclame_gene_lengths):
            if (float(read.reference_length) / aclame_gene_lengths[read.reference_name]) > 0.5:
                gene_db_name = read.reference_name
        elif (read.reference_name in pf_gene_lengths):
            if (float(read.reference_length) / pf_gene_lengths[read.reference_name]) > 0.5:
                gene_db_name = read.reference_name
        
        if (gene_db_name != ''):
            if gene_db_name in mge_accessions:
                translated_name = mge_genes_names[gene_db_name]
            elif gene_db_name in prophage_ref_df.index: 
                translated_name = str(prophage_ref_df.loc[gene_db_name].iloc[0])
            elif gene_db_name in plasmid_ref_df.index:
                translated_name = str(plasmid_ref_df.loc[gene_db_name].iloc[0])
            # Find a better way
            elif ('Inc' in gene_db_name) or ('Rep' in gene_db_name):
                translated_name = gene_db_name
                
                
        if (translated_name != ''):
            if translated_name not in samples:
                samples[translated_name] = set()
            samples[translated_name].add(samples_names[sample_name])
                
    sam_file.close()


# In[6]:


bovine_sample_names = [samples_names['A01'], samples_names['B01'], samples_names['C01']]
mock_sample_names = [samples_names['MOCK_D01'], samples_names['MOCK_E01'], samples_names['MOCK_F01']]
loman_sample_names = [samples_names['ERR3152366'], samples_names['ERR3152367']]


# In[17]:


genes_tuples = []
for gene_name in samples:
    if (gene_name in mge_genes_types):
        gene_type = mge_genes_types[gene_name]
        if (gene_type != 'Red'):
            genes_tuples.append((gene_name, gene_type))

genes_tuples.sort(key=lambda tuple: tuple[1])


# In[ ]:


samples


# In[49]:


# Output tsv
with open('mges.csv', 'w') as csvfile:
    header = ['Gene Name', 'Gene Type',
              'fecal 2kb (TELS)', 'fecal 5kb (TELS)', 'fecal 8kb (TELS)',
              'mock 2kb (TELS)', 'mock 5kb (TELS)', 'mock 8kb (TELS)',
              'mock (GridION)', 'mock (PromethION)']
    csv_writer = csv.writer(csvfile, delimiter=',')
    csv_writer.writerow(header)
    for gene_name, samples_set in samples.items():
        if (gene_name in mge_genes_types):
            gene_type = mge_genes_types[gene_name]
            if (gene_type != 'Red'):
                line = [gene_name, gene_type]
                for e in header[2:]:
                    if (e in samples_set):
                        line.append('1')
                    else:
                        line.append('0')
                csv_writer.writerow(line)


# In[19]:


mgetick_text = []
mgetick_vals = []
label_matrix = []
current_label = genes_tuples[0][1]
label_val = float(0.0)
last_count = float(0.0)
count = float(0.0)
tickvals = []
ticklabels = []

for gene_name, gene_label in genes_tuples:
    if current_label != gene_label:
        tickvals.append((last_count + count - 1.0) / 2.0)
        ticklabels.append(current_label)

        last_count = count
        label_val += 1.0
        current_label = gene_label

    label_matrix.append([label_val])

    mgetick_text.append(gene_name)
    mgetick_vals.append(count)

    count += 1.0

bovine_matrix = []
mock_matrix = []
loman_matrix = []

represented_hierarchy = dict()
for gene_ref, gene_label in genes_tuples:
    if (gene_label not in represented_hierarchy):
        represented_hierarchy[gene_label] = 0
    represented_hierarchy[gene_label] += 1

for gene_ref, gene_label in genes_tuples:
    
    bovine_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
    for i in range(0, len(bovine_matrix_row)):
        if bovine_sample_names[i] in samples[gene_ref]:
            bovine_matrix_row[i] = 1.0
    bovine_matrix.append(bovine_matrix_row)

    mock_matrix_row = [0 for x in range(0, len(mock_sample_names))]
    for i in range(0, len(mock_matrix_row)):
        if mock_sample_names[i] in samples[gene_ref]:
            mock_matrix_row[i] = 1.0
    mock_matrix.append(mock_matrix_row)
    
    loman_matrix_row = [0 for x in range(0, len(loman_sample_names))]
    for i in range(0, len(loman_matrix_row)):
        if loman_sample_names[i] in samples[gene_ref]:
            loman_matrix_row[i] = 1.0
    loman_matrix.append(loman_matrix_row)
    

label_matrix = []
i = 1
for gene_type, n in represented_hierarchy.items():
    val = i / float(len(represented_hierarchy) + 1.0)
    for _ in range(n):
        label_matrix.append([val])
    i+=1


# In[40]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

custom_color_list = ["lightpink", "palegreen", "peachpuff", "powderblue", "plum", "lightsalmon", "burlywood", "tomato", "chocolate", "gold", "salmon", "orchid"]
n_bins = len(set([l[0] for l in label_matrix]))
custom_color = LinearSegmentedColormap.from_list('custom', custom_color_list, N=n_bins)


color_max = np.max(label_matrix)
color_min = 0

custom_color_list_on_off = ["grey", "lightsalmon"]
n_bins = 2
custom_color_on_off = LinearSegmentedColormap.from_list('custom', custom_color_list_on_off, N=n_bins)

fig, axs = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1, 3, 3, 2]}, figsize=(18, 32))
axs[0].imshow(label_matrix, aspect='auto', vmin=color_min, vmax=color_max, cmap=custom_color)
axs[1].imshow(bovine_matrix, aspect='auto', cmap=custom_color_on_off)
axs[2].imshow(mock_matrix, aspect='auto', cmap=custom_color_on_off)
axs[3].imshow(loman_matrix, aspect='auto', cmap=custom_color_on_off)

# Use the pyplot interface to change just one subplot...
plt.sca(axs[0])
plt.yticks(tickvals, ticklabels, color='black', fontsize=15)
plt.xticks([], [], color='black')

plt.sca(axs[1])
plt.xticks([0,1,2], ["fecal 2kb (TELS)", "fecal 5kb (TELS)", "fecal 8kb (TELS)"], color='black', rotation=-45, ha="left", fontsize=15)
plt.yticks([], [], color='black')

plt.sca(axs[2])
plt.xticks([0,1,2], ["mock 2kb (TELS)", "mock 5kb (TELS)", "mock 8kb (TELS)"], color='black', rotation=-45, ha="left", fontsize=15)
plt.yticks([], [], color='black')

plt.sca(axs[3])
plt.xticks([0,1], ["mock (GridION)", "mock (PromethION)"], color='black', rotation=-45, ha="left", fontsize=15)
#plt.yticks(mgetick_vals, mgetick_text, color='black')
plt.tick_params(axis='y', labelleft=False, left=False, labelright=False, labelsize=5)
plt.gcf().subplots_adjust(bottom=0.15, left=0.20)
plt.savefig('heatmap_mge.svg')


# In[ ]:





# In[ ]:





# In[70]:


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


# In[90]:


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


fig, axs = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1, 3, 3, 2]}, figsize=(20, 15))
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
plt.tick_params(axis='y', labelleft=False, labelright=True, labelsize=15)


# In[ ]:




