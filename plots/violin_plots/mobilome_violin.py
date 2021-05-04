#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

from Bio import SeqIO
import pysam
import plotly.graph_objects as go

from math import log
import os.path
import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np


# In[2]:


dataset_sizes = [0.369211314,
                 0.402193959,
                 0.363206853,
                 0.183602223,
                 0.226892212,
                 0.414064482,
                 16.506755978,
                 153.718152889]

#Get aclame, iceberg, and plasmid finder lengths for coverage
gene_lengths = {}
aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
pf_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
    gene_lengths[rec.name] = len(rec.seq)

for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
    gene_lengths[rec.name] = len(rec.seq)

for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
    gene_lengths[rec.name] = len(rec.seq)

#Input is text file with list of sam files
sbp = '/home/noyes046/moliva/projects/argmobrich/analysis_js/datasheets/tmp'
sam_files_paths = [
    sbp + '/A01_mobilome.sam',
    sbp + '/B01_mobilome.sam',
    sbp + '/C01_mobilome.sam',
    sbp + '/MOCK_D01_mobilome.sam',
    sbp + '/MOCK_E01_mobilome.sam',
    sbp + '/MOCK_F01_mobilome.sam',
    sbp + '/../loman_tmp/ERR3152366_1_mobilome.sam',
    sbp + '/../loman_tmp/ERR3152367_1_mobilome.sam'
]

sam_files = []
for path in sam_files_paths:
    sam_files.append(pysam.AlignmentFile(path, 'r'))

absolute_abundances = {}
for sam_file in sam_files:
    sample_absolute_abundance = {}
    #Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        #check coverage
        if (float(read.reference_length) / gene_lengths[read.reference_name]) > 0.5:
            if not read.reference_name in sample_absolute_abundance:
                sample_absolute_abundance[read.reference_name] = 0
            sample_absolute_abundance[read.reference_name] += 1

    filename = str(sam_file.filename, 'utf-8')
    sam_file.close()
    absolute_abundances[filename] = sample_absolute_abundance

relative_abundances = {}
i = 0
#Convert absolute abundances to relative abundance
for filename in absolute_abundances:
    sample_relative_abundance = {}
    for gene_name in absolute_abundances[filename]:
        sample_relative_abundance[gene_name] = log(absolute_abundances[filename][gene_name] / float( gene_lengths[gene_name] * dataset_sizes[i]))

    relative_abundances[filename] = sample_relative_abundance
    i +=1


# In[11]:


df = pd.DataFrame.from_dict(relative_abundances)
df.rename(columns={df.columns[0]: "fecal 2kb (TELS)", 
                   df.columns[1]: "fecal 5kb (TELS)",
                   df.columns[2]: "fecal 8kb (TELS)",
                   df.columns[3]: "mock 2kb (TELS)",
                   df.columns[4]: "mock 5kb (TELS)",
                   df.columns[5]: "mock 8kb (TELS)",
                   df.columns[6]: "mock (GridION)",
                   df.columns[7]: "mock (PromethION)"}, inplace = True)
sns.set_style("whitegrid")
sns.set_context("paper")
ax = sns.violinplot(data=df, inner='box', palette='hls')
ax.set(xlabel='Samples', ylabel='Log Relative Abundance')
ax.set_xticklabels(ax.get_xticklabels(),rotation=-30, ha="left")
plt.gcf().subplots_adjust(bottom=0.30, right=0.85)
ax.get_figure().savefig('mobilome_violin.svg')
ax.get_figure().savefig('mobilome_violin.png', dpi=500)


# In[41]:


print('mobilome')
for c in df.columns:
    print(c, df[c].median(), sep=",")
          


# In[ ]:




