#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

import subprocess
import json
import glob
import pysam
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import json
from pathlib import Path
import csv


# Get mapped and unmapped sets

def unmapped_mapped(sam_file_path):
    samfile = pysam.AlignmentFile(sam_file_path, "r")
    
    mapped   = dict()
    unmapped = dict()
    for read in samfile:
        if (read.is_unmapped):
            unmapped[read.query_name] = read.query_length
        else:
            mapped[read.query_name] = read.infer_query_length()
            
    return unmapped, mapped

# Plot histograms

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def plot_histograms(unmapped, mapped, name):
    values_list_unmapped = unmapped.values()
    values_list_mapped   =   mapped.values()

    plt.hist(reject_outliers(np.array(list(values_list_unmapped))), alpha=0.5, label='unmapped', bins=np.logspace(1,4,50))
    plt.hist(reject_outliers(np.array(list(values_list_mapped))), alpha=0.5, label='mapped', bins=np.logspace(1,4,50))
    plt.gca().set_xscale("log")
    plt.legend()
    plt.title(name)
    plt.savefig(name + '.svg')
    plt.show()
    
# Generate csvs
def generate_csv(unmapped, mapped, name):
    unmapped_out_file = open(name + '_unmapped.csv', 'w')
    mapped_out_file   =   open(name + '_mapped.csv', 'w')
    
    # header
    unmapped_out_file.write('read_name,read_length\n')
    mapped_out_file.write('read_name,read_length\n')
    
    for read, length in unmapped.items():
        unmapped_out_file.write(read + ',')
        unmapped_out_file.write(str(length) + '\n')
        
    for read, length in mapped.items():
        mapped_out_file.write(read + ',')
        mapped_out_file.write(str(length) + '\n')

file_path = '/panfs/roc/groups/11/noyes046/jsettle/argmobrich/analysis/datasheets/tmp/C01_megares_align.sam'
unmapped, mapped = unmapped_mapped(file_path)
generate_csv(unmapped, mapped, file_path[-21:-4])
plot_histograms(unmapped, mapped, file_path[-21:-4])





