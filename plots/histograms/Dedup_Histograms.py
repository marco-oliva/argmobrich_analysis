#!/usr/bin/env python
# coding: utf-8


## Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

import json
import numpy as np
import itertools
import matplotlib.pyplot as plt
import glob
import seaborn as sns


def gen_histogram_not_log(tsv_file_path, name):
    tsv_file = open(tsv_file_path, 'r')
    dups_list = []
    for line in tsv_file:
        dups_list.append(len(line.split(',')))
    data = np.array(dups_list)
    
    plt.hist(data, 100)
    plt.yscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Set Size (Number of Reads)')
    plt.title(name)
    plt.show() 

def gen_histogram(tsv_file_path, name):
    tsv_file = open(tsv_file_path, 'r')
    dups_list = []
    for line in tsv_file:
        dups_list.append(len(line.split(',')))
        if len(line.split(',')) >= 900:
            print(len(line.split(',')), line[0:100])

    # log scaled bins
    data = np.array(dups_list)
    bins = np.logspace(0, 4, 50)
    widths = (bins[1:] - bins[:-1])
    hist = np.histogram(data, bins=bins)
    hist_norm = hist[0]/1 #widths
    sns.set_context("paper")
    sns.set_palette("hls")
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
    
    plt.bar(bins[:-1], hist_norm, widths, color=[(0.4118,0.4824,0.7725)])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Set Size (Number of Reads)')
    plt.grid(axis='y', alpha=0.75)
    plt.title(name)
    plt.savefig(name + '.png', dpi=500)
    plt.savefig(name + '.svg')
    plt.show() 


for tsv in glob.glob('./*MOCK*.tsv'):
    name = ''
    if 'D' in tsv:
        name = 'Mock 2kb'
    elif 'E' in tsv:
        name = 'Mock 5kb'
    elif 'F' in tsv:
        name = 'Mock 8kb'
    else:
        name = 'Somthing went wrong'
    gen_histogram(tsv, name)
    
for tsv in glob.glob('./*1896*.tsv'):
    name = ''
    if 'A' in tsv:
        name = 'Fecal 2kb'
    elif 'B' in tsv:
        name = 'Fecal 5kb'
    elif 'C' in tsv:
        name = 'Fecal 8kb'
    else:
        name = 'Somthing went wrong'
    gen_histogram(tsv, name)



