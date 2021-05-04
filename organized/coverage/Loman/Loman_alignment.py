#!/usr/bin/env python

## Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import itertools
import matplotlib.pyplot as plt
from tqdm import tqdm
import itertools
import subprocess


# Reads Files 
reads_base_path = '/home/noyes046/moliva/projects/argmobrich/loman'
reads_files = list()
    
# ERR3152366_1
for i in range(1, 5):
    reads_files.append('ERR3152366_1.part-{:02d}.fastq'.format(i))
    
# ERR3152367_1
for i in range(1, 31):
    reads_files.append('ERR3152367_1.part-{:02d}.fastq'.format(i))
    
    

# Align reads to all references
minimap2_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/minimap2'
bwa_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/bwa'
working_base_path = '/home/noyes046/moliva/projects/argmobrich/analysis/coverage/loman/sam_files'
ref_base_path = '/home/noyes046/moliva/projects/argmobrich/analysis/coverage/references/ZymoBIOMICS.STD.refseq.v2/Genomes'
references = [
                'Bacillus_subtilis_complete_genome.fasta',
                'Cryptococcus_neoformans_draft_genome.fasta',
                'Enterococcus_faecalis_complete_genome.fasta',
                'Escherichia_coli_complete_genome.fasta',
                'Lactobacillus_fermentum_complete_genome.fasta',
                'Listeria_monocytogenes_complete_genome.fasta',
                'Pseudomonas_aeruginosa_complete_genome.fasta',
                'Saccharomyces_cerevisiae_draft_genome.fasta',
                'Salmonella_enterica_complete_genome.fasta',
                'Staphylococcus_aureus_complete_genome.fasta'
             ]

for idx, reference in enumerate(references):
    sam_files = list()
    for read_file in reads_files:
        print('[{}]: {} to {}'.format(idx, read_file, reference))
        minimap2_command = '{} -ax map-ont -t 24 {}/{} {}/{}'.format(minimap2_exe, ref_base_path, reference, reads_base_path, read_file)
        out_sam_file = '{}/{}_ato_{}_MINIMAP2.sam'.format(working_base_path, read_file[:-6], reference[:-6])
        with open(out_sam_file, 'w') as out_file:
             subprocess.run(minimap2_command.split(), stdout=out_file, shell=False)        
        
