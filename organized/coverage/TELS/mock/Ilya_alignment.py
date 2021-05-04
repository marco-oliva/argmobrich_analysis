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
import subprocess
import itertools


# Reads Files (D,E,F)
reads_base_path = '/home/noyes046/moliva/projects/argmobrich/ccs_fastqs'
reads_files = [
                'sequel-demultiplex.MOCK_D01.ccs.fastq',
                'sequel-demultiplex.MOCK_E01.ccs.fastq',
                'sequel-demultiplex.MOCK_F01.ccs.fastq'
              ]

# Align reads to all references
minimap2_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/minimap2'
bwa_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/bwa'
working_base_path = '/home/noyes046/moliva/projects/argmobrich/analysis/coverage/ilya/mock'
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
    for read_file in reads_files:
        print('[{}]: {} to {}'.format(idx, read_file, reference))
        minimap2_command = '{} -ax map-pb {}/{} {}/{}'.format(minimap2_exe, ref_base_path, reference, reads_base_path, read_file)
        with open('{}/{}_ato_{}_MINIMAP2.sam'.format(working_base_path, read_file[:-6], reference[:-6]), 'w') as out_file:
             subprocess.run(minimap2_command.split(), stdout=out_file, shell=False)
