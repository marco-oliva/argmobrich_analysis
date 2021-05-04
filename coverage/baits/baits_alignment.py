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


# Baits to fasta
baits_df = pd.read_excel(io="40168_2017_361_MOESM12_ESM.xlsx")

print('Generating baits fasta')
with open('baits.fa', "w") as output_handle:
    for index, row in baits_df.iterrows():
        record = SeqRecord(Seq(row['Sequence']), id=str(row['ProbeID']), description=str(row['Target Gene Accession']))
        SeqIO.write(record, output_handle, "fasta")


# Align baits to all references
minimap2_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/minimap2'
bwa_exe = '/home/noyes046/moliva/miniconda3/envs/argmobrich/bin/bwa'
working_base_path = '/home/noyes046/moliva/projects/argmobrich/analysis/coverage/baits'
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

run_minimap2 = True
run_bwa = False

for idx, reference in enumerate(references):
    print('[{}]: {}'.format(idx, reference))
    if (run_minimap2):
        minimap2_command = '{} -a {}/{} {}/baits.fa'.format(minimap2_exe, ref_base_path, reference, working_base_path)
        with open('{wp}/baits_ato_{r}_MINIMAP2.sam'.format(wp=working_base_path, r=reference[:-6]), 'w') as out_file:
             subprocess.run(minimap2_command.split(), stdout=out_file, shell=False)
    
    if (run_bwa):
        bwa_index_command = '{be} index {wp}/{r}'.format(be=bwa_exe, wp=working_base_path, r=reference)
        bwa_align_command = '{be} mem {wp}/{r} {wp}/baits.fa'.format(r=reference, be=bwa_exe, wp=working_base_path)
        subprocess.run(bwa_index_command.split(), shell=False)
        with open('{wp}/baits_ato_{r}_BWA.sam'.format(wp=working_base_path, r=reference[:-4]), 'w') as out_file:
            subprocess.run(bwa_align_command.split(), stdout=out_file, shell=False)





