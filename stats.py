#!/usr/bin/env python
# coding: utf-8

# Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)


from Bio import SeqIO
import numpy as np


databases = {
    'ACLAME'   : '/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta',
    'ICEBERG'  : '/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta',
    'PLASMIDS' : '/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa'
}

print('db_name, min_length, max_length, median_length, number_of_sequences')
for db_name, db_path in databases.items():
    lengths_list = []
    with open(db_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lengths_list.append(len(record.seq))
    lengths_array = np.array(lengths_list)
    print('{},{},{},{},{}'.format(db_name, np.min(lengths_array), np.max(lengths_array), np.median(lengths_array), lengths_array.size))



