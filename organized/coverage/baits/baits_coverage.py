## Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

import json
import glob
import pysam
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO


# Get Reference Genomes files info
reference_info_json = open('/home/noyes046/moliva/projects/argmobrich/analysis/coverage/references_info.json')
reference_info = json.load(reference_info_json)

print(reference_info)

# Alignment files
aligners_used = ['MINIMAP2']
alignments_base_path = '/home/noyes046/moliva/projects/argmobrich/analysis/coverage/baits'
sam_files = dict()
for organism in reference_info['organisms_list']:
    sam_files[organism] = dict()
    for als in aligners_used:
        sam_files[organism][als] = glob.glob(alignments_base_path + '/' + 'baits_ato_' +  organism + '_' + als + '.sam')

# For each alignment produce counts, if multi chromosomes then multi chromosome count
organisms_coverage = dict()
for organism in reference_info['organisms_list']:
    # Read reference in a dict<id, sequence>
    reference_file = reference_info['base_path'] + '/' + organism + reference_info['extension']
    reference = dict()
    for record in SeqIO.parse(reference_file, "fasta"):
        reference[record.id] = len(record.seq)
    print('reference: ', reference)

    # Read Sam files and update coverage
    organisms_coverage[organism] = dict()
    for als in aligners_used: 
        # Initialize tmp coverage array
        tmp_array_dict = dict()
        for record_id, record_seq_len in reference.items():
            tmp_array_dict[record_id] = np.zeros(record_seq_len) 

        # Compute coverage
        for samfile_path in sam_files[organism][als]:
            print('organisms_coverage[{}][{}]: processing {}'.format(organism, als, samfile_path))
            samfile = pysam.AlignmentFile(samfile_path)
            for read in samfile:
                if ((not read.is_unmapped) and (not read.is_secondary)):
                    # Update Hits
                    for i in range(read.reference_start, read.reference_start + read.reference_length):
                        r_id = samfile.get_reference_name(read.reference_id)
                        tmp_array_dict[r_id][i] = tmp_array_dict[r_id][i] + 1

        for ref_id in tmp_array_dict.keys():
            tmp_array_dict[ref_id] = tmp_array_dict[ref_id].tolist()

        organisms_coverage[organism][als] = tmp_array_dict

# Dump coverage arrays to file
with open('Baits_hits.json', 'w') as out_file:
    json.dump(organisms_coverage, out_file)


