#!/usr/bin/env python
# coding: utf-8

## Imports
import sys

anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

from io import StringIO
import time
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser


# Functions
def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

def download_database(organism_list, out_file_path):
    with open(out_file_path, 'w') as output_handle:
        for organism in organism_list:
            genes_list_res = REST.kegg_list(organism).read()
            genes_list = list(to_df(genes_list_res)[0])
            for gene in tqdm(genes_list):
                fasta_seq = ""
                for wait_time in [0.00001, 16, 32, 64, 128, 256, 512]:
                    try:
                        time.sleep(wait_time)
                        fasta_seq = REST.kegg_get(gene, "ntseq").read()
                    except:
                        print("Error in downloading gene:\t{}\twating:{}".format(gene, wait_time))
                        continue
                    break
                if len(fasta_seq) == 0:
                    print("Error persists, exiting")
                    exit()
                with StringIO(fasta_seq) as fasta_io:
                    records = list(SeqIO.parse(fasta_io, "fasta"))
                    record = SeqRecord(Seq(str(records[0].seq).upper()), id=records[0].id, description='')
                    SeqIO.write(record, output_handle, "fasta")


# Organisms List
mock_comunity_organisms_list = ['bsu', 'eco', 'lfe', 'lmo', 'pae',
                                'sau', 'cne', 'efa', 'sce', 'sty']
prokaryotes_default_organisms_list = ['dme', 'ath', 'sce', 'pfa', 'eco',
                                      'sty', 'hin', 'pae', 'nme', 'hpy',
                                      'rpr', 'mlo', 'bsu', 'sau', 'lla',
                                      'spn', 'cac', 'mge', 'mtu', 'ctr',
                                      'bbu', 'syn', 'aae', 'mja', 'afu',
                                      'pho', 'ape']

# Download KEGG's prokaryotes default database
download_database(prokaryotes_default_organisms_list, 'kegg_prokaryotes.fasta')

# Downloads mock comunity
download_database(prokaryotes_default_organisms_list, 'kegg_mock.fasta')




