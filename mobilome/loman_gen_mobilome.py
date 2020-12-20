from Bio import SeqIO
import pysam

import csv
import matplotlib.pyplot as plt
import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    #Get aclame, iceberg, and plasmid finder lengths for coverage
    aclame_gene_lengths = {}
    iceberg_gene_lengths = {}
    pf_gene_lengths = {}
    aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
    iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
    pf_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"


    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        aclame_gene_lengths[rec.name] = len(rec.seq)

    for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        iceberg_gene_lengths[rec.name] = len(rec.seq)

    for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
        pf_gene_lengths[rec.name] = len(rec.seq)

    reads_aligned = {}
    gene_dict = {}
    for i in range(0, int(sys.argv[2])):

        if i < 9:
            aclame_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-0" + str(i+1) + "_aclame.sam", 'r')
            iceberg_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-0" + str(i+1) + "_iceberg.sam", 'r')
            pf_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-0" + str(i+1) + "_plasmidfinder.sam", 'r')
        else:
            aclame_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-" + str(i+1) + "_aclame.sam", 'r')
            iceberg_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-" + str(i+1) + "_iceberg.sam", 'r')
            pf_sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-" + str(i+1) + "_plasmidfinder.sam", 'r')

        #Iterate through every read. Accumulate number of reads while recording read length
        for read in aclame_sam_file.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            #check coverage
            if (float(read.reference_length) / aclame_gene_lengths[read.reference_name]) > 0.5:
                if(not read.reference_name in gene_dict):
                    gene_dict[read.reference_name] = ["aclame", 1]
                else:
                    gene_dict[read.reference_name][1] += 1

                if(not read.query_name in reads_aligned):
                    reads_aligned[read.query_name] = True

        #Iterate through every read. Accumulate number of reads while recording read length
        for read in iceberg_sam_file.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            #check coverage
            if (float(read.reference_length) / iceberg_gene_lengths[read.reference_name]) > 0.5:
                if( not read.reference_name in gene_dict):
                    gene_dict[read.reference_name] = ["iceberg", 1]
                else:
                    gene_dict[read.reference_name][1] += 1

                if(not read.query_name in reads_aligned):
                    reads_aligned[read.query_name] = True

        #Iterate through every read. Accumulate number of reads while recording read length
        for read in pf_sam_file.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            #check coverage
            if (float(read.reference_length) / pf_gene_lengths[read.reference_name]) > 0.5:
                if( not read.reference_name in gene_dict):
                    gene_dict[read.reference_name] = ["pf", 1]
                else:
                    gene_dict[read.reference_name][1] += 1

                if(not read.query_name in reads_aligned):
                    reads_aligned[read.query_name] = True

    #Prepare rows of tsv
    tsv_rows = [[os.path.basename(sys.argv[1])[0:6]]]
    gene_riches = [(header, gene_dict[header][0], gene_dict[header][1]) for header in gene_dict]
    tsv_rows.append(["MGE Richness:", len(gene_riches)])
    tsv_rows.append(["Number of Reads Aligned to MGE Databases:", len(reads_aligned)])

    tsv_rows.append(["MGE Database", "MGE Header", "Num Reads"])
    for gene_count_tuple in sorted(gene_riches, key=lambda gene_count_tuple: gene_count_tuple[2], reverse=True):
        tsv_rows.append([gene_count_tuple[1], gene_count_tuple[0], gene_count_tuple[2]])

    #Write tsv
    with open(os.path.basename(sys.argv[1])[0:6] + "_mobilome.tsv", 'w') as out_tsv:
        tsv_writer = csv.writer(out_tsv, delimiter='\t')
        tsv_writer.writerows(tsv_rows)

    sys.exit(0)
