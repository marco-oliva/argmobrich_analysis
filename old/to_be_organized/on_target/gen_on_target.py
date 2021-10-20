from Bio import SeqIO
import pysam

import csv
import os.path
import sys

if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    sam_file = pysam.AlignmentFile(sys.argv[1], 'r')

    #Get megares lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    seq_dict = {}
    #write header for tsv
    tsv_rows = [["Sample header", "MEGARes header", "Coverage", "Matches (nts)", "Matches (blks)", "Ins (nts)", "Ins (blks)", \
                 "Dels(nts)", "Dels(blks)", "CIGAR length (nts)", "Length of alignment on sample",
                 "Length of alignment on MEGARes"]]

    #Iterate through every read. Record alignment info, and cigar string
    for read in sam_file.fetch():
        if read.reference_length < 0.8*megares_gene_lengths[read.reference_name]:
            continue
        if "RequiresSNPConfirmation" in read.reference_name:
            continue
 
        seq_dict[read.query_name] = True

        cigar_stats_nt = (read.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
        cigar_stats_blk = (read.get_cigar_stats())[1] #1st element is block (blk) counts of each cigar op
        tsv_row = [read.query_name, read.reference_name, float(read.reference_length) / megares_gene_lengths[read.reference_name]]
        tsv_row.extend([cigar_stats_nt[0], cigar_stats_blk[0], cigar_stats_nt[1], cigar_stats_blk[1], cigar_stats_nt[2],
                        cigar_stats_blk[2]])
        tsv_row.extend([sum(cigar_stats_nt[:3]), read.query_alignment_length, read.reference_length])

        tsv_rows.append(tsv_row)


    #Write tsv
    with open(os.path.basename(sys.argv[1]) + ".tsv", 'w') as out_tsv:
        tsv_writer = csv.writer(out_tsv, delimiter='\t')
        #first row of tsv file is raw count
        tsv_writer.writerow(["Number of on target reads: ", len(seq_dict)])
        #empty row
        tsv_writer.writerow([])
        #tsv_writer.writerows(tsv_rows)

    sys.exit(0)
