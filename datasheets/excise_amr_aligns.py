import pysam
from Bio import SeqIO

import sys


if __name__ == "__main__":

    unmapped_sequences = {}
    print("start")

    #make a sequence of 1s for each read in original fastq
    fastq_filename = sys.argv[2]
    for rec in SeqIO.parse(fastq_filename, "fastq"):
        unmapped_sequences[rec.id] = [rec.seq, [1 for x in range(0, len(rec.seq))], False]

    print("done going through original fastq")
    #get gene lengths for megares
    megares_gene_lengths = {}
    megares_reference_fasta_filename = sys.argv[4]
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    print("done getting megares gene lengths")

    #open up sam file and iterate through alignments
    sam_filename = sys.argv[1]
    samfile = pysam.AlignmentFile(sam_filename, 'r')
    for read in samfile.fetch():

        #if alignment is 1000 bps and the reference is 80% covered
        if (read.query_length - read.query_alignment_length) > 1000 and read.reference_length > 0.8*megares_gene_lengths[read.reference_name]:
            #Set to True to indicate that this sequence needs to be written to new fasta file
            unmapped_sequences[read.query_name][2] = True

            #Set all "bits" corresponding to this alignment to zero
            if read.is_reverse:
                for i in range(read.query_length - read.query_alignment_end, read.query_length - read.query_alignment_start):
                    unmapped_sequences[read.query_name][1][i] = 0
            else:
                for i in range(read.query_alignment_start, read.query_alignment_end):
                    unmapped_sequences[read.query_name][1][i] = 0

    #Now write the fasta file of unmapped regions
    headers_and_seqs = [] #stores fasta headers and sequences
    seq_str = "" #records sequence string for unmapped regions
    for item in unmapped_sequences.items():
        #Check if the sequence was marked True or not (aka if ARG aligned to it)
        if not item[1][2]:
            continue


        flag = 0 #Keeps track of if we are currently in an unmapped sequence (1) or mapped sequence (0)
        idx = 0 #Index bp in read
        start_idx = idx
        for b in item[1][1]:
            #in the middle of a mapped sequence, who cares
            if flag == 0 and b == 0:
                flag = 0

            #start of mapped sequence
            elif flag == 1 and b == 0: 
                #write out unmapped portion if it's large enough
                if len(seq_str) > 1000:
                    headers_and_seqs.append(("\n>" + item[0] + '|' + str(start_idx) + '|' + str(idx-1) + '\n', seq_str))

                #set flag and sequence string
                flag = 0
                seq_str = ""

            #start of unmapped sequence, we need to start recording bps
            elif flag == 0 and b == 1: 
                seq_str += item[1][0][idx]
                flag = 1
                start_idx = idx

            #(flag == 1 and b == 1) in the middle of unmapped sequence, record
            else: 
                seq_str += item[1][0][idx]
            idx += 1

        #At end of read, might not have written out yet
        if len(seq_str) > 1000: 
            headers_and_seqs.append(("\n>" + item[0] + '|' + str(start_idx) + '|' + str(idx-1) + '\n', seq_str))

        #Clear sequence string for next read
        seq_str = ""
            

    #Put longest regions first
    sorted_headers_and_seqs = sorted(headers_and_seqs, key=lambda header_and_seq: len(header_and_seq[1]), reverse=True)

    #Write out
    output_fasta_filename = sys.argv[3]
    with open(output_fasta_filename, 'w') as out_file:
        for h_n_s in sorted_headers_and_seqs:
            #Index 0 is header, index 1 is sequence string
            out_file.write(h_n_s[0] + h_n_s[1])
