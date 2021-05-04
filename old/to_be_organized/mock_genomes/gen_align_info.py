from Bio import SeqIO
import pysam

import csv
import sys

if __name__ == "__main__":

    #Create spreadsheet with rows for each pilot library/genome pair, with row in between pilot libs
    #columns: largest % of genome covered by one alignment; columns about that alignment: SAM flag, read header, length on reference, ins dels, matches; % of raw reads that align to genome

    lib_fastq = sys.argv[1]
    lib_prefix = sys.argv[2]
    species_names = ["bacillus_subtilis", "cryptococcus_neoformans", "enterococcus_faecalis",
                     "escherichia_coli", "lactobacillus_fermentum", "listeria_monocytogenes",
                     "pseudomonas_aeruginosa", "saccharomyces_cerevisiae", "salmonella_enterica",
                     "staphylococcus_aureus"]

    genome_alignment_sams = [pysam.AlignmentFile( species_names[0] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[1] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[2] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[3] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[4] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[5] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[6] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[7] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[8] + '.fna_' + lib_prefix + "_align.sam", 'r'),
                             pysam.AlignmentFile( species_names[9] + '.fna_' + lib_prefix + "_align.sam", 'r')]


    #First, determine number of reads in library
    num_reads = 0
    for rec in SeqIO.parse(lib_fastq, "fastq"):
        num_reads += 1
    
    species_genome_lengths = {}
    for species_name in species_names:
        for rec in SeqIO.parse(species_name + ".fna", "fasta"):
            species_genome_lengths[species_name] = len(rec)

    #Used to determine percent of reads that align to genome set
    read_dict = {}

    #Label tsv with library name, and prepare column titles
    tsv_rows = [[lib_prefix], [],
                ["Species", "Read header", "Length on reference (%)", "Length on reference (nt)", "Matches", "Ins", "Dels"]]

    #Go through alignment file for each genome and find alignment with longest length on reference
    for sam_file in genome_alignment_sams:
        species_name = str(sam_file.filename, 'utf-8').split('.')[0]
        species_genome_length = species_genome_lengths[species_name]
        
        longest_ref_length = 0
        longest_alignment = dict.fromkeys(["SAM flag", "read header", "length on reference", "percent length on reference", "ins",
                                           "dels", "matches"])
        for alignment in sam_file.fetch():
            if alignment.reference_length is None:
                continue

            if alignment.is_unmapped:
                continue

            if not alignment.query_name in read_dict:
                read_dict[alignment.query_name] = True

            #New longest alignment! Update tracked alignment
            if alignment.reference_length > longest_ref_length:
                #TODO worry about primary/supplementary?


                longest_alignment["read header"] = alignment.query_name
                longest_alignment["length on reference"] = alignment.reference_length
                longest_alignment["percent length on reference"] = alignment.reference_length / float(species_genome_length)
                longest_alignment["SAM flag"] = alignment.flag

                cigar_stats_nt = (alignment.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
                longest_alignment["matches"] = cigar_stats_nt[0]
                longest_alignment["ins"] = cigar_stats_nt[1]
                longest_alignment["dels"] = cigar_stats_nt[2]

                longest_ref_length = alignment.reference_length

        tsv_rows.append([species_name, longest_alignment["read header"], longest_alignment["percent length on reference"],
                         longest_alignment["length on reference"], longest_alignment["matches"], longest_alignment["ins"],
                         longest_alignment["dels"]])

    #Write tsv
    with open(lib_prefix + ".tsv", 'w') as out_tsv:
        tsv_writer = csv.writer(out_tsv, delimiter='\t')
        tsv_writer.writerows(tsv_rows)
        tsv_writer.writerow([])
        tsv_writer.writerow(["Percent of reads that align to mock community: ", len(read_dict)/float(num_reads)])

    print(lib_prefix + " is done")
    sys.exit(0)
