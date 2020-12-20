import pysam
from Bio import SeqIO

import argparse
import csv
import os
import sys

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Generate datasheet for ARGMOBRich")
    parser.add_argument("files", metavar='F', type=str, nargs=3)
    parser.add_argument("--no-mge", dest="no_mge", action="store_true")

    args = parser.parse_args()
    do_mge = not args.no_mge


    #Get MEGARes ontology
    ontologies = {}
    with open("/home/noyes046/jsettle/databases/megares_v1.01/annotations/megares_annotations_v1.01.csv", 'r') as annotations:
        reader = csv.reader(annotations)
        for row in reader:
            if row[0] == "header":
                continue

            ontologies[row[0]] = [row[1], row[2], row[3]]

    if do_mge:

        #initialize dictionaries and sam files for each of the MGE databases
        aclame_data = {}
        iceberg_data = {}
        pf_data = {}
        #TODO virulence finder?

        aclame_sam_filename = sys.argv[2] + "_aclame.sam"
        aclame_samfile = pysam.AlignmentFile(aclame_sam_filename, 'r')
        iceberg_sam_filename = sys.argv[2] + "_iceberg.sam" 
        iceberg_samfile = pysam.AlignmentFile(iceberg_sam_filename, 'r')
        pf_sam_filename = sys.argv[2] + "_pf.sam" 
        pf_samfile = pysam.AlignmentFile(pf_sam_filename, 'r')

        #Get gene lengths for each MGE database. Used to calculate coverage
        aclame_gene_lengths = {}
        aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
        for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
            aclame_gene_lengths[rec.name] = len(rec.seq)

        iceberg_gene_lengths = {}
        iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
        for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
            iceberg_gene_lengths[rec.name] = len(rec.seq)

        pf_gene_lengths = {}
        pf_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"
        for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
            pf_gene_lengths[rec.name] = len(rec.seq)

        #Go through each alignment to aclame
        for aclame_read in aclame_samfile.fetch():
            #Check 1) if the alignment is actually mapped
            #      2) if the alignment covers 50% of the gene
            if (not aclame_read.is_unmapped) and aclame_read.reference_length > 0.5*aclame_gene_lengths[aclame_read.reference_name]:

                #When writing the unmapped fasta files, I augment the query name with start and stop of unmapped region. 
                #This split gets that data
                split_query_name = aclame_read.query_name.split('|')
                base_query_name = split_query_name[0]
                unmapped_region_start = split_query_name[1]
                unmapped_region_stop  = split_query_name[2]

                #If this query has already been found on a different alignment, append info to the list of alignments
                if base_query_name in aclame_data:
                    #reference_end points one past the last aligned base
                    aclame_data[base_query_name].append([aclame_read.query_name, aclame_read.reference_name, "aclame",
                                                         aclame_gene_lengths[aclame_read.reference_name],
                                                         aclame_read.reference_start, aclame_read.reference_end - 1,
                                                         aclame_read.query_alignment_start, aclame_read.query_alignment_end - 1,
                                                         unmapped_region_start, unmapped_region_stop, aclame_read.cigarstring])
                #If this query has not been found, add a new list of alignments to the dictionary
                else:
                    aclame_data[base_query_name] = [[aclame_read.query_name, aclame_read.reference_name, "aclame",
                                                     aclame_gene_lengths[aclame_read.reference_name], aclame_read.reference_start,
                                                     aclame_read.reference_end - 1, aclame_read.query_alignment_start, 
                                                     aclame_read.query_alignment_end - 1, unmapped_region_start,
                                                     unmapped_region_stop, aclame_read.cigarstring]]
 
        #Go through each alignment in iceberg
        for iceberg_read in iceberg_samfile.fetch():
            #Check 1) if the alignment is actually mapped
            #      2) if the alignment covers 50% of the gene
            if (not iceberg_read.is_unmapped) and iceberg_read.reference_length > 0.5*iceberg_gene_lengths[iceberg_read.reference_name]:

                #When writing the unmapped fasta files, I augment the query name with start and stop of unmapped region. 
                #This split gets that data
                split_query_name = iceberg_read.query_name.split('|')
                base_query_name = split_query_name[0]
                unmapped_region_start = split_query_name[1]
                unmapped_region_stop  = split_query_name[2]

                #If this query has already been found on a different alignment, append info to the list of alignments
                if base_query_name in iceberg_data:
                    iceberg_data[base_query_name].append([iceberg_read.query_name, iceberg_read.reference_name, "iceberg",
                                                          iceberg_gene_lengths[iceberg_read.reference_name],
                                                          iceberg_read.reference_start, iceberg_read.reference_end - 1,
                                                          iceberg_read.query_alignment_start, iceberg_read.query_alignment_end - 1,
                                                          unmapped_region_start, unmapped_region_stop, iceberg_read.cigarstring])
                #If this query has not been found, add a new list of alignments to the dictionary
                else:
                    iceberg_data[base_query_name] = [[iceberg_read.query_name, iceberg_read.reference_name, "iceberg",
                                                      iceberg_gene_lengths[iceberg_read.reference_name],
                                                      iceberg_read.reference_start, iceberg_read.reference_end - 1,
                                                      iceberg_read.query_alignment_start, iceberg_read.query_alignment_end - 1,
                                                      unmapped_region_start, unmapped_region_stop, iceberg_read.cigarstring]]

        #Go through each alignment in plasmid finder
        for pf_read in pf_samfile.fetch():
            #Check 1) if the alignment is actually mapped
            #      2) if the alignment covers 50% of the gene
            if (not pf_read.is_unmapped) and pf_read.reference_length > 0.5*pf_gene_lengths[pf_read.reference_name]:

                #When writing the unmapped fasta files, I augment the query name with start and stop of unmapped region. 
                #This split gets that data
                split_query_name = pf_read.query_name.split('|')
                base_query_name = split_query_name[0]
                unmapped_region_start = split_query_name[1]
                unmapped_region_stop  = split_query_name[2]

                #If this query has already been found on a different alignment, append info to the list of alignments
                if base_query_name in pf_data:
                    pf_data[base_query_name].append([pf_read.query_name, pf_read.reference_name, "plasmid_finder",
                                                     pf_gene_lengths[pf_read.reference_name], pf_read.reference_start,
                                                     pf_read.reference_end - 1, pf_read.query_alignment_start,
                                                     pf_read.query_alignment_end -1, unmapped_region_start, unmapped_region_stop,
                                                     pf_read.cigarstring])
                #If this query has not been found, add a new list of alignments to the dictionary
                else:
                    pf_data[base_query_name] = [[pf_read.query_name, pf_read.reference_name, "plasmid_finder",
                                                 pf_gene_lengths[pf_read.reference_name], pf_read.reference_start,
                                                 pf_read.reference_end - 1, pf_read.query_alignment_start,
                                                 pf_read.query_alignment_end -1, unmapped_region_start, unmapped_region_stop,
                                                 pf_read.cigarstring]]


    #Get gene lengths for coverage of MEGARes genes
    megares_gene_lengths = {}
    megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    #Set up MEGARes sam file
    megares_sam_filename = sys.argv[1]
    megares_samfile = pysam.AlignmentFile(megares_sam_filename, 'r')

    #Set output file name depending on whether MGEs are included or not
    if do_mge:
        out_filename = sys.argv[3] + '/' + os.path.basename(sys.argv[2]) + "_argmobrich_colocalization.tsv"
    else:
        out_filename = sys.argv[3] + '/' + os.path.basename(sys.argv[2]) + "_argmobrich_amr.tsv"


    #Write out colocalizations to file
    with open(out_filename, 'w') as csv_out:
        csv_wr = csv.writer(csv_out, delimiter='\t', quoting=csv.QUOTE_ALL)
        if do_mge:
            #Column headers first
            headers = ["SAMPLE_TYPE", "READ_ID", "READ_LENGTH", "AMR_GN_HEADER", "AMR_CLASS", "AMR_MECH", "AMR_GROUP",
                       "AMR_LENGTH", "AMR_START", "AMR_STOP", "AMR_CIGAR", "MGE_GN_HEADER", "MGE_TYPE", "MGE_LENGTH", "MGE_START",
                       "MGE_STOP", "MGE_QUERY_START", "MGE_QUERY_STOP", "UNMAPPED_REGION_START", "UNMAPPED_REGION_STOP",
                       "MGE_OVERLAP", "MGE_CIGAR"]
            csv_wr.writerow(headers)

            #Go through all MEGARes alignments
            for megares_read in megares_samfile.fetch():
                #Skip ARGs that require SNP confirmation
                if "RequiresSNPConfirmation" in megares_read.reference_name:
                    continue
                #Only care about primary alignments
                if megares_read.is_secondary or megares_read.is_supplementary:
                    continue
                #Check 1) alignment is mapped
                #      2) alignment on query is at least 1000 bps
                #      3) alignemnt covers 80% of ARG
                if (not megares_read.is_unmapped) and (megares_read.query_length - megares_read.query_alignment_length) > 1000 \
                                                  and megares_read.reference_length > 0.8*megares_gene_lengths[megares_read.reference_name]:

                    #Used to track whether or not alignments overlaps on the read. 
                    alignment_intervals = []
                    megares_query_length = megares_read.query_length

                    #Can't remember why I had to do this, tbh
                    if( megares_query_length == 0 ):
                        megares_query_length = megares_read.infer_query_length()

                    #Code duplication below :/
                    #Also really obtuse indices for the data list... Here's what they correspond to:
                    #   0 query name
                    #   1 reference name
                    #   2 MGE database name
                    #   3 MGE length
                    #   4 start of alignment on reference
                    #   5 end of alignment on reference
                    #   6 start of alignment on query
                    #   7 end of alignment on query
                    #   8 start of unmapped region on read (unmapped = region of read not aligned to MEGARes gene)
                    #   9 end of unmapped region on read (unmapped = region of read not aligned to MEGARes gene)
                    #   10 cigar

                    #If this read aligned to a MGE in aclame
                    if megares_read.query_name in aclame_data:

                        #For all alignments this read had to aclame
                        for data in aclame_data[megares_read.query_name]:

                            #if alignment intervals are empty, nothing to check regarding overlaps
                            if not alignment_intervals:
                                #reference_end points one past the last aligned base
                                csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                 ontologies[megares_read.reference_name][0],
                                                 ontologies[megares_read.reference_name][1],
                                                 ontologies[megares_read.reference_name][2],
                                                 megares_gene_lengths[megares_read.reference_name],
                                                 megares_read.reference_start, megares_read.reference_end - 1,
                                                 megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5], data[6],
                                                 data[7], data[8], data[9], "No", data[10]])

                                #Add the query start and end pair to alignment intervals
                                alignment_intervals.append([data[6], data[7]])
                            else:

                                #Determine if this alignment overlaps or not
                                overlap = False
                                for interval in alignment_intervals:
                                    #if alignment start/end position is in one of the other alignment intervals, then it overlaps
                                    #Recall:
                                    #   6 start of alignment on query
                                    #   7 end of alignment on query
                                    if (data[6] >= interval[0] and data[6] <= interval[1]) or \
                                       (data[7] >= interval[0] and data[7] <= interval[1]):
                                        overlap = True

                                        #Expand interval if necessary
                                        if data[6] < interval[0]:
                                            interval[0] = data[6]
                                        if data[7] < interval[0]:
                                            interval[1] = data[7]
                                        break
                                
                                if overlap:
                                    #Output "yes" in overlap column
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5], data[6],
                                                     data[7], data[8], data[9], "Yes", data[10]])

                                else:
                                    #Output "no" in overlap column
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5], data[6],
                                                     data[7], data[8], data[9], "No", data[10]])

                                    #Add start/end of alignment on query to intervals to check overlap
                                    alignment_intervals.append([data[6], data[7]])

                    #If this read aligned to MGE in iceberg
                    if megares_read.query_name in iceberg_data:
                        for data in iceberg_data[megares_read.query_name]:

                            #No need to check for overlap if no other alignments considered for this read yet
                            if not alignment_intervals:
                                csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                 ontologies[megares_read.reference_name][0],
                                                 ontologies[megares_read.reference_name][1],
                                                 ontologies[megares_read.reference_name][2],
                                                 megares_gene_lengths[megares_read.reference_name], megares_read.reference_start,
                                                 megares_read.reference_end - 1, megares_read.cigarstring, data[1], data[2],
                                                 data[3], data[4], data[5], data[6], data[7], data[8], data[9], "No",
                                                 megares_read.cigarstring])

                                #Add start/end of alignment on query to intervals
                                alignment_intervals.append([data[6], data[7]])
                            else:
                                #Check for overlap for other mge
                                overlap = False
                                for interval in alignment_intervals:
                                    #if alignment start/end position is in one of the other alignment intervals, then it overlaps
                                    #Recall:
                                    #   6 start of alignment on query
                                    #   7 end of alignment on query
                                    if (data[6] >= interval[0] and data[6] <= interval[1]) or \
                                       (data[7] >= interval[0] and data[7] <= interval[1]):
                                        overlap = True
                                        if data[6] < interval[0]:
                                            interval[0] = data[6]
                                        if data[7] < interval[0]:
                                            interval[1] = data[7]
                                        break
                                
                                if overlap:
                                    # Write "yes" in overlap column
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5],
                                                     data[6], data[7], data[8], data[9], "Yes", megares_read.cigarstring])

                                else:
                                    # Write "no" in overlap column
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5],
                                                     data[6], data[7], data[8], data[9], "No", megares_read.cigarstring])

                                    #Add start/end of alignment on query to intervals to check for overlap
                                    alignment_intervals.append([data[6], data[7]])

                    #If colocalized with plasmid finder MGE
                    if megares_read.query_name in pf_data:
                        for data in pf_data[megares_read.query_name]:

                            #No other alignments considered, no need to check for overlap
                            if not alignment_intervals:
                                csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                 ontologies[megares_read.reference_name][0],
                                                 ontologies[megares_read.reference_name][1],
                                                 ontologies[megares_read.reference_name][2],
                                                 megares_gene_lengths[megares_read.reference_name], megares_read.reference_start,
                                                 megares_read.reference_end - 1, megares_read.cigarstring, data[1], data[2],
                                                 data[3], data[4], data[5], data[6], data[7], data[8], data[9], "No",
                                                 megares_read.cigarstring])

                                #Check this interval in the future
                                alignment_intervals.append([data[6], data[7]])
                            else:
                                #Have to check for overlap
                                overlap = False
                                for interval in alignment_intervals:
                                    #if alignment start/end position is in one of the other alignment intervals, then it overlaps
                                    #Recall:
                                    #   6 start of alignment on query
                                    #   7 end of alignment on query
                                    if (data[6] >= interval[0] and data[6] <= interval[1]) or \
                                       (data[7] >= interval[0] and data[7] <= interval[1]):
                                        overlap = True
                                        if data[6] < interval[0]:
                                            interval[0] = data[6]
                                        if data[7] < interval[0]:
                                            interval[1] = data[7]
                                        break
                                
                                if overlap:
                                    #Write "yes" in overlap
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5], data[6],
                                                     data[7], data[8], data[9], "Yes", megares_read.cigarstring])

                                else:
                                    #Write "no" in overlap
                                    csv_wr.writerow(["", data[0], megares_query_length, megares_read.reference_name,
                                                     ontologies[megares_read.reference_name][0],
                                                     ontologies[megares_read.reference_name][1],
                                                     ontologies[megares_read.reference_name][2],
                                                     megares_gene_lengths[megares_read.reference_name],
                                                     megares_read.reference_start, megares_read.reference_end - 1,
                                                     megares_read.cigarstring, data[1], data[2], data[3], data[4], data[5], data[6],
                                                     data[7], data[8], data[9], "No", megares_read.cigarstring])

                                    #Look at this alignment interval next time for overlap
                                    alignment_intervals.append([data[6], data[7]])


        else:
            #Without MGEs...not currently used as we do analysis in the resistome folder now
            headers = ["SAMPLE_TYPE", "READ_ID", "READ_LENGTH", "AMR_GN_HEADER", "AMR_CLASS", "AMR_MECH", "AMR_GROUP", "AMR_LENGTH"
                       "AMR_START", "AMR_STOP", "AMR_CIGAR"]
            csv_wr.writerow(headers)

            for megares_read in megares_samfile.fetch():
                if "RequiresSNPConfirmation" in megares_read.reference_name:
                    continue
                if megares_read.is_secondary or megares_read.is_supplementary:
                    continue
                if (not megares_read.is_unmapped) and (megares_read.query_length - megares_read.query_alignment_length) > 1000 \
                                                  and megares_read.reference_length > 0.8*megares_gene_lengths[megares_read.reference_name]:
                    megares_query_length = megares_read.query_length

                    #Forget exactly why this is necessary
                    if( megares_read.query_length == 0 ):
                        megares_query_length = megares_read.infer_query_length()

                    csv_wr.writerow(["", megares_read.query_name, megares_query_length, megares_read.reference_name,
                                     ontologies[megares_read.reference_name][0], ontologies[megares_read.reference_name][1],
                                     ontologies[megares_read.reference_name][2], megares_gene_lengths[megares_read.reference_name],
                                     megares_read.reference_start, megares_read.reference_end - 1, megares_read.cigarstring])

    #Give output file read/write permissions
    os.chmod(out_filename, 0o660)
