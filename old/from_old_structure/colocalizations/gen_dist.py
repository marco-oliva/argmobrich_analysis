import pysam

import csv
import matplotlib.pyplot as plt
import os.path
import sys

if __name__ == "__main__":
    #Go through colocalization results, 
    results_filename = sys.argv[1]
    distance_info = {}
    with open(results_filename, 'r') as results_csv:
        print(results_filename)
        results_reader = csv.reader(results_csv, delimiter='\t')
        for row in results_reader:
            #Skip column names
            if row[0] == "SAMPLE_TYPE":
                continue

            #group = row[6]
            #if group == "TETO":
            read_id_with_unmapped = row[1].split('|')
            #Append megares and mge cigar strings to handle cases where same read has multiple colocalizations
            #unique_read_id = read_id_with_unmapped[0] + "||" + row[10] +  "||" + row[21]
            read_id = read_id_with_unmapped[0] #+ "||" + row[10] +  "||" + row[21]
            read_dict = {}
            #distance_info[read_id]["unmapped region start"]     = read_id_with_unmapped[1]
            #distance_info[read_id]["unmapped region stop"]      = read_id_with_unmapped[2]
            read_dict["read id"]                        = read_id
            read_dict["read length"]                    = int(row[2])
            read_dict["reference megares length"]       = int(row[7])
            read_dict["reference megares align start"]  = int(row[8])
            read_dict["reference megares align stop"]   = int(row[9])
            read_dict["megares cigar"]                  = row[10]
            read_dict["reference mge name"]             = row[11]
            read_dict["reference mge length"]           = int(row[13])
            read_dict["reference mge align start"]      = int(row[14])
            read_dict["reference mge align stop"]       = int(row[15])
            read_dict["read mge align start"]           = int(row[16])
            read_dict["read mge align stop"]            = int(row[17])
            read_dict["unmapped region start"]          = int(row[18])
            read_dict["unmapped region stop"]           = int(row[19])
            read_dict["mge cigar"]                      = row[21]
            read_dict["primary alignment"]              = True

            if read_id in distance_info:
                distance_info[read_id].append(read_dict)
            else:
                distance_info[read_id] = [read_dict]


    megares_align_file = pysam.AlignmentFile(sys.argv[2], 'r')
    for megares_alignment in megares_align_file.fetch():
        megares_query_name_with_unmapped = megares_alignment.query_name.split('|')
        megares_query_name = megares_query_name_with_unmapped[0]
        if megares_query_name in distance_info:
            read_dicts = distance_info[megares_query_name]
            for read_dict in read_dicts:
                if megares_alignment.cigarstring == read_dict["megares cigar"]:
                    
                    #Only care about primary alignments
                    if megares_alignment.is_secondary or megares_alignment.is_supplementary:
                        read_dict["primary alignment"] = False

                    cigar_stats_nt = (megares_alignment.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
                    read_dict["reference megares ins"]        = cigar_stats_nt[1]
                    read_dict["reference megares dels"]       = cigar_stats_nt[2]
                    read_dict["reference megares name"]       = megares_alignment.reference_name
                    if megares_alignment.is_reverse:
                        read_dict["read megares align start"]   = megares_alignment.query_length - megares_alignment.query_alignment_end
                        read_dict["read megares align stop"]    = megares_alignment.query_length - megares_alignment.query_alignment_start
                    else:
                        read_dict["read megares align start"]   = megares_alignment.query_alignment_start
                        read_dict["read megares align stop"]    = megares_alignment.query_alignment_end

                    #Determine nt distance
                    if read_dict["read megares align start"] > read_dict["read mge align stop"] + read_dict["unmapped region start"]:
                        read_dict["distance"] = read_dict["read megares align start"] - read_dict["read mge align stop"] \
                                                                                      - read_dict["unmapped region start"]
                    else:
                        read_dict["distance"] = read_dict["read mge align start"] + read_dict["unmapped region start"] \
                                                                                  - read_dict["read megares align stop"]


    aclame_align_file = pysam.AlignmentFile(sys.argv[3], 'r')
    for aclame_alignment in aclame_align_file.fetch():
        if aclame_alignment.is_unmapped:
            continue

        aclame_query_name_with_unmapped = aclame_alignment.query_name.split('|')
        aclame_query_name = aclame_query_name_with_unmapped[0]
        if aclame_query_name in distance_info:
            read_dicts = distance_info[aclame_query_name]
            for read_dict in read_dicts:
                if aclame_alignment.cigarstring == read_dict["mge cigar"]:
                    
                    #Only care about primary alignments
                    if aclame_alignment.is_secondary or aclame_alignment.is_supplementary:
                        read_dict["primary alignment"] = False

                    cigar_stats_nt = (aclame_alignment.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
                    read_dict["reference mge ins"]          = cigar_stats_nt[1]
                    read_dict["reference mge dels"]         = cigar_stats_nt[2]

    iceberg_align_file = pysam.AlignmentFile(sys.argv[4], 'r')
    for iceberg_alignment in iceberg_align_file.fetch():
        if iceberg_alignment.is_unmapped:
            continue

        iceberg_query_name_with_unmapped = iceberg_alignment.query_name.split('|')
        iceberg_query_name = iceberg_query_name_with_unmapped[0]
        if iceberg_query_name in distance_info:
            read_dicts = distance_info[iceberg_query_name]
            for read_dict in read_dicts:
                if iceberg_alignment.cigarstring == read_dict["mge cigar"]:
                    
                    #Only care about primary alignments
                    if iceberg_alignment.is_secondary or iceberg_alignment.is_supplementary:
                        read_dict["primary alignment"] = False

                    cigar_stats_nt = (iceberg_alignment.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
                    read_dict["reference mge ins"]          = cigar_stats_nt[1]
                    read_dict["reference mge dels"]         = cigar_stats_nt[2]

    pf_align_file = pysam.AlignmentFile(sys.argv[5], 'r')
    for pf_alignment in pf_align_file.fetch():
        if pf_alignment.is_unmapped:
            continue

        pf_query_name_with_unmapped = pf_alignment.query_name.split('|')
        pf_query_name = pf_query_name_with_unmapped[0]
        if pf_query_name in distance_info:
            read_dicts = distance_info[pf_query_name]
            for read_dict in read_dicts:
                if pf_alignment.cigarstring == read_dict["mge cigar"]:
                    
                    #Only care about primary alignments
                    if pf_alignment.is_secondary or pf_alignment.is_supplementary:
                        read_dict["primary alignment"] = False

                    cigar_stats_nt = (pf_alignment.get_cigar_stats())[0] #0th element is nucleotide (nt) counts of each cigar op
                    read_dict["reference mge ins"]          = cigar_stats_nt[1]
                    read_dict["reference mge dels"]         = cigar_stats_nt[2]



    #Output in blocks by the mge/plasmid
    #sorted_read_ids = sorted(distance_info, key = lambda read_id: read_id.split("||")[0])
    distance_info_list = []
    for read_id in distance_info:
        for read_dict in distance_info[read_id]:
            if not read_dict["primary alignment"]:
                continue
            if not "read megares align start" in read_dict or not "reference mge ins" in read_dict:
                continue
            if read_dict["read megares align start"] >= 0:
                distance_info_list.append(read_dict)
    sorted_dist_list = sorted(distance_info_list, key = lambda read_dict: read_dict["read id"])


    #Write tsv for group colocalizations
    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_distance.tsv", 'w') as dist_tsv:
        dist_writer = csv.writer(dist_tsv, delimiter='\t')
        #dist_writer.writerow([os.path.splitext(os.path.basename(sys.argv[1]))[0][:6] + " colocalization distances for TETO alignments"])
        dist_writer.writerow([os.path.splitext(os.path.basename(sys.argv[1]))[0][:6] + " colocalization distances all alignments"])
        dist_writer.writerow([])
        dist_writer.writerow(["Total number:", len(sorted_dist_list)])
        dist_writer.writerow([])
        dist_writer.writerow(["Read id", "MEGARes name", "MGE name", "Distance between genes", "Read length", "Read AMR start",
                              "Read AMR stop", "MEGARes align start", "MEGARes align stop", "MEGARes align ins",
                              "MEGARes align dels", "MEGARes sequence length", "Read MGE start", "Read MGE stop",
                              "MGE align start", "MGE align stop", "MGE align ins", "MGE align dels", "MGE sequence length"])
        for read_dict in sorted_dist_list:
            dist_writer.writerow([read_dict["read id"], read_dict["reference megares name"], read_dict["reference mge name"],
                                    #read_dict["unmapped region start"], read_dict["unmapped region stop"],
                                    read_dict["distance"],
                                    read_dict["read length"],
                                    read_dict["read megares align start"], read_dict["read megares align stop"],
                                    read_dict["reference megares align start"], read_dict["reference megares align stop"],
                                    read_dict["reference megares ins"], read_dict["reference megares dels"],
                                    read_dict["reference megares length"],
                                    read_dict["read mge align start"] + read_dict["unmapped region start"],
                                    read_dict["read mge align stop"] + read_dict["unmapped region start"],
                                    read_dict["reference mge align start"], read_dict["reference mge align stop"],
                                    read_dict["reference mge ins"], read_dict["reference mge dels"],
                                    read_dict["reference mge length"]])


    #Distance vs num colocalizations
    distances     = [read_dict["distance"] for read_dict in sorted_dist_list if "distance" in read_dict]

    if distances:
        fig, ax = plt.subplots(tight_layout=True)
        ax.set_xlabel("Absolute distance between AMR genes and MGEs (nt)")
        ax.set_ylabel("Number of colocalizations")
        ax.hist(distances, bins="fd", range=(min(distances), max(distances)))

        basename = os.path.splitext(os.path.basename(sys.argv[1]))[0]
        if "MOCK" in basename:
            ax.set_title(basename[:7] + " distances vs num colocalizations")
            plt.savefig(basename[:7] + "_dist_hist.png")
        else:
            ax.set_title(basename[:2] + " distances vs num colocalizations")
            plt.savefig( basename[:2] + "_dist_hist.png")

    #Coverage vs num colocalizations
    gene_coverages = [(read_dict["reference mge align stop"] - read_dict["reference mge align start"]) / read_dict["reference mge length"] \
                        for read_dict in sorted_dist_list if "reference mge align start" in read_dict]
    if gene_coverages:
        fig, ax = plt.subplots(tight_layout=True)
        ax.set_title("MGE coverage vs num colocalizations")
        ax.set_xlabel("Coverage of MGE alignments (%)")
        ax.set_ylabel("Number of colocalizations")
        ax.hist(gene_coverages, bins="fd", range=(min(gene_coverages), max(gene_coverages)))

        plt.savefig(os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_coverage_dist_hist.png")
