from Bio import SeqIO
import pysam

from matplotlib_venn import venn2

import csv
import matplotlib.pyplot as plt
import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

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

    gene_lengths = { "aclame": aclame_gene_lengths,
                     "iceberg": iceberg_gene_lengths,
                     "pf": pf_gene_lengths}

    a01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(), }
    a01_aclame_sam_file = pysam.AlignmentFile(sys.argv[1]+"_aclame.sam", 'r')
    a01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[1]+"_iceberg.sam", 'r')
    a01_pf_sam_file = pysam.AlignmentFile(sys.argv[1]+"_plasmidfinder.sam", 'r')
    a01_sam_files = { "aclame": a01_aclame_sam_file,
                      "iceberg": a01_iceberg_sam_file,
                      "pf": a01_pf_sam_file}

    b01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(),}
    b01_aclame_sam_file = pysam.AlignmentFile(sys.argv[2]+"_aclame.sam", 'r')
    b01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[2]+"_iceberg.sam", 'r')
    b01_pf_sam_file = pysam.AlignmentFile(sys.argv[2]+"_plasmidfinder.sam", 'r')
    b01_sam_files = { "aclame": b01_aclame_sam_file,
                      "iceberg": b01_iceberg_sam_file,
                      "pf": b01_pf_sam_file}

    c01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(),}
    c01_aclame_sam_file = pysam.AlignmentFile(sys.argv[3]+"_aclame.sam", 'r')
    c01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[3]+"_iceberg.sam", 'r')
    c01_pf_sam_file = pysam.AlignmentFile(sys.argv[3]+"_plasmidfinder.sam", 'r')
    c01_sam_files = { "aclame": c01_aclame_sam_file,
                      "iceberg": c01_iceberg_sam_file,
                      "pf": c01_pf_sam_file}

    d01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(),}
    d01_aclame_sam_file = pysam.AlignmentFile(sys.argv[4]+"_aclame.sam", 'r')
    d01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[4]+"_iceberg.sam", 'r')
    d01_pf_sam_file = pysam.AlignmentFile(sys.argv[4]+"_plasmidfinder.sam", 'r')
    d01_sam_files = { "aclame": d01_aclame_sam_file,
                      "iceberg": d01_iceberg_sam_file,
                      "pf": d01_pf_sam_file}

    e01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(),}
    e01_aclame_sam_file = pysam.AlignmentFile(sys.argv[5]+"_aclame.sam", 'r')
    e01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[5]+"_iceberg.sam", 'r')
    e01_pf_sam_file = pysam.AlignmentFile(sys.argv[5]+"_plasmidfinder.sam", 'r')
    e01_sam_files = { "aclame": e01_aclame_sam_file,
                      "iceberg": e01_iceberg_sam_file,
                      "pf": e01_pf_sam_file}

    f01_aligns = { "aclame": set(),
                   "iceberg": set(),
                   "pf": set(),}
    f01_aclame_sam_file = pysam.AlignmentFile(sys.argv[6]+"_aclame.sam", 'r')
    f01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[6]+"_iceberg.sam", 'r')
    f01_pf_sam_file = pysam.AlignmentFile(sys.argv[6]+"_plasmidfinder.sam", 'r')
    f01_sam_files = { "aclame": f01_aclame_sam_file,
                      "iceberg": f01_iceberg_sam_file,
                      "pf": f01_pf_sam_file,}

    align_dicts = [a01_aligns, b01_aligns, c01_aligns, d01_aligns, e01_aligns, f01_aligns]
    sam_file_dicts = [a01_sam_files, b01_sam_files, c01_sam_files, d01_sam_files, e01_sam_files, f01_sam_files]

    feces_richnesses = {}
    mock_richnesses = {}

    for i in range(0,len(align_dicts)):
        align_dict = align_dicts[i]
        sam_file_dict = sam_file_dicts[i]
        
        for key in sam_file_dict:
            sam_file = sam_file_dict[key]
            align_set = align_dict[key]

            for read in sam_file.fetch():
                if read.is_unmapped:
                    continue
                if (float(read.reference_length) / gene_lengths[key][read.reference_name]) < 0.5:
                    continue

                align_dict[key].add(read.reference_name)

                if key == "aclame":
                    if i < len(align_dicts)/2:
                        if not read.reference_name in feces_richnesses:
                            feces_richnesses[read.reference_name] = 1
                        else:
                            feces_richnesses[read.reference_name] += 1
                    else:
                        if not read.reference_name in mock_richnesses:
                            mock_richnesses[read.reference_name] = 1
                        else:
                            mock_richnesses[read.reference_name] += 1
            sam_file.close()


    plt.rc('axes', titlesize=24)
    plt.rc('figure', titlesize=32)
    #Venn diag for:
    #mock v feces aclame
    fig = plt.figure()
    ax = fig.subplots()
    aclame_venn = venn2([(a01_aligns["aclame"] | b01_aligns["aclame"] | c01_aligns["aclame"]),
                              (d01_aligns["aclame"] | e01_aligns["aclame"] | f01_aligns["aclame"])],
                              set_labels=["Feces", "Mock"]) 

    ax.set_title("Aclame")
    plt.savefig("aclame_venn.png")

    #mock v feces iceberg
    fig = plt.figure()
    ax = fig.subplots()
    iceberg_venn = venn2([(a01_aligns["iceberg"] | b01_aligns["iceberg"] | c01_aligns["iceberg"]),
                          (d01_aligns["iceberg"] | e01_aligns["iceberg"] | f01_aligns["iceberg"])],
                          set_labels=["Feces", "Mock"]) 
    ax.set_title("Iceberg")
    plt.savefig("iceberg_venn.png")

    #mock v feces plasmid finder
    fig = plt.figure()
    ax = fig.subplots()
    pf_venn = venn2([(a01_aligns["pf"] | b01_aligns["pf"] | c01_aligns["pf"]),
                     (d01_aligns["pf"] | e01_aligns["pf"] | f01_aligns["pf"])],
                     set_labels=["Feces", "Mock"]) 

    ax.set_title("Plasmid Finder")
    plt.savefig("pf_venn.png")

    feces_only = (a01_aligns["aclame"] | b01_aligns["aclame"] | c01_aligns["aclame"]) - \
                 (d01_aligns["aclame"] | e01_aligns["aclame"] | f01_aligns["aclame"])

    mock_only = (d01_aligns["aclame"] | e01_aligns["aclame"] | f01_aligns["aclame"]) - \
                (a01_aligns["aclame"] | b01_aligns["aclame"] | c01_aligns["aclame"])

    
    feces_only = sorted(feces_only, key=lambda ref:feces_richnesses[ref], reverse=True)
    mock_only = sorted(mock_only, key=lambda ref:mock_richnesses[ref], reverse=True)






















































    a01_aclame_sam_file = pysam.AlignmentFile(sys.argv[1]+"_aclame.sam", 'r')
    a01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[1]+"_iceberg.sam", 'r')
    a01_pf_sam_file = pysam.AlignmentFile(sys.argv[1]+"_plasmidfinder.sam", 'r')
    a01_sam_files = { "aclame": a01_aclame_sam_file,
                      "iceberg": a01_iceberg_sam_file,
                      "pf": a01_pf_sam_file}

    b01_aclame_sam_file = pysam.AlignmentFile(sys.argv[2]+"_aclame.sam", 'r')
    b01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[2]+"_iceberg.sam", 'r')
    b01_pf_sam_file = pysam.AlignmentFile(sys.argv[2]+"_plasmidfinder.sam", 'r')
    b01_sam_files = { "aclame": b01_aclame_sam_file,
                      "iceberg": b01_iceberg_sam_file,
                      "pf": b01_pf_sam_file}

    c01_aclame_sam_file = pysam.AlignmentFile(sys.argv[3]+"_aclame.sam", 'r')
    c01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[3]+"_iceberg.sam", 'r')
    c01_pf_sam_file = pysam.AlignmentFile(sys.argv[3]+"_plasmidfinder.sam", 'r')
    c01_sam_files = { "aclame": c01_aclame_sam_file,
                      "iceberg": c01_iceberg_sam_file,
                      "pf": c01_pf_sam_file}

    d01_aclame_sam_file = pysam.AlignmentFile(sys.argv[4]+"_aclame.sam", 'r')
    d01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[4]+"_iceberg.sam", 'r')
    d01_pf_sam_file = pysam.AlignmentFile(sys.argv[4]+"_plasmidfinder.sam", 'r')
    d01_sam_files = { "aclame": d01_aclame_sam_file,
                      "iceberg": d01_iceberg_sam_file,
                      "pf": d01_pf_sam_file}

    e01_aclame_sam_file = pysam.AlignmentFile(sys.argv[5]+"_aclame.sam", 'r')
    e01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[5]+"_iceberg.sam", 'r')
    e01_pf_sam_file = pysam.AlignmentFile(sys.argv[5]+"_plasmidfinder.sam", 'r')
    e01_sam_files = { "aclame": e01_aclame_sam_file,
                      "iceberg": e01_iceberg_sam_file,
                      "pf": e01_pf_sam_file}

    f01_aclame_sam_file = pysam.AlignmentFile(sys.argv[6]+"_aclame.sam", 'r')
    f01_iceberg_sam_file = pysam.AlignmentFile(sys.argv[6]+"_iceberg.sam", 'r')
    f01_pf_sam_file = pysam.AlignmentFile(sys.argv[6]+"_plasmidfinder.sam", 'r')
    f01_sam_files = { "aclame": f01_aclame_sam_file,
                      "iceberg": f01_iceberg_sam_file,
                      "pf": f01_pf_sam_file,}

    sam_file_dicts = [a01_sam_files, b01_sam_files, c01_sam_files, d01_sam_files, e01_sam_files, f01_sam_files]

    mock_only_reads = [set(), set(), set()]
    for i in range(3, len(sam_file_dicts)):
        sam_file_dict = sam_file_dicts[i]
        mock_only_read_set = mock_only_reads[i-3]
        
        for key in sam_file_dict:
            sam_file = sam_file_dict[key]

            for read in sam_file.fetch():
                if read.is_unmapped:
                    continue
                if (float(read.reference_length) / gene_lengths[key][read.reference_name]) < 0.5:
                    continue
                mock_only_read_set.add(read.query_name)


    mock_fastqs = ["../../ccs_fastqs/sequel-demultiplex.MOCK_D01.ccs.fastq", "../../ccs_fastqs/sequel-demultiplex.MOCK_E01.ccs.fastq",
                   "../../ccs_fastqs/sequel-demultiplex.MOCK_F01.ccs.fastq"]

    out_mock_fastqs = ["MOCK_D01_unique_mge_reads.fastq", "MOCK_E01_unique_mge_reads.fastq",
                       "MOCK_F01_unique_mge_reads.fastq"]


    for i in range(0, len(mock_fastqs)):
        mock_only_read_set = mock_only_reads[i]
        records = []
        for rec in SeqIO.parse(mock_fastqs[i], "fastq"):
            if rec.name in mock_only_read_set:
                records.append(rec)
        with open(out_mock_fastqs[i], 'w') as out_fastq:
            SeqIO.write(records, out_fastq, "fastq")

    sys.exit(0)
