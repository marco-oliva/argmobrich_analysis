from Bio import SeqIO
import pysam

from matplotlib_venn import venn2, venn3

import csv
import matplotlib.pyplot as plt
import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    #Get megares lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    #Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    ontology_filename = "/home/noyes046/jsettle/argmobrich/MEGARESONTOLOGY.tsv"
    with open(ontology_filename, 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv, delimiter='\t')
        for row in ontology_reader:
            #Skip column names
            if row[0] == "header":
                continue

            #FIll in our dict
            megares_ontology[row[0]] = { "class"        : row[1],
                                         "mechanism"    : row[2],
                                         "group"        : row[3]
                                       }

    a01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    a01_sam_file = pysam.AlignmentFile(sys.argv[1], 'r')

    b01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    b01_sam_file = pysam.AlignmentFile(sys.argv[2], 'r')

    c01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    c01_sam_file = pysam.AlignmentFile(sys.argv[3], 'r')

    d01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    d01_sam_file = pysam.AlignmentFile(sys.argv[4], 'r')

    e01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    e01_sam_file = pysam.AlignmentFile(sys.argv[5], 'r')

    f01_aligns = { "genes": set(),
                   "groups": set(),
                   "mechanisms": set(),
                   "classes": set()}
    f01_sam_file = pysam.AlignmentFile(sys.argv[6], 'r')

    align_dicts = [a01_aligns, b01_aligns, c01_aligns, d01_aligns, e01_aligns, f01_aligns]
    sam_files = [a01_sam_file, b01_sam_file, c01_sam_file, d01_sam_file, e01_sam_file, f01_sam_file]

    for i in range(0,len(align_dicts)):
        align_dict = align_dicts[i]
        sam_file = sam_files[i]

        for read in sam_file.fetch():
            #Check SNP confirmation and coverage
            if "RequiresSNPConfirmation" in read.reference_name:
                continue
            if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) < 0.8:
                continue

            classname = megares_ontology[read.reference_name]["class"]
            group = megares_ontology[read.reference_name]["group"]
            mech = megares_ontology[read.reference_name]["mechanism"]

            align_dict["genes"].add(read.reference_name)
            align_dict["groups"].add(group)
            align_dict["mechanisms"].add(mech)
            align_dict["classes"].add(classname)


    plt.rc('axes', titlesize=24)
    plt.rc('figure', titlesize=32)
    #Venn diag for:
    #3 feces
    fig, axs = plt.subplots(2,2, figsize=(15,15))
    feces_genes_venn = venn3([a01_aligns["genes"], b01_aligns["genes"], c01_aligns["genes"]], set_labels=["A01", "B01", "C01"], 
                             ax=axs[0][0])
    feces_groups_venn = venn3([a01_aligns["groups"], b01_aligns["groups"], c01_aligns["groups"]], set_labels=["A01", "B01", "C01"],
                              ax=axs[0][1])
    feces_mechanisms_venn = venn3([a01_aligns["mechanisms"], b01_aligns["mechanisms"], c01_aligns["mechanisms"]],
                                   set_labels=["A01", "B01", "C01"], ax=axs[1][0])
    feces_classes_venn = venn3([a01_aligns["classes"], b01_aligns["classes"], c01_aligns["classes"]], set_labels=["A01", "B01", "C01"],
                                ax=axs[1][1])

    fig.suptitle("Feces MEGARes Ontology")
    axs[0][0].set_title("Genes")
    axs[0][1].set_title("Groups")
    axs[1][0].set_title("Mechanisms")
    axs[1][1].set_title("Classes")
    plt.savefig("feces_venn.png")

    #3 mock
    fig, axs = plt.subplots(2,2, figsize=(15,15))
    mock_genes_venn = venn3([d01_aligns["genes"], e01_aligns["genes"], f01_aligns["genes"]], set_labels=["D01", "E01", "F01"], 
                             ax=axs[0][0])
    mock_groups_venn = venn3([d01_aligns["groups"], e01_aligns["groups"], f01_aligns["groups"]], set_labels=["D01", "E01", "F01"],
                              ax=axs[0][1])
    mock_mechanisms_venn = venn3([d01_aligns["mechanisms"], e01_aligns["mechanisms"], f01_aligns["mechanisms"]],
                                   set_labels=["D01", "E01", "F01"], ax=axs[1][0])
    mock_classes_venn = venn3([d01_aligns["classes"], e01_aligns["classes"], f01_aligns["classes"]], set_labels=["D01", "E01", "F01"],
                                ax=axs[1][1])

    fig.suptitle("Mock MEGARes Ontology")
    axs[0][0].set_title("Genes")
    axs[0][1].set_title("Groups")
    axs[1][0].set_title("Mechanisms")
    axs[1][1].set_title("Classes")
    plt.savefig("mock_venn.png")

    #mock vs feces
    fig, axs = plt.subplots(2,2, figsize=(15,15))
    genes_venn = venn2([(a01_aligns["genes"] | b01_aligns["genes"] | c01_aligns["genes"]),
                        (d01_aligns["genes"] | e01_aligns["genes"] | f01_aligns["genes"])],
                        set_labels=["Mock", "Feces"], ax=axs[0][0])
    groups_venn = venn2([(a01_aligns["groups"] | b01_aligns["groups"] | c01_aligns["groups"]),
                        (d01_aligns["groups"] | e01_aligns["groups"] | f01_aligns["groups"])],
                        set_labels=["Mock", "Feces"], ax=axs[0][1])
    mechanisms_venn = venn2([(a01_aligns["mechanisms"] | b01_aligns["mechanisms"] | c01_aligns["mechanisms"]),
                        (d01_aligns["mechanisms"] | e01_aligns["mechanisms"] | f01_aligns["mechanisms"])],
                        set_labels=["Mock", "Feces"], ax=axs[1][0])
    classes_venn = venn2([(a01_aligns["classes"] | b01_aligns["classes"] | c01_aligns["classes"]),
                        (d01_aligns["classes"] | e01_aligns["classes"] | f01_aligns["classes"])],
                        set_labels=["Mock", "Feces"], ax=axs[1][1])

    fig.suptitle("Mock vs. Feces MEGARes Ontology")
    axs[0][0].set_title("Genes")
    axs[0][1].set_title("Groups")
    axs[1][0].set_title("Mechanisms")
    axs[1][1].set_title("Classes")
    plt.savefig("mvf_venn.png")


    sys.exit(0)
