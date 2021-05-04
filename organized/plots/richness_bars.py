from Bio import SeqIO
import pysam

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

    #Plot title on first line of text file
    #Subsequent lines are "dataset name, alignment file"
    params = []
    with open(sys.argv[1], 'r') as input_handle:
        plot_title = input_handle.readline()[:-1]
        for line in input_handle:
            line_params = line.split(',')
            descriptor = line_params[0]
            sam_filename = line_params[1][:-1]
            params.append((descriptor, sam_filename))

    cls = set()
    richnesses = {}
    for param in params:
        #param[0] is dataset name
        #param[1] is alignmen file

        sam_file = pysam.AlignmentFile(param[1], 'r')
        gene_dict = {}
        class_dict = {}
        mech_dict = {}
        group_dict = {}
        #Iterate through every read. Accumulate number of reads while recording read length
        for read in sam_file.fetch():
            if "RequiresSNPConfirmation" in read.reference_name:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
                
            #check coverage
            if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:
                classname = megares_ontology[read.reference_name]["class"]
                mech = megares_ontology[read.reference_name]["mechanism"]
                group = megares_ontology[read.reference_name]["group"]

                #update gene dict
                if( not read.reference_name in gene_dict):
                    gene_dict[read.reference_name] = 1
                else:
                    gene_dict[read.reference_name] += 1

                #update class dict
                if( not classname in class_dict):
                    class_dict[classname] = 1
                    cls.add(classname)
                else:
                    class_dict[classname] += 1

                #update mechanism dict
                if( not mech in mech_dict):
                    mech_dict[mech] = 1
                else:
                    mech_dict[mech] += 1

                #update group dict
                if( not group in group_dict):
                    group_dict[group] = 1
                else:
                    group_dict[group] += 1

        sam_file.close()
        richnesses[param[0]] = { "class" : class_dict }

    #Use unique set of classes to make bar graph labels
    labels = list(cls)

    #Make sure that each class dict has all classes we want to label by inserting zeros where necessary
    for param in params:
        class_dict = richnesses[param[0]]["class"]
        for label in labels:
            if label not in class_dict:
                class_dict[label] = 0

    #Put classes with highest abundance first
    labels.sort(key=lambda label: richnesses[params[0][0]]["class"][label], reverse=True)

    #Each bar position will have 3 bars, one for 2kb, 5kb, and 8kb samples
    bar_width = 0.6
    bar_pos = list(range(0,2*len(labels),2))
    grouped_bar_pos = [[pos - bar_width/2.0 for pos in bar_pos], 
                       [pos + bar_width/2.0 for pos in bar_pos],
                       [pos + 1.5*bar_width for pos in bar_pos]]
    print(bar_pos)

    #Easier to use matplotlib here
    plt.rcParams.update({'figure.autolayout': True})
    fig, ax = plt.subplots()

    #Add each grouped bar to the plot
    i = 0
    counts = []
    for param in params:
        counts = [richnesses[param[0]]["class"][label] for label in labels]
        ax.bar(grouped_bar_pos[i], counts, bar_width, label=param[0])
        i += 1

    #Label plot
    ax.set_xticks(bar_pos)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    plt.savefig(plot_title + ".png")

    sys.exit(0)
