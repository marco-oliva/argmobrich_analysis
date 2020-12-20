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

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    if sys.argv[2] == "1":
        sam_file = pysam.AlignmentFile(sys.argv[1] + "_megares_mapped.sam", 'r')

        for read in sam_file.fetch():
            if "RequiresSNPConfirmation" in read.reference_name:
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

    else:
        for i in range(1, int(sys.argv[2])+1):
            if i < 10:
                sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-0" + str(i) + "_megares_mapped.sam", 'r')
            else:
                sam_file = pysam.AlignmentFile(sys.argv[1] + ".part-" + str(i) + "_megares_mapped.sam", 'r')

            for read in sam_file.fetch():
                if "RequiresSNPConfirmation" in read.reference_name:
                    continue
                    
                #check coverage
                if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:
                    #split header to get ontology
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

    #Prepare rows of diversity tsv
    tsv_rows = [["MEGARes Gene Header", "Num Reads", "", "Group", "Num Reads", "", "Mechanism", "Num Reads", "", "Class", "Num Reads"]]
    for gene in sorted(gene_dict, key=lambda gene: gene_dict[gene], reverse=True):
        tsv_rows.append([gene, gene_dict[gene]])

    i=1
    for group in sorted(group_dict, key=lambda group: group_dict[group], reverse=True):
        tsv_rows[i].extend(["", group, group_dict[group]])
        i+=1

    i=1
    for mech in sorted(mech_dict, key=lambda mech: mech_dict[mech], reverse=True):
        tsv_rows[i].extend(["", mech, mech_dict[mech]])
        i+=1

    i=1
    for class_name in sorted(class_dict, key=lambda class_name: class_dict[class_name], reverse=True):
        tsv_rows[i].extend(["", class_name, class_dict[class_name]])
        i+=1

    #Write diversity tsv
    with open(os.path.basename(sys.argv[1]) + "_amr_diversity.tsv", 'w') as out_tsv:
        tsv_writer = csv.writer(out_tsv, delimiter='\t')
        tsv_writer.writerows(tsv_rows)

    #Prepare rows of richness tsv
    tsv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    tsv_rows.append([len(gene_dict), len(class_dict), len(mech_dict), len(group_dict)])

    #Write richness tsv
    with open(os.path.basename(sys.argv[1]) + "_amr_richness.tsv", 'w') as out_tsv:
        tsv_writer = csv.writer(out_tsv, delimiter='\t')
        tsv_writer.writerows(tsv_rows)

    sys.exit(0)
