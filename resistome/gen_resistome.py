from Bio import SeqIO
import pysam
import argparse

import csv

GLOBAL_AMR_THRESHOLD = 0.8

GLOBAL_MEGARES_ONTOLOGY_PATH = "/blue/boucher/marco.oliva/data/MEGARes/V2/megares_modified_annotations_v2.00.csv"
GLOBAL_MEGARES_SEQS_PATH = "/blue/boucher/marco.oliva/data/MEGARes/V2/megares_full_database_v2.00.fasta"


def main():
    parser = argparse.ArgumentParser(description='Compute resistome')
    parser.add_argument('-s', help='Alignment file', dest='sam_file', required=True)
    parser.add_argument('-o', help='Output Prefix', dest='out_prefix', required=True)
    args = parser.parse_args()

    # Get megares lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = GLOBAL_MEGARES_SEQS_PATH
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    sam_file = pysam.AlignmentFile(args.sam_file, 'r')

    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    ontology_filename = GLOBAL_MEGARES_ONTOLOGY_PATH
    with open(ontology_filename, 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            # FIll in our dict
            megares_ontology[row[0]] = {"class": row[1],
                                        "mechanism": row[2],
                                        "group": row[3]
                                        }

    gene_dict = {}
    class_dict = {}
    mech_dict = {}
    group_dict = {}

    # Iterate through every read. Accumulate number of reads while recording read length
    for read in sam_file.fetch():
        if read.is_unmapped:
            continue
        if "RequiresSNPConfirmation" in read.reference_name:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        # check coverage
        if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > GLOBAL_AMR_THRESHOLD:
            classname = megares_ontology[read.reference_name]["class"]
            mech = megares_ontology[read.reference_name]["mechanism"]
            group = megares_ontology[read.reference_name]["group"]

            # update gene dict
            if (not read.reference_name in gene_dict):
                gene_dict[read.reference_name] = 1
            else:
                gene_dict[read.reference_name] += 1

            # update class dict
            if (not classname in class_dict):
                class_dict[classname] = 1
            else:
                class_dict[classname] += 1

            # update mechanism dict
            if (not mech in mech_dict):
                mech_dict[mech] = 1
            else:
                mech_dict[mech] += 1

            # update group dict
            if (not group in group_dict):
                group_dict[group] = 1
            else:
                group_dict[group] += 1

    # Prepare rows of diversity csv
    csv_rows = [
        ["MEGARes Gene Header", "Num Reads", "Group", "Num Reads", "Mechanism", "Num Reads", "Class",
         "Num Reads"]]
    for gene in sorted(gene_dict, key=lambda gene: gene_dict[gene], reverse=True):
        csv_rows.append([gene, gene_dict[gene]])

    i = 1
    for group in sorted(group_dict, key=lambda group: group_dict[group], reverse=True):
        csv_rows[i].extend([group, group_dict[group]])
        i += 1

    i = 1
    for mech in sorted(mech_dict, key=lambda mech: mech_dict[mech], reverse=True):
        csv_rows[i].extend([mech, mech_dict[mech]])
        i += 1

    i = 1
    for class_name in sorted(class_dict, key=lambda class_name: class_dict[class_name], reverse=True):
        csv_rows[i].extend([class_name, class_dict[class_name]])
        i += 1

    # Write diversity tsv
    with open(args.out_prefix + "_amr_diversity.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

    # Prepare rows of richness tsv
    csv_rows = [["Gene Richness", "Class Richness", "Mechanism Richness", "Group Richness"]]
    csv_rows.append([len(gene_dict), len(class_dict), len(mech_dict), len(group_dict)])

    # Write richness tsv
    with open(args.out_prefix + "_amr_richness.csv", 'w') as out_csv:
        csv_writer = csv.writer(out_csv)
        csv_writer.writerows(csv_rows)

if __name__ == "__main__":
    main()

