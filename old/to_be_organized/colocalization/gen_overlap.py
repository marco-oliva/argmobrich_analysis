import csv
import itertools
import os.path
import sys

if __name__ == "__main__":
    #Create ontology dictionary from MEGARes ontology file
    #megares_ontology = {}
    #ontology_filename = "/home/noyes046/jsettle/argmobrich/MEGARESONTOLOGY.tsv"
    #with open(ontology_filename, 'r') as ontology_tsv:
        #ontology_reader = csv.reader(ontology_tsv, delimiter='\t')
        #for row in ontology_reader:
            ##Skip column names
            #if row[0] == "header":
                #continue
#
            ##FIll in our dict
            #megares_ontology[row[0]] = { "class"        : row[1],
                                         #"mechanism"    : row[2],
                                         #"group"        : row[3]
                                       #}



    #Go through colocalization results, looking for overlaps
    results_filename = sys.argv[1]
    overlap_dict ={}
    tmp_overlaps = []
    with open(results_filename, 'r') as results_csv:
        results_reader = csv.reader(results_csv, delimiter='\t')
        for row in results_reader:
            #Skip column names
            if row[0] == "SAMPLE_TYPE":
                continue

            if row[18] == "No":
                #If there were overlaps, record them
                if len(tmp_overlaps) > 1:
                    overlaps_processed = []
                    for overlap in itertools.product(tmp_overlaps, repeat=2):
                        #Not interested in alignments to same read
                        if overlap[0] == overlap[1]:
                            continue

                        #(A,B) = (B,A) for this purpose
                        if (tuple(sorted(overlap)) in overlaps_processed):
                            continue

                        if overlap in overlap_dict:
                            overlap_dict[overlap] += 1
                        else:
                            overlap_dict[overlap] = 1
                        overlaps_processed.append(tuple(sorted(overlap)))

                tmp_overlaps = [row[11]]
            else: #(row[16] == "Yes")
                tmp_overlaps.append(row[11])

    sorted_overlaps = sorted(overlap_dict, key = lambda overlap: overlap_dict[overlap], reverse=True)

    #Write tsv for overlap counts
    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_overlaps.tsv", 'w') as coloc_tsv:
        coloc_writer = csv.writer(coloc_tsv, delimiter='\t')
        coloc_writer.writerow([os.path.splitext(os.path.basename(sys.argv[1]))[0][:6] + " overlaps"])
        coloc_writer.writerow([])
        coloc_writer.writerow(["Overlap Pair", "Occurrences"])
        for overlap in sorted_overlaps:
            coloc_writer.writerow([overlap, overlap_dict[overlap]]) 

