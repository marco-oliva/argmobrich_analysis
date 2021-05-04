import csv
import os.path
import sys

if __name__ == "__main__":
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



    #Go through colocalization results, looking at those whose AMR gene is in the given groups
    total_colocs = 0
    results_filename = sys.argv[1]
    #coloc_dicts = {"TETQ" : {}, "TETW" : {}, "TETO" : {}, "NORB" : {}}
    coloc_dicts = {}
    with open(results_filename, 'r') as results_csv:
        results_reader = csv.reader(results_csv, delimiter='\t')
        for row in results_reader:
            #Skip column names
            if row[0] == "SAMPLE_TYPE":
                continue

            total_colocs += 1
            group =  megares_ontology[row[3]]["group"]
            if group in coloc_dicts:
                if row[11] in coloc_dicts[group]:
                    coloc_dicts[group][row[11]] += 1
                    coloc_dicts[group]["Total"] += 1
                else:
                    coloc_dicts[group][row[11]] = 1
                    coloc_dicts[group]["Total"] += 1
            else:
                coloc_dicts[group] = {}
                coloc_dicts[group][row[11]] = 1
                coloc_dicts[group]["Total"] = 1

    sorted_groups = sorted(coloc_dicts, key = lambda coloc_dict: coloc_dicts[coloc_dict]["Total"], reverse=True)

    #Write tsv for group colocalizations
    with open(os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_diversity.tsv", 'w') as coloc_tsv:
        coloc_writer = csv.writer(coloc_tsv, delimiter='\t')
        coloc_writer.writerow([os.path.splitext(os.path.basename(sys.argv[1]))[0][:6] + " colocalization results by MEGARes group"])
        coloc_writer.writerow(["Total colocalizations: " + str(total_colocs)])
        coloc_writer.writerow([])
        for group in sorted_groups:
            coloc_writer.writerow([group + " colocalizations"])
            mobilome_hits = coloc_dicts[group].keys()
            mobilome_hits = sorted(mobilome_hits, key= lambda mobilome_name: coloc_dicts[group][mobilome_name], reverse=True)
            coloc_writer.writerow(mobilome_hits)

            coloc_counts = sorted(coloc_dicts[group].values(), reverse=True)
            coloc_writer.writerow(coloc_counts)
            coloc_writer.writerow([])

