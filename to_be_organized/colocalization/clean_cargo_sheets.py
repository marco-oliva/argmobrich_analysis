
import csv
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


    in_csv = sys.argv[1]
    out_csv = sys.argv[2]

    data_rows = {}
    with open(in_csv, 'r') as in_csv_handle:
        in_reader = csv.reader(in_csv_handle, delimiter='\t')
        for row in in_reader:
            split_id = row[0].split(":::")
            read_name = split_id[0]
            megares_group = megares_ontology[split_id[1]]["group"]
            if not read_name in data_rows:
                data_rows[read_name] = [read_name, megares_group, split_id[2]]

            data_rows[read_name].extend([row[1] + " [" + row[6] + ", " + row[7] + ']'])

    sorted_read_names = sorted(data_rows, key = lambda read_name: data_rows[read_name][1])
    with open(out_csv, 'w') as out_csv_handle:
        out_writer = csv.writer(out_csv_handle, delimiter='\t')

        header = ["Read Name", "MEGARes Group", "MGE Header", "List of potential cargo genes"]
        out_writer.writerow(header)
        for read_name in sorted_read_names:
            out_writer.writerow(data_rows[read_name])

    sys.exit(0)
