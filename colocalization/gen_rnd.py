import csv
import os.path
import sys

#Essentially just a way to define equality/uniqueness
class Colocalization:
    distance_cutoff = 250

    #def __init__(self, umi = "", amr_group = "", mge_header = "", distance = -1):
    def __init__(self, amr_group = "", mge_header = "", distance = -1):
        #self.umi = umi
        self.amr_group = amr_group
        self.mge_header = mge_header
        self.distance = distance

    def __str__(self):
        #return self.umi + ":::" + self.amr_group + ":::" + self.mge_header + ":::" + str(self.distance)
        return self.amr_group + ":::" + self.mge_header + ":::" + str(self.distance)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        #return self.umi == other.umi and \
        return self.amr_group == other.amr_group and \
               self.mge_header == other.mge_header and \
               abs(self.distance - other.distance) <= Colocalization.distance_cutoff

    def __lt__(self, other):
        #if self.umi < other.umi:
            #return True
        if self.amr_group < other.amr_group:
            return True
        elif self.mge_header < other.mge_header:
            return True
        else:
            return self.distance < other.distance

    def __le__(self, other):
        #if self.umi <= other.umi:
            #return True
        if self.amr_group <= other.amr_group:
            return True
        elif self.mge_header <= other.mge_header:
            return True
        else:
            return self.distance <= other.distance


    def __gt__(self, other):
        #if self.umi > other.umi:
            #return True
        if self.amr_group > other.amr_group:
            return True
        elif self.mge_header > other.mge_header:
            return True
        else:
            return self.distance > other.distance


    def __ge__(self, other):
        #if self.umi >= other.umi:
            #return True
        if self.amr_group >= other.amr_group:
            return True
        elif self.mge_header >= other.mge_header:
            return True
        else:
            return self.distance >= other.distance



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

    #Get read to umi dictionary
    #read_umi_dict = {}
    #read_umi_file = "/home/noyes046/jsettle/argmobrich/analysis/umi/read_umi_convert.tsv"
    #with open(read_umi_file, 'r') as tsv_file:
        #tsv_reader = csv.reader(tsv_file, delimiter='\t')
        #for row in tsv_reader:
            #read_umi_dict[row[0]] = row[1]

    #Go through colocalization dist sheet, and find unique colocalizations
    dists_filename = sys.argv[1]
    colocs = []
    coloc_counts = {}
    with open(dists_filename, 'r') as dists_csv:
        dists_reader = csv.reader(dists_csv, delimiter='\t')
        for row in dists_reader:
            #Skip rows without coloc
            if not row or row[0][0] != 'm':
                continue

            #umi = read_umi_dict[row[0]]
            amr_group =  megares_ontology[row[1]]["group"]
            mge_header = row[2]
            distance = int(row[3])

            #coloc = Colocalization(umi, amr_group, mge_header, distance)
            coloc = Colocalization(amr_group, mge_header, distance)
            ref_coloc = Colocalization()
            is_new_coloc = True
            for unique_coloc in colocs:
                if unique_coloc == coloc:
                    ref_coloc = unique_coloc
                    is_new_coloc = False 
                    break

            if is_new_coloc:
                colocs.append(coloc)
                coloc_counts[coloc] = 1
            else:
                coloc_counts[ref_coloc] += 1

    #sorted_colocs = sorted(colocs)

    #unique_colocs = {}
    #current_coloc = Colocalization()
    #for coloc in sorted_colocs:
        #if current_coloc == coloc:
            #unique_colocs[current_coloc] += 1

        #else:
            #current_coloc = coloc
            #unique_colocs[current_coloc] = 1
            #if current_coloc.amr_group == "TETO" and "97796" in current_coloc.mge_header:
                #print(str(current_coloc))

    #print("===================================================================")

    #sorted_unique_colocs = sorted(unique_colocs, key = lambda coloc: unique_colocs[coloc], reverse=True)
    sorted_unique_colocs = sorted(colocs, key = lambda coloc: coloc_counts[coloc], reverse=True)
    #sorted_unique_colocs = sorted(colocs, key = lambda coloc:coloc.umi)
    #Write tsv for group colocalizations
    #umis = set()
    #for coloc in sorted_unique_colocs:
        #umis.add(coloc.umi)

    name = os.path.splitext(os.path.basename(sys.argv[1]))[0]
    #with open(name + "_unique_colocs_umi.tsv", 'w') as coloc_tsv:
    with open(name + "_unique_colocs.tsv", 'w') as coloc_tsv:
        coloc_writer = csv.writer(coloc_tsv, delimiter='\t')
        if "MOCK" in name:
            coloc_writer.writerow([name[:8] + " colocalization diversity/richness"])
        else:
            coloc_writer.writerow([name[:3] + " colocalization diversity/richness"])
        coloc_writer.writerow([])
        coloc_writer.writerow(["Number of unique colocalizations:", len(sorted_unique_colocs)])
        coloc_writer.writerow([])
        #coloc_writer.writerow(["Number of UMIs:", len(umis)])
        #coloc_writer.writerow([])
        coloc_writer.writerow(["MEGARes group", "MGE header", "Distance (within " + str(Colocalization.distance_cutoff) + " nts)", "Occurences"])
        #coloc_writer.writerow(["UMI", "MEGARes group", "MGE header", "Occurences"])
        for coloc in sorted_unique_colocs:
            coloc_writer.writerow([coloc.amr_group, coloc.mge_header, coloc.distance, coloc_counts[coloc]])
            #coloc_writer.writerow([coloc.umi, coloc.amr_group, coloc.mge_header, coloc_counts[coloc]])

