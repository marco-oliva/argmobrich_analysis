from Bio import SeqIO
import igraph
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
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

    hierarchy_dict = {}
    class_samples = {}
    mech_samples = {}
    group_samples = {}

    #Create ontology dictionary from MEGARes ontology file
    megares_ontology = {}
    ontology_filename = "/home/noyes046/jsettle/argmobrich/MEGARESONTOLOGY.tsv"
    with open(ontology_filename, 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv, delimiter='\t')
        for row in ontology_reader:
            #Skip column names
            if row[0] == "header":
                continue


            cl =  row[1]
            mech = row[2]
            group = row[3]

            #Set up hiearachy dict. This will be our tree structure
            if not cl in hierarchy_dict:
                hierarchy_dict[cl] = {}

            if not mech in hierarchy_dict[cl]:
                hierarchy_dict[cl][mech] = []

            if not group in hierarchy_dict[cl][mech]:
                hierarchy_dict[cl][mech].append(group)

            #Make sure each label of each category is represented so we can color nodes appropriately
            class_samples[cl] = set()
            mech_samples[mech] = set()
            group_samples[group] = set()

            #FIll in our dict
            megares_ontology[row[0]] = { "class"        : cl,
                                         "mech"         : mech,
                                         "group"        : group
                                       }

    #Input is text file with list of sam files
    sam_files = []
    with open(sys.argv[1], 'r') as input_handle:
        for line in input_handle:
            sam_files.append(pysam.AlignmentFile(line[:-1], 'r'))

    for sam_file in sam_files:
        #Iterate through every read. Accumulate number of reads while recording read length
        for read in sam_file.fetch():
            if "RequiresSNPConfirmation" in read.reference_name:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
                
            #check coverage
            if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:
                cl = megares_ontology[read.reference_name]["class"]
                mech = megares_ontology[read.reference_name]["mech"]
                group = megares_ontology[read.reference_name]["group"]
            
                filename = str(sam_file.filename, 'utf-8')

                class_samples[cl].add(filename)
                mech_samples[mech].add(filename)
                group_samples[group].add(filename)

        sam_file.close()

    represented_hierarchy = {}
    for cl in class_samples:
        if len(class_samples[cl]) > 0:
            represented_hierarchy[cl] = {}
            for mech in hierarchy_dict[cl]:
                if len(mech_samples[mech]) > 0:
                    represented_hierarchy[cl][mech] = []

    #TODO don't hardcode names of files
    cat1 = "../datasheets/tmp/A01_megares_mapped_reads.sam"
    cat2 = "../datasheets/tmp/B01_megares_mapped_reads.sam"
    cat3 = "../datasheets/tmp/C01_megares_mapped_reads.sam"

    cat4 = "../datasheets/tmp/MOCK_D01_megares_mapped_reads.sam"
    cat5 = "../datasheets/tmp/MOCK_E01_megares_mapped_reads.sam"
    cat6 = "../datasheets/tmp/MOCK_F01_megares_mapped_reads.sam"

    #Categories for x axis are just samples
    bovine_sample_names = [cat1, cat2, cat3]
    mock_sample_names = [cat4, cat5, cat6]

    #Categories for y axis are mechanisms 
    mechtick_vals = []
    mechtick_text = []
    classtick_vals = []
    classtick_text = []
    left = 0
    right = 0
    mech_names = []
    for cl in represented_hierarchy:
        left = len(mech_names)
        right = left + len(represented_hierarchy[cl])
        classtick_pos = (right+left-1)/2.0
        classtick_vals.append(classtick_pos)
        classtick_text.append(cl)
        i = left+1
        for mech in represented_hierarchy[cl]:
            current_tick_pos = (left + i - 1)/2.0
            mechtick_vals.append(current_tick_pos)
            mechtick_text.append(mech)

            mech_names.append(mech)
            i += 1
            left += 1

    num_mechs = len(mech_names)

    #Just want 2 colors for presence/absence
    color_scheme = [[0.0, 'rgb(198, 198, 198)'], 
                    [0.5, 'rgb(198, 198, 198)'],
                    [0.5, 'rgb(194, 59, 34)'],
                    [1.0, 'rgb(194, 59, 34)']]

    label_heatmap_matrix = []
    i = 0
    for cl in represented_hierarchy:
        val = i / float(len(mech_names))
        for mech in represented_hierarchy[cl]:
            label_heatmap_matrix.append([val])
        i+=1

    j = 1
    bovine_presence_matrix = []
    for cl in represented_hierarchy:
        val = j / float(len(mech_names)+1.0)
        for mech in represented_hierarchy[cl]:
            bovine_presence_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
            for i in range(0, len(bovine_presence_matrix_row)):
                if bovine_sample_names[i] in mech_samples[mech]:
                    bovine_presence_matrix_row[i] = val
            bovine_presence_matrix.append(bovine_presence_matrix_row)
        j+=1

    j = 1
    mock_presence_matrix = []
    for cl in represented_hierarchy:
        val = j / float(len(mech_names)+1.0)
        for mech in represented_hierarchy[cl]:
            mock_presence_matrix_row = [0 for x in range(0, len(mock_sample_names))]
            for i in range(0, len(mock_presence_matrix_row)):
                if mock_sample_names[i] in mech_samples[mech]:
                    mock_presence_matrix_row[i] = val
            mock_presence_matrix.append(mock_presence_matrix_row)
        j+=1

    fig = make_subplots(rows=1, cols=7, specs = [[{}, {"colspan":3},None,None,{"colspan":3, "secondary_y":True},None,None]])
    custom_color = ["lightpink", "palegreen", "peachpuff", "powderblue", "plum", "lightsalmon", "burlywood", "tomato", "chocolate", "gold", "salmon", "orchid"]
    _custom_color = ["grey", "lightpink", "palegreen", "peachpuff", "powderblue", "plum", "lightsalmon", "burlywood", "tomato", "chocolate", "gold", "salmon", "orchid"]

    fig.add_trace(go.Heatmap(z = label_heatmap_matrix,
                             colorscale = custom_color,
                             showscale=False),
                             row=1, col=1)

    fig.add_trace(go.Heatmap(z = bovine_presence_matrix,
                             colorscale = _custom_color,
                             showscale=False,
                             xgap=5,
                             ygap=3),
                             row=1, col=2)
    
    fig.add_trace(go.Heatmap(z = mock_presence_matrix,
                             colorscale = _custom_color,
                             showscale=False,
                             xgap=5,
                             ygap=3),
                             secondary_y=True, row=1, col=5)

    fig.update_layout(title="Binary heatmap showing presence or absence (grey) <br> of MEGARes mechanisms",
                      titlefont_size=32,
                      title_x=0.5,
                      width=2000,
                      height=1400,
                      showlegend=True)

    fig.update_xaxes(row=1,col=1,tickmode="array", tickvals=[], ticktext=[])
    fig.update_yaxes(row=1,col=1,tickmode="array", tickvals=classtick_vals, ticktext=classtick_text, tickfont_size=36)


    fig.update_xaxes(row=1,col=2,tickmode="array", tickvals=[0,1,2], ticktext=["bovine_2kb", "bovine_5kb", "bovine_8kb"],tickfont_size=24)
    fig.update_yaxes(row=1,col=2,tickmode="array", tickvals=[], ticktext=[])

    fig.update_xaxes(row=1,col=5,tickmode="array", tickvals=[0,1,2], ticktext=["mock_2kb", "mock_5kb", "mock_8kb"],tickfont_size=24)
    fig.update_yaxes(row=1,col=5,tickmode="array", tickvals=[], ticktext=[], secondary_y=False)
    fig.update_yaxes(row=1,col=5,tickmode="array", tickvals=mechtick_vals, ticktext=mechtick_text, secondary_y=True, tickfont_size=24)
    fig.write_image("megares_heatmap.svg")

    sys.exit(0)
