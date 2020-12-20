from Bio import SeqIO
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pysam

import csv
import matplotlib.pyplot as plt
import os.path
import sys


if __name__ == "__main__":

    if(len(sys.argv) < 2):
        sys.exit(1)

    #Get aclame, iceberg, and plasmid finder lengths for coverage
    aclame_gene_lengths = {}
    #iceberg_gene_lengths = {}
    pf_gene_lengths = {}
    aclame_reference_fasta_filename = "/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
    #iceberg_reference_fasta_filename = "/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
    pf_reference_fasta_filename = "/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

    for rec in SeqIO.parse(aclame_reference_fasta_filename, "fasta"):
        aclame_gene_lengths[rec.name] = len(rec.seq)

    #for rec in SeqIO.parse(iceberg_reference_fasta_filename, "fasta"):
        #iceberg_gene_lengths[rec.name] = len(rec.seq)

    for rec in SeqIO.parse(pf_reference_fasta_filename, "fasta"):
        pf_gene_lengths[rec.name] = len(rec.seq)

    #Get informative names for aclame from excel
    plasmid_ref_df = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=0, header=1, usecols="A,F", index_col=0)
    plasmid_ref_df_2 = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=1, header=0, usecols="A,F", index_col=0)
    #print(plasmid_ref_df)
    plasmid_ref_df = plasmid_ref_df.append(plasmid_ref_df_2)
    #print(plasmid_ref_df)
    prophage_ref_df = pd.read_excel("aclame_genes_prophages_0.4.xls", header=1, usecols="A,F", index_col=0)


    #Input is text file with list of sam files
    sam_file_prefixes = []
    with open(sys.argv[1], 'r') as input_handle:
        for line in input_handle:
            sam_file_prefixes.append(line[:-1])

    aclame_samples = {}
    #iceberg_samples = {}
    pf_samples = {}

    for sam_file_prefix in sam_file_prefixes:
        aclame_sam_file = pysam.AlignmentFile(sam_file_prefix + "_aclame.sam", 'r')
        #iceberg_sam_file = pysam.AlignmentFile(sam_file_prefix + "_iceberg.sam", 'r')
        pf_sam_file = pysam.AlignmentFile(sam_file_prefix + "_pf.sam", 'r')
        #Iterate through every read. Accumulate number of reads while recording read length
        for read in aclame_sam_file.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            #check coverage
            if (float(read.reference_length) / aclame_gene_lengths[read.reference_name]) > 0.5:
                if not read.reference_name in aclame_samples:
                    aclame_samples[read.reference_name] = set()
                aclame_samples[read.reference_name].add(sam_file_prefix)

        aclame_sam_file.close()

        #Iterate through every read. Accumulate number of reads while recording read length
        #for read in iceberg_sam_file.fetch():
            #if read.is_unmapped or read.is_secondary or read.is_supplementary:
                #continue
            #check coverage
            #if (float(read.reference_length) / iceberg_gene_lengths[read.reference_name]) > 0.5:
                #if not read.reference_name in iceberg_samples:
                    #iceberg_samples[read.reference_name] = set()
                #iceberg_samples[read.reference_name].add(sam_file_prefix)

        #iceberg_sam_file.close()

        #print(iceberg_samples)

        #Iterate through every read. Accumulate number of reads while recording read length
        for read in pf_sam_file.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            #check coverage
            if (float(read.reference_length) / pf_gene_lengths[read.reference_name]) > 0.5:
                if not read.reference_name in pf_samples:
                    pf_samples[read.reference_name] = set()
                pf_samples[read.reference_name].add(sam_file_prefix)

        pf_sam_file.close()

    #TODO don't hardcode names of files
    cat1 = "../datasheets/tmp/A01_unmapped"
    cat2 = "../datasheets/tmp/B01_unmapped"
    cat3 = "../datasheets/tmp/C01_unmapped"

    cat4 = "../datasheets/tmp/MOCK_D01_unmapped"
    cat5 = "../datasheets/tmp/MOCK_E01_unmapped"
    cat6 = "../datasheets/tmp/MOCK_F01_unmapped"

    #Categories for x axis are just samples
    bovine_sample_names = [cat1, cat2, cat3]
    mock_sample_names = [cat4, cat5, cat6]

    aclame_tuples = []
    for aclame_ref in aclame_samples:
        if "proph" in aclame_ref:
            if aclame_ref in prophage_ref_df.index:
                aclame_tuples.append((aclame_ref, prophage_ref_df.loc[aclame_ref].iloc[0]))
            #else:
                #print(aclame_ref)
        else:
            if aclame_ref in plasmid_ref_df.index:
                aclame_tuples.append((aclame_ref, plasmid_ref_df.loc[aclame_ref].iloc[0]))
            #else:
                #print(aclame_ref)

    aclame_tuples.sort(key=lambda tuple: tuple[1])
    #Categories for y axis are refs 
    #aclame_names = sorted([ref_name for ref_name in aclame_samples])
    #iceberg_names = [ref_name for ref_name in iceberg_samples]
    pf_names = sorted([ref_name for ref_name in pf_samples])


    mgetick_text = []
    mgetick_vals = []
    label_matrix = []
    current_label = aclame_tuples[0][1]
    label_val = float(0.0)
    last_count = float(0.0)
    count = float(0.0)
    tickvals = []
    ticklabels = []
    for aclame_tuple in aclame_tuples:
        aclame_ref = aclame_tuple[0]
        aclame_label = aclame_tuple[1]

        if current_label != aclame_label:
            tickvals.append((last_count + count - 1.0) / 2.0)
            ticklabels.append(current_label)

            last_count = count
            label_val += 1.0
            current_label = aclame_label

        label_matrix.append([label_val])

        mgetick_text.append(aclame_ref)
        mgetick_vals.append(count)

        count += 1.0

    tickvals.append((last_count + count - 1.0) / 2.0)
    ticklabels.append(current_label)

    for ref_name in pf_samples:
        last_count = count
        count += 1.0
        tickvals.append((last_count + count - 1.0) / 2.0)
        ticklabels.append(ref_name)
        

    bovine_matrix = []
    mock_matrix = []
    #y_axis_labels = []
    #y_axis_label_vals = []
    for aclame_tuple in aclame_tuples:
        aclame_ref = aclame_tuple[0]
        aclame_label = aclame_tuple[1]

        #y_axis_labels.append(aclame_label)
        #y_axis_label_vals.append(len(y_axis_labels) - 1.0)

        aclame_bovine_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
        for i in range(0, len(aclame_bovine_matrix_row)):
            if bovine_sample_names[i] in aclame_samples[aclame_ref]:
                aclame_bovine_matrix_row[i] = 1.0
        bovine_matrix.append(aclame_bovine_matrix_row)

        aclame_mock_matrix_row = [0 for x in range(0, len(mock_sample_names))]
        for i in range(0, len(aclame_mock_matrix_row)):
            if mock_sample_names[i] in aclame_samples[aclame_ref]:
                aclame_mock_matrix_row[i] = 1.0
        mock_matrix.append(aclame_mock_matrix_row)

    #for row in aclame_bovine_matrix:
        #print(row)

    #iceberg_bovine_matrix = []
    #iceberg_mock_matrix = []
    #for iceberg_ref in iceberg_names:
        #iceberg_bovine_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
        #for i in range(0, len(iceberg_bovine_matrix_row)):
            #if bovine_sample_names[i] in iceberg_samples[iceberg_ref]:
                #iceberg_bovine_matrix_row[i] = 1.0
        #iceberg_bovine_matrix.append(iceberg_bovine_matrix_row)

        #iceberg_mock_matrix_row = [0 for x in range(0, len(mock_sample_names))]
        #for i in range(0, len(iceberg_mock_matrix_row)):
            #if mock_sample_names[i] in iceberg_samples[iceberg_ref]:
                #iceberg_mock_matrix_row[i] = 1.0
        #iceberg_mock_matrix.append(iceberg_mock_matrix_row)

    #bovine_matrix = []
    #mock_matrix = []
    for pf_ref in pf_names:
        #y_axis_labels.append(pf_ref)
        #y_axis_label_vals.append(len(y_axis_labels) - 1.0)
        pf_bovine_matrix_row = [0 for x in range(0, len(bovine_sample_names))]
        for i in range(0, len(pf_bovine_matrix_row)):
            if bovine_sample_names[i] in pf_samples[pf_ref]:
                pf_bovine_matrix_row[i] = 1.0
        bovine_matrix.append(pf_bovine_matrix_row)

        pf_mock_matrix_row = [0 for x in range(0, len(mock_sample_names))]
        for i in range(0, len(pf_mock_matrix_row)):
            if mock_sample_names[i] in pf_samples[pf_ref]:
                pf_mock_matrix_row[i] = 1.0
        mock_matrix.append(pf_mock_matrix_row)

        label_val += 1.0
        label_matrix.append([label_val])

    #print(label_matrix)
    label_matrix = [ [row[0] / float(label_matrix[-1][0])] for row in label_matrix]
    #print(label_matrix)

    #print(len(label_matrix), len(bovine_matrix), len(mock_matrix))

    #Just want 2 colors for presence/absence
    aclame_color_scheme = [[0.0, 'rgb(198, 198, 198)'], 
                           [0.5, 'rgb(198, 198, 198)'],
                           [0.5, 'rgb(194, 59, 34)'],
                           [1.0, 'rgb(194, 59, 34)']]

    fig = make_subplots(rows=1, cols=7, specs = [[{}, {"colspan":3},None,None,{"colspan":3, "secondary_y":True},None,None]])

    fig.add_trace(go.Heatmap(z = label_matrix,
                             colorscale = "spectral",
                             showscale=False),
                             row=1, col=1)


    fig.add_trace(go.Heatmap(z = bovine_matrix,
                             colorscale = aclame_color_scheme,
                             showscale=False,
                             xgap=5,
                             ygap=2),
                             row=1, col=2)

    fig.add_trace(go.Heatmap(z = mock_matrix,
                             colorscale = aclame_color_scheme,
                             showscale=False,
                             xgap=5,
                             ygap=2),
                             secondary_y=True, row=1, col=5)
    
    #fig.add_trace(go.Heatmap(z = iceberg_bovine_matrix,
                             #colorscale = color_scheme,
                             #showscale=False),
                             #row=2, col=1)
    
    #fig.add_trace(go.Heatmap(z = iceberg_mock_matrix,
                             #colorscale = color_scheme,
                             #showscale=False),
                             #row=2, col=2)
    
    #fig.add_trace(go.Heatmap(z = pf_bovine_matrix,
                             #colorscale = pf_color_scheme,
                             #showscale=False,
                             #xgap=10,
                             #ygap=2), row=4, col=1)
   
    #fig.add_trace(go.Heatmap(z = pf_mock_matrix,
                             #colorscale = none_pf_color_scheme,
                             #showscale=False,
                             #xgap=10,
                             #ygap=2), row=4, col=2)

    fig.update_layout(title="Binary heatmap showing presence (red) or <br> absence (grey) of mobile genetic elements",
                      titlefont_size=40,
                      title_x=0.5,
                      title_y=0.986,
                      title_yanchor="top",
                      width=3000,
                      height=2800)

    #fig.update_yaxes(title="Aclame", row=1, col=1)
    #fig.update_yaxes(title="Plasmid Finder", row=4, col=1)
    #fig.layout.yaxis.title.font.size = 50
    #fig.layout.yaxis3.title.font.size = 50

    #fig.update_xaxes(row=1,col=1,tickmode="array", tickvals=[], ticktext=[])
    #fig.update_yaxes(row=1,col=1,tickmode="array", tickvals=[], ticktext=[])
    #fig.update_xaxes(row=1,col=2,tickmode="array", tickvals=[], ticktext=[])

    fig.update_xaxes(row=1,col=1,tickmode="array", tickvals=[], ticktext=[])
    fig.update_yaxes(row=1,col=1,tickmode="array", tickvals=tickvals, ticktext=ticklabels, tickfont_size=32)
    fig.update_xaxes(row=1,col=2,tickmode="array", tickvals=[0,1,2], ticktext=["bovine_2kb", "bovine_5kb", "bovine_8kb"],tickfont_size=48)
    fig.update_yaxes(row=1,col=2,tickmode="array", tickvals=[], ticktext=[])
    fig.update_xaxes(row=1,col=5,tickmode="array", tickvals=[0,1,2], ticktext=["mock_2kb", "mock_5kb", "mock_8kb"],tickfont_size=48)
    fig.update_yaxes(row=1,col=5,tickmode="array", tickvals=[], ticktext=[], secondary_y=False)
    print(mgetick_vals)
    print(mgetick_text)
    fig.update_yaxes(row=1,col=5,tickmode="array", tickvals=mgetick_vals, ticktext=mgetick_text, secondary_y=True, tickfont_size=28)

    #fig.update_yaxes(row=1,col=5,tickmode="array", tickvals=[], ticktext=[])
    fig.write_image("mge_heatmap.png")
    fig.write_image("mge_heatmap.svg")

    sys.exit(0)
