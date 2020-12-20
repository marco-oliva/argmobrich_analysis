from Bio import SeqIO
import pandas as pd
import pysam

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.colors import unlabel_rgb, label_rgb

import csv
import sys



if __name__ == "__main__":

    #Get megares lengths for coverage
    megares_gene_lengths = {}
    megares_reference_fasta_filename = "/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
    for rec in SeqIO.parse(megares_reference_fasta_filename, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)


    #Read in metadata about plasmids and prophages in aclame
    plasmid_ref_df = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=0, header=1, usecols="A,F", index_col=0)
    plasmid_ref_df_2 = pd.read_excel("aclame_genes_plasmids_0.4.xls", sheet_name=1, header=0, usecols="A,F", index_col=0)
    plasmid_ref_df = plasmid_ref_df.append(plasmid_ref_df_2)
    prophage_ref_df = pd.read_excel("aclame_genes_prophages_0.4.xls", header=1, usecols="A,F", index_col=0)

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

            #Make sure each label of each category is represented so we can color nodes appropriately
            #class_samples[cl] = set()
            #mech_samples[mech] = set()
            #group_samples[group] = set()

            #FIll in our dict
            megares_ontology[row[0]] = { "class"        : cl,
                                         "mech"         : mech,
                                         "group"        : group
                                       }

    #Input is text file with "sam filename,tsv filename" on each line
    sam_files_bov = []
    tsv_files_bov = []
    sam_files_mock = []
    tsv_files_mock = []
    i = 0
    with open(sys.argv[1], 'r') as input_handle:
        for line in input_handle:
            line_params = line.split(',')
            if i > 2:
                sam_files_bov.append(pysam.AlignmentFile(line_params[0], 'r'))
                tsv_files_bov.append(line_params[1][:-1])
            else:
                sam_files_mock.append(pysam.AlignmentFile(line_params[0], 'r'))
                tsv_files_mock.append(line_params[1][:-1])
            i+=1

    #These two dictionaries will track stats on ARGs and their colocalizations in order to plot
    info_dict_bov = {}
    info_dict_mock = {}
    for sam_file in sam_files_bov:
        for read in sam_file.fetch():
            #Skip alignments whose MEGARes gene requires SNP confirmation or are not primary 
            if "RequiresSNPConfirmation" in read.reference_name:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
                
            #check coverage
            if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:

                #Get info from ontology
                cl = megares_ontology[read.reference_name]["class"]
                mech = megares_ontology[read.reference_name]["mech"]
                group = megares_ontology[read.reference_name]["group"]
        
                #Initialize mechanism dictionary. Will track richness, MGE colocalizations, and its MEGARes class
                if not mech in info_dict_bov:
                    info_dict_bov[mech] = {"richness" : 0, "colocs" : 0, "mges" : {}, "class" : cl}

                #Update mechanism richness
                info_dict_bov[mech]["richness"] += 1
                #filename = str(sam_file.filename, 'utf-8')

                #class_samples[cl].add(filename)
                #mech_samples[mech].add(filename)
                #group_samples[group].add(filename)

        sam_file.close()

    for sam_file in sam_files_mock:
        for read in sam_file.fetch():
            #Skip alignments whose MEGARes gene requires SNP confirmation or are not primary 
            if "RequiresSNPConfirmation" in read.reference_name:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
                
            #check coverage
            if (float(read.reference_length) / megares_gene_lengths[read.reference_name]) > 0.8:

                #Get info from ontology
                cl = megares_ontology[read.reference_name]["class"]
                mech = megares_ontology[read.reference_name]["mech"]
                group = megares_ontology[read.reference_name]["group"]
          
                #Initialize mechanism dictionary. Will track richness, MGE colocalizations, and its MEGARes class
                if not mech in info_dict_mock:
                    info_dict_mock[mech] = {"richness" : 0, "colocs" : 0, "mges" : {}, "class" : cl}

                info_dict_mock[mech]["richness"] += 1
                #filename = str(sam_file.filename, 'utf-8')

                #class_samples[cl].add(filename)
                #mech_samples[mech].add(filename)
                #group_samples[group].add(filename)

        sam_file.close()

 
        
    hierarchy_dict_bov = {} #Will represent the subset of the MEGARes hierarchy present. Used for labeling on chart
    total_richness_bov = 0 #Used as normalization factor later

    #Go through colocalization results, only looking the ARGs selected above
    for tsv_file in tsv_files_bov:
        with open(tsv_file, 'r') as results_tsv:
            results_reader = csv.reader(results_tsv, delimiter='\t')
            for row in results_reader:
                #Skip column names
                if row[0] == "SAMPLE_TYPE":
                    continue

                #Use ontology
                cl =  megares_ontology[row[3]]["class"]
                mech =  megares_ontology[row[3]]["mech"]
                group =  megares_ontology[row[3]]["group"]

                if mech not in info_dict_bov:
                    #TODO datasheet contains secondary megares alignments. 12/7: I believe this was fixed
                    continue

                #Set up hiearachy dict
                if not cl in hierarchy_dict_bov:
                    hierarchy_dict_bov[cl] = {}

                if not mech in hierarchy_dict_bov[cl]:
                    total_richness_bov += info_dict_bov[mech]["richness"]
                    hierarchy_dict_bov[cl][mech] = []

                if not group in hierarchy_dict_bov[cl][mech]:
                    hierarchy_dict_bov[cl][mech].append(group)

                #Count total colocalizations
                info_dict_bov[mech]["colocs"] += 1

                #Record MGE and increase the (ARG, MGE) colocalization richness
                mge_ref_name = row[11]
                if not mge_ref_name in info_dict_bov[mech]["mges"]:
                    info_dict_bov[mech]["mges"][mge_ref_name] = 0
                info_dict_bov[mech]["mges"][mge_ref_name] += 1

    hierarchy_dict_mock = {} #Will represent the subset of the MEGARes hierarchy present. Used for labeling on chart
    total_richness_mock = 0 #Used as normalization factor later
    #Go through colocalization results, only looking the ARGs selected above
    for tsv_file in tsv_files_mock:
        with open(tsv_file, 'r') as results_tsv:
            results_reader = csv.reader(results_tsv, delimiter='\t')
            for row in results_reader:
                #Skip column names
                if row[0] == "SAMPLE_TYPE":
                    continue

                #Use ontology
                cl =  megares_ontology[row[3]]["class"]
                mech =  megares_ontology[row[3]]["mech"]
                group =  megares_ontology[row[3]]["group"]

                if mech not in info_dict_mock:
                    #TODO datasheet contains secondary megares alignments. 12/7: again, I believe this was fixed
                    continue

                #Set up hiearachy dict
                if not cl in hierarchy_dict_mock:
                    hierarchy_dict_mock[cl] = {}

                if not mech in hierarchy_dict_mock[cl]:
                    total_richness_mock += info_dict_mock[mech]["richness"]
                    hierarchy_dict_mock[cl][mech] = []

                if not group in hierarchy_dict_mock[cl][mech]:
                    hierarchy_dict_mock[cl][mech].append(group)

                #Count total colocalizations
                info_dict_mock[mech]["colocs"] += 1

                #Record MGE and increase the (ARG, MGE) colocalization numbers
                mge_ref_name = row[11]
                if not mge_ref_name in info_dict_mock[mech]["mges"]:
                    info_dict_mock[mech]["mges"][mge_ref_name] = 0
                info_dict_mock[mech]["mges"][mge_ref_name] += 1

    #Setup colors. Want color of MGE bars to be similar color, different shade of the mech color
    mech_colors = ["rgb(255, 154, 162)", "rgb(255, 183, 178)", "rgb(255, 218, 193)", "rgb(226, 240, 203)", "rgb(181, 234, 215)",
                   "rgb(199, 206, 234)", "rgb(133, 222, 119)"]
    #mge_colors =["rgb(120, 162, 204)", "rgb(136, 174, 208)", "rgb(150, 185, 208)", "rgb(164, 195, 210)", "rgb(174, 203, 214)",
                 #"rgb(174, 203, 214)", "rgb(218, 227, 217)", "rgb(212, 181, 157)", "rgb(206, 156, 111)", "rgb(87, 45, 21)",
                 #"rgb(182, 77, 58)", "rgb(183, 198, 139)", "rgb(244, 240, 203)", "rgb(222, 210, 158)", "rgb(179, 165, 128)",
                 #"rgb(162, 149, 116)", "rgb(245, 198, 145)", "rgb(254, 216, 178)", "rgb(254, 233, 206)", "rgb(191, 191, 191)",
                 #"rgb(172, 172, 172)", "rgb(50, 50, 50)", "rgb(224, 187, 228)", "rgb(149, 125, 173)", "rgb(210, 145, 188)",
                 #"rgb(254, 200, 216)", "rgb(255, 223, 211)", "rgb(255, 154, 162)", "rgb(255, 183, 178)", "rgb(255, 218, 193)",
                 #"rgb(226, 240, 203)", "rgb(181, 234, 215)", "rgb(199, 206, 234)", "rgb(133, 222, 119)"]

    #Map names of mechs and mges to color for plotting
    mech_color_dict = {}
    mge_color_dict = {}

    mechs_to_plot = 0 #Determines number of rows in the graph
    for mech in info_dict_bov:
        #Normalize richnesses
        info_dict_bov[mech]["richness"] = info_dict_bov[mech]["richness"] / total_richness_bov

        #Only consider mechanisms that have a richness above cutoff value
        if info_dict_bov[mech]["richness"] > 0.05:
            mechs_to_plot += 1

            #Get a new mechanism color if needed
            if not mech in mech_color_dict:
                mech_color_dict[mech] = mech_colors.pop()

            for mge in info_dict_bov[mech]["mges"]:
                #Normalize MGE's colocalization richness 
                mge_rel_richness = info_dict_bov[mech]["mges"][mge] / float(info_dict_bov[mech]["colocs"])

                #Ignore MGE if it doesn't meet richness criterion
                if mge_rel_richness > 0.05:
                    if not mge in mge_color_dict:
                        mge_color_dict[mge] = mech_color_dict[mech]


    for mech in info_dict_mock:
        #Normalize richness
        info_dict_mock[mech]["richness"] = info_dict_mock[mech]["richness"] / total_richness_mock

        #Only consider mechanisms that have a richness above cutoff value
        if info_dict_mock[mech]["richness"] > 0.05:
            mechs_to_plot += 1

            #Get a new mechanism color if necessary
            if not mech in mech_color_dict:
                mech_color_dict[mech] = mech_colors.pop()

            for mge in info_dict_mock[mech]["mges"]:
                #Normalize MGE's colocalization richness
                mge_rel_richness = info_dict_mock[mech]["mges"][mge] / float(info_dict_mock[mech]["colocs"])

                #Ignore MGE if doesn't meet richness criterion
                if mge_rel_richness > 0.05:
                    if not mge in mge_color_dict:
                        mge_color_dict[mge] = mech_color_dict[mech]

    #specs = []
    #for x in range(mechs_to_plot):
        #if x == 3:
            #specs.append([{}, {}])
        #specs.append([{}, {"secondary_y":True}])

    #Use hierarchy dict to go class by class, making bar graphs for each mech and its colocalizations
    fig = make_subplots(rows=mechs_to_plot+1, cols=2, vertical_spacing=0.005, #specs = specs,
                        subplot_titles=("MEGARes Mechanism Relative Richness", "Colocalized MGE Relative Richness"))
    for i in fig["layout"]["annotations"]:
        i["font"]["size"] = 50
    #fig = make_subplots(rows=1, cols=5)

    row = 1
    for cl in hierarchy_dict_bov:
        #mech_vals = []
        #mech_labels = []

        for mech in hierarchy_dict_bov[cl]:
            #mech_vals.append(info_dict[mech]["richness"])
            #mech_labels.append(mech)
            #Skip mechanisms below cutoff
            if info_dict_bov[mech]["richness"] < 0.05:
                continue

            mge_vals = []
            mge_labels = []
            for mge in info_dict_bov[mech]["mges"]:
                mge_rel_richness = info_dict_bov[mech]["mges"][mge] / float(info_dict_bov[mech]["colocs"])
                #Skip MGEs below cutof
                if mge_rel_richness > 0.05:
                    mge_vals.append(info_dict_bov[mech]["mges"][mge] / float(info_dict_bov[mech]["colocs"]))
            
                    #Append relevant name of plasmid/prophage instead of awful database accession label
                    if mge in plasmid_ref_df.index:
                        mge_labels.append(str(plasmid_ref_df.loc[mge].iloc[0]))
                    if mge in prophage_ref_df.index:
                        mge_labels.append(str(prophage_ref_df.loc[mge].iloc[0]))

            #Try to display mechanisms with long names nicely by inserting line breaks
            mech_name = mech
            if mech_name.count(' ') > 2:
                idx1 = mech_name.find(' ')
                idx2 = mech_name.find(' ', idx1+1)
                mech_name = mech_name[:idx2] + "<br>" + mech_name[idx2:]

            if mech_name.count(' ') > 5:
                idx3 = mech_name.find(' ', idx2+1)
                idx4 = mech_name.find(' ', idx3+1)
                idx5 = mech_name.find(' ', idx4+1)
                mech_name = mech_name[:idx5] + "<br>" + mech_name[idx5:]

            #mge_bar_colors = []
            #for mge__ in mge_labels:
                #mge_bar_colors.append(mge_color_dict[mge__])

            #Done with this mech, add trace (bar graphs) for mech and mges
            fig.add_bar(y=[mech_name], x=[info_dict_bov[mech]["richness"]], marker_color=mech_color_dict[mech], orientation='h', showlegend=False, row=row, col=1)
            fig.add_bar(y=mge_labels, x=mge_vals, marker_color=mech_color_dict[mech], orientation='h', showlegend=False, row=row, col=2)#, secondary_y=True)


            #Configure the graphs
            fig.layout["xaxis" + str(2*row - 1)]["range"]=[0,1]
            fig.layout["xaxis" + str(2*row)]["range"]=[0,1]

            fig.layout["xaxis" + str(2*row - 1)]["tickfont"]["size"]=48
            fig.layout["xaxis" + str(2*row)]["tickfont"]["size"]=48
            fig.layout["yaxis" + str(2*row - 1)]["tickfont"]["size"]=48
            fig.layout["yaxis" + str(2*row)]["tickfont"]["size"]=38
            
            row += 1

            #If we are done, insert a title here. 
            if row > 3:
                fig.layout["xaxis" + str(2*row - 3)]["title"]="Mock"
                fig.layout["xaxis" + str(2*row - 3)]["title"]["font"]["size"]=56
                fig.layout["xaxis" + str(2*row - 2)]["title"]="Mock"
                fig.layout["xaxis" + str(2*row - 2)]["title"]["font"]["size"]=56


    row += 1
    for cl in hierarchy_dict_mock:
        #mech_vals = []
        #mech_labels = []

        for mech in hierarchy_dict_mock[cl]:
            #mech_vals.append(info_dict[mech]["richness"])
            #mech_labels.append(mech)
            if info_dict_mock[mech]["richness"] < 0.05:
                continue

            mge_vals = []
            mge_labels = []
            for mge in info_dict_mock[mech]["mges"]:
                mge_rel_richness = info_dict_mock[mech]["mges"][mge] / float(info_dict_mock[mech]["colocs"])
                if mge_rel_richness > 0.05:
                    mge_vals.append(info_dict_mock[mech]["mges"][mge] / float(info_dict_mock[mech]["colocs"]))

                    #Get appropriate name instead of ugly/meaningless aclame accession label
                    if mge in plasmid_ref_df.index:
                        mge_labels.append(str(plasmid_ref_df.loc[mge].iloc[0]))
                    if mge in prophage_ref_df.index:
                        mge_labels.append(str(prophage_ref_df.loc[mge].iloc[0]))

            #Try to deisplay mechanism names properly
            mech_name = mech
            if mech_name.count(' ') > 2:
                idx1 = mech_name.find(' ')
                idx2 = mech_name.find(' ', idx1+1)
                mech_name = mech_name[:idx2] + "<br>" + mech_name[idx2:]

            if mech_name.count(' ') > 5:
                idx3 = mech_name.find(' ', idx2+1)
                idx4 = mech_name.find(' ', idx3+1)
                idx5 = mech_name.find(' ', idx4+1)
                mech_name = mech_name[:idx5] + "<br>" + mech_name[idx5:]

            #mge_bar_colors = []
            #for mge__ in mge_labels:
                #mge_bar_colors.append(mge_color_dict[mge__])

            #print(mge_bar_colors)
            #Done with this mech, add trace for mech and mges
            fig.add_bar(y=[mech_name], x=[info_dict_mock[mech]["richness"]], marker_color=mech_color_dict[mech], orientation='h', showlegend=False, row=row, col=1)
            fig.add_bar(y=mge_labels, x=mge_vals, marker_color=mech_color_dict[mech], orientation='h', showlegend=False, row=row, col=2)#, secondary_y=True)

            fig.layout["xaxis" + str(2*row - 1)]["range"]=[0,1]
            fig.layout["xaxis" + str(2*row)]["range"]=[0,1]

            fig.layout["xaxis" + str(2*row - 1)]["tickfont"]["size"]=48
            fig.layout["xaxis" + str(2*row)]["tickfont"]["size"]=48
            fig.layout["yaxis" + str(2*row - 1)]["tickfont"]["size"]=48
            fig.layout["yaxis" + str(2*row)]["tickfont"]["size"]=30
            
            row += 1

            #If we are done, add title 
            if row > 8:
                fig.layout["xaxis" + str(2*row - 3)]["title"]="Bovine"
                fig.layout["xaxis" + str(2*row - 3)]["title"]["font"]["size"]=56
                fig.layout["xaxis" + str(2*row - 2)]["title"]="Bovine"
                fig.layout["xaxis" + str(2*row - 2)]["title"]["font"]["size"]=56
    
    for i in range(0, mechs_to_plot):
        if i == 3: #Skip row with no graphs
            continue
        fig.update_xaxes(row=i,col=1,tickmode="array", tickvals=[0.0, 0.2, 0.4, 0.8, 1.0], ticktext=["","","","",""])
        fig.update_xaxes(row=i,col=2,tickmode="array", tickvals=[0.0, 0.2, 0.4, 0.8, 1.0], ticktext=["","","","",""])

    fig.update_layout(width=4000,height=3600)
    fig.write_image("coloc_bar.png")
    fig.write_image("coloc_bar.svg") #Ilya wants svgs for his manipulation of figures

    sys.exit(0)
