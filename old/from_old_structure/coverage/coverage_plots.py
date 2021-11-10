#!/usr/bin/env python
# coding: utf-8

## Imports
import sys
anaconda_path = '/home/noyes046/moliva/miniconda3/envs/argmobrich_3.7/lib/python3.7/site-packages'
if anaconda_path not in sys.path:
    sys.path.insert(1, anaconda_path)

import json
import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.patches import ConnectionPatch

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


# Functions
def calcRegionBounds(bool_array):
    '''
    Returns the lower and upper bounds of contiguous regions.

    Parameters
    ==========
    bool_array	1-D Binary numpy array
    '''
    assert(bool_array.dtype == 'bool')
    idx = np.diff(np.r_[0, bool_array, 0]).nonzero()[0]
    assert(len(idx)%2 == 0)
    return np.reshape(idx, (-1,2))


# Json files
print('Open:\tBaits')
baits_json_file_handler = open('/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/coverage/hits/Baits_hits.json')
baits_json_file = json.load(baits_json_file_handler)

print('Open:\tIlya')
ilya_json_file_handler = open('/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/coverage/hits/Ilya_hits.json')
ilya_json_file = json.load(ilya_json_file_handler)

print('Open:\tLoman')
loman_json_file_handler = open('/panfs/roc/groups/11/noyes046/moliva/projects/argmobrich/analysis/coverage/hits/Loman_hits.json')
loman_json_file = json.load(loman_json_file_handler)


# Extract baits total coverage for each organism
baits_coverage = dict()
for organism in baits_json_file.keys():
    print('Baits:\t{}'.format(organism))
    for ref_seq in baits_json_file[organism]['MINIMAP2'].keys():
        if organism not in baits_coverage:
            baits_coverage[organism] = list()
        baits_coverage[organism].extend(baits_json_file[organism]['MINIMAP2'][ref_seq])
    baits_coverage[organism] = np.array(baits_coverage[organism])    
    
    ba = np.zeros(len(baits_coverage[organism]), dtype=bool)
    for i in np.nonzero(baits_coverage[organism]):
        ba[i] = 1
    baits_coverage[organism] = ba


# Extract Ilya total coverage for each organism
ilya_coverage = dict()
for organism in ilya_json_file.keys():
    print('Ilya:\t{}'.format(organism))
    for l in ilya_json_file[organism].keys():
        for ref_seq in ilya_json_file[organism][l].keys():
            if organism not in ilya_coverage:
                ilya_coverage[organism] = dict()
            if l not in ilya_coverage[organism]:
                ilya_coverage[organism][l] = list()
            ilya_coverage[organism][l].extend(ilya_json_file[organism][l][ref_seq])
        ilya_coverage[organism][l] = np.array(ilya_coverage[organism][l])


# Extract Loman total coverage for each organism
loman_coverage = dict()
for organism in loman_json_file.keys():
    print('Loman:\t{}'.format(organism))
    for acc in loman_json_file[organism].keys():
        for ref_seq in loman_json_file[organism][acc].keys():
            if organism not in loman_coverage:
                loman_coverage[organism] = dict()
            if acc not in loman_coverage[organism]:
                loman_coverage[organism][acc] = list()
            loman_coverage[organism][acc].extend(loman_json_file[organism][acc][ref_seq])
        loman_coverage[organism][acc] = np.array(loman_coverage[organism][acc])


############################################################
# Select organism and regions

selected_points = dict()
selected_points['Staphylococcus_aureus_complete_genome'] = [25000 + 3000, 475000 + 3000, 2050000 + 6000, 2730000 - 5000]
selected_points['Salmonella_enterica_complete_genome'] = [95000 + 1000, 1238000 - 4000, 1560000 - 8000, 3773000]
selected_points['Saccharomyces_cerevisiae_draft_genome'] = [520000, 3220000, 5960000, 11280000 - 10000]
selected_points['Pseudomonas_aeruginosa_complete_genome'] = [440000, 1185000, 2860000, 5280000]
selected_points['Listeria_monocytogenes_complete_genome'] = [180000, 1335000, 1868000, 2630000]
selected_points['Lactobacillus_fermentum_complete_genome'] = [100000, 410000, 1360000, 1890000]
selected_points['Escherichia_coli_complete_genome'] = [170000, 2205000, 3610000, 3800000]
selected_points['Enterococcus_faecalis_complete_genome'] = [155000, 685000, 1620000, 2170000]
selected_points['Cryptococcus_neoformans_draft_genome'] = [29050000]
selected_points['Bacillus_subtilis_complete_genome'] = [80000, 632500, 2900000, 3689000]

margins = dict()
margins['Staphylococcus_aureus_complete_genome'] = [5000, 5000, 10000, 10000] 
margins['Salmonella_enterica_complete_genome'] = [13000, 16000, 24000, 13000]
margins['Saccharomyces_cerevisiae_draft_genome'] = [10000, 20000, 20000, 60000]
margins['Pseudomonas_aeruginosa_complete_genome'] = [30000, 30000, 20000, 30000]
margins['Listeria_monocytogenes_complete_genome'] = [20000, 20000, 60000, 30000]
margins['Lactobacillus_fermentum_complete_genome'] = [15000, 15000, 15000, 15000]
margins['Escherichia_coli_complete_genome'] = [30000, 30000, 30000, 30000]
margins['Enterococcus_faecalis_complete_genome'] = [30000, 30000, 5000, 5000]
margins['Cryptococcus_neoformans_draft_genome'] = [100000]
margins['Bacillus_subtilis_complete_genome'] = [10000, 10000, 70000, 20000]

############################################################


for organism in baits_json_file.keys():
    print('Plotting: {}'.format(organism))
    bait_regions_bounds = calcRegionBounds(baits_coverage[organism])
    # Create main container with size of 20x10
    fig = plt.figure(figsize=(40, 10))
    plt.subplots_adjust(bottom = 0, left = 0, top = 1., right = 1)

    # Zooms, one for each bait region
    num_of_zooms = len(selected_points[organism])
    sub_plots = list()
    for idx,point in enumerate(selected_points[organism]):
        zoom_id = idx + 1
        sub_n = fig.add_subplot(2, num_of_zooms, zoom_id)
        sub_n.plot(loman_coverage[organism]['ERR3152366_1'], '-', color='#cfd8dc', label='Mock (GridION)')
        #sub_n.fill_between(range(0, len(loman_coverage[organism]['ERR3152366_1'])), 0, loman_coverage[organism]['ERR3152366_1'], facecolor='#cfd8dc')
        sub_n.plot(ilya_coverage[organism]['D'], '-', label='Mock 2kb (TELS)')
        sub_n.plot(ilya_coverage[organism]['E'], '-', label='Mock 5kb (TELS)')
        sub_n.plot(ilya_coverage[organism]['F'], '-', label='Mock 8kb (TELS)')
        sub_n.set_xlim(point - margins[organism][idx], point + margins[organism][idx])
        sub_n.set_ylim(bottom=1)
        sub_n.set_yscale('log')
        sub_plots.append(sub_n)
        for region_sub_it in bait_regions_bounds:
            sub_n.axvspan(region_sub_it[0], region_sub_it[1], facecolor='r', alpha=0.2) 

    # Create bottom axes, a combination of remaining cells
    sub_b = fig.add_subplot(2, num_of_zooms,(num_of_zooms + 1, 2*num_of_zooms))
    sub_b.plot(loman_coverage[organism]['ERR3152366_1'], '-', label='Mock (GridION)', color='#cfd8dc')
    #sub_b.fill_between(range(0, len(loman_coverage[organism]['ERR3152366_1'])), 0, loman_coverage[organism]['ERR3152366_1'], facecolor='#cfd8dc')
    sub_b.plot(ilya_coverage[organism]['D'], '-', label='Mock 2kb (TELS)')
    sub_b.plot(ilya_coverage[organism]['E'], '-', label='Mock 5kb (TELS)')
    sub_b.plot(ilya_coverage[organism]['F'], '-', label='Mock 8kb (TELS)')
    sub_b.set_yscale('log')
    sub_b.set_ylim(bottom=1)
    sub_b.legend(loc='upper right')
    for region in bait_regions_bounds:
        sub_b.axvspan(region[0], region[1], facecolor='r', alpha=0.2)    

    # Create Connection patches 
    for idx,point in enumerate(selected_points[organism]):
        xy_l = (point - margins[organism][idx], 1)
        con_pl = ConnectionPatch(xyA=xy_l, xyB=xy_l, coordsA=sub_plots[idx].transData, coordsB=sub_b.transData, color = 'black')
        xy_r = (point + margins[organism][idx], 1)
        con_pr = ConnectionPatch(xyA=xy_r, xyB=xy_r, coordsA=sub_plots[idx].transData, coordsB=sub_b.transData, color = 'black')
        fig.add_artist(con_pl)
        fig.add_artist(con_pr)

    plt.savefig('zoomed_plots/{}.svg'.format(organism))
    plt.savefig('zoomed_plots/{}.png'.format(organism))
