#!/bin/sh

readonly colocalization=/home/noyes046/jsettle/argmobrich/analysis/colocalization

a01=${colocalization}/A01_unmapped_argmobrich_colocalization_distance.tsv
b01=${colocalization}/B01_unmapped_argmobrich_colocalization_distance.tsv
c01=${colocalization}/C01_unmapped_argmobrich_colocalization_distance.tsv

d01=${colocalization}/MOCK_D01_unmapped_argmobrich_colocalization_distance.tsv
e01=${colocalization}/MOCK_E01_unmapped_argmobrich_colocalization_distance.tsv
f01=${colocalization}/MOCK_F01_unmapped_argmobrich_colocalization_distance.tsv

readonly pilot_colocalization_distances=(${a01} ${b01} ${c01} ${d01} ${e01} ${f01})

for distance_tsv_file in ${pilot_colocalization_distances[@]}; do
    python gen_rnd.py "${distance_tsv_file}"
done
