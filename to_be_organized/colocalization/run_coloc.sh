#!/bin/sh
readonly data_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets/data

a01=${data_dir}/A01_unmapped_argmobrich_colocalization.tsv
b01=${data_dir}/B01_unmapped_argmobrich_colocalization.tsv
c01=${data_dir}/C01_unmapped_argmobrich_colocalization.tsv

d01=${data_dir}/MOCK_D01_unmapped_argmobrich_colocalization.tsv
e01=${data_dir}/MOCK_E01_unmapped_argmobrich_colocalization.tsv
f01=${data_dir}/MOCK_F01_unmapped_argmobrich_colocalization.tsv

readonly pilot_colocalization_results=(${a01} ${b01} ${c01} ${e01} ${d01} ${f01})

for result_tsv_file in ${pilot_colocalization_results[@]}; do
    python gen_coloc.py "${result_tsv_file}"
done
