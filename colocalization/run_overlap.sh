#!/bin/sh

readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly bovine=${argmobrich}/bovine_feces
readonly mock=${argmobrich}/mock

a01_00=${bovine}/A01_00_unmapped_argmobrich_colocalization.csv
b01_00=${bovine}/B01_00_unmapped_argmobrich_colocalization.csv
b01_11=${bovine}/B01_11_unmapped_argmobrich_colocalization.csv

a01_11=${mock}/A01_11_unmapped_argmobrich_colocalization.csv
a01_22=${mock}/A01_22_unmapped_argmobrich_colocalization.csv
b01_22=${mock}/B01_22_unmapped_argmobrich_colocalization.csv

readonly pilot_colocalization_results=(${a01_00} ${b01_00} ${b01_11} ${a01_11} ${a01_22} ${b01_22})

for result_csv_file in ${pilot_colocalization_results[@]}; do
    python gen_overlap.py "${result_csv_file}"
done
