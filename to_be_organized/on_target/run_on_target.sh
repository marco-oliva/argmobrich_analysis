#!/bin/sh

readonly tmp_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets/tmp

readonly a01=${tmp_dir}/A01_megares_mapped_reads.sam
readonly b01=${tmp_dir}/B01_megares_mapped_reads.sam
readonly c01=${tmp_dir}/C01_megares_mapped_reads.sam

readonly d01=${tmp_dir}/MOCK_D01_megares_mapped_reads.sam
readonly e01=${tmp_dir}/MOCK_E01_megares_mapped_reads.sam
readonly f01=${tmp_dir}/MOCK_F01_megares_mapped_reads.sam

readonly pilot_megares_sams=(${a01} ${b01} ${c01} ${d01} ${e01} ${f01})

for sam_file in ${pilot_megares_sams[@]}; do
    python gen_on_target.py ${sam_file}
done
