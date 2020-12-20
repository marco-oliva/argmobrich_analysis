#!/bin/sh

readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly ccs_dir=${argmobrich}/analysis/deduplication

readonly a01=${ccs_dir}/deduplicated_sequel-demultiplex.1896_A01.ccs.fastq
readonly b01=${ccs_dir}/deduplicated_sequel-demultiplex.1896_B01.ccs.fastq
readonly c01=${ccs_dir}/deduplicated_sequel-demultiplex.1896_C01.ccs.fastq

readonly d01=${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_D01.ccs.fastq
readonly e01=${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq
readonly f01=${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_F01.ccs.fastq

readonly pilot_samples=(${a01} ${b01} ${c01} ${d01} ${e01} ${f01})

for sample_file in ${pilot_samples[@]}; do
    python gen_summary.py ${sample_file}
done
