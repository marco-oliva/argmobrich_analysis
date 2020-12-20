#!/bin/sh

readonly ccs_dir=/home/noyes046/jsettle/argmobrich/ccs_fastqs

d01_pre=D01
d01_fq="${ccs_dir}/sequel-demultiplex.MOCK_D01.ccs.fastq"

e01_pre=E01
e01_fq="${ccs_dir}/sequel-demultiplex.MOCK_E01.ccs.fastq"

f01_pre=F01
f01_fq="${ccs_dir}/sequel-demultiplex.MOCK_F01.ccs.fastq"


readonly inputs=(${d01_fq} ${d01_pre} ${e01_fq} ${e01_pre} ${f01_fq} ${f01_pre})
readonly indices=(0 2 4)

for i in ${indices[@]}; do
    j=$(($i+1))
    python gen_align_info.py "${inputs[$i]}" "${inputs[$j]}"
done
