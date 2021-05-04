#!/bin/sh

readonly argmobrich=/home/noyes046/jsettle/argmobrich
readonly ccs_fastqs=${argmobrich}/ccs_fastqs
readonly colocalization=${argmobrich}/analysis/colocalization

readonly a01_tsv=${colocalization}/A01_unmapped_argmobrich_colocalization_distance.tsv
readonly a01_fasta=${ccs_fastqs}/sequel-demultiplex.1896_A01.ccs.fasta
readonly a01_out=${colocalization}/cargo/A01

readonly b01_tsv=${colocalization}/B01_unmapped_argmobrich_colocalization_distance.tsv
readonly b01_fasta=${ccs_fastqs}/sequel-demultiplex.1896_B01.ccs.fasta
readonly b01_out=${colocalization}/cargo/B01

readonly c01_tsv=${colocalization}/C01_unmapped_argmobrich_colocalization_distance.tsv
readonly c01_fasta=${ccs_fastqs}/sequel-demultiplex.1896_C01.ccs.fasta
readonly c01_out=${colocalization}/cargo/C01

readonly d01_tsv=${colocalization}/MOCK_D01_unmapped_argmobrich_colocalization_distance.tsv
readonly d01_fasta=${ccs_fastqs}/sequel-demultiplex.MOCK_D01.ccs.fasta
readonly d01_out=${colocalization}/cargo/MOCK_D01

readonly e01_tsv=${colocalization}/MOCK_E01_unmapped_argmobrich_colocalization_distance.tsv
readonly e01_fasta=${ccs_fastqs}/sequel-demultiplex.MOCK_E01.ccs.fasta
readonly e01_out=${colocalization}/cargo/MOCK_E01

readonly f01_tsv=${colocalization}/MOCK_F01_unmapped_argmobrich_colocalization_distance.tsv
readonly f01_fasta=${ccs_fastqs}/sequel-demultiplex.MOCK_F01.ccs.fasta
readonly f01_out=${colocalization}/cargo/MOCK_F01

readonly files=(${a01_tsv} ${a01_fasta} ${a01_out} ${b01_tsv} ${b01_fasta} ${b01_out} ${c01_tsv} ${c01_fasta} ${c01_out} ${d01_tsv} ${d01_fasta} ${d01_out} ${e01_tsv} ${e01_fasta} ${e01_out} ${f01_tsv} ${f01_fasta} ${f01_out})
readonly indices=(0 3 6 9 12 15)

for i in ${indices[@]}; do
    j=$(($i+1))
    k=$(($i+2))
    python make_cargo_fastas.py "${files[$i]}" "${files[$j]}" "${files[$k]}"
done
