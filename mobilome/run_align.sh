#!/bin/sh

ACLAME="/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
ICEBERG="/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
PLASMIDS="/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

ccs_dir="/home/noyes046/jsettle/argmobrich/analysis/deduplication"
a01_pre="A01"
a01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.1896_A01.ccs.fastq"
b01_pre="B01"
b01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.1896_B01.ccs.fastq"
c01_pre="C01"
c01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.1896_C01.ccs.fastq"

d01_pre="MOCK_D01"
d01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_D01.ccs.fastq"
e01_pre="MOCK_E01"
e01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq"
f01_pre="MOCK_F01"
f01_fq="${ccs_dir}/deduplicated_sequel-demultiplex.MOCK_F01.ccs.fastq"

files=(${a01_pre} ${a01_fq} ${b01_pre} ${b01_fq} ${c01_pre} ${c01_fq} ${d01_pre} ${d01_fq} ${e01_pre} ${e01_fq} ${f01_pre} ${f01_fq})
indices=(0 2 4 6 8 10)

minimap2="/home/noyes046/shared/tools/minimap2/minimap2"
for i in ${indices[@]}; do
    j=$(($i+1))
    ${minimap2} -ax map-pb $ACLAME "${files[$j]}" > "${files[$i]}_aclame.sam"
    ${minimap2} -ax map-pb $ICEBERG "${files[$j]}" > "${files[$i]}_iceberg.sam"
    ${minimap2} -ax map-pb $PLASMIDS "${files[$j]}" > "${files[$i]}_plasmidfinder.sam"
done

