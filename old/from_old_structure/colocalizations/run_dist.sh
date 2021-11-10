#!/bin/sh

readonly data_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets/data
readonly tmp_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets/tmp

a01=${data_dir}/A01_unmapped_argmobrich_colocalization.tsv
a01_megares=${tmp_dir}/A01_megares_mapped_reads.sam
a01_aclame=${tmp_dir}/A01_unmapped_aclame.sam
a01_iceberg=${tmp_dir}/A01_unmapped_iceberg.sam
a01_pf=${tmp_dir}/A01_unmapped_pf.sam

b01=${data_dir}/B01_unmapped_argmobrich_colocalization.tsv
b01_megares=${tmp_dir}/B01_megares_mapped_reads.sam
b01_aclame=${tmp_dir}/B01_unmapped_aclame.sam
b01_iceberg=${tmp_dir}/B01_unmapped_iceberg.sam
b01_pf=${tmp_dir}/B01_unmapped_pf.sam

c01=${data_dir}/C01_unmapped_argmobrich_colocalization.tsv
c01_megares=${tmp_dir}/C01_megares_mapped_reads.sam
c01_aclame=${tmp_dir}/C01_unmapped_aclame.sam
c01_iceberg=${tmp_dir}/C01_unmapped_iceberg.sam
c01_pf=${tmp_dir}/C01_unmapped_pf.sam

d01=${data_dir}/MOCK_D01_unmapped_argmobrich_colocalization.tsv
d01_megares=${tmp_dir}/MOCK_D01_megares_mapped_reads.sam
d01_aclame=${tmp_dir}/MOCK_D01_unmapped_aclame.sam
d01_iceberg=${tmp_dir}/MOCK_D01_unmapped_iceberg.sam
d01_pf=${tmp_dir}/MOCK_D01_unmapped_pf.sam

e01=${data_dir}/MOCK_E01_unmapped_argmobrich_colocalization.tsv
e01_megares=${tmp_dir}/MOCK_E01_megares_mapped_reads.sam
e01_aclame=${tmp_dir}/MOCK_E01_unmapped_aclame.sam
e01_iceberg=${tmp_dir}/MOCK_E01_unmapped_iceberg.sam
e01_pf=${tmp_dir}/MOCK_E01_unmapped_pf.sam

f01=${data_dir}/MOCK_F01_unmapped_argmobrich_colocalization.tsv
f01_megares=${tmp_dir}/MOCK_F01_megares_mapped_reads.sam
f01_aclame=${tmp_dir}/MOCK_F01_unmapped_aclame.sam
f01_iceberg=${tmp_dir}/MOCK_F01_unmapped_iceberg.sam
f01_pf=${tmp_dir}/MOCK_F01_unmapped_pf.sam

readonly files=(${a01} ${a01_megares} ${a01_aclame} ${a01_iceberg} ${a01_pf} ${b01} ${b01_megares} ${b01_aclame} ${b01_iceberg} ${b01_pf} ${c01} ${c01_megares} ${c01_aclame} ${c01_iceberg} ${c01_pf} ${d01} ${d01_megares} ${d01_aclame} ${d01_iceberg} ${d01_pf} ${e01} ${e01_megares} ${e01_aclame} ${e01_iceberg} ${e01_pf} ${f01} ${f01_megares} ${f01_aclame} ${f01_iceberg} ${f01_pf})
#readonly indices=(0 3 6 9 12 15)
readonly indices=(0 5 10 15 20 25)

for i in ${indices[@]}; do
    j=$(($i+1))
    k=$(($i+2))
    l=$(($i+3))
    m=$(($i+4))
    echo python gen_dist.py "${files[$i]}" "${files[$j]}" "${files[$k]}" "${files[$l]}" "${files[$m]}"
    python gen_dist.py "${files[$i]}" "${files[$j]}" "${files[$k]}" "${files[$l]}" "${files[$m]}"
    #python gen_dist.py "${files[$i]}" "${files[$j]}" "${files[$k]}"
done
