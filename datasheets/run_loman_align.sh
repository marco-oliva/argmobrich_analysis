#!/bin/sh

MEGARES="/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
ACLAME="/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
ICEBERG="/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
PLASMIDS="/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

argmobrich="/home/noyes046/jsettle/argmobrich"
excise_script="${argmobrich}/analysis/datasheets/excise_amr_aligns.py"
minimap2="/home/noyes046/shared/tools/minimap2/minimap2"

dir_pre="${argmobrich}/analysis/datasheets"
data_out="${dir_pre}/loman_data"
tmp_out="${dir_pre}/loman_tmp"


module load samtools

#Align to megares
${minimap2} -ax map-ont $MEGARES "$1" > "${tmp_out}/$2_megares_align.sam"

#Create SAM file of only mapped reads (a flag of 4 means unmapped, and -F excludes those flags)
samtools view -h -F 4 "${tmp_out}/$2_megares_align.sam" > "${tmp_out}/$2_megares_mapped.sam"
set -x
python ${excise_script} "${tmp_out}/$2_megares_mapped.sam" "$1" "${tmp_out}/$2_megares_unmapped_sequences.fasta" $MEGARES

#Align sequences identified in previous step to various MGE databases
${minimap2} -ax map-ont $ACLAME "${tmp_out}/$2_megares_unmapped_sequences.fasta" > "${tmp_out}/$2_unmapped_aclame.sam"
${minimap2} -ax map-ont $ICEBERG "${tmp_out}/$2_megares_unmapped_sequences.fasta" > "${tmp_out}/$2_unmapped_iceberg.sam"
${minimap2} -ax map-ont $PLASMIDS "${tmp_out}/$2_megares_unmapped_sequences.fasta" > "${tmp_out}/$2_unmapped_pf.sam"
