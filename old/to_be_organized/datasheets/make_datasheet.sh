#!/bin/sh

#Database files
megares="/home/noyes046/shared/databases/megares_v1.01/megares_database_v1.01.fasta"
aclame="/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
iceberg="/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
plasmids="/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

#Define path where the deduplicated data is
argmobrich="/home/noyes046/jsettle/argmobrich"
datadir="${argmobrich}/analysis/deduplication"

#Each of these pairs is an output file prefix and the corresponding fastq data file
b01_pre="B01"
b01_fq="${datadir}/deduplicated_sequel-demultiplex.1896_B01.ccs.fastq"

c01_pre="C01"
c01_fq="${datadir}/deduplicated_sequel-demultiplex.1896_C01.ccs.fastq"

f01_pre="MOCK_F01"
f01_fq="${datadir}/deduplicated_sequel-demultiplex.MOCK_F01.ccs.fastq"

a01_pre="A01"
a01_fq="${datadir}/deduplicated_sequel-demultiplex.1896_A01.ccs.fastq"

d01_pre="MOCK_D01"
d01_fq="${datadir}/deduplicated_sequel-demultiplex.MOCK_D01.ccs.fastq"

e01_pre="MOCK_E01"
e01_fq="${datadir}/deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq"

#Array of inputs to python script
inputs=(${b01_pre} ${b01_fq} ${c01_pre} ${c01_fq} ${f01_pre} ${f01_fq} ${a01_pre} ${a01_fq} ${d01_pre} ${d01_fq} ${e01_pre} ${e01_fq})

#Indices to pick out pairs of inputs to give to python script, e.g. inputs 0 and 1 are used for B01 data
indices=(0 2 4 6 8 10)

#Define which aligner to use
minimap2="/home/noyes046/shared/tools/minimap2/minimap2"

#Prefix to this directory
dir_pre="${argmobrich}/analysis/datasheets"

#Script to take out portions of read that align to ARG
excise_script="${dir_pre}/excise_amr_aligns.py"

#Script to make colocalization tsvs
datasheet_script="${dir_pre}/gen_datasheet.py"

#Directories to dump final and intermidate data to, respectively
data_out="${dir_pre}/data"
tmp_out="${dir_pre}/tmp"

#Make sure samtools is available
module load samtools

#Loop over all indexes, one for each input fastq
#i represents prefix
#j represents fastq file name
for i in ${indices[@]}; do
    j=$(($i+1))

    #Align to megares
    ${minimap2} -ax map-pb $megares "${inputs[$j]}" > "${tmp_out}/${inputs[$i]}_megares_align.sam"

    #Create SAM file of only mapped reads (a flag of 4 means unmapped, and -F excludes those flags)
    samtools view -h -F 4 "${tmp_out}/${inputs[$i]}_megares_align.sam" > "${tmp_out}/${inputs[$i]}_megares_mapped_reads.sam"

    #TODO this keeps getting killed when run on queue: is there a way to lessen the burden on memory maybe?
    #Create fasta file of sequences to examine for colocalization
    python ${excise_script} "${tmp_out}/${inputs[$i]}_megares_mapped_reads.sam" "${inputs[$j]}" \
        "${tmp_out}/${inputs[$i]}_megares_unmapped_sequences.fasta" $megares

    #Align sequences identified in previous step to various MGE databases
    ${minimap2} -ax map-pb $aclame "${tmp_out}/${inputs[$i]}_megares_unmapped_sequences.fasta" > "${tmp_out}/${inputs[$i]}_unmapped_aclame.sam"
    ${minimap2} -ax map-pb $iceberg "${tmp_out}/${inputs[$i]}_megares_unmapped_sequences.fasta" > "${tmp_out}/${inputs[$i]}_unmapped_iceberg.sam"
    ${minimap2} -ax map-pb $plasmids "${tmp_out}/${inputs[$i]}_megares_unmapped_sequences.fasta" > "${tmp_out}/${inputs[$i]}_unmapped_pf.sam"

    #Create csv with MGEs
    python ${datasheet_script} "${tmp_out}/${inputs[$i]}_megares_mapped_reads.sam" "${tmp_out}/${inputs[$i]}_unmapped" "${data_out}"

    #Create csv without MGEs, kinda outdated because of resistome folder
    #python ${datasheet_script} "${tmp_out}/${inputs[$i]}_megares_mapped_reads.sam" "${data_out}/${inputs[$i]}" --no-mge
done

exit 0
