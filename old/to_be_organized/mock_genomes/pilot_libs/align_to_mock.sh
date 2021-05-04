#!/bin/sh

readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly mock_genomes=${argmobrich}/analysis/mock_genomes
readonly ccs_dir=${argmobrich}/ccs_fastqs

readonly b_sub=${mock_genomes}/bacillus_subtilis.fna
readonly c_neo=${mock_genomes}/cryptococcus_neoformans.fna
readonly e_fae=${mock_genomes}/enterococcus_faecalis.fna
readonly e_col=${mock_genomes}/escherichia_coli.fna
readonly l_fer=${mock_genomes}/lactobacillus_fermentum.fna
readonly l_mon=${mock_genomes}/listeria_monocytogenes.fna
readonly p_aer=${mock_genomes}/pseudomonas_aeruginosa.fna
readonly s_cer=${mock_genomes}/saccharomyces_cerevisiae.fna
readonly s_ent=${mock_genomes}/salmonella_enterica.fna
readonly s_aur=${mock_genomes}/staphylococcus_aureus.fna

readonly genome_fastas=(${b_sub} ${c_neo} ${e_fae} ${e_col} ${l_fer} ${l_mon} ${p_aer} ${s_cer} ${s_ent} ${s_aur})

readonly minimap2=${user_dir}/tools/minimap2/minimap2

for genome in ${genome_fastas[@]}; do
    ${minimap2} -ax map-pb ${genome} ${ccs_dir}/sequel-demultiplex.MOCK_D01.ccs.fastq > "${genome}_D01_align.sam"
    ${minimap2} -ax map-pb ${genome} ${ccs_dir}/sequel-demultiplex.MOCK_E01.ccs.fastq > "${genome}_E01_align.sam"
    ${minimap2} -ax map-pb ${genome} ${ccs_dir}/sequel-demultiplex.MOCK_F01.ccs.fastq > "${genome}_F01_align.sam"
done
