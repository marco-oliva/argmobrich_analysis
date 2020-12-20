#!/bin/sh

#Helpful shortcuts for directories
readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly loman=${argmobrich}/loman
readonly mock_genomes=${argmobrich}/analysis/mock_genomes

#All genome fastas
#readonly b_sub=${mock_genomes}/bacillus_subtilis.fna
readonly c_neo=${mock_genomes}/cryptococcus_neoformans.fna
readonly e_fae=${mock_genomes}/enterococcus_faecalis.fna
readonly e_col=${mock_genomes}/escherichia_coli.fna
readonly l_fer=${mock_genomes}/lactobacillus_fermentum.fna
readonly l_mon=${mock_genomes}/listeria_monocytogenes.fna
readonly p_aer=${mock_genomes}/pseudomonas_aeruginosa.fna
readonly s_cer=${mock_genomes}/saccharomyces_cerevisiae.fna
readonly s_ent=${mock_genomes}/salmonella_enterica.fna
readonly s_aur=${mock_genomes}/staphylococcus_aureus.fna

#readonly genome_fastas=(${b_sub} ${c_neo} ${e_fae} ${e_col} ${l_fer} ${l_mon} ${p_aer} ${s_cer} ${s_ent} ${s_aur})
readonly genome_fastas=(${c_neo} ${e_fae} ${e_col} ${l_fer} ${l_mon} ${p_aer} ${s_cer} ${s_ent} ${s_aur})

#Relevant parameters for Loman files
readonly err3152366_file=${loman}/ERR3152366_1
readonly err3152366_parts=4
readonly err3152366=(${err3152366_file} ${err3152366_parts})

readonly err3152367_file=${loman}/ERR3152367_1
readonly err3152367_parts=30
readonly err3152367=(${err3152367_file} ${err3152367_parts})

readonly loman_samples=(${err3152366[@]} ${err3152367[@]})
readonly indices=(0 2)

readonly minimap2=${user_dir}/tools/minimap2/minimap2


#align all genomes with all loman files
for genome in ${genome_fastas[@]}; do
    for i in ${indices[@]}; do
        declare -i j=$(($i+1))
        declare -i part_num=1
        while [ ${part_num} -le ${loman_samples[$j]} ]
        do
            if [ ${part_num} -lt 10 ]
            then
                #echo "$(basename -- ${genome})_$(basename -- ${loman_samples[$i]}).part-0${part_num}_align.sam"
                ${minimap2} -ax map-ont ${genome} "${loman_samples[$i]}.part-0${part_num}.fastq" > "$(basename -- ${genome})_$(basename -- ${loman_samples[$i]}).part-0${part_num}_align.sam"
            else
                #echo "$(basename -- ${genome})_$(basename -- ${loman_samples[$i]}).part-${part_num}_align.sam"
                ${minimap2} -ax map-ont ${genome} "${loman_samples[$i]}.part-${part_num}.fastq"  > "$(basename -- ${genome})_$(basename -- ${loman_samples[$i]}).part-${part_num}_align.sam"
            fi
            ((part_num=part_num+1))
        done
    done
done
