#!/bin/bash


a01_fasta=sequel-demultiplex.1896_A01.ccs

for i in {0..199}; do
    fasta_name="${a01_fasta}_${i}.fasta"
    psl_name="./psls/${a01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
