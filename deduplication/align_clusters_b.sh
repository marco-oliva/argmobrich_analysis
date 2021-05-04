#!/bin/bash

b01_fasta=sequel-demultiplex.1896_B01.ccs

for i in {0..199}; do
    fasta_name="${b01_fasta}_${i}.fasta"
    psl_name="./psls/${b01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
