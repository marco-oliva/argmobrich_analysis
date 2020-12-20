#!/bin/bash

c01_fasta=sequel-demultiplex.1896_C01.ccs

for i in {0..199}; do
    fasta_name="${c01_fasta}_${i}.fasta"
    psl_name="./psls/${c01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
