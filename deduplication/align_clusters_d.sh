#!/bin/bash

d01_fasta=sequel-demultiplex.MOCK_D01.ccs

for i in {0..199}; do
    fasta_name="${d01_fasta}_${i}.fasta"
    psl_name="./psls/${d01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
