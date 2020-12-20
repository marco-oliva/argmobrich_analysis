#!/bin/bash

f01_fasta=sequel-demultiplex.MOCK_F01.ccs

for i in {0..199}; do
    fasta_name="${f01_fasta}_${i}.fasta"
    psl_name="./psls/${f01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
