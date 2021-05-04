#!/bin/bash

e01_fasta=sequel-demultiplex.MOCK_E01.ccs

for i in {0..199}; do
    fasta_name="${e01_fasta}_${i}.fasta"
    psl_name="./psls/${e01_fasta}_${i}.psl"
    blat ${fasta_name} ${fasta_name} ${psl_name}
done
