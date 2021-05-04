ACLAME="/home/noyes046/shared/databases/aclame/aclame_genes_all_0.4.fasta"
ICEBERG="/home/noyes046/shared/databases/ice_berg/ICEberg_seq.fasta"
PLASMIDS="/home/noyes046/shared/databases/plasmid_finder/plasmids_combined.fsa"

#Align to mobilome
minimap2="/home/noyes046/shared/tools/minimap2/minimap2"

indices=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
for i in "${indices[@]}"; do
    ${minimap2} -ax map-ont $ACLAME "/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-${i}.fastq" > "ERR3152367_1.part-${i}_aclame.sam"
    ${minimap2} -ax map-ont $ICEBERG "/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-${i}.fastq" > "ERR3152367_1.part-${i}_iceberg.sam"
    ${minimap2} -ax map-ont $PLASMIDS "/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-${i}.fastq" > "ERR3152367_1.part-${i}_plasmidfinder.sam"
done
