#!/bin/bash
#SBATCH --job-name=amges
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=288:00:00
#SBATCH --output=amges_%A-%a.out    # Standard output and error log
#SBATCH --array=1-16                # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/mobilome/sam_files"

MGES_BASE="/blue/boucher/marco.oliva/data/MGE_DBs"
KEGG_DB="${MGES_BASE}/KEGG/kegg_prokaryotes.fasta"
ACLAME_DB="${MGES_BASE}/ACLAME/aclame_genes_all_0.4.fasta"
ICEBERG_DB="${MGES_BASE}/ICEBERG/ICEberg_seq.fasta"
PLASMIDS_DB="${MGES_BASE}/PLASMIDS/plasmids_combined.fsa"

PROFILER="/usr/bin/time --verbose"

module load python
module load samtools
module load minimap

file_list=(${READS_DIR}/*.fasta.gz)

# Working on i-th file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
FILE_NAME=${file_list[${FILE_NUM}]}
BASE_NAME=$(basename ${FILE_NAME} .fasta.gz)

echo "Working on ${FILE_NUM}: ${FILE_NAME}"

MINIMAP_FLAGS="-ax map-pb"
mkdir -p "${OUT_DIR_BASE}/KEGG"
mkdir -p "${OUT_DIR_BASE}/ACLAME"
mkdir -p "${OUT_DIR_BASE}/ICEBERG"
mkdir -p "${OUT_DIR_BASE}/PLASMIDS"

# KEGG
if [ ! -f "${OUT_DIR_BASE}/KEGG/${BASE_NAME}_ato_KEGG.sam" ]
then
  minimap2 ${MINIMAP_FLAGS} ${KEGG_DB} ${FILE_NAME} > ${OUT_DIR_BASE}/KEGG/${BASE_NAME}_ato_KEGG.sam
fi

# ACLAME
if [ ! -f "${OUT_DIR_BASE}/ACLAME/${BASE_NAME}_ato_ACLAME.sam" ]
then
  minimap2 ${MINIMAP_FLAGS} ${ACLAME_DB} ${FILE_NAME} > ${OUT_DIR_BASE}/ACLAME/${BASE_NAME}_ato_ACLAME.sam
fi

# ICEBERG
if [ ! -f "${OUT_DIR_BASE}/ICEBERG/${BASE_NAME}_ato_ICEBERG.sam" ]
then
  minimap2 ${MINIMAP_FLAGS} ${ICEBERG_DB} ${FILE_NAME} > ${OUT_DIR_BASE}/ICEBERG/${BASE_NAME}_ato_ICEBERG.sam
fi

# PLASMIDS
if [ ! -f "${OUT_DIR_BASE}/PLASMIDS/${BASE_NAME}_ato_PLASMIDS.sam" ]
then
  minimap2 ${MINIMAP_FLAGS} ${PLASMIDS_DB} ${FILE_NAME} > ${OUT_DIR_BASE}/PLASMIDS/${BASE_NAME}_ato_PLASMIDS.sam
fi