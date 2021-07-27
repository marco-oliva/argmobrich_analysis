#!/bin/bash
#SBATCH --job-name=fcoloc
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=288:00:00
#SBATCH --output=fcoloc_%A-%a.out    # Standard output and error log
#SBATCH --array=1-97                # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated"
MGE_ALIGNMENTS_BASE_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/data_analysis/mobilome/sam_files"
ARG_ALIGNMENTS_BASE_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/data_analysis/resistome/sam_files"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026/data_analysis/colocalizations"
SCRIPT="/blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/colocalizations/find_colocalizations.py"

MGES_BASE="/blue/boucher/marco.oliva/data/MGE_DBs"
KEGG_DB="${MGES_BASE}/KEGG/kegg_prokaryotes.fasta"
ACLAME_DB="${MGES_BASE}/ACLAME/aclame_genes_all_0.4.fasta"
ICEBERG_DB="${MGES_BASE}/ICEBERG/ICEberg_seq.fasta"
PLASMIDS_DB="${MGES_BASE}/PLASMIDS/plasmids_combined.fsa"

PROFILER="/usr/bin/time --verbose"

module load python
module load samtools
module load minimap

file_list=(${READS_DIR}/*.fastq.gz)

# Working on i-th file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
FILE_NAME=${file_list[${FILE_NUM}]}
BASE_NAME=$(basename ${FILE_NAME} .fastq.gz)

if [[ $BASE_NAME =~ "V2" ]]
then
  SKIP_BEGIN=3
  SKIP_END=66
else
  SKIP_BEGIN=0
  SKIP_END=0
fi


MEGARES_ALN_FILE="${ARG_ALIGNMENTS_BASE_DIR}/MEGARes/${BASE_NAME}_ato_megaresv2.sam"
ACLAME_ALN_FILE="${MGE_ALIGNMENTS_BASE_DIR}/ACLAME/${BASE_NAME}_ato_ACLAME.sam"
ICEBERG_ALN_FILE="${MGE_ALIGNMENTS_BASE_DIR}/ICEBERG/${BASE_NAME}_ato_ICEBERG.sam"
PLASMIDS_ALN_FILE="${MGE_ALIGNMENTS_BASE_DIR}/PLASMIDS/${BASE_NAME}_ato_PLASMIDS.sam"
KEGG_ALN_FILE="${MGE_ALIGNMENTS_BASE_DIR}/KEGG/${BASE_NAME}_ato_KEGG.sam"


echo "Working on ${FILE_NUM}: ${FILE_NAME}"

set -x

# Find colocalizations
python ${SCRIPT} -k ${KEGG_ALN_FILE} \
                 -p ${PLASMIDS_ALN_FILE} \
                 -i ${ICEBERG_ALN_FILE} \
                 -a ${ACLAME_ALN_FILE} \
                 -m ${MEGARES_ALN_FILE} \
                 -r ${FILE_NAME} > ${OUT_DIR_BASE}/${BASE_NAME}_colocalizations.csv \
                 -e ${SKIP_END} \
                 -b ${SKIP_BEGIN}

set +x