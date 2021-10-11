#!/bin/bash
#SBATCH --job-name=fcoloc
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=288:00:00
#SBATCH --output=fcoloc_%A-%a.out    # Standard output and error log
#SBATCH --array=1-97                # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026/data_analysis/colocalizations/richness"
IN_BASE_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/data_analysis/colocalizations"
SCRIPT="/blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/colocalizations/colocalization_richness.py"

PROFILER="/usr/bin/time --verbose"

module load python
module load samtools
module load minimap

file_list=(${READS_DIR}/*.fastq.gz)

# Working on i-th file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
FILE_NAME=${file_list[${FILE_NUM}]}
BASE_NAME=$(basename ${FILE_NAME} .fastq.gz)

mkdir -p ${OUT_DIR_BASE}

echo "Working on ${FILE_NUM}: ${FILE_NAME}"

set -x

# Find colocalizations
python ${SCRIPT} -c "${IN_BASE_DIR}/${BASE_NAME}_colocalizations.csv" > "${OUT_DIR_BASE}/${BASE_NAME}_colocalizations_richness.csv"

set +x