#!/bin/bash
#SBATCH --job-name=mobilome
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --output=%A-%a_mobilome.out    # Standard output and error log
#SBATCH --array=1-16                 # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta"
SCRIPT="python /blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/mobilome/gen_mobilome.py"
SAM_FILES_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/mobilome/sam_files"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/mobilome"

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
${PROFILER} ${SCRIPT} -o ${OUT_DIR_BASE}/${BASE_NAME} \
                      -a  ${SAM_FILES_BASE}/ACLAME/${BASE_NAME}_ato_ACLAME.sam \
                      -i ${SAM_FILES_BASE}/ICEBERG/${BASE_NAME}_ato_ICEBERG.sam \
                      -p ${SAM_FILES_BASE}/PLASMIDS/${BASE_NAME}_ato_PLASMIDS.sam