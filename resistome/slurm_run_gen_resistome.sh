#!/bin/bash
#SBATCH --job-name=resistome
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=288:00:00
#SBATCH --output=%A-%a_resistome.out    # Standard output and error log
#SBATCH --array=1-16                   # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta"
SAM_FILES_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/resistome/sam_files/MEGARes"
OUT_BASE_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/resistome"
SCRIPT="python /blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/resistome/gen_resistome.py"

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
${PROFILER} ${SCRIPT} -s ${SAM_FILES_BASE}/${BASE_NAME}_ato_megaresv2.sam \
                      -o ${OUT_BASE_DIR}/${BASE_NAME}
