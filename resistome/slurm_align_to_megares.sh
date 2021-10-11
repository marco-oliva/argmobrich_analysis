#!/bin/bash
#SBATCH --job-name=amegares
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=288:00:00
#SBATCH --output=amegares_%A-%a.out    # Standard output and error log
#SBATCH --array=1-16                   # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026_assemblies/Assemblies_Fasta/resistome/sam_files/MEGARes/"
DATABASE="/blue/boucher/marco.oliva/data/MEGARes/V2/megares_full_database_v2.00.fasta"

PROFILER="/usr/bin/time --verbose"

module load python
module load samtools
module load minimap

file_list=(${READS_DIR}/*.fasta.gz)

# Working on i-th file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
FILE_NAME=${file_list[${FILE_NUM}]}
BASE_NAME=$(basename ${FILE_NAME} .fasta.gz)

MINIMAP_FLAGS="-ax map-pb"

echo "Working on ${FILE_NUM}: ${FILE_NAME}"
mkdir -p ${OUT_DIR_BASE}
if [ ! -f "${OUT_DIR_BASE}/${BASE_NAME}_ato_megaresv2.sam" ]
then
  minimap2 ${MINIMAP_FLAGS} ${DATABASE} ${FILE_NAME} > ${OUT_DIR_BASE}/${BASE_NAME}_ato_megaresv2.sam
fi