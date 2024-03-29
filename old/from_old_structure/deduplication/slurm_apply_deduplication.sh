#!/bin/bash
#SBATCH --job-name=adedup
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=144:00:00
#SBATCH --output=adedup_%A-%a.out    # Standard output and error log
#SBATCH --array=1-4                 # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads"
SCRIPT="/blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/deduplication/deduplicate.py"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated"
PROFILER="/usr/bin/time --verbose"

module load python

#file_list=(${READS_DIR}/*.fastq.gz)
mapfile -t file_list < /blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/deduplication/missed_files_done.txt

# Dedup on read file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
file_name=${file_list[${FILE_NUM}]}
base_name=$(basename ${file_name} .fastq.gz)

tsv_file="${OUT_DIR_BASE}/${base_name}/pls_files/duplicates.tsv"
if [ ! -f "/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated/deduplicated_${base_name}.fastq.gz" ]; then
    echo "Running on /blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated/deduplicated_${base_name}.fastq.gz"
    python ${SCRIPT} -r "${file_name}" -d ${tsv_file} > "/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated/deduplicated_${base_name}.fastq.gz"
fi