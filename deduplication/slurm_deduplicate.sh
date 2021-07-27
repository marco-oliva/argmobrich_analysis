#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=100gb
#SBATCH --time=288:00:00
#SBATCH --output=dedup_%A-%a.out    # Standard output and error log
#SBATCH --array=1-2                 # Array range

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
READS_DIR="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads"
SCRIPT="/blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/deduplication/find_duplicates.py"
OUT_DIR_BASE="/blue/boucher/marco.oliva/data/Noyes_Project_026/Reads/Deduplicated"
PROFILER="/usr/bin/time --verbose"

module load python
module load blat

#file_list=(${READS_DIR}/*.fastq.gz)
mapfile -t file_list < /blue/boucher/marco.oliva/projects/remote/argmobrich_analysis/deduplication/missed_files.txt


# Dedup on read file
FILE_NUM=$(( $SLURM_ARRAY_TASK_ID - 1 ))
filename=${file_list[${FILE_NUM}]}
echo "Working on ${FILE_NUM}: ${filename}"
base_name=$(basename ${filename} .fastq.gz)

out_dir="${OUT_DIR_BASE}/${base_name}"
mkdir -p ${out_dir}

python ${SCRIPT} ${filename} ${out_dir}
