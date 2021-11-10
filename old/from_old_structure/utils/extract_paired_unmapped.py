#!/usr/bin/env python3

from src.common import *
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute resistome')
    parser.add_argument('-s', help='Alignment file', dest='sam_file', required=True)
    parser.add_argument('-n', help='Library name', dest='lib_name', required=True)
    parser.add_argument('-w', help='Work directory', dest='work_dir', required=True)
    args = parser.parse_args()

    root_logger = init_logger()

    # sort sam file into bam
    root_logger.info("Sorting sam file")
    command = 'samtools sort -o {wd}/{out}_sorted.bam {in_sam}'.format(wd=args.work_dir, out=args.lib_name, in_sam=args.sam_file)
    execute_command(command)

    # R1 unmapped, R2 mapped
    root_logger.info("Get R1 unmapped, R2 mapped")
    command = 'samtools view -u -f 4 -F 264 {wd}/{out}_sorted.bam'.format(out=args.lib_name, wd=args.work_dir)
    execute_command(command, out_file_path='{wd}/{out}.r1unmap.r2map.bam'.format(wd=args.work_dir, out=args.lib_name))

    # R1 mapped, R2 unmapped
    root_logger.info("Get R1 mapped, R2 unmapped")
    command = 'samtools view -u -f 8 -F 260 {wd}/{out}_sorted.bam'.format(out=args.lib_name, wd=args.work_dir)
    execute_command(command, out_file_path='{wd}/{out}.r1map.r2unmap.bam'.format(wd=args.work_dir, out=args.lib_name))

    # R1 and R2 unmapped
    root_logger.info("Get R1 unmapped, R2 unmapped")
    command = 'samtools view -u -f 12 -F 256 {wd}/{out}_sorted.bam'.format(out=args.lib_name, wd=args.work_dir)
    execute_command(command, out_file_path='{wd}/{out}.r1unmap.r2unmap.bam'.format(wd=args.work_dir, out=args.lib_name))

    # merge
    root_logger.info("Merge sam files")
    command = 'samtools merge -u {wd}/{out}.unmapped.bam ' \
              '{wd}/{out}.r1unmap.r2map.bam ' \
              '{wd}/{out}.r1map.r2unmap.bam ' \
              '{wd}/{out}.r1unmap.r2unmap.bam'.format(out=args.lib_name, wd=args.work_dir)
    execute_command(command)

    # sort merged file
    root_logger.info("Sort merged bam file")
    command = 'samtools sort -n {wd}/{out}.unmapped.bam -o {wd}/{out}.unmapped.sorted.bam'.format(wd=args.work_dir, out=args.lib_name)
    execute_command(command)

    # extract paired ends reads
    root_logger.info("Extract paired ends reads")
    command = 'samtools fastq {wd}/{out}.unmapped.sorted.bam ' \
              '-1 {wd}/{out}.unmapped_R1.fastq.gz ' \
              '-2 {wd}/{out}.unmapped_R2.fastq.gz'.format(wd=args.work_dir, out=args.lib_name)
    execute_command(command)


if __name__ == "__main__":
    main()








