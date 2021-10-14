#!/usr/bin/env python3
import errno
import logging
import signal
import sys
import configparser
import os
import argparse
import subprocess


# ------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, out_file_path='', time_it=False, seconds=10000000):
    rootLogger = logging.getLogger()
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
        rootLogger.info("Executing: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, err = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        rootLogger.info("Error executing command line")
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        rootLogger.info("Command exceeded timeout")
        return None
    if output and out_file_path != '':
        rootLogger.info('Writing output to {}'.format(out_file_path))
        with open(out_file_path, "wb") as out_file:
            out_file.write(output)
    if err:
        err = err.decode("utf-8")
        rootLogger.info("\n" + err)
    return output


# mkdir -p
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise  # nop


def align_to_megares(config, fastq_file_path, out_dir):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    megares_path = config['DATABASE']['MEGARES']

    mkdir_p(out_dir)
    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_MEGARES']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=megares_path,
        i_file=fastq_file_path)
    execute_command(align_command, out_file_path=out_file)


def align_to_kegg(config, fastq_file_path, out_dir):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    kegg_path = config['DATABASE']['KEGG']

    mkdir_p(out_dir)
    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_KEGG']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=kegg_path,
        i_file=fastq_file_path)
    execute_command(align_command, out_file_path=out_file)
    return out_file


def align_to_mges(config, fastq_file_path, out_dir):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    plasmids_path = config['DATABASE']['PLASMIDS']
    aclame_path = config['DATABASE']['ACLAME']
    iceberg_path = config['DATABASE']['ICEBERG']

    mkdir_p(out_dir)

    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_PLASMIDS']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=plasmids_path,
        i_file=fastq_file_path)
    execute_command(align_command, out_file_path=out_file)

    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ACLAME']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=aclame_path,
        i_file=fastq_file_path)
    execute_command(align_command, out_file_path=out_file)

    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ICEBERG']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=iceberg_path,
        i_file=fastq_file_path)
    execute_command(align_command, out_file_path=out_file)


def gen_resistome(config, config_file_path, fastq_file_path, out_dir):
    gen_resistome_script = config['SCRIPTS']['GEN_RESISTOME']
    sam_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_MEGARES']
    out_file = out_dir + '/' + fastq_file_path

    gen_resistome_command = 'python {script} -s {sam_file} -o {out_name} -c {config_path}'.format(
        script=gen_resistome_script,
        sam_file=sam_file,
        out_name=out_file,
        config_path=config_file_path
    )
    execute_command(gen_resistome_command)


def gen_mobilome(config, config_file_path, fastq_file_path, out_dir):
    gen_mobilome_script = config['SCRIPTS']['GEN_MOBILOME']
    sam_file_plasmids = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_PLASMIDS']
    sam_file_aclame = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ACLAME']
    sam_file_iceberg = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ICEBERG']
    out_file = out_dir + '/' + fastq_file_path

    gen_mobilome_command = 'python {script} -p {sam_plasmids} -a {sam_aclame} -i {sam_iceberg} -o {out_name}  -c {config_path}'.format(
        script=gen_mobilome_script,
        sam_plasmids=sam_file_plasmids,
        sam_aclame=sam_file_aclame,
        sam_iceberg=sam_file_iceberg,
        out_name=out_file,
        config_path=config_file_path
    )
    execute_command(gen_mobilome_command)


def gen_colocalizations(config, config_file_path, fastq_file_path, out_dir):
    gen_colocalizations_script = config['SCRIPTS']['FIND_COLOCALIZATIONS']
    sam_file_plasmids = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_PLASMIDS']
    sam_file_aclame = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ACLAME']
    sam_file_iceberg = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_ICEBERG']
    sam_file_kegg = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_KEGG']
    sam_file_megares = out_dir + '/' + fastq_file_path + config['EXTENSION']['A_TO_MEGARES']
    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['COLOCALIZATIONS']
    skip_begin = 0
    skip_end = 0
    if "V2" in fastq_file_path:
        skip_begin = config['MISC']['V2_SKIP_BEGIN']
        skip_end = config['MISC']['V2_SKIP_END']

    find_colocalizations_command = 'python {script} -p {sam_plasmids} -a {sam_aclame} -i {sam_iceberg} -k {sam_kegg}' \
                                   ' -m {sam_megares} -r {reads} -e {skip_end} -b {skip_begin} -c {config_file}'.format(
        script=gen_colocalizations_script,
        sam_plasmids=sam_file_plasmids,
        sam_aclame=sam_file_aclame,
        sam_iceberg=sam_file_iceberg,
        sam_kegg=sam_file_kegg,
        sam_megares=sam_file_megares,
        skip_begin=skip_begin,
        reads=fastq_file_path,
        skip_end=skip_end,
        config_file=config_file_path
    )
    execute_command(find_colocalizations_command, out_file_path=out_file)


def gen_colocalizations_richness(config, config_file_path, fastq_file_path, out_dir):
    gen_colocalizations_richness_script = config['SCRIPTS']['COLOCALIZATIONS_RICHNESS']
    colocalizations_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['COLOCALIZATIONS']
    out_file = out_dir + '/' + fastq_file_path + config['EXTENSION']['COLOCALIZATIONS_RICHNESS']

    colocalizations_richness_command = 'python {script} -i {coloc_file} -c {config_path}'.format(
        script=gen_colocalizations_richness_script,
        coloc_file=colocalizations_file,
        config_path=config_file_path
    )
    execute_command(colocalizations_richness_command, out_file_path=out_file)


def main():
    parser = argparse.ArgumentParser(description='Colocalizations Pipeline')
    parser.add_argument('-i', help='Input file', dest='input_path', default='')
    parser.add_argument('-o', help='Output directory', dest='output_dir_path', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', default='./config.ini')
    args = parser.parse_args()

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    config = configparser.ConfigParser()
    config.read(args.config_path)

    # Single input file
    if config['PIPELINE_STEPS']['ALIGN_TO_MEGARES']:
        root_logger.info("Aligning to Megares")
        align_to_megares(config, args.input_path, args.output_dir_path)

    if config['PIPELINE_STEPS']['ALIGN_TO_MGES']:
        root_logger.info("Aligning to MGEs")
        align_to_mges(config, args.input_path, args.output_dir_path)

    if config['PIPELINE_STEPS']['ALIGN_TO_KEGG']:
        root_logger.info("Aligning to KEGG")
        align_to_kegg(config, args.input_path, args.output_dir_path)

    if config['PIPELINE_STEPS']['COMPUTE_RESISTOME']:
        root_logger.info("Generating Resistome")
        gen_resistome(config, args.config_path, args.input_path, args.output_dir_path)

    if config['PIPELINE_STEPS']['COMPUTE_MOBILOME']:
        root_logger.info("Generating Mobilome")
        gen_mobilome(config, args.config_path, args.input_path, args.output_dir_path)

    if config['PIPELINE_STEPS']['COMPUTE_COLOCALIZATIONS']:
        root_logger.info("Generating Colocalizations")
        gen_colocalizations(config, args.config_path, args.input_path, args.output_dir_path)
        gen_colocalizations_richness(config, args.config_path, args.input_path, args.output_dir_path)


if __name__ == "__main__":
    main()
