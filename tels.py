#!/usr/bin/env python3
from Bio import SeqIO
import pysam
import statistics

from src.common import *

def find_duplicates(config, TELS_statistcs):
    import shutil
    if shutil.which('blat') is None:
        logging.getLogger().info("blat missing, needed for deduplication")
        exit(1)
    tmp_dir = config['OUTPUT']['OUT_DIR'] + '/tmp_files'
    mkdir_p(tmp_dir)
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['DUPLICATES']
    command = 'python {script} -r {in_file} -o {out_dir} -t {threads} -n {clusters}'.format(
        script=config['SCRIPTS']['FIND_DUPLICATES'],
        in_file=config['INPUT']['INPUT_FILE'],
        out_dir=tmp_dir,
        clusters=config['MISC']['DEDUP_CLUSTERS'],
        threads=config['MISC']['HELPER_THREADS'])
    execute_command(command, out_file_path=out_file)

def deduplicate(config, TELS_statistcs):
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['DEDUPLICATED']
    command = 'python {script} -d {duplicates_csv} -r {in_file}'.format(
        script=config['SCRIPTS']['DEDUPLICATE'],
        duplicates_csv=config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['DUPLICATES'],
        in_file=config['INPUT']['INPUT_FILE'])
    execute_command(command, out_file_path=out_file)
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(out_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(out_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    TELS_statistcs['READS_AFTER_DEDUPLICATION'] = sum(1 for record in SeqIO.parse(open(out_file, 'r'), read_file_type(out_file)))
    TELS_statistcs['READS_AFTER_DEDUPLICATION_PERC'] = (float(TELS_statistcs['READS_AFTER_DEDUPLICATION']) / float(TELS_statistcs['READS_BEFORE_DEDUPLICATION'])) * 100


def align_to_megares(config, TELS_statistcs):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    megares_path = config['DATABASE']['MEGARES']
    megares_gene_lengths = dict() # This should not be here
    for rec in SeqIO.parse(megares_path, "fasta"):
        megares_gene_lengths[rec.name] = len(rec.seq)

    mkdir_p(config['OUTPUT']['OUT_DIR'])
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_MEGARES']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=megares_path,
        i_file=config['INPUT']['INPUT_FILE'])
    execute_command(align_command, out_file_path=out_file)


def align_to_kegg(config, TELS_statistcs):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    kegg_path = config['DATABASE']['KEGG']

    mkdir_p(config['OUTPUT']['OUT_DIR'])
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_KEGG']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=kegg_path,
        i_file=config['INPUT']['INPUT_FILE'])
    execute_command(align_command, out_file_path=out_file)
    return out_file


def align_to_mges(config, TELS_statistcs):
    aligner_exe = config['TOOLS']['ALIGNER']
    aligner_flags = config['TOOLS']['ALIGNER_FLAGS']
    aligner_flags = aligner_flags + config['TOOLS']['ALINGER_PB_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_ONT_OPTION']
    aligner_flags = aligner_flags + config['TOOLS']['ALIGNER_HIFI_OPTION']

    plasmids_path = config['DATABASE']['PLASMIDS']
    aclame_path = config['DATABASE']['ACLAME']
    iceberg_path = config['DATABASE']['ICEBERG']

    mkdir_p(config['OUTPUT']['OUT_DIR'])

    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_PLASMIDS']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=plasmids_path,
        i_file=config['INPUT']['INPUT_FILE'])
    execute_command(align_command, out_file_path=out_file)

    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ACLAME']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=aclame_path,
        i_file=config['INPUT']['INPUT_FILE'])
    execute_command(align_command, out_file_path=out_file)

    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ICEBERG']
    align_command = '{exe} {flags} {db} {i_file}'.format(
        exe=aligner_exe,
        flags=aligner_flags,
        db=iceberg_path,
        i_file=config['INPUT']['INPUT_FILE'])
    execute_command(align_command, out_file_path=out_file)


def gen_resistome(config, TELS_statistcs):
    gen_resistome_script = config['SCRIPTS']['GEN_RESISTOME']
    sam_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_MEGARES']
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT']

    gen_resistome_command = 'python {script} -s {sam_file} -o {out_name} -c {config_path}'.format(
        script=gen_resistome_script,
        sam_file=sam_file,
        out_name=out_file,
        config_path=config['MISC']['CONFIG_FILE']
    )
    execute_command(gen_resistome_command)


def gen_mobilome(config, TELS_statistcs):
    gen_mobilome_script = config['SCRIPTS']['GEN_MOBILOME']
    sam_file_plasmids = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + \
                        config['EXTENSION']['A_TO_PLASMIDS']
    sam_file_aclame = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ACLAME']
    sam_file_iceberg = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ICEBERG']
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT']

    gen_mobilome_command = 'python {script} -p {sam_plasmids} -a {sam_aclame} -i {sam_iceberg} -o {out_name}  -c {config_path}'.format(
        script=gen_mobilome_script,
        sam_plasmids=sam_file_plasmids,
        sam_aclame=sam_file_aclame,
        sam_iceberg=sam_file_iceberg,
        out_name=out_file,
        config_path=config['MISC']['CONFIG_FILE']
    )
    execute_command(gen_mobilome_command)


def gen_colocalizations(config, TELS_statistcs):
    gen_colocalizations_script = config['SCRIPTS']['FIND_COLOCALIZATIONS']
    sam_file_plasmids = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + \
                        config['EXTENSION']['A_TO_PLASMIDS']
    sam_file_aclame = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ACLAME']
    sam_file_iceberg = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_ICEBERG']
    sam_file_kegg = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_KEGG']
    sam_file_megares = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['A_TO_MEGARES']
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['COLOCALIZATIONS']
    skip_begin = 0
    skip_end = 0
    if "V2" in config['INPUT']['INPUT_FILE_NAME_EXT']:
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
        reads=config['INPUT']['INPUT_FILE'],
        skip_end=skip_end,
        config_file=config['MISC']['CONFIG_FILE']
    )
    execute_command(find_colocalizations_command, out_file_path=out_file)


def gen_colocalizations_richness(config, TELS_statistcs):
    gen_colocalizations_richness_script = config['SCRIPTS']['COLOCALIZATIONS_RICHNESS']
    colocalizations_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + \
                           config['EXTENSION']['COLOCALIZATIONS']
    out_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION'][
        'COLOCALIZATIONS_RICHNESS']

    colocalizations_richness_command = 'python {script} -i {coloc_file} -c {config_path}'.format(
        script=gen_colocalizations_richness_script,
        coloc_file=colocalizations_file,
        config_path=config['MISC']['CONFIG_FILE']
    )
    execute_command(colocalizations_richness_command, out_file_path=out_file)

    arg_read_lengths = list()
    TELS_statistcs['N_ARG_ON_TARGET_READS'] = 0
    genes_list_csvs = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['GENES_LIST']
    with open(genes_list_csvs, 'r') as genes_list_csv_handler:
        reader = csv.reader(genes_list_csv_handler)
        header = next(reader)
        for row in reader:
            read_name  = row[0]
            amr_genes  = row[1].split(';')
            mge_genes  = row[2].split(';')
            kegg_genes = row[3].split(';')

            if (len(amr_genes) > 0):
                arg_read_lengths.append(TELS_statistcs['READ_LENGTHS'][read_name])
    if len(arg_read_lengths) > 0:
        TELS_statistcs['N_ARG_ON_TARGET_READS'] = str(len(arg_read_lengths))
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_MEAN'] = str(statistics.mean(arg_read_lengths))
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_RANGE'] = str((min(arg_read_lengths), max(arg_read_lengths)))
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_STD_DEV'] = str(statistics.stdev(arg_read_lengths))
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_VARIANCE'] = str(statistics.variance(arg_read_lengths))
        from scipy.stats import kurtosis
        from scipy.stats import skew
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_SKEW'] = str(skew(arg_read_lengths))
        TELS_statistcs['ARG_ON_TARGET_READ_LENGTH_KURTOSIS'] = str(kurtosis(arg_read_lengths))

def print_statistics(config, TELS_statistics):
    with open('{}/{}_stats.csv'.format(config['OUTPUT']['OUT_DIR'], config['INPUT']['INPUT_FILE_NAME_EXT']), 'w') as stats_csv_file:
        stats_writer = csv.writer(stats_csv_file)
        header = ['Metric', 'Value']
        stats_writer.writerow(header)

        for stat_name, stat_value in TELS_statistics.items():
            if (stat_name != 'READ_LENGTHS'):
                line = [stat_name, stat_value]
                stats_writer.writerow(line)


def main():
    parser = argparse.ArgumentParser(description='Colocalizations Pipeline')
    parser.add_argument('-i', help='Input file', dest='input_path', default='')
    parser.add_argument('-o', help='Output directory', dest='output_dir_path', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', default='./config.ini')
    parser.add_argument('-t', help='Helper threads', dest='helper_threads', default='1')
    args = parser.parse_args()

    root_logger = init_logger()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    if config['SCRIPTS']['BASE_PATH'] == '':
        config['SCRIPTS']['BASE_PATH'] = os.path.dirname(os.path.abspath(__file__))

    config['SCRIPTS']['FIND_DUPLICATES'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['FIND_DUPLICATES'])
    config['SCRIPTS']['DEDUPLICATE'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['DEDUPLICATE'])
    config['SCRIPTS']['GEN_MOBILOME'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['GEN_MOBILOME'])
    config['SCRIPTS']['GEN_RESISTOME'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['GEN_RESISTOME'])
    config['SCRIPTS']['FIND_COLOCALIZATIONS'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['FIND_COLOCALIZATIONS'])
    config['SCRIPTS']['COLOCALIZATIONS_RICHNESS'] = os.path.join(config['SCRIPTS']['BASE_PATH'], config['SCRIPTS']['COLOCALIZATIONS_RICHNESS'])

    config['INPUT'] = dict()
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.input_path)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.input_path))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])
    config['OUTPUT'] = dict()
    config['OUTPUT']['OUT_DIR'] = os.path.abspath(args.output_dir_path)
    config['MISC']['HELPER_THREADS'] = args.helper_threads
    config['MISC']['CONFIG_FILE'] = os.path.abspath(args.config_path)

    TELS_statistcs = dict()
    TELS_statistcs['READ_LENGTHS'] = dict()
    TELS_statistcs['READS_BEFORE_DEDUPLICATION'] = 0
    TELS_statistcs['READS_AFTER_DEDUPLICATION'] = 0

    with open(config['INPUT']['INPUT_FILE'], 'r') as reads_file_handle:
        for record in SeqIO.parse(reads_file_handle, read_file_type(config['INPUT']['INPUT_FILE'])):
            read_len = len(record.seq)
            TELS_statistcs['READS_BEFORE_DEDUPLICATION'] += 1
            TELS_statistcs['READ_LENGTHS'][record.name] = read_len
        root_logger.info("Reads befgore deduplication {}".format(TELS_statistcs['READS_BEFORE_DEDUPLICATION']))

    if config['PIPELINE_STEPS']['DEDUPLICATE'] in ['True', 'true']:
        root_logger.info("Deduplicating: Finding duplicates")
        find_duplicates(config, TELS_statistcs)
        root_logger.info("Deduplicating: Filtering duplicates")
        deduplicate(config, TELS_statistcs)
    else:
        deduped_file = config['OUTPUT']['OUT_DIR'] + '/' + config['INPUT']['INPUT_FILE_NAME_EXT'] + config['EXTENSION']['DEDUPLICATED']
        if (os.path.isfile(deduped_file)):
            TELS_statistcs['READS_AFTER_DEDUPLICATION'] = sum(1 for record in SeqIO.parse(open(deduped_file, 'r'), read_file_type(deduped_file)))
            TELS_statistcs['READS_AFTER_DEDUPLICATION_PERC'] = (float(TELS_statistcs['READS_AFTER_DEDUPLICATION']) / float(TELS_statistcs['READS_BEFORE_DEDUPLICATION'])) * 100
            config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(deduped_file)
            config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
            config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(deduped_file))
            config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    if config['PIPELINE_STEPS']['ALIGN_TO_MEGARES'] in ['True', 'true']:
        root_logger.info("Aligning to Megares")
        align_to_megares(config, TELS_statistcs)

    if config['PIPELINE_STEPS']['ALIGN_TO_MGES'] in ['True', 'true']:
        root_logger.info("Aligning to MGEs")
        align_to_mges(config, TELS_statistcs)

    if config['PIPELINE_STEPS']['ALIGN_TO_KEGG'] in ['True', 'true']:
        root_logger.info("Aligning to KEGG")
        align_to_kegg(config, TELS_statistcs)

    if config['PIPELINE_STEPS']['COMPUTE_RESISTOME'] in ['True', 'true']:
        root_logger.info("Generating Resistome")
        gen_resistome(config, TELS_statistcs)

    if config['PIPELINE_STEPS']['COMPUTE_MOBILOME'] in ['True', 'true']:
        root_logger.info("Generating Mobilome")
        gen_mobilome(config, TELS_statistcs)

    if config['PIPELINE_STEPS']['COMPUTE_COLOCALIZATIONS'] in ['True', 'true']:
        root_logger.info("Generating Colocalizations")
        gen_colocalizations(config, TELS_statistcs)
        gen_colocalizations_richness(config, TELS_statistcs)

    print_statistics(config, TELS_statistcs)


if __name__ == "__main__":
    main()