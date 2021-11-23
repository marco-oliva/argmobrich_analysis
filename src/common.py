import errno
import logging
import signal
import sys
import configparser
import os
import argparse
import subprocess
import gzip
import csv
import numpy as np


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

def remove(path):
    if os.path.isdir(path):
        try:
            os.rmdir(path)
        except OSError as e:  ## if failed, report it back to the user ##
            print("Error: %s - %s." % (e.filename, e.strerror))
    elif os.path.isfile(path):
        try:
            os.remove(path)
        except OSError as e:  ## if failed, report it back to the user ##
            print("Error: %s - %s." % (e.filename, e.strerror))

def init_logger():
    root_logger = logging.getLogger(__name__)
    root_logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(f'%(asctime)s - %(script_name) - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    return root_logger

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def reject_outliers(data, m = 2.):
    data = np.array(data)
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

# Input: config
def read_megares_ontology(config):
    # Create ontology dictionary from MEGARes ontology file
    megares_ontology = dict()
    hierarchy_dict = dict()
    with open(config['DATABASE']['MEGARES_ONTOLOGY'], 'r') as ontology_tsv:
        ontology_reader = csv.reader(ontology_tsv)
        for row in ontology_reader:
            # Skip column names
            if row[0] == "header":
                continue

            cl = row[1]
            mech = row[2]
            group = row[3]

            # Set up hiearachy dict. This will be our tree structure
            if not cl in hierarchy_dict:
                hierarchy_dict[cl] = {}

            if not mech in hierarchy_dict[cl]:
                hierarchy_dict[cl][mech] = []

            if not group in hierarchy_dict[cl][mech]:
                hierarchy_dict[cl][mech].append(group)

            # FIll in our dict
            megares_ontology[row[0]] = {"class": cl,
                                        "mech": mech,
                                        "group": group
                                        }
    return megares_ontology, hierarchy_dict
