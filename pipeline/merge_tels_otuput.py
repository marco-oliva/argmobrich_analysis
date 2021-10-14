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