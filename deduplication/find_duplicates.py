import signal
import numpy as np
import os.path
import gzip
import subprocess
import errno
import csv
import argparse
from multiprocessing import Pool

from Bio import SeqIO, SearchIO
from sklearn.cluster import KMeans


# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, time_it=False, seconds=1000000):
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
        print('Executing: {}'.format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        return None
    if output:
        output = output.decode('utf-8')
    if err:
        err = err.decode('utf-8')
    return output

# mkdir -p
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

def run_blat(input_file, out_pls_file):
    blat_command = 'blat {fasta_name} {fasta_name} {psl_name}'.format(fasta_name=input_file,psl_name=out_pls_file)
    execute_command(blat_command)
    return out_pls_file

def main():
    parser = argparse.ArgumentParser(description='Cluster reads based on read length')
    parser.add_argument('-r', help='Reads File (FASTA/FASTQ)', type=str, dest='reads_file', required=True)
    parser.add_argument('-o', help='Output directory', type=str, dest='out_dir', required=True)
    parser.add_argument('-n', help='Number of clusters to generate', type=int, dest='num_clusters', default=200)
    parser.add_argument('-t', help='Number of threads', type=int, dest='num_threads', default=1)
    args = parser.parse_args()

    reads_file = args.reads_file
    output_dir = args.out_dir
    num_clusters = args.num_clusters

    fastq_extensions = ['fastq', 'fq', 'FASTQ', 'FQ', 'fq.gz', 'fastq.gz', 'FQ.gz', 'FASTQ.gz']
    fasta_extensions = ['fasta', 'fa', 'FASTA', 'FA', 'fa.gz', 'fasta.gz', 'FA.gz', 'FASTA.gz']

    file_type = ''
    for p_ext in fastq_extensions:
        if (reads_file.endswith(p_ext)):
            file_type = 'fastq'
    for p_ext in fasta_extensions:
        if (reads_file.endswith(p_ext)):
            file_type = 'fasta'

    if file_type == '':
        print('File has to be either fastq or fasta')
        exit()

    # Iterate through every read. Accumulate number of reads while recording read length
    read_lengths = ([], [])

    if (reads_file.endswith('.gz')):
        print('Opening gzipped file')
        file_handler = gzip.open(reads_file, 'rt')
    else:
        print('Opening uncompressed file')
        file_handler = open(reads_file, 'rt')

    for record in SeqIO.parse(file_handler, file_type):
        read_len = len(record.seq)
        read_lengths[0].append(read_len)
        read_lengths[1].append(record)

    X = np.array(read_lengths[0])

    if (len(read_lengths[0]) < num_clusters):
        num_clusters = int(len(read_lengths[0]) / 10)
    print("Num of cluster used: {}".format(num_clusters))

    kmeans = KMeans(n_clusters=num_clusters).fit(X.reshape((X.shape[0], 1)))

    fasta_files = list()
    for i in range(0, num_clusters):
        sequences = []
        for j in range(0, kmeans.labels_.shape[0]):
            if kmeans.labels_[j] == i:
                sequences.append(read_lengths[1][j])

        output_filename = output_dir + '/' + os.path.splitext(os.path.basename(reads_file))[0] + '_' + str(i) + '.fasta'
        fasta_files.append(output_filename)
        with open(output_filename, 'w') as output_fasta:
            SeqIO.write(sequences, output_fasta, 'fasta')

    # run blat on each file, using multiple threads
    mkdir_p(output_dir + '/pls_files')
    input_tuples = list()
    for i in range(0, num_clusters):
        # here blat was called
        out_pls = output_dir + '/pls_files/' + str(i) + '.pls'
        fasta_filename = output_dir + '/' + os.path.splitext(os.path.basename(reads_file))[0] + '_' + str(i) + '.fasta'
        input_tuples.append((fasta_filename, out_pls))

    with Pool(args.num_threads) as threads_pool:
        pls_files = threads_pool.starmap(run_blat, input_tuples)

    # find duplicates
    tsv_name = output_dir + '/pls_files/duplicates.tsv'
    for i in range(0, num_clusters):
        duplication_sets = []  # container of sets, where each set represents a clique of CCS that closely align
        # each alignment in psl: if qualifications met, those two belong to the same set
        cur_pls_file = pls_files[i]
        print('Reading in ' + cur_pls_file, 'blat-psl')
        for qresult in SearchIO.parse(cur_pls_file, 'blat-psl'):
            for hit in qresult:
                for hsp in hit.hsps:
                    # print( abs(hsp.hit_end - hsp.hit_start) / float(hit.seq_len) )
                    if sum(hsp.hit_span_all) >= 0.9 * hit.seq_len and sum(
                            hsp.query_span_all) >= 0.9 * qresult.seq_len:
                        if not duplication_sets:
                            duplication_sets.append(set())
                            duplication_sets[0].add(qresult.id)
                            duplication_sets[0].add(hit.id)
                        else:
                            alignment_is_in_set = False
                            for dup_set in duplication_sets:
                                if (qresult.id in dup_set):
                                    dup_set.add(hit.id)
                                    alignment_is_in_set = True
                                    break
                                elif (hit.id in dup_set):
                                    dup_set.add(qresult.id)
                                    alignment_is_in_set = True
                                    break
                            if not alignment_is_in_set:
                                duplication_sets.append(set())
                                duplication_sets[-1].add(qresult.id)
                                duplication_sets[-1].add(hit.id)

        non_singleton_dup_sets = [dup_set for dup_set in duplication_sets if len(dup_set) > 1]
        # append dup set to tsv
        print('Writing to ' + tsv_name)
        with open(tsv_name, 'a') as tsv_handle:
            tsv_writer = csv.writer(tsv_handle)
            tsv_writer.writerows(sorted(non_singleton_dup_sets, key=lambda dup_set: len(dup_set), reverse=True))

if __name__ == "__main__":
    main()


