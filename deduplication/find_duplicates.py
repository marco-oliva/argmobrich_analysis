from Bio import SearchIO

import csv
import os
import sys

if __name__ == "__main__":

    basename = sys.argv[1]
    tsv_name = os.path.basename(sys.argv[1]) + ".tsv"

    for i in range(0,200):
        
        duplication_sets = [] #container of sets, where each set represents a clique of CCS that closely align
        #each alignment in psl: if qualifications met, those two belong to the same set
        print("Reading in " + basename + '_' + str(i) + ".psl", "blat-psl")
        for qresult in SearchIO.parse(basename + '_' + str(i) + ".psl", "blat-psl"):
            for hit in qresult:
                for hsp in hit.hsps:
                    #print( abs(hsp.hit_end - hsp.hit_start) / float(hit.seq_len) )
                    if sum(hsp.hit_span_all) >= 0.9*hit.seq_len and sum(hsp.query_span_all) >= 0.9*qresult.seq_len:
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
        #append dup set to tsv
        print("Writing to " + tsv_name)
        with open(tsv_name, 'a') as tsv_handle:
            tsv_writer = csv.writer(tsv_handle)
            tsv_writer.writerows(sorted(non_singleton_dup_sets, key = lambda dup_set: len(dup_set), reverse=True))
                

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
                    if sum(hsp.hit_span_all) >= 0.9 * hit.seq_len and sum(hsp.query_span_all) >= 0.9 * qresult.seq_len:
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