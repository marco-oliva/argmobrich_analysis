from Bio import SeqIO

import matplotlib.pyplot as plt

import csv
import os
import sys
import gzip

if __name__ == "__main__":
    fastq_with_dups = sys.argv[1]
    sets_csv = sys.argv[2]

    sets = []
    with open(sets_csv, 'r') as sets_handle:
        sets_csv_reader = csv.reader(sets_handle)
        for row in sets_csv_reader:
            curr_set_elements = []
            for read_id in row:
                curr_set_elements.append(read_id)

            #curr_set = frozenset(curr_set_elements)
            sets.append(set(curr_set_elements))
            #set_deduplicated[curr_set] = False

    #Sets are not disjoint :(
    are_sets_disjoint = False
    itr = 0
    while(not are_sets_disjoint):
        print("iteration: " + str(itr))
        print(len(sets))
        itr += 1
        are_sets_disjoint = True
        combined_sets = []
        used_set_indices = {}
        for i in range(0, len(sets)):
            used_set_indices[i] = False

        for i in range(0, len(sets)):
            if not used_set_indices[i]:
                combined_set = sets[i]
                used_set_indices[i] = True
            else:
                #print("continue to next from " + str(i))
                continue

            for j in range(i+1, len(sets)):
                if not used_set_indices[j]:
                    if not len(combined_set & sets[j]) == 0:
                        #print("merging " + str(i) + ", " + str(j))
                        combined_set = combined_set | sets[j]
                        are_sets_disjoint = False
                        used_set_indices[j] = True
                else:
                    continue

            #print("appending: " + str(i))
            combined_sets.append(combined_set)

        print(len(combined_sets))
        sets = combined_sets
                

    set_deduplicated = {}
    frozen_sets = []
    print(len(sets))
    for curr_set in sets:
        frozen_set = frozenset(curr_set)
        frozen_sets.append(frozen_set)
        set_deduplicated[frozen_set] = False

    dedup_records = []
    set_sizes = []
    num_singletons = 0
    for record in SeqIO.parse(fastq_with_dups, "fastq"):
        record_in_set = False
        for dup_set in frozen_sets:
            if record.id in dup_set:
                record_in_set = True
                curr_record_set = dup_set
                break

        if record_in_set:
            if not set_deduplicated[curr_record_set]:
                dedup_records.append(record)
                set_sizes.append(len(curr_record_set))
                set_deduplicated[curr_record_set] = True
        else: #Singletons are not in the set csv
            num_singletons += 1
            dedup_records.append(record)

    with gzip.open("deduplicated_" + os.path.splitext(os.path.basename(fastq_with_dups))[0] + ".fastq.gz", 'w') as out_handle:
        SeqIO.write(dedup_records, out_handle, "fastq")

    #Plot histogram
    fig, ax = plt.subplots(tight_layout=True)
    ax.set_title("deduplicated_" + os.path.splitext(os.path.basename(fastq_with_dups))[0] + " (singletons: " + str(num_singletons) + ")")
    ax.set_xlabel("Set size (number of reads)")
    ax.set_ylabel("Frequency")
    bin_tuple = ax.hist(set_sizes, bins=75)
    plt.savefig("deduplicated_" + os.path.splitext(os.path.basename(fastq_with_dups))[0] + ".svg")


