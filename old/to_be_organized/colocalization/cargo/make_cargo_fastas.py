from Bio import SeqIO

import csv
import sys

if __name__ == "__main__":

    ccs_info = {}
    dists_filename = sys.argv[1]
    with open(dists_filename, 'r') as dists_csv:
        dists_reader = csv.reader(dists_csv, delimiter='\t')
        for row in dists_reader:
            #Skip rows without coloc
            if not row or row[0][0] != 'm':
                continue

            ccs_id      = row[0]
            amr_header  = row[1]
            amr_start   = row[4]
            amr_end     = row[5]
            mge_header  = row[2]
            mge_start   = row[11]
            mge_end     = row[12]

            if ccs_id not in ccs_info:
                ccs_info[ccs_id] = []

            current_ccs_dict = { "amr_header"   : amr_header,
                                 "amr_start"    : int(amr_start),
                                 "amr_end"      : int(amr_end),
                                 "mge_header"   : mge_header,
                                 "mge_start"    : int(mge_start),
                                 "mge_end"      : int(mge_end),
                               }
            ccs_info[ccs_id].append(current_ccs_dict)

    
    fasta_filename = sys.argv[2]
    for rec in SeqIO.parse(fasta_filename, "fasta"):
        id_for_file = rec.id.replace("/", "_")
        output_fasta_filename = sys.argv[3] + '/' + id_for_file + "_cargo_genes.fasta"
        if rec.id not in ccs_info:
            continue

        #for each colocalization "present" on the record
        for coloc in ccs_info[rec.id]:
            #find portions of record that are not amr gene or mge
            
            cargo_seqs = []
            #figure out whether amr gene or mge comes first on ccs
            if coloc["amr_start"] < coloc["mge_start"]:
                first_stop = coloc["amr_start"]

                second_start = coloc["amr_end"]
                second_stop  = coloc["mge_start"]

                third_start = coloc["mge_end"]
            else:
                first_stop = coloc["mge_start"]

                second_start = coloc["mge_end"]
                second_stop  = coloc["amr_start"]

                third_start = coloc["amr_end"]

            #is there an unaligned sequence at the beginning of the record?
            if first_stop > 14:
                cargo_seqs.append(rec.seq[0:first_stop])

            #is there an unaligned sequence between the amr gene and the mge?
            if second_stop - second_start > 14:
                cargo_seqs.append(rec.seq[second_start:second_stop])

            #is there an unaligned sequence at the end of the record?
            if len(rec.seq) - third_start > 14:
                cargo_seqs.append(rec.seq[third_start:])

            #append those sequences to record's fasta with appropriate headers
            with open(output_fasta_filename, 'a') as out:
                header_base = '>' + rec.id + ":::" + coloc["amr_header"] + ":::" + coloc["mge_header"]
                i = 0
                for seq in cargo_seqs:
                    out.write(header_base + '_' + str(i) + '\n')
                    out.write(str(seq) + '\n')
                    i += 1
    sys.exit(0)
