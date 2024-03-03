import os
import sys
from utils import *
import time

import itertools

if __name__ == "__main__":
    # generate all combination based on idx_dir
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} <idx_dir> <out_prefix>")
        sys.exit(1)

    _, idx_dir, out_prefix = sys.argv

    if idx_dir[-1] != "/":
        idx_dir += "/"
    if out_prefix[-1] == '_':
        out_prefix = out_prefix[:-1]

    rseq_dicts = []
    with open(idx_dir + "comps_nme.txt", "r") as fd_ref:
        for line in fd_ref:
            ref_file = line.strip()
            rseq_dicts.append({})
            handle = parse_reads(ref_file)
            if handle == None:
                raise Exception(f"{ref_file}, file format un-recognized")
            for record in handle:
                rseq_dicts[-1][record.id] = str(record.seq)
        fd_ref.close()


    idx_arrs = read_idx_file(idx_dir + "comps_idx.txt")

    lens_per_components = [[_ for _ in range(len(column))] for column in idx_arrs]

    out_idx = out_prefix + "_comb_idx.txt"
    out_seq = out_prefix + "_comb_seq.fasta"
    Create(out_idx)
    Create(out_seq)
    with open(out_idx, "w") as idx_fd:
        with open(out_seq, "w") as seq_fd:
            i = 0
            for comb in itertools.product(*lens_per_components):
                comb_names = []
                comb_seqs = []
                for idx, jdx in enumerate(comb):
                    sid = idx_arrs[idx][jdx]
                    comb_names.append(sid)
                
                # twice the sequence, cyclic
                for _ in range(2):
                    for idx, jdx in enumerate(comb):
                        sid = idx_arrs[idx][jdx]
                        seq = rseq_dicts[idx][sid]
                        comb_seqs.append(seq)

                seq_str = "".join(comb_seqs)
                names_str=",".join(comb_names)
                idx_fd.write(f"{names_str}\n")
                seq_fd.write(f">{i}\n{seq_str}\n")
                i += 1
            seq_fd.close()
        idx_fd.close()