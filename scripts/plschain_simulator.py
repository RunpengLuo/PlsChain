import sys
import os
from utils import *
import random
import itertools

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(f"{sys.argv[0]} <idx_dir> <out_prefix> <num_reads> <error_rate> t1:t2:t3:... <seed: 0>")
        print(f"{sys.argv[0]} <idx_dir> <out_prefix> <num_reads_per_comb> <error_rate> all <seed: 0>")
        # error rate in percentage
        sys.exit(1)
    
    if len(sys.argv) == 7:
        random.seed(int(sys.argv[6]))
    else:
        random.seed(0)

    _, idx_dir, out_prefix = sys.argv[:3]
    num_reads, err_rate = [int(sys.argv[3]), float(sys.argv[4])]
    if sys.argv[5].startswith("all"):
        type_ratios = None
    else:
        type_ratios = list(map(float, sys.argv[5].split(":")))
        sum_ratios = sum(type_ratios)

    # [ins, del, mis, ns]
    sel_weights = [0.25, 0.25, 0.25, 0.25]
    alpha_prob = 0.25 # probability for selecting a \in {A,C,T,G}

    if idx_dir[-1] != "/":
        idx_dir += "/"

    rseq_dicts = []
    with open(idx_dir + "comps_nme.txt", "r") as rname_fd:
        ri = 0
        for line in rname_fd:
            rname = line.strip()
            data = get_file(rname, True)
            rseq_dicts.append({})
            for (sid, seq) in data:
                rseq_dicts[ri][sid] = seq
            ri += 1
        rname_fd.close()
    print([sx.keys() for sx in rseq_dicts])

    idx_arrs = read_idx_file(idx_dir + "comps_idx.txt")

    lens_per_components = [[_ for _ in range(len(column))] for column in idx_arrs]
    comb_list = list(itertools.product(*lens_per_components))
    total_combs = len(comb_list)

    if type_ratios != None:
        num_types = len(type_ratios)
    else:
        # select all types
        num_types = total_combs

    part_len = round(total_combs / num_types)

    alphas = ['A', 'C', 'T', 'G']

    out_idx = out_prefix + "_read_idx.txt"
    out_seq = out_prefix + "_read_seq.fasta"
    Create(out_idx)
    Create(out_seq)

    idx_fd = open(out_idx, "w")
    seq_fd = open(out_seq, "w")

    read_glb_idx = 0
    for part_idx in range(num_types):
        rand_idx = random.randint(part_idx * part_len, (part_idx + 1) * part_len - 1)
        # selected combo
        if type_ratios != None:
            part_reads = round((num_reads * type_ratios[part_idx]) / sum_ratios)
        else:
            part_reads = num_reads
        comb_sids = []
        comb_seqs = []
        raw_seq = ""

        for idx, jdx in enumerate(comb_list[rand_idx]):
            sid = idx_arrs[idx][jdx]
            comb_sids.append(sid)
            raw_seq += rseq_dicts[idx][sid]

        num_shift = 0
        len_raw_seq = len(raw_seq)
        num_err = round((len_raw_seq * err_rate) / 100)
        for i in range(read_glb_idx, read_glb_idx+part_reads):
            # err_stat = {"I": 0, "D": 0, "M": 0, "N": 0}
            j = 0
            seq_with_err = list(raw_seq)
            # randomly inserting single nucleotide errors
            while j < num_err:
                j += 1
                rand_loc = random.randint(0, len(seq_with_err)-1)
                choice = random.choices(['I', 'D', 'M', 'N'], weights=sel_weights, k=1)[0]
                # err_stat[choice] += 1
                if choice == 'N':
                    seq_with_err[rand_loc] = 'N'
                elif choice == 'D':
                    seq_with_err.pop(rand_loc)
                else:
                    rand_char = alphas[random.randint(0,3)]
                    if choice == 'I':
                        seq_with_err.insert(rand_loc+1, rand_char)
                    else:
                        # M
                        seq_with_err[rand_loc] = rand_char
            
            # randomly determine a cyclic shift position.
            perform_shift = random.randint(0, 1)
            if perform_shift == 1:
                # randomly selected pivot at RAND[1, N-1] for length-N sequence.
                shift_position = random.randint(1, len(seq_with_err)-2)
                seq_with_err = seq_with_err[shift_position:]+seq_with_err[:shift_position]
                num_shift += 1
            else:
                pass
            
            seq_str="".join(seq_with_err)
            seq_fd.write(f">{i}\n{seq_str}\n")
        idx_fd.write(",".join(comb_sids) + f",{read_glb_idx},{read_glb_idx+part_reads-1}\n")
        read_glb_idx += part_reads
        print(f"Randomly generate {part_reads} for plasmid type: {comb_sids}, num cyclic shift: {num_shift}")
    print("total simulated types: ", num_types)
    idx_fd.close()
    seq_fd.close()