import sys
import os
from utils import *
import random
import itertools

sel_weights = [0.25, 0.25, 0.25, 0.25]
alpha_prob = 0.25 # probability for selecting a \in {A,C,T,G}
alphas = ['A', 'C', 'T', 'G']

def parse_conf_file(conf_file: str):
    mode = ""
    general_opts = {"idx_dir": "", "out_prefix": "", "error_rate": 0.0, "rand_seed": 0}
    mode_opts = {"ratios": None, "num_reads": 0, "num_reads_per_comb": 0, "total_group_file": "", "filter_ratio": 0}

    with open(conf_file, "r") as fd:
        for line in fd:
            if line.startswith("#") or line.startswith("\n"):
                pass
            elif line.startswith("mode="):
                mode = line.strip()[len("mode="):]
                if mode not in {"sub_sampling", "all_sampling", "real_sampling"}:
                    print("invalid mode")
                    sys.exit(1)
            elif line.startswith("idx_dir="):
                general_opts["idx_dir"]=line.strip()[len("idx_dir="):]
                if general_opts["idx_dir"][-1] == "/":
                    general_opts["idx_dir"] = general_opts["idx_dir"][:-1]
            elif line.startswith("out_prefix="):
                general_opts["out_prefix"]=line.strip()[len("out_prefix="):]
                if general_opts["out_prefix"].endswith("_"):
                    general_opts["out_prefix"] = general_opts["out_prefix"][:-1]
            elif line.startswith("error_rate="):
                general_opts["error_rate"] = float(line.strip()[len("error_rate="):])
            elif line.startswith("rand_seed="):
                general_opts["rand_seed"] = int(line.strip()[len("rand_seed="):])
            elif line.startswith("ratios="):
                mode_opts["ratios"] = line.strip()[len("ratios="):]
            elif line.startswith("num_reads="):
                mode_opts["num_reads"] = int(line.strip()[len("num_reads="):])
            elif line.startswith("num_reads_per_comb="):
                mode_opts["num_reads_per_comb"] = int(line.strip()[len("num_reads_per_comb="):])
            elif line.startswith("total_group_file="):
                mode_opts["total_group_file"] = line.strip()[len("total_group_file="):]
            elif line.startswith("filter_ratio="):
                mode_opts["filter_ratio"] = float(line.strip()[len("filter_ratio="):])
            else:
                pass
        fd.close()
    
    print("Mode: ", mode)
    print(general_opts.items())
    print(mode_opts.items())
    return mode, general_opts, mode_opts

def get_rseq_dict(cnme_file: str):
    rseq_dicts = []
    with open(cnme_file, "r") as rname_fd:
        ri = 0
        for line in rname_fd:
            rname = line.strip()
            data = get_file(rname, True)
            rseq_dicts.append({})
            for (sid, seq) in data:
                rseq_dicts[ri][sid] = seq
            ri += 1
        rname_fd.close()
    return rseq_dicts

def sim_reads(seq_fd, err_rate: float, idx_arrs: list, rseq_dicts: list, read_glb_idx, part_reads, combs=None, comb_tuple=None):
    if comb_tuple == None:
        comb_sids = []
        comb_seqs = []
        for idx, jdx in enumerate(combs):
            sid = idx_arrs[idx][jdx]
            comb_sids.append(sid)
            comb_seqs.append(rseq_dicts[idx][sid])
        raw_seq = "".join(comb_seqs)
    else:
        (comb_sids, raw_seq) = comb_tuple

    num_shift = 0
    len_raw_seq = len(raw_seq)
    num_err = round(len_raw_seq * err_rate)
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
    # print(f"Randomly generate {part_reads} reads for plasmid type: {comb_sids}, num cyclic shift: {num_shift}")
    return comb_sids

def sub_sampling(general_opts: dict, mode_opts: dict):
    type_ratios = list(map(float, mode_opts["ratios"].split(":")))
    sum_ratios = sum(type_ratios)

    rseq_dicts = get_rseq_dict(general_opts["idx_dir"] + "/comps_nme.txt")
    idx_arrs = read_idx_file(general_opts["idx_dir"] + "/comps_idx.txt")

    lens_per_components = [[_ for _ in range(len(column))] for column in idx_arrs]
    comb_list = list(itertools.product(*lens_per_components))
    num_types = len(type_ratios)

    part_len = round(len(comb_list) / num_types)
    print("total simulated types: ", num_types)

    out_idx = general_opts["out_prefix"] + "_read_idx.txt"
    out_seq = general_opts["out_prefix"] + "_read_seq.fasta"
    Create(out_idx)
    Create(out_seq)

    idx_fd = open(out_idx, "w")
    seq_fd = open(out_seq, "w")
    read_glb_idx = 0
    for part_idx in range(num_types):
        rand_idx = random.randint(part_idx * part_len, (part_idx + 1) * part_len - 1)
        # selected combo
        part_reads = round((mode_opts["num_reads"] * type_ratios[part_idx]) / sum_ratios)
        comb_sids = sim_reads(seq_fd, general_opts["error_rate"], idx_arrs, rseq_dicts, read_glb_idx, part_reads, combs=comb_list[rand_idx])
        idx_fd.write(",".join(comb_sids) + f",{read_glb_idx},{read_glb_idx+part_reads-1}\n")
        read_glb_idx += part_reads
    idx_fd.close()
    seq_fd.close()
    return


def all_sampling(general_opts: dict, mode_opts: dict):
    rseq_dicts = get_rseq_dict(general_opts["idx_dir"] + "/comps_nme.txt")
    idx_arrs = read_idx_file(general_opts["idx_dir"] + "/comps_idx.txt")

    lens_per_components = [[_ for _ in range(len(column))] for column in idx_arrs]
    comb_list = list(itertools.product(*lens_per_components))
    total_combs = len(comb_list)
    num_types = total_combs

    print("total simulated types: ", num_types)

    out_idx = general_opts["out_prefix"] + "_read_idx.txt"
    out_seq = general_opts["out_prefix"] + "_read_seq.fasta"
    Create(out_idx)
    Create(out_seq)

    idx_fd = open(out_idx, "w")
    seq_fd = open(out_seq, "w")
    read_glb_idx = 0
    for part_idx in range(num_types):
        part_reads = mode_opts["num_reads_per_comb"]
        comb_sids = sim_reads(seq_fd, general_opts["error_rate"], idx_arrs, rseq_dicts, read_glb_idx, part_reads, combs=comb_list[part_idx])
        idx_fd.write(",".join(comb_sids) + f",{read_glb_idx},{read_glb_idx+part_reads-1}\n")
        read_glb_idx += part_reads
    print("total simulated types: ", num_types)
    idx_fd.close()
    seq_fd.close()
    return

def real_sampling(general_opts: dict, mode_opts: dict):
    rseq_dicts = get_rseq_dict(general_opts["idx_dir"] + "/comps_nme.txt")
    idx_arrs = read_idx_file(general_opts["idx_dir"] + "/comps_idx.txt")

    lens_per_components = [[_ for _ in range(len(column))] for column in idx_arrs]
    comb_list = list(itertools.product(*lens_per_components))
    total_combs = len(comb_list)

    all_comb_dict = {}
    for part_idx in range(total_combs):
        combs = comb_list[part_idx]
        comb_sids = []
        comb_seqs = []
        for idx, jdx in enumerate(combs):
            sid = idx_arrs[idx][jdx]
            comb_sids.append(sid)
            comb_seqs.append(rseq_dicts[idx][sid])
        all_comb_dict[",".join(comb_sids)] = "".join(comb_seqs)

    ds_comb_dict = {}
    with open(mode_opts["total_group_file"], "r") as fd:
        for line in fd:
            splited = line.strip().split(",")
            name, count = ",".join(splited[:-1]), int(splited[-1])
            if name == "fail" or name == "contamination" or '*' in name:
                continue
            ds_comb_dict[name] = count
        fd.close()
    
    shuffled_ds_combs = list(ds_comb_dict.keys())
    random.shuffle(shuffled_ds_combs)

    num_filtered = round(len(shuffled_ds_combs) * mode_opts["filter_ratio"])

    out_idx = general_opts["out_prefix"] + "_read_idx.txt"
    out_seq = general_opts["out_prefix"] + "_read_seq.fasta"
    Create(out_idx)
    Create(out_seq)

    idx_fd = open(out_idx, "w")
    seq_fd = open(out_seq, "w")
    read_glb_idx = 0
    # filter first #num_filtered combs
    for combs in shuffled_ds_combs[num_filtered:]:
        part_reads = ds_comb_dict[combs]
        comb_sids = sim_reads(seq_fd, general_opts["error_rate"], None, None, read_glb_idx, part_reads, comb_tuple=(combs, all_comb_dict[combs]))
        idx_fd.write(combs + f",{read_glb_idx},{read_glb_idx+part_reads-1}\n")
        read_glb_idx += part_reads
    
    out_fil = general_opts["out_prefix"] + "_fil_idx.txt"
    out_inc = general_opts["out_prefix"] + "_inc_idx.txt"
    Create(out_fil)
    Create(out_inc)
    fil_fd = open(out_fil, "w")
    inc_fd = open(out_inc, "w")

    for combs in shuffled_ds_combs[:num_filtered]:
        fil_fd.write(combs + f",{ds_comb_dict[combs]}\n")

    for combs in shuffled_ds_combs[num_filtered:]:
        inc_fd.write(combs + f",{ds_comb_dict[combs]}\n")

    idx_fd.close()
    seq_fd.close()

    fil_fd.close()
    inc_fd.close()
    return



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"{sys.argv[0]} <sim_conf.txt>")
        sys.exit(1)
    
    _, conf_file = sys.argv
    mode, general_opts, mode_opts = parse_conf_file(conf_file)

    random.seed(general_opts["rand_seed"])

    {"sub_sampling": sub_sampling, 
     "all_sampling": all_sampling, 
     "real_sampling": real_sampling}[mode](general_opts, mode_opts)
