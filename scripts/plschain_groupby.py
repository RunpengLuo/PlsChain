import sys
import os
from utils import *

if __name__ == "__main__":
    if len(sys.argv) != 3: 
        print(f"{sys.argv[0]} <query_dir> <idx_dir>")
        sys.exit(1)
    _, out_dir, idx_dir = sys.argv
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]    
    qry_total = out_dir + "/qry_total.csv"
    idx_nmes = idx_dir + "/comps_nme.txt"

    comp_dicts = []
    comp_names = []
    # group by components
    with open(idx_nmes, "r") as fd_ref:
        for line in fd_ref:
            ref_file = line.strip()
            comp_dicts.append({})
            comp_names.append(ref_file)
            handle = parse_reads(ref_file)
            if handle == None:
                raise Exception(f"{ref_file}, file format un-recognized")
            for record in handle:
                comp_dicts[-1][record.id] = 0
        fd_ref.close()

    # group by combinations
    gb_dict = {}
    with open(qry_total, "r") as fd_total:
        for line in fd_total:
            splited = line.strip().split(",")
            rid, names = splited[0], ",".join(splited[1:])
            if names not in gb_dict:
                gb_dict[names] = 0
            gb_dict[names] += 1
            if names == "fail":
                continue
            # for idx, name in enumerate(names.split(",")):
            #     comp_dicts[idx][name] += 1
        fd_total.close()
    # with open(out_dir + "/qry_total_comp.group.csv", "w") as fd_cgroup:
    #     for (comp_name, comp_dict) in zip(comp_names, comp_dicts):
    #         fd_cgroup.write(f"{comp_name}\n")
    #         for comp_idx, comp_count in comp_dict.items():
    #             fd_cgroup.write(f"{comp_idx},{comp_count}\n")
    #         fd_cgroup.write("\n")
    #     fd_cgroup.close()

    with open(out_dir + "/qry_total.group.csv", "w") as fd_group:
        for name, count in sorted(gb_dict.items(), key=lambda tp: tp[1], reverse=True):
            fd_group.write(f"{name},{count}\n")
        fd_group.close()

    print("Done")
