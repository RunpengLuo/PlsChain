import sys
import os

def read_idx_file(idx_file: str):
    idx_arrs = None
    with open(idx_file, "r") as idx_fd:
        num_comp = int(idx_fd.readline().strip())
        idx_arrs = [[] for _ in range(num_comp)]
        glb_idx = 0 # index between components
        for line in idx_fd:
            if line == "\n":
                glb_idx += 1
                continue
            
            sid = line.strip()
            idx_arrs[glb_idx].append(sid)
        idx_fd.close()
    return idx_arrs

def display(out_dir):
    [tota, full, part, cont, uncl] = [0,0,0,0,0]
    with open(out_dir + "/qry_total.csv", 'r') as fd:
        for line in fd:
            tota += 1
            splited = line.strip().split(",")
            rid, names = splited[0], splited[1:]
            name_str = ",".join(names)
            if "contamination" in name_str:
                cont += 1
            elif "fail" in name_str:
                uncl += 1
            elif "*" in name_str:
                part += 1
            else:
                full += 1
        fd.close()
    
    print(f"PlsChain\t{tota}\t{full}\t{part}\t{cont}\t{uncl}")
    
    [tota, full, part, cont, uncl] = [0,0,0,0,0]
    with open(out_dir + "/qry_total.fuzzy.csv", 'r') as fd:
        for line in fd:
            tota += 1
            splited = line.strip().split(",")
            rid, names = splited[0], splited[1:]
            name_str = ",".join(names)
            if "contamination" in name_str:
                cont += 1
            elif "fail" in name_str:
                uncl += 1
            elif "*" in name_str:
                part += 1
            else:
                full += 1
        fd.close()
    print(f"PlsChain-Fuzzy\t{tota}\t{full}\t{part}\t{cont}\t{uncl}")
    return

if __name__ == "__main__":
    if len(sys.argv) != 3: 
        print(f"{sys.argv[0]} <query_dir> <index_dir>")
        sys.exit(1)
    _, out_dir, idx_dir = sys.argv
    
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    if idx_dir[-1] == '/':
        idx_dir = idx_dir[:-1]

    qry_total = out_dir + "/qry_total.csv"
    idx_arrs = read_idx_file(idx_dir + "/comps_idx.txt")

    # Step 1, parse + fuzzy match (if set) + group by
    with open(out_dir + "/qry_total.fuzzy.csv", 'w') as fd_fuzzy:
        gb_dict = {}
        gb_fuzzy_dict = {}
        with open(qry_total, "r") as fd_total:
            for line in fd_total:
                splited = line.strip().split(",")
                rid, names = splited[0], splited[1:]

                name = ",".join(names)
                if name not in gb_dict:
                    gb_dict[name] = 0
                gb_dict[name] += 1

                fnames = list(names)
                if '*' in fnames:
                    for i, v in enumerate(list(fnames)):
                        if v == '*':
                            if len(idx_arrs[i]) == 1:
                                # fuzzy inference
                                fnames[i] = idx_arrs[i][0]

                fname = ",".join(fnames)
                if fname not in gb_fuzzy_dict:
                    gb_fuzzy_dict[fname] = 0
                gb_fuzzy_dict[fname] += 1
                fd_fuzzy.write(f"{rid},{fname}\n")

            fd_total.close()
        fd_fuzzy.close()

    display(out_dir)

    with open(out_dir + "/qry_total.group.csv", "w") as fd_group:
        fd_group.write(f"total_count,{sum(gb_dict.values())}\n")
        for name, count in sorted(gb_dict.items(), key=lambda tp: tp[1], reverse=True):
            fd_group.write(f"{name},{count}\n")
        fd_group.close()
    
    with open(out_dir + "/qry_total.group.fuzzy.csv", "w") as fd_group:
        fd_group.write(f"total_count,{sum(gb_fuzzy_dict.values())}\n")
        for name, count in sorted(gb_fuzzy_dict.items(), key=lambda tp: tp[1], reverse=True):
            fd_group.write(f"{name},{count}\n")
        fd_group.close()

    print("Done")
