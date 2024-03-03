import os

import gzip
from Bio import SeqIO

def System(command: str):
    return os.system(command)

def Create(file: str):
    return System("echo "" > " + file)

def Create_Dir(dir: str):
    if not os.path.exists(dir):
        # create new directory
        os.makedirs(dir)
        print("Directory: ", dir, "created.")
        return True
    else:
        print(f"Directory: {dir} exists")
        return False

def parse_reads(read_file: str):
    if read_file.endswith("fastq"):
        return SeqIO.parse(open(read_file, "r"), "fastq")
    elif read_file.endswith("fasta"):
        return SeqIO.parse(open(read_file, "r"), "fasta")
    elif read_file.endswith("fastq.gz"):
        return SeqIO.parse(gzip.open(read_file, "rt"), "fastq")
    elif read_file.endswith("fasta.gz"):
        return SeqIO.parse(gzip.open(read_file, "rt"), "fasta")
    else:
        return None

def line_counter(file: str):
    c = 0
    with open(file, "r") as fd:
        for line in fd:
            if line.startswith(">") or line.startswith("@"):
                c += 1
        fd.close()
    return c

def get_file(file: str, is_fasta: bool):
    if is_fasta:
        return get_fasta(file)
    else:
        return get_fastq(file)

def get_fasta(fasta_file: str):
    res = []
    with open(fasta_file, "r") as fa_fd:
        sid = ""
        seq = ""
        for line in fa_fd:
            if line.startswith(">"):
                # process previous entry
                if sid != "":
                    res.append((sid, seq))
                sid = line.strip().replace(' ', '\t').split("\t")[0][1:]
                seq = ""
            else:
                seq += line.strip()
        fa_fd.close()
        if sid != "":
            res.append((sid, seq))
    return res

def get_fastq(fastq_file: str):
    reads = []
    counter = 0
    with open(fastq_file, "r") as fd:
        prev_id = ""
        prev_seq = ""
        i = 0
        for l in fd:
            if i == 0:
                assert l.startswith("@")
                counter += 1

                if prev_id != "":
                    reads.append((prev_id, prev_seq))
                
                # read id, get first non-space id, without @ symbol
                prev_id = l.strip().replace(' ', '\t').split("\t")[0][1:]
                prev_seq = ""
            elif i == 1:
                prev_seq += l.strip()
            elif i == 2:
                if not l == "+\n":
                    # multiline sed
                    prev_seq += l.strip()
                    continue
            i = (i+1) % 4
        fd.close()
        if prev_id != "":
            reads.append((prev_id, prev_seq))

    print("total reads: ", counter)
    return reads

def get_reads_todict(fastq_file: str):
    reads = {}
    counter = 0
    with open(fastq_file, "r") as fd:
        prev_id = ""
        prev_seq = ""
        i = 0
        for l in fd:
            if i == 0:
                assert l.startswith("@")
                counter += 1

                if prev_id != "":
                    reads[prev_id] = prev_seq
                
                # read id, get first non-space id
                prev_id = l.strip()
                prev_seq = ""
            elif i == 1:
                prev_seq += l.strip()
            elif i == 2:
                if not l == "+\n":
                    # multiline sed
                    prev_seq += l.strip()
                    continue
            i = (i+1) % 4
        fd.close()
        if prev_id != "":
            reads[prev_id] = prev_seq

    print("total reads: ", counter)
    return reads

def reformat_fasta(in_file: str, out_file: str):
    os.system("echo "" > " + out_file)
    with open(in_file, "r") as in_fd:
        with open(out_file, "w") as out_fd:
            sid = ""
            seq = ""
            for line in in_fd:
                if line.startswith(">"):
                    # process previous entry
                    if sid != "":
                        out_fd.write(sid + "\n")
                        out_fd.write(seq + "\n")
                    sid = line.strip()
                    
                    seq = ""
                else:
                    seq += line.strip().upper()
            out_fd.close()
            if sid != "":
                out_fd.write(sid + "\n")
                out_fd.write(seq + "\n")
        in_fd.close()

def store_reads(reads: dict, only_id=True):
    # store basic read informations
    System("echo "" > reads.txt")
    with open("reads.txt", "w") as fd:
        for rid, rseq in reads.items():
            if only_id:
                fd.write(f"{rid}\n")
            else:
                fd.write(f"{rid}:{rseq}\n")
        fd.close()

def read_idx_file(idx_file: str):
    idx_arrs = None
    # read idx_file
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

def proc_minimap_res(pfile: str):
    mapq_dict = {}
    mrate_dict_acc = {}
    num_mapped = 0
    with open(pfile, "r") as paf_fd:
        for line in paf_fd:
            splited = line.strip().split("\t")
            seg_no = splited[0]
            mapq = int(splited[11])
            match_base = int(splited[9]) # matching base
            map_base = int(splited[10]) # mapping base
            if seg_no not in mapq_dict:
                mapq_dict[seg_no] = mapq
                mrate_dict_acc[seg_no] = [(match_base, map_base)]
                num_mapped += 1
            else:
                mapq_dict[seg_no] = min(mapq, mapq_dict[seg_no])
                mrate_dict_acc[seg_no].append((match_base, map_base))
        paf_fd.close()
    mrate_dict = {}
    for seg_no, mrates in mrate_dict_acc.items():
        match_bases = sum(mb1 for (mb1, _) in mrates)
        map_bases = sum(mb2 for (_, mb2) in mrates)
        mrate_dict[seg_no] = match_bases / map_bases
    
    return mapq_dict, mrate_dict, num_mapped

# get query - ref mapping, if multiple mapping on same query, 
def proc_minimap_identity(pfile: str):
    aln_mapping = {}
    with open(pfile, "r") as paf_fd:
        for line in paf_fd:
            splited = line.strip().split("\t")
            seg_no = splited[0]
            ref_no = splited[5]
            if seg_no not in aln_mapping:
                aln_mapping[seg_no] = set()
            aln_mapping[seg_no].add(int(ref_no))
        paf_fd.close()

    return aln_mapping