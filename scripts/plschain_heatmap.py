import sys
import os
import utils
import itertools

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


"""
heatmap plot a 46 * (6*23) for a group of datasets, and map average per
cell with color weight
"""
def get_labels(file: str):
    handle = utils.parse_reads(file)
    if handle == None:
        raise Exception(f"{file}, file format un-recognized")
    return [rec.id for rec in handle]

def save_to_csv(matrix, row_header: list, col_header: list, out_file: str) -> None:
    with open(out_file, 'w') as fd:
        ch = ",".join(f"({col[0]}&{col[1]})" for col in col_header)
        fd.write(f",{ch}\n")
        for i, row_name in enumerate(row_header):
            to_list = [row_name]
            for j in range(len(col_header)):
                to_list.append(str(matrix[i][j]))
            fd.write(",".join(to_list) + "\n")
        fd.close()
    return None

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(f"{sys.argv[0]} <dir1> [<dir2> ...] <promoter.fasta> <peptide.fasta> <terminator.fasta> <out_prefix>")
        sys.exit(1)

    dirs = sys.argv[1:-4]
    [promoter_file, peptide_file, terminator_file, out_prefix] = sys.argv[-4:]

    if out_prefix.endswith("_"):
        out_prefix = out_prefix[:-1]

    labels_promoter = get_labels(promoter_file) # 23
    labels_peptide = get_labels(peptide_file) # 46
    labels_terminator = get_labels(terminator_file) # 6

    prod_promoter_terminator = list(itertools.product(labels_promoter, labels_terminator))

    # column index mapping
    rev_cidx = {}
    for i, (promoter, terminator) in enumerate(prod_promoter_terminator):
        rev_cidx[(promoter, terminator)] = i
    
    # row index mapping
    rev_ridx = {}
    for i, peptide in enumerate(labels_peptide):
        rev_ridx[peptide] = i

    matrices = []
    matrices_opt = []

    print("dimension", len(labels_peptide), len(prod_promoter_terminator))
    total_comb = len(labels_peptide) * len(prod_promoter_terminator)
    print("dataset,#unclassified,(%),#contamination,(%),#partial-classified+missing necessary components,(%),#total,#reads-heatmap.csv,(%),#read-heatmap_opt.csv,(%),combinations(%),combinations(opt, %)")
    for dir in dirs:
        if dir.endswith('/'): 
            dir = dir[:-1]
        qry_total = f"{dir}/qry_total.csv"
        assert os.path.exists(qry_total)

        matrix = np.zeros((len(labels_peptide), len(prod_promoter_terminator)), dtype=np.int64)
        matrix_opt = np.zeros((len(labels_peptide), len(prod_promoter_terminator)), dtype=np.int64)
        num_succ = 0
        num_fail = 0
        num_contam = 0
        num_necessary_missed = 0
        num_optional_missed = 0
        num_succ_partial = 0
        num_total = 0

        count_terminator = dict.fromkeys(labels_terminator, 0)
        with open(qry_total, 'r') as fd:
            for line in fd:
                num_total += 1
                splited = line.strip().split(",")
                rid, entries = splited[0], splited[1:]
                name =  ",".join(entries)
                if name == "fail":
                    num_fail += 1
                    continue
                if name == "contamination":
                    num_contam += 1
                    continue
                promoter  = entries[1]
                peptide   = entries[2]
                terminator = entries[5]
                if promoter == '*' or peptide == '*' or terminator == '*':
                    num_necessary_missed += 1
                    continue
                count_terminator[terminator] += 1

                matrix_opt[rev_ridx[peptide], rev_cidx[(promoter, terminator)]] += 1
                num_succ_partial += 1
                if not any(name == '*' for name in entries):
                    matrix[rev_ridx[peptide], rev_cidx[(promoter, terminator)]] += 1
                    num_succ += 1

            fd.close()

        matrices.append(matrix)
        matrices_opt.append(matrix_opt)
     
        save_to_csv(matrix, labels_peptide, prod_promoter_terminator, f"{dir}/heatmap.csv")
        save_to_csv(matrix_opt, labels_peptide, prod_promoter_terminator, f"{dir}/heatmap_opt.csv")

        to_display = [
            f"{num_fail}", f"{round(100*num_fail/num_total,2)}%",
            f"{num_contam}", f"{round(100*num_contam/num_total,2)}%",
            f"{num_necessary_missed}", f"{round(100*num_necessary_missed/num_total,2)}%",
            f"{num_total}",
            f"{num_succ}", f"{round(100*num_succ/num_total,2)}%",
            f"{num_succ_partial}", f"{round(100*num_succ_partial/num_total,2)}%"
            f"{round(100*np.count_nonzero(matrix)/total_comb, 2)}%",
            f"{round(100*np.count_nonzero(matrix_opt)/total_comb, 2)}%",
        ]

        print(f"{dir.split('/')[-1]}," + ",".join(to_display))
    
    sys.exit(0)
    num_mat = len(dirs)
    div_round = lambda x: x // num_mat

    labels_yaxis = labels_peptide
    labels_xaxis = [f"({col[0]},{col[1]})" for col in prod_promoter_terminator]

    # -------------------------------------
    matrix_ave = np.zeros((len(labels_yaxis), len(labels_xaxis)), dtype=np.int64)
    for mat in matrices:
        matrix_ave = matrix_ave + mat
    matrix_ave = div_round(matrix_ave)
    save_to_csv(matrix_ave, labels_yaxis, prod_promoter_terminator, f"{out_prefix}_heatmap.csv")

    # -------------------------------------
    matrix_ave_opt = np.zeros((len(labels_yaxis), len(labels_xaxis)), dtype=np.int64)
    for mat in matrices_opt:
        matrix_ave_opt = matrix_ave_opt + mat
    matrix_ave_opt = div_round(matrix_ave_opt)
    save_to_csv(matrix_ave_opt, labels_yaxis, prod_promoter_terminator, f"{out_prefix}_heatmap_opt.csv")

    row_sub = 5
    col_sub = 15

    cmaps = [plt.cm.viridis, plt.cm.plasma, plt.cm.inferno, 
             plt.cm.Blues, plt.cm.Greens, plt.cm.Reds]
    for i, cmap in enumerate(cmaps):
        cmap.set_under(color='white')

        # -------------------------------------
        fig = plt.figure(figsize=(36,24))
        sns.heatmap(data=matrix_ave, xticklabels=labels_xaxis, yticklabels=labels_yaxis, cmap=cmap, vmin=1)
        plt.savefig(f"{out_prefix}_heatmap_c{i}.png", dpi=500)
        plt.close(fig)

        # -------------------------------------
        fig = plt.figure(figsize=(24,16))
        sns.heatmap(data=matrix_ave[:row_sub,:col_sub], xticklabels=labels_xaxis[:col_sub], yticklabels=labels_yaxis[:row_sub], cmap=cmap, vmin=1)
        plt.savefig(f"{out_prefix}_heatmap_c{i}_portion.png", dpi=500)
        plt.close(fig)

        # -------------------------------------
        fig = plt.figure(figsize=(36,24))
        sns.heatmap(data=matrix_ave_opt, xticklabels=labels_xaxis, yticklabels=labels_yaxis, cmap=cmap, vmin=1)
        plt.savefig(f"{out_prefix}_heatmap_opt_c{i}.png", dpi=500)
        plt.close(fig)

        # -------------------------------------
        fig = plt.figure(figsize=(24,16))
        sns.heatmap(data=matrix_ave_opt[:row_sub,:col_sub], xticklabels=labels_xaxis[:col_sub], yticklabels=labels_yaxis[:row_sub], cmap=cmap, vmin=1)
        plt.savefig(f"{out_prefix}_heatmap_opt_c{i}_portion.png", dpi=500)
        plt.close(fig)









