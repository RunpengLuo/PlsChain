import sys
import os
import importlib

has_pandas = True

try:
    import pandas as pd
except ModuleNotFoundError: # python v >= 3.6
    # Error handling
    has_pandas = False

if __name__ == "__main__":
    if len(sys.argv) != 2: 
        print(f"{sys.argv[0]} <query_dir>")
        sys.exit(1)
    _, out_dir = sys.argv
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]    
    qry_total = out_dir + "/qry_total.csv"

    if has_pandas:
        df_total = pd.read_csv(qry_total, header=None)
        counts = df_total.groupby(list(df_total.columns)[1:])[0].count() # group by combination type
        counts.to_csv(out_dir + "/qry_total.group.csv", header=False)
    else:
        print("Pandas does not exists, try traditional dict groupby")
        gb_dict = {}
        with open(qry_total, "r") as fd_total:
            for line in fd_total:
                splited = line.strip().split(",")
                rid, names = splited[0], ",".join(splited[1:])
                if names not in gb_dict:
                    gb_dict[names] = 0
                gb_dict[names] += 1
            fd_total.close()
        with open(out_dir + "/qry_total.group.csv", "w") as fd_group:
            for name, count in gb_dict.items():
                fd_group.write(f"{name},{count}\n")
            fd_group.close()
    print("Done")
