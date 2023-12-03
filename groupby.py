import sys
import os
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 2: 
        print(f"{sys.argv[0]} <query_dir>")
        sys.exit(1)
    _, out_dir = sys.argv
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]

    
    qry_total = out_dir + "/qry_total.csv"
    df_total = pd.read_csv(qry_total, header=None)
    counts = df_total.groupby(list(df_total.columns)[1:])[0].count() # group by combination type
    counts.to_csv(out_dir + "/qry_total.group.csv", header=False)
    print("Done")
