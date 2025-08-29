#! /usr/bin/env python

import csv
import gzip
import matplotlib.pyplot as plt
from upsetplot import plot
import argparse
from pathlib import Path
import pandas as pd

def get_default_output(input_file, plottype='upset'):
    fp = Path(input_file)
    base_path = fp.with_name(fp.stem).with_suffix(f'.{plottype}.png')
    print(base_path.name)
    return f"{base_path.name}"

def main():
    p = argparse.ArgumentParser()

    p.add_argument('ranktable', type=str, nargs='+', help='The pangenome ranktable(s) for processing.')
    p.add_argument('-f','--filter',type=float,default=0.0,help="Remove all hashvals from the plot <= filter value")
    p.add_argument('-o','--output',type=str,nargs='?',help='Output plot filename (defaulting to first input filename + .png)')

    args = p.parse_args()
    assert len(args.ranktable) > 1
    sets = []
    set_names = []

    if args.output is None:
        args.output = get_default_output(args.ranktable[0])

    for i,table in enumerate(args.ranktable, start=1):
        filepath = Path(table)
        set_names.append(f'{filepath.stem}_hashvals')
        hashval_set = set()
        
        open_func = gzip.open if filepath.suffix == '.gz' else open
        
        with open_func(table, mode='rt', newline='') as fp:
            reader = csv.DictReader(fp)
            
            for row in reader:
                freq = float(row['freq'])

                if freq <= args.filter:
                    break
                
                hashval = row['hashval']
                hashval_set.add(hashval)

        sets.append(hashval_set)

    print(f"Found {len(sets)} sets.")
    for name, s in zip(set_names, sets):
        print(f"{name}: {len(s)} hashvals")

    all_elems = list(set().union(*sets))
    df = pd.DataFrame([[e in st for st in sets] for e in all_elems], columns = set_names)

    df_up = df.groupby(set_names).size()

    plot(df_up, orientation='horizontal')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to: {args.output}")


if __name__ == '__main__':
    main()
