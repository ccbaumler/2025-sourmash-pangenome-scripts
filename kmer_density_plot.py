#! /usr/bin/env python

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import seaborn as sns
import sys
import argparse
import math

sys.stdout.reconfigure(encoding='utf-8')

def main():
    p = argparse.ArgumentParser(description="Count frequency bins using Polars.")
    p.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[1, 0.9, 0.1, 0.0],
        help="List of thresholds to create bins (e.g. --thresholds 0.9 0.1 0.0)"
    )
    p.add_argument(
        "--ranktable",
        type=str,
        default=None,
        help="Path to input CSV file"
    )
    p.add_argument(
        "--remainder",
        action='store_true',
        help="Keep all values <= the lowest value in --thresholds (<=0.1)"
    )
    args = p.parse_args()

    df = pl.read_csv(args.ranktable)
    #df = df.cast({"hashval": pl.String})
    thresholds = sorted(args.thresholds, reverse=True)
    ranges = list(zip(thresholds[:-1], thresholds[1:]))

    expr = None
    for upper, lower in ranges:
        label = f"{upper} > freq > {lower}"
        condition = (pl.col("freq") <= upper) & (pl.col("freq") > lower)
        if expr is None:
            expr = pl.when(condition).then(pl.lit(label))
        else:
            expr = expr.when(condition).then(pl.lit(label))

    if args.remainder:
        expr = expr.otherwise(pl.lit(f"{thresholds[-1]} >="))
        df = df.with_columns(expr.alias("freq_bin"))
    else:
        df = df.with_columns(expr.alias("freq_bin"))
        df = df.filter(pl.col("freq_bin").is_not_null())

    print(df)
    bin_counts = (
        df.group_by("freq_bin")
        .len()
        .sort("freq_bin", descending=True)
    )
    
    print(bin_counts)

    # --- Plot distributions ---

#    sns.set(style="whitegrid")
#    fig, ax = plt.subplots(figsize=(6,6))

    sns.displot(df, x="hashval", hue="freq_bin", kind="kde", fill=True, bw_adjust=.25)

    #ax.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.3, 1.0), title="Frequency Bins")

#    ax.set_aspect("equal")
#    ax.axis("off")
    plt.title("Distribution of pangenomic elements")
    plt.show()

    plt.clf()

    sns.displot(df, x="hashval", hue="freq_bin", kde=True, fill=True, bins=1000)
    plt.title("Distribution of pangenomic elements")
    plt.show()

    plt.clf()

    sns.displot(df, x="hashval", hue="freq_bin", kind='ecdf')
    plt.title("Distribution of pangenomic elements")
    plt.show()



if __name__ == "__main__":
    main()
