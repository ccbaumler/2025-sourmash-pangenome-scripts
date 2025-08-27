#! /usr/bin/env python

import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main():
    p = argparse.ArgumentParser(description="Sourmash hash distribution using Polars")
    p.add_argument(
        "--ranktable",
        type=str,
        default=None,
        help="Path to input pangenome CSV file"
        )
    p.add_argument(
        "--filter",
        type=float,
        default=0.0,
        help="Remove all hashvals from the plot >= filter value"
        )
    args=p.parse_args()

    df = pl.read_csv(args.ranktable)
    df = df.filter(pl.col("freq") > args.filter)

    plt.figure(figsize=(8, 6))
    sns.histplot(
        df['freq'],
        bins=[x / 100 for x in range(0, 105, 5)],  # 5% bin width
        color='skyblue',
        edgecolor='black'
    )
    
    plt.xlabel('Frequency (Normalized Percent)')
    plt.ylabel('Count')
    plt.title('Histogram of Frequency (5% bins)')
    plt.xticks(
    [x / 100 for x in range(0, 105, 5)],     # tick positions (e.g., 0.0, 0.05, ..., 1.0)
    [f"{x}%" for x in range(0, 105, 5)]      # tick labels (e.g., "0%", "5%", ..., "100%")
    )
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
