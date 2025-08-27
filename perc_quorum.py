#! /usr/bin/env python

import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')

def get_default_output(input_file):
    fp = Path(input_file)
    base_path = fp.with_name(fp.stem).with_suffix('.quorum.png')
    print(base_path.name)
    return f"{base_path.name}"

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
    p.add_argument(
        "--xbins",
        type=int,
        default=100,
        help="Number of bins for x-axis normalization (E.g. 5, 20, 100)"
    )
    # From https://gitlab.ub.uni-bielefeld.de/gi/pangrowth/-/blob/main/scripts/plot_hist.py?ref_type=heads
    p.add_argument('--norm-y', choices=['multiplicity', 'percentage', 'both', 'none'], 
                        default="none", 
                        help='Normalize the y-axis values: '
                         '"multiplicity" multiply each bar by the number of genomes it appears (unique items are multiplied by one, items appearing in two genomes are multiplied by two, and so on..), '
                         '"percentage" normalizes by converting counts to a percentage of the total, '
                         '"both" applies both multiplicity and percentage normalizations in sequence, '
                         '"none" applies no normalization.')
    p.add_argument(
        "--output",
        type=str,
        nargs='?',
        help='Output plot filename (defaulting to ranktable filename + .png',
    )
    args=p.parse_args()

    bin_width = 1.0 / args.xbins

    if args.output is None:
        args.output = get_default_output(args.ranktable)

    df = pl.read_csv(args.ranktable)

    if args.norm_y == 'multiplicity':
        print(df)
        df_count = (
            df
            .group_by('abund')
            .len()
            .sort("abund")
            .with_columns(
                (pl.col("abund") * pl.col("len")).alias("c"),
                (pl.col('abund') / pl.col('abund').max()).round(2).alias('freq'),
            )
        )
        print(df_count)
        (pl.col('c') / pl.col('c').sum()).alias('y_per')
    if args.norm_y == 'percentage':
        print(df)
        df_count = (
            df
            .group_by('abund')
            .len()
            .sort("abund")
            .with_columns(
                (pl.col('abund') / pl.col('abund').max()).round(2).alias('freq'),
                (pl.col('len') / pl.col('len').sum()).alias('y_perc')
            )
        )
        print(df_count)
    if args.norm_y == 'both':
        print(df)
        df_count = (
            df
            .group_by('abund')
            .len()
            .sort("abund")
            .with_columns(
                (pl.col("abund") * pl.col("len")).alias("c"),
            (
                ((pl.col('abund') / pl.col('abund').max()) / bin_width)
                .round() * bin_width
            ).alias("freq_bin")
           )
        )
        print(df_count)
        df_norm = (
            df_count
            .with_columns(
                (pl.col('c') / pl.col('c').sum()).alias('y_perc')
            )
        )
        df_filter = (
            df_norm
            .group_by("freq_bin")  # Group by the rounded freq
            .agg(
                pl.sum("y_perc").alias("norm_y")              # Count the rows in each group
            )
            .sort("freq_bin")     # Sort by freq_bin (optional, for clean output)
            .filter(pl.col("freq_bin") >= args.filter)
        )
        print(df_filter)
    else:
        df_filter = (
            df.with_columns(
                pl.col("freq").round(2).alias("freq_bin")  # Round and create new column
            )
            .group_by("freq_bin")  # Group by the rounded freq
            .len()              # Count the rows in each group
            .sort("freq_bin")     # Sort by freq_bin (optional, for clean output)
            .filter(pl.col("freq_bin") >= args.filter)
        )
    
    print(df_filter)
    
    
    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df_filter,
        x="freq_bin",  # X = position
        y='len' if not args.norm_y else 'norm_y',
        marker='o'
    )
    
    ticks = df_filter["freq_bin"].to_list()
    tick_percents = [int(x * 100) for x in ticks]
    
    labels = [
        f"{p}%" if p % 5 == 0 else ""
        for p in tick_percents
    ]
    
    for i, (p, label) in enumerate(zip(tick_percents, labels)):
        if label != "":
            continue
    
        has_close_label = any(
            abs(p - other_p) <= 3 and labels[j] != ""
            for j, other_p in enumerate(tick_percents) if j != i
        )
    
        if not has_close_label:
            labels[i] = f"{p}%"
    
    plt.xticks(ticks=ticks, labels=labels)
   
    plt.xlabel('Percentage Quorum')
    plt.ylabel('Hash frequency (Est. Core)')
    plt.title('Line Plot of Frequency (90% - 100%)')
    plt.tight_layout()
    plt.savefig(args.output, dpi=100)

if __name__ == "__main__":
    main()
