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
        "ranktable",
        type=str,
        nargs="+",
        help="Path(s) to input pangenome CSV file(s)"
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
    p.add_argument('--norm-y', choices=['multiplicity', 'percentage', 'both', 'none'],
        default="none",
        help='Normalization method for the y-axis.')
    p.add_argument(
        "--output",
        type=str,
        nargs='?',
        help='Output plot filename (defaulting to first input filename + .png)',
    )
    args = p.parse_args()

    bin_width = 1.0 / args.xbins

    if args.output is None:
        args.output = get_default_output(args.ranktable[0])

    all_df_list = []

    for table_path in args.ranktable:
        df = pl.read_csv(table_path)
        filename = Path(table_path).stem
        filename_with_max_abund = f"{filename} (Genome count: {df['abund'].max()})"

        if args.norm_y == 'multiplicity':
            df_count = (
                df
                .group_by('abund')
                .len()
                .sort('abund')
                .with_columns([
                    (pl.col('abund') * pl.col('len')).alias('c'),
                    (pl.col('abund') / pl.col('abund').max()).round(2).alias('freq_bin')
                ])
                .with_columns(
                    (pl.col('c') / pl.col('c').sum()).alias('norm_y')
                )
            )
            df_filter = (
                df_count
                .group_by('freq_bin')
                .agg(pl.sum('norm_y').alias('norm_y'))
                .sort('freq_bin')
                .filter(pl.col('freq_bin') >= args.filter)
            )

        elif args.norm_y == 'percentage':
            df_count = (
                df
                .group_by('abund')
                .len()
                .sort('abund')
                .with_columns([
                    (pl.col('abund') / pl.col('abund').max()).round(2).alias('freq_bin'),
                    (pl.col('len') / pl.col('len').sum()).alias('norm_y')
                ])
            )
            df_filter = (
                df_count
                .group_by('freq_bin')
                .agg(pl.sum('norm_y').alias('norm_y'))
                .sort('freq_bin')
                .filter(pl.col('freq_bin') >= args.filter)
            )

        elif args.norm_y == 'both':
            df_count = (
                df
                .group_by('abund')
                .len()
                .sort('abund')
                .with_columns([
                    (pl.col('abund') * pl.col('len')).alias('c'),
                    (
                        ((pl.col('abund') / pl.col('abund').max()) / bin_width)
                        .round() * bin_width
                    ).alias('freq_bin')
                ])
            )
            df_norm = (
                df_count
                .with_columns(
                    (pl.col('c') / pl.col('c').sum()).alias('norm_y')
                )
            )
            df_filter = (
                df_norm
                .group_by('freq_bin')
                .agg(pl.sum('norm_y').alias('norm_y'))
                .sort('freq_bin')
                .filter(pl.col('freq_bin') >= args.filter)
            )

        else:  # norm_y == 'none'
            df_filter = (
                df
                .with_columns(
                    (pl.col('abund') / pl.col('abund').max()).round(2).alias('freq_bin')
                )
                .group_by('freq_bin')
                .len()
                .sort('freq_bin')
                .filter(pl.col('freq_bin') >= args.filter)
                .rename({'len': 'norm_y'})
            )

        # this adds a label for seaborn hue and the legend
        df_filter = df_filter.with_columns(
            pl.lit(filename_with_max_abund).alias('label')
        )
        all_df_list.append(df_filter)

    combined_df = pl.concat(all_df_list)

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=combined_df,
        x="freq_bin",
        y="norm_y",
        hue="label",
        marker='o',
        alpha=.5
    )

    ticks = combined_df["freq_bin"].unique()
    tick_percents = [int(x * 100) for x in ticks]
    labels = [f"{p}%" if p % 5 == 0 else "" for p in tick_percents]

    for i, (p, label) in enumerate(zip(tick_percents, labels)):
        if label != "":
            continue
        has_close_label = any(
            abs(p - other_p) <= 3 and labels[j] != ""
            for j, other_p in enumerate(tick_percents) if j != i
        )
        if not has_close_label:
            labels[i] = f"{p}%"

    plt.legend(title="Pangenome Ranktable")
    plt.xticks(ticks=ticks, labels=labels)
    plt.xlabel('Percentage Quorum')
    plt.ylabel('Hash frequency')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to: {args.output}")

if __name__ == "__main__":
    main()
