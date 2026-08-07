#! /usr/bin/env python 

import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse
from pathlib import Path
from scipy.signal import find_peaks

sys.stdout.reconfigure(encoding='utf-8') 

def get_default_output(input_file):
    fp = Path(input_file)
    base_path = fp.with_name(fp.stem).with_suffix('.quorum.png')
    print(base_path.name)
    return f"{base_path.name}"

def truncate_front(text, max_len=20):
    if len(text) <= max_len:
        return text
    return f"...{text[-max_len:]}"

def truncate_middle(text, max_len=20):
    if len(text) <= max_len:
        return text
    # Calculate how many characters to keep on each side of the "..."
    front_len = (max_len - 3) // 2 + (max_len - 3) % 2
    back_len = (max_len - 3) // 2
    return f"{text[:front_len]}...{text[-back_len:]}"


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
    p.add_argument("--abund",
                   type=int,
                   default=1,
                   help="Number of genomes that must exist within the pangenome. I.e. the max_abund of the ranktable.")
    p.add_argument("--peaks", type=float, default=0.01,
                   help="Min peak prominence relative to neighbors to be highlighted.")
    p.add_argument("--text-len", type=int, default=20,
                   help="Max length a legend text label may be")
    args = p.parse_args()

    bin_width = 1.0 / args.xbins

    if args.output is None:
        args.output = get_default_output(args.ranktable[0])

    all_df_list = []

    for table_path in args.ranktable:
        df = pl.read_csv(table_path)
        max_abund = df['max_abund'].max()
        if max_abund < args.abund:
            continue

        filename = Path(table_path).stem
        filename_with_max_abund = f"{filename} (Genome count: {max_abund})"

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
        #df_filter = df_filter.with_columns(
        #    pl.lit(filename_with_max_abund).alias('label')
        #)
        y_vals = df_filter['norm_y'].to_numpy()
        
        # 'prominence' measures height of peak relative to surrounding low points
        peak_indices, _ = find_peaks(
            y_vals, 
            prominence=args.peaks,
            distance=2  # Optional: minimum spacing (in index steps) between peaks
        )

        # Flag peak rows in Polars DataFrame
        is_peak = [i in peak_indices for i in range(len(y_vals))]
        df_filter = df_filter.with_columns(
            pl.Series('is_peak', is_peak),
            pl.lit(filename_with_max_abund).alias('label')
        )

        all_df_list.append(df_filter)

    combined_df = pl.concat(all_df_list)
    print(combined_df)
    combined_df = pl.concat(all_df_list)
    
    # Convert Polars DataFrame to Pandas for Seaborn compatibility
    combined_pd = combined_df.to_pandas()

    plt.figure(figsize=(10, 6))

    # Main line plot
    sns.lineplot(
        data=combined_pd,
        x="freq_bin",
        y="norm_y",
        hue="label",
        marker='o',
        alpha=0.5
    )

    # -------------------------------------------------------------
    # PLOT & ANNOTATE PEAKS
    # -------------------------------------------------------------
    peaks_pd = combined_pd[combined_pd['is_peak']]

    if not peaks_pd.empty:
        # Overlay larger highlighted dots on detected peaks
        sns.scatterplot(
            data=peaks_pd,
            x="freq_bin",
            y="norm_y",
            hue="label",
            s=120,            # Marker size
            edgecolor="red",  # Red ring around peak markers
            linewidth=1.5,
            legend=False      # Avoid duplicate legend entries
        )

        # Add text labels displaying the peak y-value above the peak
        for _, row in peaks_pd.iterrows():
            plt.annotate(
                f"{row['norm_y']:.2f}",
                (row['freq_bin'], row['norm_y']),
                textcoords="offset points",
                xytext=(0, 8),
                ha='center',
                fontsize=8,
                fontweight='bold'
            )

    ticks = combined_pd["freq_bin"].unique()
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

    ax = plt.gca()
    handles, legend_labels = ax.get_legend_handles_labels()

    #truncated_labels = [truncate_middle(lbl, max_len=55) for lbl in legend_labels]
    truncated_labels = [truncate_front(lbl, max_len=args.text_len) for lbl in legend_labels]

    plt.legend(
        handles, 
        truncated_labels,
        title="Pangenome Ranktable",
        loc='upper center',
        shadow=True,
        fontsize=8,
        title_fontsize=10
    )

    plt.xticks(ticks=ticks, labels=labels)
    plt.xlabel('Percentage Quorum')
    plt.ylabel('Hash frequency')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to: {args.output}")
if __name__ == "__main__":
    main()
