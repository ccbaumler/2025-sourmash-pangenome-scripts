#! /usr/bin/env python

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import seaborn as sns
import sys
import argparse
import math
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8')

def get_default_output(input_file):
    fp = Path(input_file)
    base_path = fp.with_name(fp.stem).with_suffix('.circles.png')
    print(base_path.name)
    return f"{base_path.name}"

def main():
    p = argparse.ArgumentParser(description="Count frequency bins using Polars.")
    p.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[1, 0.95, 0.1, 0.0],
        help="List of thresholds to create bins (e.g. --thresholds 0.95 0.1 0.0)"
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
    p.add_argument(
        "-o",
        "--output",
        type=str,
        nargs='?',
        help='Output plot filename (default to ranktable filename + .png'
    )
    p.add_argument(
        "--title",
        type=str,
        default="Concentric Circles of pangenomic elements",
        help="Add a custom title to the plot"
    )
    args = p.parse_args()

    if args.output is None:
        args.output = get_default_output(args.ranktable)

    df = pl.read_csv(args.ranktable)
    thresholds = sorted(args.thresholds, reverse=True)
    ranges = list(zip(thresholds[:-1], thresholds[1:]))

    # create freq_bin labels
    bins = [
        f"{upper} > freq > {lower}"
        for upper, lower in ranges
    ]
    
    if args.remainder:
        bins.append(f"{thresholds[-1]} >=")
    
    # build freq_bin assignment
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
    genome_count = df.select(pl.col("max_abund")).to_series().first()
    
    # count observed bins
    bin_counts = (
        df.group_by("freq_bin")
          .len()
          .rename({"len": "count"})
    )
    
    # create all bins with zero counts
    all_bins = pl.DataFrame({
        "freq_bin": bins,
        "count": [0] * len(bins)
    })
    
    # fill in observed counts
    bin_counts = (
        all_bins.join(bin_counts, on="freq_bin", how="left")
        .with_columns(
            pl.coalesce(["count_right", "count"])
            .alias("count")
        )
        .select(["freq_bin", "count"])
    )
    
    # remove remainder bin if requested
    if not args.remainder:
        bin_counts = bin_counts.filter(
            pl.col("freq_bin") != f"{thresholds[-1]} >="
        )
    
    print(bin_counts)
    # --- Plot concentric circles ---

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(6,6))

    center = (0, 0)

    # Extract bins and counts as lists
    bins = bin_counts["freq_bin"].to_list()
    counts = bin_counts["count"].to_list()

    radii = [math.sqrt(c) for c in counts]

    cmap = plt.get_cmap("viridis", len(bins))
    patches = []
    white_circle = plt.Circle(
        (0, 0),
        radius=sum(radii)*0.01,
        color="white" if counts[0] == 0 else cmap(0),
        zorder=10
    )
    ax.add_artist(white_circle)

    prev_radius = sum(radii)*0.01

    #create a set of wedges centered at origin
    for i, (label, count, radius) in enumerate(zip(bins, counts, radii)):
        color = "white" if count == 0 else cmap(i)
        prev_radius += radius
        wedge = Wedge(center,
                      prev_radius,
                      0,
                      360,
                      width=prev_radius,
                      facecolor=color,
                      edgecolor="black" if count == 0 else None,
                      zorder=-i)
        ax.add_patch(wedge)
        wedge.set_label(f"{label} -- {count} hashes")
        patches.append(wedge)


    legend = ax.legend(handles=patches, loc='lower left', title=f"Frequency Bins for {genome_count} genomes")
    legend.get_title().set_fontweight("bold")

    ax.set_xlim(-sum(radii)*1.1, sum(radii)*1.1)
    ax.set_ylim(-sum(radii)*1.1, sum(radii)*1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(args.output, dpi=100)

if __name__ == "__main__":
    main()
