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
    args = p.parse_args()

    if args.output is None:
        args.output = get_default_output(args.ranktable)

    df = pl.read_csv(args.ranktable)
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

    # --- Plot concentric circles ---

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(6,6))

    center = (0, 0)

    # Extract bins and counts as lists
    bins = bin_counts["freq_bin"].to_list()
    counts = bin_counts["len"].to_list()

    radii = [math.sqrt(c) for c in counts]

    # What seaborn colormap is best?
    cmap = plt.cm.get_cmap("viridis", len(bins))
    patches = []
    white_circle = plt.Circle((0, 0), radius=sum(radii)*0.01, color='white', zorder=10)
    ax.add_artist(white_circle)

    prev_radius = sum(radii)*0.01

    #create a set of wedges centered at origin
    for i, (label, count, radius) in enumerate(zip(bins, counts, radii)):
        color = cmap(i)
        prev_radius += radius
        wedge = Wedge(center, prev_radius, 0, 360, width=prev_radius, facecolor=color, zorder=-i)
        ax.add_patch(wedge)
        wedge.set_label(f"{label} -- {count} hashes")
        patches.append(wedge)

    ax.legend(handles=patches, loc='lower left', title="Frequency Bins")

    ax.set_xlim(-sum(radii)*1.1, sum(radii)*1.1)
    ax.set_ylim(-sum(radii)*1.1, sum(radii)*1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.title("Concentric Circles of pangenomic elements")
    plt.tight_layout()
    plt.savefig(args.output, dpi=100)

if __name__ == "__main__":
    main()
