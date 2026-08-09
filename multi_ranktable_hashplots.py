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
        nargs="+",
        required=True,
        default=None,
        help="Path to input CSV file"
    )
    p.add_argument(
        "--remainder",
        action='store_true',
        help="Keep all values <= the lowest value in --thresholds (<=0.1)"
    )
    args = p.parse_args()

    all_dfs = []

    for filepath in args.ranktable:   # pass list of ranktable csvs to --ranktable-files
        df = pl.read_csv(filepath)
        df = df.with_columns(pl.lit(filepath).alias("source"))
        all_dfs.append(df)
    combined_df = pl.concat(all_dfs)
    # --- Plot distributions ---

#    sns.set(style="whitegrid")
#    fig, ax = plt.subplots(figsize=(6,6))

    sns.displot(combined_df, x="hashval", hue="source", kind="kde", fill=True, bw_adjust=.25)

    #ax.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.3, 1.0), title="Frequency Bins")

#    ax.set_aspect("equal")
#    ax.axis("off")
    plt.title("Distribution of pangenomic elements")
    plt.show()

if __name__ == "__main__":
    main()
