import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse
sys.stdout.reconfigure(encoding='utf-8')

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
        default=0.9999999990,
        help="Remove all hashvals from the plot >= filter value"
        )
    args=p.parse_args()

    df = pl.read_csv(args.ranktable)
    
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
        y='len',
        marker='o'
    )
    
    plt.xticks(
        ticks=df_filter["freq_bin"],
        labels=[
            f"{int(x * 100)}%" if int(x * 100) % 2 == 0 else ""
            for x in df_filter["freq_bin"]
        ]
    )
    
    plt.xlabel('Percentage Quorum')
    plt.ylabel('Hash frequency (Est. Core)')
    plt.title('Line Plot of Frequency (90% - 100%)')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
