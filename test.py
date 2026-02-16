import polars as pl
import sys

sys.stdout.reconfigure(encoding='utf-8')


df = pl.DataFrame({"freq": [1.0, 0.95, 0.92, 0.85, 0.5, 0.2, 0.05]})
thresholds = [1, 0.9, 0.1, 0.0]
thresholds = sorted(thresholds, reverse=True)
ranges = list(zip(thresholds[:-1], thresholds[1:]))

expr = None
for upper, lower in ranges:
    label = f"{upper} > freq > {lower}"
    condition = (pl.col("freq") <= upper) & (pl.col("freq") > lower)
    if expr is None:
        expr = pl.when(condition).then(pl.lit(label))
    else:
        expr = expr.when(condition).then(pl.lit(label))

expr = expr.otherwise(pl.lit(f"<= {thresholds[-1]}"))

df = df.with_columns(expr.alias("freq_bin"))

bin_counts = df.group_by("freq_bin").len()
print(bin_counts)

