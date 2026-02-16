#! /usr/bin/env python

import polars as pl
import numpy as np
from scipy.stats import linregress

# Example data: n = total tokens, V = vocab size
data = [
    {"n": 1000, "V": 314},
    {"n": 5000, "V": 918},
    {"n": 10000, "V": 1325},
    {"n": 50000, "V": 3150},
    {"n": 100000, "V": 5000},
]

df = pl.DataFrame(data)

# Log-transform
df = df.with_columns([
    pl.col("n").log10().alias("log_n"),
    pl.col("V").log10().alias("log_V")
])

# Linear regression in log-log space
slope, intercept, r_value, _, _ = linregress(df["log_n"], df["log_V"])

beta = slope
k = 10 ** intercept  # because logs are base-10

print(f"Estimated β: {beta:.4f}")
print(f"Estimated k: {k:.2f}")

