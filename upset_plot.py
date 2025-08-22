import polars as pl
import pandas as pd
from upsetplot import plot
import matplotlib.pyplot as plt

# Define sets
set1 = {0, 1, 2, 3, 4, 5}
set2 = {3, 4, 5, 6, 10}
set3 = {0, 5, 6, 7, 8, 9}

set_names = ['set1', 'set2', 'set3']
all_elems = sorted(set1 | set2 | set3)

# Build Polars DataFrame showing membership
df = pl.DataFrame([
    [e in set1, e in set2, e in set3]
    for e in all_elems
], columns=set_names)

# Group by membership and count
df_up = df.groupby(set_names).len().rename({"len": "count"})

# Convert to Pandas for upsetplot
df_pd = df_up.to_pandas()

# Set MultiIndex from boolean columns
df_pd.set_index(set_names, inplace=True)

# Plot
plot(df_pd["count"], orientation='horizontal')
plt.tight_layout()
plt.show()
