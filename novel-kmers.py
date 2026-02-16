#! /usr/bin/env python

import polars as pl
from itertools import accumulate
import argparse
import matplotlib.pyplot as plt
import os
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.stats import describe
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import csv
import numpy as np
import random

GLOBAL_COLUMN_ARRAYS = None
GLOBAL_MASK = None
GLOBAL_VERBOSE = False
GLOBAL_PERM_NOVEL = None
GLOBAL_PERM_COUNTS = None


def init_worker(column_arrays, mask, verbose):
    global GLOBAL_COLUMN_ARRAYS, GLOBAL_MASK
    GLOBAL_COLUMN_ARRAYS = column_arrays
    GLOBAL_MASK = mask
    GLOBAL_VERBOSE = verbose

def init_saturation_worker(permutation_novel_counts, permutation_counts):
    global GLOBAL_PERM_NOVEL, GLOBAL_PERM_COUNTS
    GLOBAL_PERM_NOVEL = permutation_novel_counts
    GLOBAL_PERM_COUNTS = permutation_counts

def run_permutation(i):
    rng = random.Random(i)

    shuffled_indices = list(range(len(GLOBAL_COLUMN_ARRAYS)))
    rng.shuffle(shuffled_indices)

    shuffled_columns = [GLOBAL_COLUMN_ARRAYS[j] for j in shuffled_indices]
    shuffled_df = pl.DataFrame(shuffled_columns)

    if GLOBAL_VERBOSE:
        print(f"[Worker {os.getpid()}] Running permutation {i}")
    return compute_novel_counts(shuffled_df, GLOBAL_MASK, GLOBAL_VERBOSE)

def run_saturation(i):
    print(f"[Worker {os.getpid()}] Starting run_saturation({i})", flush=True)

    novel = GLOBAL_PERM_NOVEL[i]
    counts = GLOBAL_PERM_COUNTS[i]

    result = pangenome_saturation(novel, counts)

    print(f"[Worker {os.getpid()}] Finished run_saturation({i})", flush=True)
    return result

def process_permutation(perm, do_sort=False):
    values = [v[1] for v in perm.values()]
    if do_sort:
        values = sorted(values, reverse=True)
    acc = list(accumulate(values))
    return acc

def compute_novel_counts(df, mask, verbose=False):
    novel_counts = {}
    tots = 0
    for col in df.columns:
        col_values = df.get_column(col)
        tots += col_values.sum()
        count = col_values.filter(mask).sum()
        novel_counts[col] = [tots, count]
        if verbose: print(f"{col}: Total {tots} and novel {count} true values")
        mask = mask & ~col_values  # exclude True value
    return novel_counts

def vocab_saturation(n, V_max, alpha):
    return V_max * (1 - np.exp(-alpha * n))

def object_fun(p, x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    y_hat = p[0] * x**p[1] + p[2]
    J = np.sum((y - y_hat)**2)
    return J

def fitting(x_val, params):
    return params[0] * x_val**params[1] + params[2]

def plot_totals(data, lineage, path, ext, h=5, w=12, line_width=1, marker='none',
                yaxis='kmers', scaled=1000, permutations=0, alpha=0.3, plot_original=True, saturation=None):
    fpath = os.path.join(path + ext)
    plt.figure(figsize=(w,h))

    y_data = data if yaxis == "hashes" else [i * scaled for i in data]
    x = list(range(len(data)))

    if permutations:
        if yaxis == "kmers":
            permutations = [[i * scaled for i in p] for p in permutations]

        permutations_np = np.array(permutations)
        min_vals = permutations_np.min(axis=0)
        max_vals = permutations_np.max(axis=0)
        mean_vals = permutations_np.mean(axis=0)

        x = list(range(len(data)))
        plt.fill_between(x, min_vals, max_vals, color='gray', alpha=alpha, label=f'{len(permutations)} permutations range')
        plt.plot(x, mean_vals, linestyle='--', color='black', label='Permutation mean')

    if plot_original:
        plt.plot(x, y_data, linewidth=line_width, marker=marker, label='Original')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel('New Genomes')

    if yaxis == "hashes":
        plt.ylabel('Total novel hashes')
        plt.title(f'Total novel pangenome hashes per genome from Lineage: {lineage}')
    else:
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
        plt.ylabel('Total novel k-mers')
        plt.title(f'Total novel pangenome k-mers per genome from Lineage: {lineage}')

    if saturation:
        plt.text(
            0.95, 0.05,  # X, Y in axes fraction coords
            f"Mean: {saturation:,.0f}",
            horizontalalignment='right',
            verticalalignment='bottom',
            transform=plt.gca().transAxes,
            fontsize=10,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )

    plt.legend()
    plt.savefig(fpath, bbox_inches="tight", transparent=True)


def pangenome_saturation(novel_counts, total_counts, debug=True):
    if debug:
        print(f"[{os.getpid()}] Starting pangenome_saturation", flush=True)
        print(f"  novel_counts sample: {list(novel_counts.items())[:5]}", flush=True)
        print(f"  total_counts sample: {total_counts[:5]}", flush=True)
    df = pl.DataFrame(novel_counts).transpose(include_header=False)
    df = df.rename({"column_0": "n", "column_1": "V"})

    df = df.with_columns([
        pl.col("n").cast(pl.Float64).log10().alias("log_n"),
        pl.col("V").cast(pl.Float64).log10().alias("log_V"),
    ])
    df = df.drop_nulls().drop_nans()
    df = df.filter((~pl.col("log_n").is_infinite()) & (~pl.col("log_V").is_infinite()))
    print(df)
    slope, intercept, r_value, _, _ = linregress(df["log_n"], df["log_V"])
    
    beta = slope
    k = 10 ** intercept  # because logs are base-10
    n = df.select(pl.max('n'))
    V_n = max(total_counts)

#    print(f"Estimated β: {beta:.4f}")
#    print(f"Estimated k: {k:.2f}")
#    print(k*n**beta)
#    print(-(np.log(max(total_counts)/k)/np.log(len(total_counts))))

    n_array = df["n"].to_numpy()
    V_array = df["V"].to_numpy()
    
    params, _ = curve_fit(vocab_saturation, n_array, V_array, p0=[max(V_array)*1.5, 1e-5], bounds=(0, np.inf))
    
    V_max_fit, alpha_fit = params
    
#    print(f"Estimated V_max: {V_max_fit:.2f}")
#    print(f"Estimated alpha: {alpha_fit:.8f}")

    x = pl.DataFrame({
            "nth.genome": list(range(1, len(total_counts) + 1)),
            "base.pairs": total_counts
        })
    
    x_norm = x['base.pairs'] / x['base.pairs'].max()
    x_vals = x['nth.genome'].to_numpy() if isinstance(x['nth.genome'], pl.Series) else np.asarray(x['nth.genome'])
    y_vals = x_norm.to_numpy() if isinstance(x_norm, pl.Series) else np.asarray(x_norm)

    initial_params = [0, 0, 0]
    bounds = [(-100, 100), (-100, 100), (-100, 100)]  # add bound for p[3] too
    result = minimize(object_fun, initial_params, args=(x_vals, y_vals), method='L-BFGS-B', bounds=bounds)
    
#    print(result)
    
    z = x['base.pairs'].max()
    m = z / 1e9
    params = result.x
    
#    print(x['base.pairs'].min())
#    print(x['base.pairs'].max())
#    print(params[0] * z)
#    print(params[1] * z)
#    print(params[2] * z)
   
    n = x['nth.genome'].max()
#    print("Predicted increase in pangenome size from last genome", z * (fitting(n, params) - fitting(n - 1, params)))
#    print("Increase from first genomes to second genome", z * (fitting(2, params) - fitting(1, params)))

    a = result.x[0]    # fitted p[0]
    gamma = result.x[1]  # fitted p[1]
    c = result.x[2]    # fitted p[2]
    z = x['base.pairs'].max()  # original max base pair count (for scaling back)
    
    # Define the fitted function
    def f(n):
        return a * n**gamma + c
    
    # Set the minimum change threshold (e.g., unicity of 3?)
    threshold_bp = 3
    threshold_norm = threshold_bp / z  # normalize since model was fit on normalized data
    
    max_n = 100_000  # or whatever is reasonable
    n = 2
    while n < max_n:
        delta = f(n) - f(n - 1)
        if delta * z < threshold_bp:
            break
        n += 1
    
    if n >= max_n:
        print(f"[{os.getpid()}] WARNING: Saturation loop reached max_n={max_n} without convergence", flush=True)
        return -1  # or a fallback value   

    print(f"Estimated number of genomes needed for pangenome saturation: {n}")
    return n

def plot_novels(data, lineage, path, ext, h=5, w=12, line_width=1, marker='none', yaxis='kmers', scaled=1000):
    fpath = os.path.join(path + ext)
    plt.figure(figsize=(w,h))
    if yaxis=="hashes":
        plt.plot(range(len(data)), data, linewidth=line_width, marker=marker) #marker='o' for circles
    elif yaxis=="kmers":
        plt.plot(range(len(data)), [i * scaled for i in data], linewidth=line_width, marker=marker)
    else:
        print("Invalid value for yaxis. Defaulting to 'kmers' and 'scaled=1000'.")
        plt.plot(range(len(data)), [i * scaled for i in data], linewidth=line_width, marker=marker)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel('New Genomes')
    if yaxis=="hashes":
        plt.ylabel('Novel hashes')
        plt.title(f'Novel pangenome hashes per genome for Lineage: {lineage}')
    elif yaxis=="kmers":
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
        plt.ylabel('Novel k-mers')
        plt.title(f'Novel pangenome k-mers per genome for Lineage: {lineage}')
    else:
        #print("Invalid value for yaxis. Defaulting to 'kmers'.")
        plt.ylabel('Novel k-mers')
        plt.title(f'Novel pangenome k-mers per genome for Lineage: {lineage}')
    plt.savefig(fpath, bbox_inches="tight", transparent=True)

def main():
    p = argparse.ArgumentParser(prog='plot_diff', description='Plot the novel k-mers/hashes of a sourmash pangenome.')

    p.add_argument('plot', choices=['novel','total'], help='Choose between "novel" or "total" plot type')
    p.add_argument('data', metavar='CSV_DATA', help='This is the processed CSV data from the pangenome CSV file')
    p.add_argument('--n-permutations', type=int)
    p.add_argument('--verbose', action='store_true')
    p.add_argument('-l', '--lineage', metavar='LINEAGE_STRING', help='This is the lineage name to isolate from the dataset')
    p.add_argument('-o', '--output', help='The output path for the plot image file')
    p.add_argument('-H', '--height', type=int, default=5, help='The plot image height')
    p.add_argument('-w', '--width', type=int, default=12, help='The plot image width')
    p.add_argument('-lw', '--line-width', type=int, default=1, help='The line width of the plotted lines')
    p.add_argument('-m', '--marker', type=str, default='_', help='The marker decoration to place on top of each line (e.g. "o", "_")')
    p.add_argument('-y', '--yaxis', type=str, default='kmers', choices=['kmers','kmer','hashes','hash'], help='Choose between "hashes" or "kmers" for the y-axis')
    p.add_argument('-s', '--scaled', type=int, default=1000, help='The scaled value of the sourmash database')
    p.add_argument('--n-cols', type=int)
    p.add_argument('--sorted', action='store_true')
    p.add_argument('--num-threads', type=int)

    args = p.parse_args()

    print(f"Loading {args.data} as polars dataframe...")

    if args.n_cols:
        print(f'    Loading only {args.n_cols} columns...')
        with open(args.data, 'r') as f:
            header = next(csv.reader(f))
        filtered = [col for col in header if col != 'hashval']
        selected_cols = filtered[:args.n_cols]

        df_lazy = pl.read_csv(args.data, columns=selected_cols).lazy()
    else:
        df_lazy = pl.scan_csv(args.data).drop("hashval")

    print(df_lazy.explain(optimized=True))

    boolean_scan = df_lazy.with_columns([
        pl.col(col).cast(pl.Boolean) for col in df_lazy.collect_schema().names()
    ])

    df = boolean_scan.collect()

    novel_counts = {}
    tots = 0

    mask = pl.Series([True] * df.height)

    # === ORIGINAL ORDER ===
    novel_counts = compute_novel_counts(df, mask, verbose=False)

    # === PERMUTATIONS ===
    permutation_novel_counts = []
    columns = df.columns

#    if args.n_permutations > 0:
#        for i in range(args.n_permutations):
#            shuffled_cols = columns[:]
#            random.shuffle(shuffled_cols)
#            shuffled_df = df.select(shuffled_cols)
#            counts = compute_novel_counts(shuffled_df, mask)
#            permutation_novel_counts.append(counts)
#            if (i + 1) % 10 == 0 or i == 0:
#                print(f"Completed {i + 1} permutations")
#
    column_array = [df[col].to_list() for col in columns]

    print('multiprocessing')
    with ProcessPoolExecutor(
        max_workers=args.num_threads,
        initializer=init_worker,
        initargs=(column_array, mask, args.verbose)
    ) as executor:

        futures = [executor.submit(run_permutation, i) for i in range(args.n_permutations)]
    
        permutation_novel_counts = []
        for idx, future in enumerate(as_completed(futures), 1):
            result = future.result()
            permutation_novel_counts.append(result)
    
            if idx % 10 == 0 or idx == 1:
                print(f"Completed {idx} permutations")

    root, ext = os.path.splitext(args.output)

    if args.plot == 'novel':
        print(pl.Series("Novel counts:", novel_counts))
        plot_novels(novel_counts, path=root, ext=ext, lineage=args.lineage)
    elif args.plot == 'total':
        #total_counts = list(accumulate(sorted([v[1] for v in novel_counts.values()], reverse=True))) if args.sorted else list(accumulate(v[1] for v in novel_counts.values()))

        # === Prepare original data ===
        original_values = [v[1] for v in novel_counts.values()]
        if args.sorted:
            original_values = sorted(original_values, reverse=True)
        total_counts = list(accumulate(original_values))

        with ProcessPoolExecutor(max_workers=args.num_threads) as executor:
            futures = [
                executor.submit(process_permutation, perm, args.sorted)
                for perm in permutation_novel_counts
            ]
        
            permutation_counts = []
            for i, future in enumerate(as_completed(futures), 1):
                acc = future.result()
                permutation_counts.append(acc)
        
                if args.verbose and (i == 1 or i % 10 == 0):
                    print(f"[Main] Processed {i} permutations for accumulate")
        
        print(pl.Series("Total counts:", total_counts))
        print(pl.Series("permutations counts:", permutation_counts))

        per_np = np.array(permutation_counts)
        print(per_np.shape)

        alone = pangenome_saturation(novel_counts, total_counts)

        saturation_counts = []
        for i in range(len(permutation_counts)):
            n = pangenome_saturation(permutation_novel_counts[i], permutation_counts[i])
            saturation_counts.append(n)
    #    with ProcessPoolExecutor(
    #        max_workers=args.num_threads,
    #        initializer=init_saturation_worker,
    #        initargs=(permutation_novel_counts, permutation_counts)
    #    ) as executor:
    #        futures = [executor.submit(run_saturation, i) for i in range(len(permutation_counts))]
    #    
    #        saturation_counts = []
    #        for i, future in enumerate(as_completed(futures), 1):
    #            result = future.result()
    #            saturation_counts.append(result)
    #    
    #            if args.verbose and (i == 1 or i % 10 == 0):
    #                print(f"[Main] Processed {i} permutations for saturation")

        result = describe(saturation_counts)
        print(f"""The pangenome saturation metric:
                  {alone}

                  "Count:", {result.nobs}
                  "Min:", {result.minmax[0]}
                  "Max:", {result.minmax[1]}
                  "Mean:", {result.mean}
                  "Variance:", {result.variance}
                  "Skewness:", {result.skewness}
                  "Kurtosis:", result.kurtosis

               """
              )

        plot_totals(total_counts,
                    path=root,
                    ext=ext,
                    yaxis=args.yaxis,
                    lineage=args.lineage,
                    permutations=permutation_counts,
                    saturation=result.mean,
                   )
    else:
        sys.exit("Choose a plot argument type!")

if __name__ == '__main__':
    main()
