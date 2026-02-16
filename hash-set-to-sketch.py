#!/usr/bin/env python

import argparse
import gzip
import os
import csv
from pathlib import Path

import sourmash
from sourmash import sourmash_args
import sourmash_utils

def hash_generator(filehandle, block_size=10000):
    block = []
    line_count = 0
    for line in filehandle:
        block.append(line)
        line_count += 1
        if line_count % (block_size * 10) == 0:
            print(f"Processed {line_count:,} lines...")
        if len(block) == block_size:
            yield block
            block = []
    if block:
        yield block
        print(f"Processed {line_count:,} lines... (final block)")

def main():
    parser = argparse.ArgumentParser(description="Filter and hash input data into a sourmash signature.")
    parser.add_argument('input', help='Input gzip-compressed file with hash values (tab-separated).')
    parser.add_argument('-o', '--output', required=True, help='Output .sig file path.')
    parser.add_argument('-b', '--filter-by', type=int, help='Column index (0-based) to filter by.')
    parser.add_argument('-f', '--filter-amount', type=int, help='Value to filter for in selected column. Set to greater than the value given.')
    parser.add_argument('--abund', action='store_true', help='Track and store abundance values.')

    sourmash_utils.add_standard_minhash_args(parser)
    args = parser.parse_args()

    if os.path.exists(args.output):
        raise FileExistsError(f"Output file '{args.output}' already exists. Aborting.")

    print(f"Creating base MinHash with parameters from arguments...")
    base_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"Base MinHash created: {base_mh}")

    if args.abund:
        base_mh.track_abundance = True
        abundance_counts = {}
    else:
        abundance_counts = None

    print(f"Reading from '{args.input}'...")

    open_func = gzip.open if args.input.endswith('.gz') else open

    with open_func(args.input, 'rt', newline='') as fp:
        for block in hash_generator(fp):
            fields = csv.reader(block, delimiter=',')
            for line in fields:
                if not line or not line[0].isdigit():
                    continue  # skip invalid lines

                # Filter
                if args.filter_by is not None and args.filter_amount is not None:
                    if not line[args.filter_by].isdigit() or int(line[args.filter_by]) <= args.filter_amount:
                        continue

                hashval = int(line[0])
                if args.abund:
                    abundance_counts[hashval] = abundance_counts.get(hashval, 0) + 1
                else:
                    base_mh.add_hash(hashval)

    # If abundance mode, apply counts
    if args.abund:
        print("Setting abundances in MinHash...")
        mh = base_mh.copy_and_clear()
        mh.track_abundance = True
        mh.set_abundances(abundance_counts)
    else:
        mh = base_mh

    print(f'Collected {len(mh)} hashes')

    path = Path(args.input)
    if path.suffix == '.gz':
        sig_name = path.with_suffix('').stem  # Remove .gz then .csv
    else:
        sig_name = path.stem
    signature = sourmash.SourmashSignature(mh, name=sig_name)

    print(f"Saving signature to '{args.output}'...")
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        save_sigs.add(signature)

    print("Done.")

if __name__ == '__main__':
    main()
