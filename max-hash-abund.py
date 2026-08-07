#! /usr/bin/env python

import sourmash
import argparse
import csv
from collections import Counter

p = argparse.ArgumentParser()
p.add_argument("signature")
p.add_argument("--taxa-file")
p.add_argument("-r", "--rank")
p.add_argument("-t","--target")
args = p.parse_args()

if args.signature == "-":
    sigs = sourmash.load_file_as_signatures("-", ksize=31)
else:
    sigs = sourmash.load_file_as_signatures(args.signature)

for sig in sigs:
    mh = sig.minhash
    print("track_abundance:", mh.track_abundance)
    hashval = max(mh.hashes)
    abund = mh.hashes[hashval]
    print(abund)
    print(sig.name)
    print("max hash:", hashval)
    print("abundance:", abund)
    max_abund = max(mh.hashes.values())
    print("Max abundance:", max_abund)

counts = Counter()

with open(args.taxa_file, newline="") as fp:
    reader = csv.DictReader(fp)

    for row in reader:
        c = row[args.rank]
        if c == args.target:
            counts[c] += 1
print(f"Total {args.target} genomes found:", counts[args.target])
