#! /usr/bin/env python

import argparse
import sys

import sourmash
from sourmash import sourmash_args
#import sourmash_utils

def main():
    p = argparse.ArgumentParser()

    p.add_argument('db')
    p.add_argument(
        "--picklist",
        metavar="FILE:COL:MODE",
        help=(
            "Picklist in the form 'filename:column_name:include|exclude'. "
            "Example: idents.txt:ident:include"
        )
    )
    p.add_argument("--moltype", default="DNA", type=str, help="Path to output file.")
    p.add_argument("--ksize", default=31, type=int, help="Path to output file.")
    p.add_argument("--scaled", default=1000, type=int, help="Path to output file.")

    args = p.parse_args()

    print(f"Loading sketches from: {args.db}")
    db = sourmash.load_file_as_index(args.db)
    db = db.select(moltype=args.moltype, ksize=args.ksize, scaled=args.scaled)
    print(f"Loaded {len(db)} sketches from {args.db}")

    # --- Load picklist ---
    pl = sourmash_args.load_picklist(args)

    if pl:
        db = db.select(picklist=pl)
        print(f"Subset contains {len(db)} sketches after applying picklist.")
    return db

if __name__ == "__main__":
    main()
