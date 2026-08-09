#! /usr/bin/env python

import argparse
import csv
from sourmash.picklist import SignaturePicklist

import sourmash_utils
from sourmash import sourmash_args
from search_example import search_tax_file
from sourmash_plugin_tables import Command_Hash_Tables


def main():
    p = argparse.ArgumentParser()

    p.add_argument('ranktable', help="Input csv containing classified hashes for specified organism")
    p.add_argument('sketches', nargs="+", help="Input file with sketches to process")
    p.add_argument('-o', '--output', required=True, help="Output CSV/Parquet file")
    p.add_argument('-v', '--verbose', action='store_true', help="Please flood my terminal with output. Thx.")

    p.add_argument('--pattern')
    p.add_argument('--taxonomy-csv')
    p.add_argument('-d','--database')
    p.add_argument('-c','--count', action='store_true')
    p.add_argument('--force', action='store_true')

    p.add_argument('-po','--picklist-output')
    sourmash_utils.add_standard_minhash_args(p)

    args = p.parse_args()

    moltype = sourmash_args.calculate_moltype(args)
    header, matches = search_tax_file(args.taxonomy_csv, "s", args.pattern, args.count)

    with open(args.picklist_output, 'w', newline='') as fp:
        writer = csv.writer(fp)
        writer.writerow(header.split(','))
        for item in matches:
            writer.writerow(item.split(','))

    picklist = SignaturePicklist.from_picklist_args(f"{args.picklist_output}:ident:ident")
    picklist.load()

    idx = sourmash_args.load_file_as_index(args.database, yield_all_files=args.force)

    idx = idx.select(ksize=args.ksize, moltype=moltype, picklist=picklist)
    
    #manifest = idx.manifest

    sub_idx = idx.select(picklist=picklist)

    args.sketches = sub_idx.signatures()
    Command_Hash_Tables.main()

    for ss in sub_idx.signatures():
        print(ss)

if __name__ == '__main__':
    main()
