#! /usr/bin/env python

import argparse
import sys
import csv
from sourmash.picklist import SignaturePicklist
from sourmash import sourmash_args

def search_tax_file(filename, rank_initial: str, rank_name: str, count: bool = False) -> dict:
    """
    Find passport number and series by interating through text records
    :param filename:csv filename path
    :param rank_initial: Initial of tax rank (E.g. species == s)
    :param rank_name: Name of tax at rank(E.g. 'Escherichia Coli')
    :return:
    """
    pattern = rank_initial + "__" + rank_name
    matches = []

    with open(filename, 'r', encoding='utf_8_sig') as csvfile:
        if count:
            sys.exit(f"{pattern}: {sum(1 for line in csvfile if pattern in line)}")
        sniffer = csv.Sniffer()
        has_header = sniffer.has_header(csvfile.read(2048))
        csvfile.seek(0)
        if has_header:
            header = csvfile.readline().strip('\n')
        for line in csvfile:
            if pattern in line.strip():
                matches.append(line.strip())

    if matches:
        return header, matches
    else:
        sys.exit("Pattern not found in database")

def main():
    p = argparse.ArgumentParser()

    p.add_argument('pattern')
    p.add_argument('-i','--input')
    p.add_argument('-c','--count', action='store_true')

    p.add_argument('-po','--picklist-output') 

    args = p.parse_args()

    header, result = search_tax_file(args.input, "s", args.pattern, args.count)
    print(result)
    print(len(result))

    header, matches = search_tax_file(args.input, "s", args.pattern, args.count)

    with open(args.picklist_output, 'w', newline='') as fp:
        writer = csv.writer(fp)
        writer.writerow(header.split(','))
        for item in matches:
            writer.writerow(item.split(','))

    picklist = SignaturePicklist.from_picklist_args("my_picklist.csv:ident:ident")
    picklist.load()

    idx = sourmash_args.load_file_as_index(filename, yield_all_files=args.force)

    idx = idx.select(ksize=args.ksize, moltype=moltype, picklist=picklist)


if __name__ == '__main__':
    main()
