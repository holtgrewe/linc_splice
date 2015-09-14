#!/usr/bin/env python3

import argparse
import sys
import os.path

import gffutils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-name', type=str, default='transcripts.db')
    parser.add_argument('--gencode-gtf', type=str, required=True)

    args = parser.parse_args()

    db = gffutils.create_db(
            args.gencode_gtf, dbfn=args.db_name, force=True, keep_order=True,
            merge_strategy='merge', sort_attribute_values=True,
            disable_infer_genes=True, disable_infer_transcripts=True)


if __name__ == '__main__':
    sys.exit(main())
