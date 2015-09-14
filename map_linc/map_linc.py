#!/usr/bin/env python3

import argparse
import sys
import os.path

import gffutils
import pysam

def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-name', type=str, required=True)
    parser.add_argument('--alignment-bam', type=str, required=True)
    args = parser.parse_args()

    # load database
    db = gffutils.FeatureDB('test.db', keep_order=True)

    # open BAM file
    sam_file = pysam.AlignmentFile(args.alignment_bam, 'r')
    print(sam_file, file=sys.stderr)


if __name__ == '__main__':
    sys.exit(main())
