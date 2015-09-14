#!/usr/bin/env python3

import argparse
import sys
import os.path

import gffutils
import pysam


class TranscriptIDInfo:
    """Stores transcript ID information"""

    @staticmethod
    def parse(tx_meta):
        """Parse transcript ID information from FASTA meta line"""
        return TranscriptIDInfo(*tx_meta.split('|'))

    def __init__(self, ensembl_tx_id, ensembl_gene_id, *args):
        """Initialize object with the given IDs"""
        self.ensembl_tx_id = ensembl_tx_id
        self.ensembl_gene_id = ensembl_gene_id
        self.rest = args


def chunk_by_query(ali_file):
    """Read in the pysam AlignmentFile ali_file chunk-wise

    Return iterator of lists of pysam alignment records.
    """
    result = []
    for x in ali_file:
        if result and x.query_name != result[0].query_name:
            yield result
            result = [x]
        else:
            result.append(x)
    if result:
        yield result



class App:
    """Code and state for the mapping application"""

    def __init__(self, args):
        self.args = args
        # load database and open BAM file
        self.db = gffutils.FeatureDB(args.db_name, keep_order=True)
        # open BAM file and iterate over it
        self.sam_file = pysam.AlignmentFile(args.alignment_bam, 'r')

    def run(self):
        """Kick off the processing

        The SAM/BAM file is read in chunk-wise (with matching query names
        used for grouping).
        """
        for chunk in chunk_by_query(self.sam_file):
            print('STARTING CHUNK', file=sys.stderr)
            if not chunk:
                continue  # ignore empty chunks
            self.handle_chunk(chunk)

    def handle_chunk(self, chunk):
        """Handle one chunk of lincRNA transcript alignments
        
        For each chunk, each match is considered in handle_match.
        """
        # get information for lincRNA transcript in chunk
        first = chunk[0]  # shortcut
        id_info = TranscriptIDInfo.parse(first.query_name)
        linc_tx = self.db[id_info.ensembl_tx_id]
        if 'lincRNA' not in linc_tx['gene_type']:
            return  # only consider of is lincRNA transcript
        for match in chunk:
            self.handle_match(match)

    def handle_match(self, match):
        """Handle one match of a lincRNA against the genome

        For each match, look at all overlapping exons and consider them
        as candidate lincRNA-to-coding gene interactions.
        """
        # look for exons overlapping with the lincRNA match
        match_strand = ('-' if match.flag & 16 else '+')
        region = (self.sam_file.getrname(match.reference_id),
                  match.pos, match.reference_end)
        print('Querying for exons...', file=sys.stderr)
        for exon in self.db.region(region=region, featuretype=['exon']):
            if match_strand == exon.strand:
                continue  # does not match signal we are looking for
            print('----\nMATCH\n{}\n----\nEXON\n{}\n----\nIN TRANSCRIPT\n{}\n'.format(match, exon, self.db[exon.attributes['transcript_id'][0]]),
                  file=sys.stderr)


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-name', type=str, default='transcripts.db')
    parser.add_argument('--alignment-bam', type=str, required=True)
    args = parser.parse_args()

    App(args).run()


if __name__ == '__main__':
    sys.exit(main())
