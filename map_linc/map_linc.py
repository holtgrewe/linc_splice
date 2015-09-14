#!/usr/bin/env python3

import argparse
import gzip
import shlex
import sys
import os.path

from data import TranscriptInfo

import GTF
import pickle
import pysam
import tabix


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


GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end',
              'score', 'strand', 'frame']


class GTFFeature:
    """Represents a feature from a GTF file"""

    @staticmethod
    def parse(line=None, arr=None):
        if not arr:
            arr = line.rstrip().split('\t')
        attributes = {}
        for x in arr[-1].split(';'):
            tokens = shlex.split(x)
            if not tokens:
                pass
            elif len(tokens) == 1:
                attributes[tokens[0]] = True
            elif len(tokens) == 2:
                attributes[tokens[0]] = tokens[1]
            else:
                attributes[tokens[0]] = tokens[1:]
        arr[-1] = attributes
        return GTFFeature(*arr)

    def __init__(self, seqname, source, feature, start, end, score,
                 strand, frame, attrs):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = int(start) - 1
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attrs = attrs

    def __str__(self):
        return 'GTFFeature({})'.format(', '.join(map(repr, map(str, [
            self.seqname, self.source, self.feature, self.start + 1,
            self.end, self.score, self.strand, self.frame,
            self.attrs]))))


class App:
    """Code and state for the mapping application"""

    def __init__(self, args):
        self.args = args
        # parse out TransciptInfos
        print('Loading transcripts...', file=sys.stderr)
        self.tx_infos = self._parse_tx_infos(args.gencode_gtf)
        self.tx_info_by_id = dict([(info.transcript_id, info) for info in self.tx_infos])
        # open tabix file
        print('Opening tabix file...', file=sys.stderr)
        self.tabix = tabix.open(args.gencode_gtf)
        # open BAM file and iterate over it
        print('Opening BAM file...', file=sys.stderr)
        self.sam_file = pysam.AlignmentFile(args.alignment_bam, 'r')

    def _parse_tx_infos(self, gtf_path):
        """Parse transcript infos from GTF file or load from cache

        In case of successful loading from GTF, result will be cached.
        """
        if os.path.exists('_tx_cache.bin'):
            with open('_tx_cache.bin', 'rb') as f:
                return pickle.load(f)
        result = []
        with gzip.open(gtf_path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 1000 == 0:
                    print('processed {}'.format(i), file=sys.stderr)
                if line.startswith('#'):
                    continue
                if line.split('\t', 3)[2] != 'transcript':
                    continue
                record = GTFFeature.parse(line)
                if record.feature != 'transcript':
                    continue
                result.append(
                    TranscriptInfo(record.attrs['gene_id'],
                                   record.attrs['transcript_id'],
                                   record.attrs['transcript_type']))
        with open('_tx_cache.bin', 'wb') as g:
            pickle.dump(result, g)
            print(len(result), file=sys.stderr)
        return result


    def run(self):
        """Kick off the processing

        The SAM/BAM file is read in chunk-wise (with matching query names
        used for grouping).
        """
        for chunk in chunk_by_query(self.sam_file):
            #print('STARTING CHUNK', file=sys.stderr)
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
        linc_tx = self.tx_info_by_id[id_info.ensembl_tx_id]
        if 'lincRNA' not in linc_tx.transcript_type:
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
        #print('Querying for exons...', file=sys.stderr)
        try:
            for arr in self.tabix.query(*region):
                record = GTFFeature.parse(arr=arr)
                if record.feature != 'exon':
                    continue  # we look for overlapping transcripts
                if record.attrs['transcript_type'] != 'protein_coding':
                    continue  # we are only interested in these
                print(record)
        except tabix.TabixError as e:
            pass  # swallow, probably some unplaced region
        #for exon in self.db.region(region=region, featuretype=['exon']):
        #    transcript = self.db[exon.attributes['transcript_id'][0]]
        #    if match_strand == exon.strand:
        #        continue  # does not match signal we are looking for
        #    if transcript.attributes['transcript_type'] != 'protein_coding':
        #        continue  # does not code for a protein
        #    print('----\nMATCH\n{}\n----\nEXON\n{}\n----\nIN TRANSCRIPT\n{}\n'.format(match, exon, transcript),
        #          file=sys.stderr)


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--gencode-gtf', type=str, required=True)
    parser.add_argument('--alignment-bam', type=str, required=True)
    args = parser.parse_args()

    App(args).run()


if __name__ == '__main__':
    sys.exit(main())
