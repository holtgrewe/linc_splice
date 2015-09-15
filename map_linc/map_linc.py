#!/usr/bin/env python3

import argparse
import gzip
import sys
import os.path

from data import TranscriptInfo, Region, GTFFeature

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

CIGAR_MAP = dict((x, i) for i, x in enumerate('MIDNSHP=X'))


def proc_cigar(cigar_s):
    result = []
    curr = []
    for x in cigar_s:
        if x.isdigit():
            curr.append(x)
        else:
            result.append((CIGAR_MAP[x], int(''.join(curr))))
            curr = []
    return tuple(result)


def expand_xa_tag(record, ali_file):
    result = [(record, int(dict(record.tags)['NM']))]
    tags = dict(record.tags)
    if not tags.get('XA'):
        return result
    for xa in tags['XA'].split(';'):
        if not xa:
            continue  # ignore last
        chrom, pos, cigar, nm = xa.split(',')
        strand = pos[0]
        pos = pos[1:]
        a = pysam.AlignedSegment()
        a.query_name = record.query_name
        a.query_sequence = record.query_sequence
        a.flag = (0 if strand == '+' else 16)
        a.reference_id = ali_file.gettid(chrom)
        a.reference_start = int(pos) - 1
        if strand == '-':
            a.cigar = list(reversed(proc_cigar(cigar)))
        else:
            a.cigar = proc_cigar(cigar)
        result.append((a, nm))
    return result


def chunk_by_query(ali_file, expand_xa=False):
    """Read in the pysam AlignmentFile ali_file chunk-wise

    Return iterator of lists of pysam alignment records.

    Actually, (record, nm) pairs as there is a problem with tags and
    Python 3 with pysam.
    """
    result = []
    for x in ali_file:
        if result and x.query_name != result[0][0].query_name:
            yield result
            if expand_xa:
                result = expand_xa_tag(x, ali_file)
            else:
                result = [(x, int(dict(x.tags)['NM']))]
        else:
            if expand_xa:
                result += expand_xa_tag(x, ali_file)
            else:
                result.append((x, int(dict(x.tags)['NM'])))
    if result:
        yield result


class OutputRecord:

    def __init__(self, match_region, exon, linc_tx, linc_length,
                 linc_ali_length, errors, window_5, window_3, classes=[]):
        #: Region of the match
        self.match_region = match_region
        #: GTFFeature of the target exon
        self.exon = exon
        #: GTFFeature of the lincRNA transcript
        self.linc_tx = linc_tx
        #: length of the lincRNA transcript
        self.linc_length = linc_length
        #: length of alignment in lincRNA
        self.linc_ali_length = linc_ali_length
        #: number of alignment errors
        self.errors = errors
        #: Window relative to 5' end
        self.window_5 = window_5
        #: Window relative to 3' end
        self.window_3 = window_3
        #: Classes
        self.classes = list(classes)

    def to_list(self):
        result = [
            self.linc_tx.transcript_id,
            self.linc_tx.gene_id,
            self.linc_tx.seqname,
            self.linc_tx.start + 1,
            self.linc_tx.end,
            self.linc_length,
            self.exon.attrs['transcript_id'],
            self.exon.attrs['gene_id'],
            self.exon.attrs['exon_id'],
            self.exon.attrs['exon_number'],
            self.exon.attrs['level'],
            self.exon.seqname,
            self.exon.start + 1,
            self.exon.end,
            self.match_region.seqname,
            self.match_region.start + 1,
            self.match_region.end,
            self.linc_ali_length,
            '{:.1f}'.format(
                100 * self.linc_ali_length / self.linc_length),
            '{:.1f}'.format(self.errors / self.linc_ali_length),
        ]
        result += list(self.window_5.to_tuple())
        result += list(self.window_3.to_tuple())
        result += ['&'.join(sorted(set(self.classes)))
                   if self.classes else '.']
        return result

    def to_tsv(self):
        return '\t'.join(map(str, self.to_list()))

    @staticmethod
    def get_header_fields():
        return ['linc_tx_id',
                'linc_gene_id',
                'linc_chrom',
                'linc_start',
                'linc_end',
                'linc_length',
                'target_tx_id',
                'target_gene_id',
                'target_exon_id',
                'target_exon_no',
                'target_tx_level',
                'target_exon_chrom',
                'target_exon_start',
                'target_exon_end',
                'match_chrom',
                'match_start',
                'match_end',
                'match_linc_length',
                'match_linc_coverage',
                'match_linc_erate',
                'exon_5_prime_window_start',
                'exon_5_prime_window_end',
                'exon_3_prime_window_start',
                'exon_3_prime_window_end',
                'classes',
                ]


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
                                   record.attrs['transcript_type'],
                                   record.seqname,
                                   record.start,
                                   record.end))
        with open('_tx_cache.bin', 'wb') as g:
            pickle.dump(result, g)
            print(len(result), file=sys.stderr)
        return result


    def run(self):
        """Kick off the processing

        The SAM/BAM file is read in chunk-wise (with matching query names
        used for grouping).
        """
        print('#' + '\t'.join(OutputRecord.get_header_fields()),
              file=self.args.output_tsv)
        for chunk in chunk_by_query(self.sam_file, expand_xa=True):
            #print('STARTING CHUNK', file=sys.stderr)
            if not chunk:
                continue  # ignore empty chunks
            self.handle_chunk(chunk)

    def handle_chunk(self, chunk):
        """Handle one chunk of lincRNA transcript alignments
        
        For each chunk, each match is considered in handle_match.
        """
        # get information for lincRNA transcript in chunk
        first = chunk[0][0]  # shortcut
        id_info = TranscriptIDInfo.parse(first.query_name)
        linc_tx = self.tx_info_by_id[id_info.ensembl_tx_id]
        if 'lincRNA' not in linc_tx.transcript_type:
            return  # only consider of is lincRNA transcript
        for match, nm in chunk:
            self.handle_match(linc_tx, match, nm)

    def handle_match(self, linc_tx, match, nm):
        """Handle one match of a lincRNA against the genome

        For each match, look at all overlapping exons and consider them
        as candidate lincRNA-to-coding gene interactions.
        """
        # look for exons overlapping with the lincRNA match
        match_strand = ('-' if match.flag & 16 else '+')
        region = Region(self.sam_file.getrname(match.reference_id),
                        match.pos, match.reference_end)
        #print('Querying for exons...', file=sys.stderr)
        try:
            first = True
            for arr in self.tabix.query(*region.to_tuple()):
                exon = GTFFeature.parse(arr=arr)
                if exon.feature != 'exon':
                    continue  # we look for overlapping transcripts
                if exon.attrs['transcript_type'] != 'protein_coding':
                    continue  # we are only interested in these
                if match_strand == exon.strand:
                    continue  # must be on different strands
                overlap_type = self.classify_overlap(
                    [region.start, region.end],
                    [exon.start, exon.end])
                if first:
                    # print('MATCH', match, file=sys.stderr)
                    first = False
                # print('TARGET', exon, file=sys.stderr)
                window_5 = exon.get_5_prime_window(region)
                window_3 = exon.get_3_prime_window(region)
                classes = self.compute_classes(exon, window_5, window_3)

                def ali_length(cigar):
                    OP = 'MIDNSHP=X'
                    return sum(num for op, num in cigar if OP[op] in 'MI=X')

                # generate and write out OutputRecord
                out = OutputRecord(
                        region, exon, linc_tx, len(match.query_sequence),
                        ali_length(match.cigar), int(nm),
                        window_5, window_3, classes)
                print(out.to_tsv(), file=self.args.output_tsv)
        except tabix.TabixError as e:
            pass  # swallow, probably some unplaced region

    def compute_classes(self, exon, window_5, window_3):
        classes = []
        if window_5.is_overlapping:
            classes.append('retain_intron_5')
        if window_3.is_overlapping:
            classes.append('retain_intron_3')
        if window_5.is_overlapping and window_3.is_overlapping:
            classes.append('skip_exon')
        return classes

    def classify_overlap(self, region1, region2):
        (s1, e1) = region1
        (s2, e2) = region2
        
        if s1 < s2 and e1 >= s2 and e1 <= e2:
            return "lincRNA partially overlaps exon's upstream"
        elif e1 > e2 and s1 >= s2 and s1 <= e2:
            return "lincRNA partially overlaps exon's downstream"
        elif s1 <= s2 and e1>=e2:
            return "exon is flanked by lincrna"
        elif s1 >= s2 and e1<=e2:       
            return "exon flanks the lincrna"


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--gencode-gtf', type=str, required=True)
    parser.add_argument('--alignment-bam', type=str, required=True)
    parser.add_argument('--output-tsv', type=argparse.FileType('wt'),
                        default=sys.stdout)
    args = parser.parse_args()

    App(args).run()


if __name__ == '__main__':
    sys.exit(main())
