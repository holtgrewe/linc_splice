#!/usr/bin/env python3

import shlex


class TranscriptInfo:
    """Information about a transcript"""

    def __init__(self, gene_id, transcript_id, transcript_type,
                 seqname, start, end):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type
        self.seqname = seqname
        self.start = start
        self.end = end

    def region_str(self, show_len=False):
        length = ''
        if show_len:
            length = '(len:{})'.format(self.end - self.start)
        return '{}:{}-{}{}'.format(self.seqname, self.start + 1,
                                   self.end, length)


class Region:

    def __init__(self, seqname, start, end):
        self.seqname = seqname
        self.start = start
        self.end = end

    def length(self):
        return self.end - self.start

    def to_tuple(self):
        return (self.seqname, self.start, self.end)

    def region_str(self, show_len=False):
        length = ''
        if show_len:
            length = '(len:{})'.format(self.end - self.start)
        return '{}:{}-{}{}'.format(self.seqname, self.start + 1,
                                   self.end, length)


GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end',
              'score', 'strand', 'frame']


class Window:
    """Window around an end of a GTFFeature"""

    def __init__(self, side, left, right):
        self.side = side
        self.left = left
        self.right = right

    def to_tuple(self):
        """Return (left, right)"""
        return (self.left, self.right)

    @property
    def is_overlapping(self):
        if self.left == 0 and self.right == 0:
            return False
        else:
            return (self.left > 0) != (self.right > 0)

    def __str__(self):
        return '({})[{}, {}]'.format(self.side, self.left, self.right)


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

    def length(self):
        return self.end - self.start

    def get_5_prime_window(self, region):
        """Return Window relative to 5' end of self"""
        return Window('5_prime', region.start - self.start,
                      region.end - self.start)

    def get_3_prime_window(self, region):
        """Return Window relative to 3' end of self"""
        return Window('3_prime', region.start - self.end,
                      region.end - self.end)

    def region_str(self, show_len=False):
        length = ''
        if show_len:
            length = '(len:{})'.format(self.end - self.start)
        return '{}:{}-{}{}'.format(self.seqname, self.start + 1,
                                   self.end, length)

    def __str__(self):
        return 'GTFFeature({})'.format(', '.join(map(repr, map(str, [
            self.seqname, self.source, self.feature, self.start + 1,
            self.end, self.score, self.strand, self.frame,
            self.attrs]))))


