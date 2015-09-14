#!/usr/bin/env python3

class TranscriptInfo:
    """Information about a transcript"""

    def __init__(self, gene_id, transcript_id, transcript_type):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type
