# Downloaded Files
* LincRNA - fasta: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
* LincRNA - GTF: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz

# File Formats

## Mapping LINC to Genome

### LINC Information

1. LINC ENSEMBL transcript ID
2. LINC ENSEMBL gene ID of transcript
3. match, chromosome
4. match, begin position on chromosome
5. match, end position on chromosome

### Target Information

6. Target transcript ID or NA
7. Target gene ID or NA
10. match, chromosome
11. match, begin position on chromosome
12. match, end position on chromosome

### Candidate Events

12. comma-separated list of 1-based number of exon that is possibly skipped
13. comma-separated list of 1-based number of intron where bases might be retained or "0" for none
