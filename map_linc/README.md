# LINC mapping scripts

## Setup environment

```
$ virtualenv -p python3 _venv
$ source _venv/bin/activate
$ pip install -r requirements.txt
```

## Build database

```
$ ./init_db.py --gencode-gtf ../data/gencode.v19.annotation.gtf.gz
```

## Mapping lincRNA to gene model

```
$ ./map_linc.py --alignment-bam ../align_linc/alignments.bam
```
