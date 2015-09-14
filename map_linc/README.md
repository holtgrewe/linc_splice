# LINC mapping scripts

## Setup environment

```
$ virtualenv -p python3 _venv
$ source _venv/bin/activate
$ pip install -r requirements.txt
```

## Build database

```
$ ./init_db.py --gencode-gtf ../gencode.v19.annotation.gtf.gz
```
