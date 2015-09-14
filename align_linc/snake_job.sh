#!/bin/bash

#$ -j y
#$ -o /scratch/biocore/mholtgr/projects/2015-09-14_hackathon/linc_splice/align_linc/sge_log

DIR=/scratch/biocore/mholtgr/projects/2015-09-14_hackathon/linc_splice/align_linc

module load Python/3.4.3-foss-2015a
#source ${HOME}/local/venv/bin/activate
source /scratch/biocore/mholtgr/venv/bin/activate

set -x

cd ${DIR}
mkdir -p sge_log

ls /home/mholtgr/Development/cubi_wrappers/cubi/snake/ngs_mapping_config.json
ls /home/mholtgr/local/venv/bin/python3
sleep 5  # automounter, yay!

snakemake \
    --cluster "qsub -V -j yes -o ${PWD}/sge_log" \
    -j \
    -k \
    -p
