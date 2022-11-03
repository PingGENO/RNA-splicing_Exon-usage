#!/bin/bash
#$ -cwd -N exon -P {group}.prjc -q short.qc
#$ -j y -b n 
#$ -pe shmem 12

# Make the script stop on any error.
set -uxeo pipefail

echo "job started on" `date +"%T %Y-%m-%d"`
module load R/3.6.2-foss-2019b

Rscript splicing_CD4.r fcount.DEXSeq.no.r.txt

echo "job finished on" `date +"%T %Y-%m-%d"`
