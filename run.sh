#!/bin/bash
#$ -cwd -N CD4_C1 -P jknight.prjc -q short.qc
#$ -j y -b n 
#$ -pe shmem 12
#$ -o /users/jknight/kwz374/R/AS.patinets.splicing
#$ -e /users/jknight/kwz374/R/AS.patinets.splicing

echo "job started on" `date +"%T %Y-%m-%d"`
module load R/3.6.2-foss-2019b

Rscript splicing_CD4.r fcount.DEXSeq_AS.patients_CD4_Cohort1.and.2.txt

echo "job finished on" `date +"%T %Y-%m-%d"`
