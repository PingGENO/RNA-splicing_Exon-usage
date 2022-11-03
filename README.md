# RNA_splicing-Exon_usage

## DEXseq
### 1. Generate gtf file
#### Download the GTF file
**Link** : ftp://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

`dexseq_prepare_annotation2` was downloaded from [here](https://pypi.org/project/dexseq_prepare_annotation2/#description)

```
python python_scripts/dexseq_prepare_annotation2.py -r no -f Homo_sapiens.GRCh37.87_DEXSeq.gtf  Homo_sapiens.GRCh37.87.gtf.gz Homo_sapiens.GRCh37.87_DEXSeq.gff.
```

### 2. Generate exonic_part count file
```
GTF=Homo_sapiens.GRCh37.87_DEXSeq.gtf
/apps/htseq/subread/bin/featureCounts -p -f -O -s 2 -F GTF -a $GTF \
-t exonic_part -o fcount.DEXSeq.no.r.txt \
*.bam.
```

### 3. Differential exon usage analysis via DEXseq
```
qusb run.sh
```



