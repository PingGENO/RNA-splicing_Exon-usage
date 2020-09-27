# RNA_splicing

## DEXseq
### 1. Generate the exon count file
#### Download the GTF file
**Link** : ftp://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
$ python python_scripts/dexseq_prepare_annotation2.py -r no -f Homo_sapiens.GRCh37.87_DEXSeq.gtf  Homo_sapiens.GRCh37.87.gtf.gz Homo_sapiens.GRCh37.87_DEXSeq.gff

### 2. Generate the featureCounts file
$ GTF=/well/jknight/ping/gtfs/Homo_sapiens.GRCh38.84_DEXSeq.counts.gtf
$ /apps/htseq/subread/bin/featureCounts -p -f -O -s 2 -F GTF -a $GTF \
-t exonic_part -o fcount.DEXSeq.no.r.txt \
*.bam

### 3. Differential exon usage analysis via DEXseq
qusb 



