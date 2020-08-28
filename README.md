# RNA_splicing


# 1. Generate the annotation file 
$ python DEXSeq.splicing/dexseq_prepare_annotation2_r_no.py -f Homo_sapiens.GRCh38.84_DEXSeq.counts.gtf /Users/pingzhang/lec03_Analysis/RNA-seq_analysis/refs/Homo_sapiens.GRCh38.84.gtf.gz Homo_sapiens.GRCh38.84_DEXSeq.counts.gff

# 2. Generate the featureCounts file
$ GTF=/well/jknight/ping/gtfs/Homo_sapiens.GRCh38.84_DEXSeq.counts.gtf
$ /apps/htseq/subread/bin/featureCounts -p -f -O -s 2 -F GTF -a $GTF \
-t exonic_part -o fcount.DEXSeq.no.r.txt \
/well/jknight/ping/RNA-seqs/bam_hisat2_macrophage_hg38/*merge.bam


