
# This program can be run as
# Rscript --vanilla xxxx.r input.txt out.txt
# Produces a table with AS genes. on the standard output.
# To install the requirements run the program with the 'install` parameter.

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}


library(dplyr)
library(data.table)
library(DEXSeq)

# Take a fcount file and convert it to dcounts for dexseq
message("Reading and adding Exon IDs for DEXSeq")

 #--------------Cohort1+2--------#
 dcounts <- as.data.frame(fread(args[1]))
 dcounts = dcounts[,c(1,7:57)]
 
 # clean the columns
 colClean <- function(x){ colnames(x) <- gsub("/well/jknight/ping/Ankylosing.spondylitis.RNAseq/CD8.bam.files/", "", colnames(x)); x } 
 dcounts = colClean(dcounts)
 colClean <- function(x){ colnames(x) <- gsub(".nodup.bam", "", colnames(x)); x } 
 dcounts = colClean(dcounts)
 
 
 colData <- data.frame(fread("CD8_Cases_CD8_Controls.colnames.txt"))
 #
 colData2 <- colData[colData$type %in% colnames(dcounts),]
 #
 sampleData = colData2[ order(match(colData2$type, colnames(dcounts))), ]
 
 sampleData$condition = as.factor(sampleData$condition)
 sampleData$condition <- relevel(sampleData$condition, ref = "controls")
 
 

#-----------------------------------------------#
# #keep genes  with counts >=10 in more than %30 samples
x = dcounts %>%
 group_by(Geneid) %>%
 summarise_all(funs(sum))
x2 <- x[rowSums(x[,-1] >=10) >= as.integer(ncol(x[,-1])*0.3)+1,]
dcounts = dcounts[dcounts$Geneid %in% x2$Geneid,]


#####
id <- as.character(dcounts[,1])
n <- id
split(n,id) <- lapply(split(n ,id), seq_along )
rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
dcounts <- dcounts[,2:ncol(dcounts)]

dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name 

## get genes and exon names out
splitted <- strsplit(rownames(dcounts), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)


## parse the flattened file
aggregates <- read.delim("/well/jknight/ping/gtfs/Homo_sapiens.GRCh37.87_DEXSeq.gtf", stringsAsFactors = FALSE, 
                          header = FALSE)
colnames(aggregates) <- c("chr", "source", "class", "start", 
                          "end", "ex", "strand", "ex2", "attr")
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
aggregates <- aggregates[which(aggregates$class == "exonic_part"),] #########
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
                          aggregates$attr)
# trim the gene_ids to 255 chars in order to match with featurecounts
longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
#warning(paste0(longIDs, 
#             " aggregate geneIDs were found truncated in featureCounts output"), 
#       call. = FALSE)
aggregates$gene_id <- substr(aggregates$gene_id,1,255)

transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
                    aggregates$attr)
transcripts <- strsplit(transcripts, "\\+")
exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", # exonic_part_number
                aggregates$attr)
exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
                                                          end = aggregates$end), strand = aggregates$strand)
names(exoninfo) <- paste(aggregates$gene_id, exonids, 
                         sep = ":E")

names(transcripts) <- names(exoninfo) 
#  if (!all(rownames(dcounts) %in% names(exoninfo))) {
#     stop("Count files do not correspond to the flattened annotation file")
#  }
matching <- match(rownames(dcounts), names(exoninfo))
stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))

######
dxd <- DEXSeqDataSet(dcounts, sampleData, design = ~sample + exon + condition:exon, exons, 
                     genesrle, exoninfo[matching], transcripts[matching])

sampleAnnotation(dxd)

# ######  the wrapper function DEXseq(), which runs the analysis shown above through a single function call.
# dxr = DEXSeq(dxd)

#---------------------------------------------
#---------------------------------------------
# only detect differences in exon usage that can be attributed to condition, independent of batch
formulaFullModel    =  ~ sample + exon + batch:exon + condition:exon
formulaReducedModel =  ~ sample + exon + batch:exon 

dxd = estimateSizeFactors( dxd )
message("1##Estimating dispersions...")
dxd = estimateDispersions( dxd, formula = formulaFullModel )
message("2##Running test...")
dxd = testForDEU( dxd, 
                  reducedModel = formulaReducedModel, 
                  fullModel = formulaFullModel )
message("3##estimate relative exon usage fold changes...")
dxd = estimateExonFoldChanges( dxd,  fitExpToVar="condition")
message("4##DEXSeqResults...")
res <- DEXSeqResults(dxd)
#---------------------------------------------
message("5##significant DEU cases...")
table( res$padj < 0.05 )
#message("generate HTML report...")
#DEXSeqHTML( res, FDR=0.05, color=c("#FF000080", "#0000FF80") )
#---------------------------------------------#

message("6##Save results")
results = as.data.frame(res)
write.csv(results, "DEXSeq.results.CD8_cohort1.2.csv")
  

